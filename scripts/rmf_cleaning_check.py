#!/usr/bin/env python3
"""Is the divergence-cleaning potential eating the answer?

This is the gate on every long RMF run. `twoFluidSimplifiedRMFBC` writes the
applied transverse field into the ghost state at every boundary face, against
whatever interior field the run has produced. The jump in the normal component
of B that this creates is what the cleaning potential psi (component 17) exists
to carry away, and on the shipped edge-driven deck it is fed faster than it can
carry: rms|psi|/rms|E| passes 12 within 210 steps and is still climbing. psi
enters the induction equation in the slot E occupies, so that is not a
bookkeeping quantity growing off to one side. See docs/known-issues.md section
13.

Two things are measured, and the second is the one that matters:

  * psi/E over time in one run. A number that grows without bound is the
    symptom.
  * the difference between two runs that differ ONLY in DIVB_SPEED. That is the
    disease: if the answer moves when the cleaning speed moves, the answer is
    partly the cleaning scheme. On the shipped deck B_x differs by 92% of its
    own rms within 210 steps.

Usage:

    # one run: is psi growing?
    python3 scripts/rmf_cleaning_check.py run_dir/

    # two runs differing only in DIVB_SPEED: does the answer depend on it?
    python3 scripts/rmf_cleaning_check.py cleaning_on/ cleaning_off/

Run this BEFORE spending days on a long run, on the same deck at a few hundred
steps. It costs about a minute and it decides whether the long run is worth
starting.
"""

import argparse
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'scripts'))
sys.path.insert(0, os.path.join(REPO, 'test'))

try:
    import numpy as np
except ImportError:
    np = None

# How far apart two runs may drift before the answer is the scheme's rather than
# the plasma's. 5% of rms|B| is generous: it is already far above discretisation
# error for two runs of identical resolution on identical input.
TOLERANCE = 0.05


def frames_in(directory):
    """[(index, path)] for the .vtu files in a run directory, in time order."""
    names = [f for f in os.listdir(directory) if f.endswith('.vtu')]
    if not names:
        raise SystemExit(f'{directory}: no .vtu frames. Did the run finish?')
    return sorted(((int(n.rsplit('_', 1)[1].split('.')[0]),
                    os.path.join(directory, n)) for n in names))


def load(path, components):
    import vtu
    import deckrun
    return deckrun.nodal(vtu.read(path), components)


def psi_growth(directory):
    """rms|psi|/rms|E| per frame. Returns [(index, rms_psi, rms_E, ratio)]."""
    import rmf_diagnostics as rd
    rows = []
    for index, path in frames_in(directory):
        _x, _y, c = load(path, (rd.E_X, rd.E_Y, rd.E_Z, rd.PHI, rd.PSI))
        e = np.sqrt(c[rd.E_X] ** 2 + c[rd.E_Y] ** 2 + c[rd.E_Z] ** 2)
        rms_e = float(np.sqrt(np.mean(e ** 2)))
        rms_psi = float(np.sqrt(np.mean(c[rd.PSI] ** 2)))
        rms_phi = float(np.sqrt(np.mean(c[rd.PHI] ** 2)))
        rows.append((index, rms_psi, rms_phi, rms_e,
                     rms_psi / rms_e if rms_e > 0 else float('nan')))
    return rows


def compare(dir_a, dir_b):
    """Per-frame difference in B between two runs. Returns [(index, dBx, dBz)]."""
    import rmf_diagnostics as rd
    fa, fb = frames_in(dir_a), frames_in(dir_b)
    common = sorted(set(i for i, _ in fa) & set(i for i, _ in fb))
    if not common:
        raise SystemExit('the two runs share no frame indices')
    pa, pb = dict(fa), dict(fb)
    rows = []
    for i in common:
        _x, _y, ca = load(pa[i], (rd.B_X, rd.B_Y, rd.B_Z))
        _x2, _y2, cb = load(pb[i], (rd.B_X, rd.B_Y, rd.B_Z))
        if ca[rd.B_X].shape != cb[rd.B_X].shape:
            raise SystemExit(f'frame {i}: the two runs have different meshes '
                             f'({ca[rd.B_X].size} against {cb[rd.B_X].size} nodes)')
        rms_bx = float(np.sqrt(np.mean(ca[rd.B_X] ** 2)))
        d_bx = float(np.sqrt(np.mean((ca[rd.B_X] - cb[rd.B_X]) ** 2)))
        d_bz = float(np.sqrt(np.mean((ca[rd.B_Z] - cb[rd.B_Z]) ** 2)))
        rows.append((i, d_bx / rms_bx if rms_bx > 0 else float('nan'),
                     d_bz / 60e-4))
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('runs', nargs='+',
                    help='one run directory, or two differing only in DIVB_SPEED')
    args = ap.parse_args(argv)
    if np is None:
        raise SystemExit('numpy is required')
    if len(args.runs) > 2:
        raise SystemExit('one or two run directories, not more')

    print(f'psi growth in {args.runs[0]}')
    print(f'  {"frame":<8}{"rms|psi|":<14}{"rms|phi|":<14}{"rms|E|":<14}{"psi/E":<10}')
    rows = psi_growth(args.runs[0])
    for index, psi, phi, e, ratio in rows:
        print(f'  {index:<8}{psi:<14.4e}{phi:<14.4e}{e:<14.4e}{ratio:<10.3f}')
    worst = max((r[4] for r in rows if r[4] == r[4]), default=float('nan'))
    print(f'\n  worst psi/E = {worst:.3f}')
    if worst > 1.0:
        print('  psi is LARGER than the electric field it shares the induction '
              'equation with.\n  Compare two cleaning speeds before trusting a '
              'long run on this deck.')
    else:
        print('  psi stays below the electric field. That is necessary, not '
              'sufficient:\n  only the two-run comparison shows whether the '
              'answer depends on the scheme.')

    if len(args.runs) == 1:
        return 0

    print(f'\nsensitivity: {args.runs[0]} against {args.runs[1]}')
    print(f'  {"frame":<8}{"rms dB_x / rms B_x":<22}{"rms dB_z / B_bias":<20}')
    diffs = compare(*args.runs)
    for index, dbx, dbz in diffs:
        print(f'  {index:<8}{dbx:<22.4f}{dbz:<20.6f}')

    worst_bx = max(d[1] for d in diffs)
    print()
    if worst_bx > TOLERANCE:
        print(f'  FAIL: B_x differs by {worst_bx:.1%} of its own rms between two '
              f'runs that\n        differ only in the divergence-cleaning speed. '
              f'Over this interval the\n        transverse field is partly a '
              f'property of the cleaning scheme, so a\n        longer run on this '
              f'deck would not measure the plasma. Use the antenna\n        deck '
              f'(see docs/known-issues.md section 13).')
        return 1
    print(f'  PASS: B_x differs by {worst_bx:.1%}, within the {TOLERANCE:.0%} '
          f'tolerance.\n        Note what this does and does not say: it covers '
          f'the interval these\n        frames span, and the divergence is linear '
          f'in time on the deck where it\n        fails, so a pass here does not '
          f'extrapolate to a run ten times longer.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
