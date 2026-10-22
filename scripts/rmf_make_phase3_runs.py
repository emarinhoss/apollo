#!/usr/bin/env python3
"""Generate the ready-to-run input folders under rmf_frc/phase3/.

Phase 3 of docs/rmf-frc-model-assessment.md is a set of runs nobody has the
compute for on a laptop. This writes each of them out as a self-contained
folder - decks plus mesh, nothing to resolve at run time - so they can be copied
to a cluster or a Colab session and started.

The decks are GENERATED rather than hand-edited, by the same `set_param` and
`apply_speedup` that `scripts/rmf_scan.py` uses and `test/test_rmf_scan.py`
checks. That matters for one specific reason: Apollo's parser reads a numeric
literal with no decimal point as an integer and then aborts with `std::bad_cast`
naming nothing (docs/known-issues.md section 12). Hand-editing a deck to
`LIGHT = 3000000` is a run that dies at setup.

    python3 scripts/rmf_make_phase3_runs.py            # write the folders
    python3 scripts/rmf_make_phase3_runs.py --check    # verify they are current

Re-run it after changing a base deck; the folders are checked in, so `--check`
is the thing to put in CI if that ever matters.
"""

import argparse
import io
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'scripts'))

import rmf_scan                      # noqa: E402
import rmf_diagnostics as rd         # noqa: E402

RMF = os.path.join(REPO, 'examples', 'unstructuredDG', 'multifluid', 'rmf_frc')
OUT = os.path.join(RMF, 'phase3')

EDGE_DECK = os.path.join(RMF, 'frc2d.pin')
EDGE_MESH = os.path.join(RMF, 'optimizedCircle2.msh')
ANT_DECK = os.path.join(RMF, 'antenna', 'frc2d.pin')
ANT_MESH = os.path.join(RMF, 'antenna', 'disc.msh')
MKDISCMESH = os.path.join(REPO, 'scripts', 'mkdiscmesh.py')

# 3 v_Te at T_e = 30 eV, with v_Te = sqrt(kT_e/m_e): the bottom of the range the
# study section 3.10 cites actually covers. See the assessment.
ME = 1.67e-27 / 1836.0
C_3VTE = 3.0 * math.sqrt(30.0 * 1.6e-19 / ME)

# Frames per RMF period. rmf_diagnostics.cycle_average refuses to call anything
# with fewer than 8 samples in a period a cycle average, and the profiles are
# noisier than the scalars, so 12 leaves margin.
FRAMES_PER_PERIOD = 12


def gauss(t):
    """A field in tesla as a filename-safe number of gauss."""
    return ('%.6g' % (t * 1e4)).replace('.', 'p') + 'G'


# Meshes generated per deck, as mkdiscmesh arguments. A mesh named here is
# written into the run folder and checked byte for byte by --check, so the
# folder stays self-contained without anyone hand-copying a mesh into it.
_DISC = ('--outer', '0.05', '--h-outer', '1.5e-3')
WINDING_MESHES = {
    'disc.msh': dict(file='disc.msh', args=_DISC + ('--ring', '0.03', '--h-inner', '9.3e-4')),
    'ring40.msh': dict(file='ring40.msh', args=_DISC + ('--ring', '0.040', '--h-inner', '9.3e-4')),
    'fine40.msh': dict(file='fine40.msh', args=_DISC + ('--ring', '0.040', '--h-inner', '6.0e-4')),
    'resolved40.msh': dict(file='resolved40.msh',
                           args=_DISC + ('--ring', '0.040', '--h-inner', '3.5e-4')),
    'rotated60.msh': dict(file='rotated60.msh',
                          args=_DISC + ('--ring', '0.03', '--h-inner', '9.3e-4', '--rotate', '60')),
}


RUNS = [
    dict(
        folder='00-divergence-gate',
        deck=EDGE_DECK, mesh=EDGE_MESH,
        periods=0.00435,          # about 210 steps
        frames=6,
        decks=[
            ('edge-cleaning-on', [('DIVB_SPEED', 1.0), ('DIVE_SPEED', 1.0)]),
            ('edge-cleaning-off', [('DIVB_SPEED', 0.0), ('DIVE_SPEED', 0.0)]),
        ],
    ),
    # The long gate. The short pair below spans 0.0040 RMF periods, which is
    # 0.17 of this deck's RISE time (3.0e-8 s against a 1.2579e-6 s period): the
    # entire measurement happens while the drive is still ramping up from zero
    # and never reaches full amplitude. A 5% tolerance measured there cannot
    # license a ten-period run, and the divergence it does show is linear -
    # fitted at 3.02 of rms B_x per period, crossing 5% at 0.0125 periods, with
    # a max residual of 2.3e-4. 0.25 periods is 10.5 x RISE and a quarter of a
    # full cycle, which is where the question can actually be asked.
    dict(
        folder='00-divergence-gate',
        deck=ANT_DECK, mesh=ANT_MESH,
        periods=0.25, frames=20,
        decks=[
            ('antenna-long-cleaning-on', [('DIVB_SPEED', 1.0), ('DIVE_SPEED', 1.0)]),
            ('antenna-long-cleaning-off', [('DIVB_SPEED', 0.0), ('DIVE_SPEED', 0.0)]),
            # The pair that actually isolates the cleaning SPEED. gamma and chi
            # are bare multiplicative factors on every term coupling phi and psi
            # to E and B, and chi also scales the charge source that generates
            # phi, so 0.0 is not a slower scheme - it is no scheme. Comparing
            # 1.0 against 0.5 keeps both non-zero, and because the Maxwell wave
            # speed is dmax(chi*c0, gamma*c0, c0) both give exactly c0 and hence
            # an identical timestep. 2.0 would double it and halve dt, which
            # would confound the comparison with a resolution change.
            ('antenna-long-cleaning-half', [('DIVB_SPEED', 0.5), ('DIVE_SPEED', 0.5)]),
        ],
    ),
    dict(
        folder='00-divergence-gate',
        deck=ANT_DECK, mesh=ANT_MESH,
        periods=0.00435, frames=6,
        # The antenna deck gets the SAME two-run test as the edge deck, not just
        # the psi/E symptom check. The whole campaign is predicated on this deck
        # being the sound one, and until both halves exist that is an assumption:
        # a deck can hold psi/E near 0.05 and still have its answer move when
        # DIVB_SPEED moves. Cheap to settle - about five minutes - and it governs
        # every long run in 01, 02 and 03.
        decks=[
            ('antenna-cleaning-on', [('DIVB_SPEED', 1.0), ('DIVE_SPEED', 1.0)]),
            ('antenna-cleaning-off', [('DIVB_SPEED', 0.0), ('DIVE_SPEED', 0.0)]),
        ],
    ),
    dict(
        folder='01-c-sensitivity',
        deck=ANT_DECK, mesh=ANT_MESH,
        periods=10.0, frames=10 * FRAMES_PER_PERIOD,
        decks=[
            ('c-as-shipped', [('LIGHT', 3.0e6)]),
            ('c-3vTe', [('LIGHT', C_3VTE)]),
        ],
    ),
    dict(
        folder='02-threshold-scan',
        deck=ANT_DECK, mesh=ANT_MESH,
        periods=5.0, frames=5 * FRAMES_PER_PERIOD,
        sweep=('Bomega', [5e-4, 10e-4, 15e-4, 20e-4, 30e-4, 50e-4, 75e-4,
                          100e-4, 150e-4]),
    ),
    # A Colab-sized pre-flight for the 20-period run below. One period is about
    # 7 h at 0.714 s/step, which fits inside a single session; the full run is
    # 139.6 h at the same rate and cannot be checkpointed, so it belongs on a
    # cluster. The probe cannot show formation - that needs many periods - but
    # it shows whether B_z on axis is moving toward reversal and whether zeta is
    # climbing, which is enough to catch a gross problem before committing six
    # days of somebody's queue.
    dict(
        folder='03-formation',
        deck=ANT_DECK, mesh=ANT_MESH,
        periods=1.0, frames=FRAMES_PER_PERIOD,
        decks=[('formation-probe', [])],
    ),
    dict(
        folder='03-formation',
        deck=ANT_DECK, mesh=ANT_MESH,
        periods=20.0, frames=20 * FRAMES_PER_PERIOD,
        decks=[('formation', [])],
    ),
    # Does the antenna instability of docs/known-issues.md 19 go away when the
    # winding is resolved? Every earlier refinement test stayed below one Debye
    # length per cell at the winding, so none could tell. This ladder crosses it,
    # with each arm's mesh generated here rather than shipped by hand. The
    # winding (COIL_R +- COIL_W/2 = 0.0335-0.0385) sits outside the shipped
    # mesh's refined region; --ring 0.040 brings it inside. lambda_D = 4.07e-4 m.
    dict(
        folder='05-winding-resolution',
        deck=ANT_DECK,
        periods=1.0, frames=2 * FRAMES_PER_PERIOD,
        decks=[
            # the shipped mesh, regenerated byte for byte: lambda_D/h 0.29 at the winding
            ('reference', [], WINDING_MESHES['disc.msh']),
            ('ring40', [], WINDING_MESHES['ring40.msh']),      # ~0.40 measured
            ('fine40', [], WINDING_MESHES['fine40.msh']),      # ~0.62 measured
        ],
    ),
    dict(
        folder='05-winding-resolution',
        deck=ANT_DECK,
        periods=0.5, frames=FRAMES_PER_PERIOD,
        decks=[
            # ~1.06 measured: the arm that actually crosses the threshold, and the
            # dearest. (4.0e-4 would measure 0.93 - mean edges run ~1.1x h-inner.)
            ('resolved40', [], WINDING_MESHES['resolved40.msh']),
            # the shipped mesh turned by 60 degrees: same cells, same dt, every
            # mesh seam moved. Peaks that rotate with it are sited by the mesh.
            ('rotated60', [], WINDING_MESHES['rotated60.msh']),
        ],
    ),
]


def render(spec):
    """[(filename, text, params)] for one run spec."""
    base = io.open(spec['deck'], encoding='utf-8').read()
    out = []
    if 'sweep' in spec:
        key, values = spec['sweep']
        cases = [('%s_%s' % (key, gauss(v)), [(key, v)]) for v in values]
    else:
        cases = spec['decks']
    for case in cases:
        name, overrides = case[0], case[1]
        mesh = case[2] if len(case) > 2 else None
        text = base
        for k, v in overrides:
            text = rmf_scan.set_param(text, k, v)
        if mesh is not None:
            text = set_gridname(text, mesh['file'])
        p = rd.deck_parameters_from_text(text)
        text = rmf_scan.set_param(text, 'TEND', spec['periods'] / p['omega'])
        text = rmf_scan.set_param(text, 'OUT', float(spec['frames']))
        out.append((name + '.pin', text, rd.deck_parameters_from_text(text), mesh))
    return out


def set_gridname(text, filename):
    """Point a deck at a different mesh file."""
    new, n = re.subn(r"^(\s*Gridname\s*=\s*')[^']*(')",
                     lambda m: m.group(1) + filename + m.group(2), text,
                     count=1, flags=re.M)
    if n != 1:
        raise SystemExit('could not find a quoted "Gridname =" line in the deck')
    return new


def make_mesh(mesh, directory):
    """Generate one mkdiscmesh mesh into `directory`; return its path."""
    path = os.path.join(directory, mesh['file'])
    done = subprocess.run([sys.executable, MKDISCMESH, '-o', path] + list(mesh['args']),
                          capture_output=True, text=True)
    if done.returncode != 0:
        raise SystemExit('mkdiscmesh failed for %s:\n%s%s'
                         % (mesh['file'], done.stdout, done.stderr))
    return path


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--check', action='store_true',
                    help='verify the checked-in folders match, and exit non-zero '
                         'if they do not')
    args = ap.parse_args(argv)

    stale, written = [], 0
    for spec in RUNS:
        folder = os.path.join(OUT, spec['folder'])
        if not args.check:
            os.makedirs(folder, exist_ok=True)
            if 'mesh' in spec:
                mesh_dst = os.path.join(folder, os.path.basename(spec['mesh']))
                if not os.path.exists(mesh_dst):
                    shutil.copy(spec['mesh'], mesh_dst)
        for name, text, params, mesh in render(spec):
            path = os.path.join(folder, name)
            if args.check:
                current = (io.open(path, encoding='utf-8').read()
                           if os.path.exists(path) else None)
                if current != text:
                    stale.append(os.path.relpath(path, REPO))
                if mesh is not None:
                    shipped = os.path.join(folder, mesh['file'])
                    with tempfile.TemporaryDirectory() as tmp:
                        with open(make_mesh(mesh, tmp), 'rb') as fh:
                            fresh = fh.read()
                    if not os.path.exists(shipped) or open(shipped, 'rb').read() != fresh:
                        stale.append(os.path.relpath(shipped, REPO))
                continue
            io.open(path, 'w', encoding='utf-8').write(text)
            written += 1
            mesh_path = make_mesh(mesh, folder) if mesh is not None else spec['mesh']
            dt, steps, secs, info = rmf_scan.estimate(params, mesh_path, ranks=4)
            print('%-24s %-26s %6.2f periods  %8d steps  %6.1f h at 4 ranks'
                  % (spec['folder'], name, params['TEND'] * params['omega'],
                     steps, secs / 3600.0))

    if args.check:
        if stale:
            print('out of date, re-run without --check:')
            for s in stale:
                print('   ', s)
            return 1
        print('phase3 run folders are up to date')
        return 0
    print('\nwrote %d decks under %s' % (written, os.path.relpath(OUT, REPO)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
