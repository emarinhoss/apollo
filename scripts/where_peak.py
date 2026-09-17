#!/usr/bin/env python3
"""Where in the domain does each component peak? Radius of the peak cell.

WHY THIS EXISTS. `run_growth.py` says WHICH component is growing; it does not
say where. That turned out to be the question that separated two completely
different failures of the same deck, which grow the same components at
indistinguishable rates:

  - a run whose peak sits at the ANTENNA WINDING is being driven, and the
    growth is the drive loading up;
  - a run whose peak has migrated to the CONDUCTING WALL is losing an
    unresolved electron sheath there - the wall cell evacuates while the
    electrons in it accelerate, and the run dies with
    "NaN after Mass matrix multiplication RHS".

Measured on a coarsened `03-formation` ladder: at 2x the shipped cell size the
peak |E_x| sits at the winding and the run survives, and at 3x it sits in a
single wall cell whose electrons go 1.2e5 -> 5.4e5 m/s over three output frames
before the run dies. The growth TABLES of those two runs look alike. The
locations do not.

    python3 scripts/where_peak.py <deck.pin> <frames...vtu>
    python3 scripts/where_peak.py <frames...vtu>          # radii, but no labels

THE DECK IS WHAT NAMES THE REGIONS. WALL_RADIUS, COIL_R, COIL_W and RAD_PLASMA
are deck values, and a tool that hard-codes one deck's geometry will label
another deck's run confidently and wrongly. Without a deck this prints the radii
- which are facts about the file - and says the labels are unavailable, rather
than inventing them.
"""

import glob
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'test'))
import numpy as np                                                   # noqa: E402
import vtu                                                           # noqa: E402

NAMES = ['e-rho', 'e-rhou', 'e-rhov', 'e-rhow', 'e-E',
         'i-rho', 'i-rhou', 'i-rhov', 'i-rhow', 'i-E',
         'E_x', 'E_y', 'E_z', 'B_x', 'B_y', 'B_z', 'phi', 'psi']
INDEX = {n: i for i, n in enumerate(NAMES)}

# The components that actually move on the RMF decks: the fields, both cleaning
# potentials, and the electron fluid. Override with --components.
DEFAULT = ['E_x', 'E_y', 'E_z', 'phi', 'psi', 'e-rho', 'e-E']

WALL_FRACTION = 0.95      # r > this * WALL_RADIUS counts as the wall ring


def frame_number(path):
    """The N in `<run>_N.vtu`, or a refusal naming the file.

    Same rule and the same reason as run_growth.py: the frames have to be
    ordered by number, not by name, or frame 10 sorts between 1 and 2.
    """
    stem = os.path.basename(path).rsplit('.', 1)[0]
    if '_' in stem:
        tail = stem.rsplit('_', 1)[1]
        if tail.isdigit():
            return int(tail)
    raise SystemExit(
        "where_peak.py: cannot tell which frame '%s' is.\n"
        "  Frames are ordered by the number the solver writes into the name:\n"
        "  <run>_0.vtu, <run>_1.vtu, ... This argument does not carry one.\n"
        "  Pass the .vtu files of one run:  where_peak.py <deck.pin> <dir>/*.vtu"
        % path)


def geometry(deck):
    """Region boundaries from the deck, or None when no deck was given."""
    if deck is None:
        return None
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import rmf_diagnostics                    # its deck reader, and its refusals
    ns = rmf_diagnostics.deck_parameters(deck)
    return {k: ns.get(k) for k in ('WALL_RADIUS', 'COIL_R', 'COIL_W', 'RAD_PLASMA')}


def label(r, geo):
    if geo is None:
        return ''
    wall, coil, width, plasma = (geo['WALL_RADIUS'], geo['COIL_R'],
                                 geo['COIL_W'], geo['RAD_PLASMA'])
    if wall and r > WALL_FRACTION * wall:
        return 'WALL'
    if coil and width and abs(r - coil) <= width:
        return 'winding'
    if plasma and r < plasma:
        return 'column'
    return 'gap'


def centroids(arrays, path):
    for key in ('Position', 'connectivity'):
        if key not in arrays:
            raise SystemExit(
                '%s carries no %s array, so the peak cannot be located.\n'
                '  Every cell-centred .vtu Apollo writes has both; a file that '
                'does not\n  is either not an Apollo frame or was truncated.'
                % (path, key))
    P = np.asarray(arrays['Position'])[:, :2]
    conn = np.asarray(arrays['connectivity']).astype(int).reshape(-1, 3)
    return P[conn].mean(axis=1)


def peak_cell(arrays, comp):
    """(value, cell) of max|q| over ALL THREE DG nodes, not just node 0.

    The solution is stored as solutiondg.N with N = node*18 + component, so a
    peak that lives on node 1 or 2 is invisible to a reader that looks only at
    the first 18 arrays.
    """
    best, cell = 0.0, -1
    for node in range(3):
        key = 'solutiondg.%d' % (node * 18 + comp)
        if key not in arrays:
            continue
        v = np.abs(np.asarray(arrays[key]))
        if not v.size:
            continue
        i = int(np.nanargmax(v))
        if v[i] > best:
            best, cell = float(v[i]), i
    return best, cell


def main(deck, frames, comps, nframes):
    geo = geometry(deck)
    paths = sorted(frames, key=frame_number)
    if geo:
        print('deck %s: wall r=%.4f, winding r=%.4f+-%.4f, column r<%.4f'
              % (os.path.basename(deck), geo['WALL_RADIUS'] or float('nan'),
                 geo['COIL_R'] or float('nan'), geo['COIL_W'] or float('nan'),
                 geo['RAD_PLASMA'] or float('nan')))
    else:
        print('no deck given, so regions cannot be named: radii only.')
        print('  pass the deck first to label them:  where_peak.py <deck.pin> *.vtu')

    at_wall = {c: 0 for c in comps}
    counted = 0
    rows = []
    for p in paths:
        a = vtu.read(p)
        cen = centroids(a, p)
        r = np.hypot(cen[:, 0], cen[:, 1])
        row = []
        for name in comps:
            value, cell = peak_cell(a, INDEX[name])
            if cell < 0:
                row.append(None)
                continue
            rr = float(r[cell])
            tag = label(rr, geo)
            if tag == 'WALL':
                at_wall[name] += 1
            row.append((value, rr, tag,
                        math.degrees(math.atan2(cen[cell][1], cen[cell][0]))))
        counted += 1
        rows.append((frame_number(p), row))

    print()
    print('%6s ' % 'frame' + ' '.join('%26s' % n for n in comps))
    for idx, row in rows[-nframes:]:
        cells = []
        for item in row:
            if item is None:
                cells.append('%26s' % '-')
            else:
                value, rr, tag, th = item
                cells.append('%11.3e r=%.4f %-7s' % (value, rr, tag))
        print('%6d ' % idx + ' '.join(cells))

    if geo:
        print()
        print('frames whose peak is in the wall ring (r > %.2f x WALL_RADIUS):'
              % WALL_FRACTION)
        for name in comps:
            print('  %-6s %3d of %3d  (%3.0f%%)'
                  % (name, at_wall[name], counted, 100.0 * at_wall[name] / max(counted, 1)))


USAGE = '''usage: where_peak.py [deck.pin] <frames...vtu>

Radius of the cell holding the peak of each component, for the last frames of a
run. Answers the question run_growth.py cannot: is the growth at the CONDUCTING
WALL, at the ANTENNA WINDING, or inside the column?

    python3 scripts/where_peak.py formation.pin results/03-formation/*.vtu

options:
  --components a,b,c   which components (default: %s)
  --frames N           how many trailing frames to print (default 6)
''' % ','.join(DEFAULT)


if __name__ == '__main__':
    args = sys.argv[1:]
    if not args or args[0] in ('-h', '--help'):
        print(USAGE, end='')
        raise SystemExit(0 if args else 2)

    comps, nframes, rest = list(DEFAULT), 6, []
    i = 0
    while i < len(args):
        if args[i] == '--components' and i + 1 < len(args):
            comps = [c.strip() for c in args[i + 1].split(',') if c.strip()]
            bad = [c for c in comps if c not in INDEX]
            if bad:
                raise SystemExit('unknown component(s): %s\nknown: %s'
                                 % (', '.join(bad), ', '.join(NAMES)))
            i += 2
        elif args[i] == '--frames' and i + 1 < len(args):
            try:
                nframes = int(args[i + 1])
            except ValueError:
                raise SystemExit('--frames wants an integer, got %r' % args[i + 1])
            i += 2
        else:
            rest.append(args[i])
            i += 1

    # A shell glob that matched nothing arrives as the literal pattern.
    expanded = []
    for a in rest:
        expanded.extend(sorted(glob.glob(a)) if any(ch in a for ch in '*?[') else [a])
    deck = None
    if expanded and expanded[0].endswith('.pin'):
        deck, expanded = expanded[0], expanded[1:]
    frames = [p for p in expanded if p.endswith('.vtu')]
    if not frames:
        raise SystemExit('no .vtu frames given.\n' + USAGE)
    main(deck, frames, comps, nframes)
