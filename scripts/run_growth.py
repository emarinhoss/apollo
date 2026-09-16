#!/usr/bin/env python3
"""Per-frame max|value| for each of the 18 components, to find what is growing.

WHY THIS EXISTS. `rmf_diagnostics.py` answers "what physics did this run
produce". When a run DIES it is the wrong question, and its aggregates can look
healthy right up to the last frame: they bin over r = 0..a and the axis, and on
the formation deck 15917 nodal points per frame fall outside that and are not
counted at all. A run that went unstable somewhere else - at the antenna
winding, at the wall - looks fine in every column.

This asks the other question: which of the eighteen components is growing, and
how fast. An instability shows up as a geometric climb frame over frame.

    python3 scripts/run_growth.py <results-dir>/*.vtu

WORKED EXAMPLE, the run this was written for. formation.pin died with
"NaN after Mass matrix multiplication RHS" at 1.73 of 20 RMF periods, ~12 hours
in, no checkpoint. The physics table showed nothing: zeta -0.0103, penetration
0.005, B_z flat on its 6 mT bias, "no layer" almost throughout - a run in which
almost nothing happened, right up to a frame that looked healthy one output
interval before the NaN. The raw magnitudes, from each field's first non-zero
frame to the last written:

    phi   1.49e-05 -> 4.44e-03   298x     <- the E-field cleaning potential
    E_y   1.02e+02 -> 1.89e+04   186x
    E_x   1.46e+02 -> 1.90e+04   131x
    psi   3.73e+01 -> 4.59e+01   1.2x     <- the B-field one, untroubled
    E_z   2.21e+02 -> 2.34e+02   1.1x     <- the physical RMF drive, flat
    e-rho 9.12e-11 -> 1.34e-10   1.5x

which names the divergence cleaning of E, and rules out the fluid and the
drive, from twenty frames that the physics diagnostics called uneventful.
"""
import glob, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                os.pardir, 'test'))
import numpy as np, vtu

NAMES = ['e-rho','e-rhou','e-rhov','e-rhow','e-E',
         'i-rho','i-rhou','i-rhov','i-rhow','i-E',
         'E_x','E_y','E_z','B_x','B_y','B_z','phi','psi']

def main(paths):
    files = sorted(paths, key=lambda p: int(os.path.basename(p).rsplit('_',1)[1].split('.')[0]))
    rows = []
    for p in files:
        a = vtu.read(p)
        peak = {}
        for k, v in a.items():
            if not k.startswith('solutiondg.'):
                continue
            comp = int(k.split('.')[1]) % 18
            m = float(np.nanmax(np.abs(np.asarray(v))))
            peak[comp] = max(peak.get(comp, 0.0), m)
        nonfinite = sum(int((~np.isfinite(np.asarray(v))).sum())
                        for k, v in a.items() if k.startswith('solutiondg.'))
        rows.append((os.path.basename(p), peak, nonfinite))

    print('%-6s %-9s ' % ('frame', 'nonfin') + ' '.join('%-11s' % n for n in NAMES))
    for i, (name, peak, nf) in enumerate(rows):
        print('%-6d %-9d ' % (i, nf) + ' '.join('%-11.3e' % peak.get(c, float('nan')) for c in range(18)))

    print()
    # BASELINE. Frame 0 is the initial condition, where the driven fields and
    # both cleaning potentials are EXACTLY zero. Ranking by last/first with
    # frame 0 as the denominator therefore scores every field that starts at
    # zero as 0 and sorts it LAST - which is precisely backwards, because a
    # field that grows from nothing is the one worth looking at.
    #
    # This is not hypothetical: on a run that died at 1.73 RMF periods, the
    # first version of this script reported e-rho (1.5x) as the fastest grower
    # and buried phi (298x), E_y (186x) and E_x (131x) at the bottom of the
    # list. The three it buried were the instability.
    #
    # So the baseline for each component is its first non-zero frame, and a
    # component that never leaves zero is reported as such rather than ranked.
    def baseline(c):
        for i, (_, peak, _) in enumerate(rows):
            if peak.get(c):
                return i, peak[c]
        return None, 0.0

    last = rows[-1][1]
    grew, never = [], []
    for c in range(18):
        i, f = baseline(c)
        if i is None:
            never.append(c)
        else:
            grew.append((last.get(c, 0.0) / f, c, i, f, last.get(c, 0.0)))
    grew.sort(reverse=True)

    print('growth, each from its own first non-zero frame:')
    for factor, c, i, f, l in grew[:8]:
        print('  %-8s frame %-2d %11.3e -> %11.3e   %8.1fx'
              % (NAMES[c], i, f, l, factor))
    if never:
        print('  never non-zero: %s' % ', '.join(NAMES[c] for c in never))
    worst = [c for _, c, _, _, _ in grew]

    print()
    print('frame-over-frame ratio for the fastest grower (%s):' % NAMES[worst[0]])
    c = worst[0]
    for i in range(1, len(rows)):
        a, b = rows[i-1][1].get(c, 0.0), rows[i][1].get(c, 0.0)
        print('  %2d -> %2d   %11.3e   %s' % (i-1, i, b, ('x%.3f' % (b/a)) if a else '-'))

if __name__ == '__main__':
    if len(sys.argv) < 2:
        raise SystemExit('usage: growth.py <frames...vtu>')
    main(sys.argv[1:])
