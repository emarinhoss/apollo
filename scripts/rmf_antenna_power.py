#!/usr/bin/env python3
"""Antenna power for the RMF example: P = -integral of E.J over the winding.

This is the number the RMF literature quotes alongside the penetration
threshold, and it is the reason for driving the plasma with a current rather
than with a prescribed edge field: a prescribed field cannot be loaded, so it
has no power. Here the antenna's current is known exactly (it is a deck
parameter) and E_z is a solved-for state variable, so the work the antenna does
on the field and plasma follows directly:

    P = - integral over the winding of E_z J_z dA        [W per metre of axial length]

with the sign such that P > 0 means the antenna is delivering power.

READ THE CYCLE AVERAGE, NOT THE INSTANTANEOUS VALUES. Most of E_z at the winding
is in quadrature with J_z, so the instantaneous P is dominated by reactive
exchange and swings by a factor of several within a fraction of a period. Only
its average over a whole RMF period is the power actually absorbed. This script
prints that average when - and only when - the frames it was given span at least
one period; with less than that it prints the instantaneous values and says so.
For scale, the plasma-free vacuum.pin case swings between about 1e6 and 7e6 W/m
while its steady absorbed power is zero, there being nothing to absorb it.

This lives in post-processing rather than in the solver because Apollo's area
integrals are not observable: WxNodalDG2dMethod reduces them into
_AgregateAreaIntegral and passes that to boundary conditions only - it is never
logged and never written to output (see src/subsolvers/wxnodaldg2dmethod.cc).

Usage:

    python3 scripts/rmf_antenna_power.py <deck.pin> <case_0.vtu> [more .vtu ...]

The deck is read for the antenna parameters, so this cannot drift out of step
with the run it is describing.
"""

import argparse
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'test'))

try:
    import numpy as np
except ImportError:
    sys.exit('numpy is required')

try:
    import vtu
except ImportError:
    sys.exit('test/vtu.py is required and was not importable')

NEQ, EZ = 18, 12          # 18-component two-fluid state; E_z is component 12


def deck_values(path):
    """Evaluate the deck's Python preamble and return the names we need.

    A .pin is Python down to the <warpx> line, so the parameters can be had by
    executing that much of it - the same thing scripts/wxinpparse.py does.
    """
    with open(path) as fh:
        text = fh.read()
    preamble = text.split('<warpx>', 1)[0]
    namespace = {}
    exec(compile(preamble, path, 'exec'), namespace)   # noqa: S102 - the deck IS Python
    missing = [k for k in ('Bomega', 'omega', 'PHASE', 'RISE', 'COIL_R', 'COIL_W', 'MU0')
               if k not in namespace]
    if missing:
        sys.exit(f'{path}: deck does not define {", ".join(missing)}; '
                 'is this an antenna deck?')
    return namespace


def current_amplitude(deck):
    """K, the same way ApRMFAntennaSrc::setup solves for it.

    B_uniform = (mu0/2) K * integral of f(s)(1 - s^2/b^2) ds, with f the
    normalised raised cosine over the winding.
    """
    rc, w = deck['COIL_R'], deck['COIL_W']
    # The decks pass WALL_RADIUS as the antenna's conductor_radius; absent it,
    # the winding is in free space and there is no image current.
    b = deck.get('WALL_RADIUS', 0.0)
    rmin, rmax = rc - 0.5 * w, rc + 0.5 * w
    s = np.linspace(rmin, rmax, 20001)
    f = 0.5 * (1.0 - np.cos(2.0 * math.pi * (s - rmin) / (rmax - rmin)))
    screen = (1.0 - s * s / (b * b)) if b > 0.0 else np.ones_like(s)
    norm = np.trapezoid(f, s) if hasattr(np, 'trapezoid') else np.trapz(f, s)
    shape = (np.trapezoid(f * screen, s) if hasattr(np, 'trapezoid')
             else np.trapz(f * screen, s))
    return 2.0 * deck['Bomega'] / (deck['MU0'] * (shape / norm)), norm


def antenna_current(deck, K, norm, x, y, t):
    """J_z at (x, y, t), the same expression ApRMFAntennaSrc::src evaluates."""
    rc, w = deck['COIL_R'], deck['COIL_W']
    rmin, rmax = rc - 0.5 * w, rc + 0.5 * w
    r = np.hypot(x, y)
    inside = (r >= rmin) & (r <= rmax)
    profile = np.where(inside,
                       0.5 * (1.0 - np.cos(2.0 * math.pi * (r - rmin) / (rmax - rmin))),
                       0.0)
    omega = 2.0 * math.pi * deck['omega']
    envelope = 1.0 - math.exp(-t / deck['RISE'])
    return (K * envelope * profile / norm
            * np.cos(np.arctan2(y, x) - omega * t - deck['PHASE']))


def frame_time(path, deck):
    """Output frames are written at TEND * index / OUT."""
    index = int(os.path.basename(path).rsplit('_', 1)[1].split('.')[0])
    return deck['TEND'] * index / deck['OUT']


def power(path, deck, K, norm):
    """-integral E_z J_z dA over the mesh, by P1 nodal quadrature."""
    data = vtu.read(path)
    points, cells = data['Position'], data['connectivity'].reshape(-1, 3)
    t = frame_time(path, deck)

    total = 0.0
    # Vertices of every cell, and the cell areas: for a linear field the exact
    # integral over a triangle is its area times the mean of the three nodes.
    p = [points[cells[:, k]] for k in range(3)]
    area = 0.5 * np.abs((p[1][:, 0] - p[0][:, 0]) * (p[2][:, 1] - p[0][:, 1])
                        - (p[2][:, 0] - p[0][:, 0]) * (p[1][:, 1] - p[0][:, 1]))
    for k in range(3):
        ez = data[f'solutiondg.{k * NEQ + EZ}']
        jz = antenna_current(deck, K, norm, p[k][:, 0], p[k][:, 1], t)
        total += np.sum(area * ez * jz) / 3.0
    return t, -total


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('deck', help='the antenna .pin the run used')
    ap.add_argument('vtu', nargs='+', help='output frames, in any order')
    args = ap.parse_args(argv)

    deck = deck_values(args.deck)
    K, norm = current_amplitude(deck)
    print(f'antenna: B_rmf = {deck["Bomega"]:.4g} T, K = {K:.6g} A/m, '
          f'r_c = {deck["COIL_R"]:.4g} m, width = {deck["COIL_W"]:.4g} m')
    print(f'{"t [s]":>12}  {"P [W/m]":>14}')
    series = []
    for path in sorted(args.vtu,
                       key=lambda f: int(os.path.basename(f).rsplit('_', 1)[1].split('.')[0])):
        t, p = power(path, deck, K, norm)
        series.append((t, p))
        print(f'{t:12.4e}  {p:14.6e}')

    period = 1.0 / deck['omega']
    span = series[-1][0] - series[0][0] if len(series) > 1 else 0.0
    print()
    if span >= period:
        # Trailing whole period only: earlier frames carry the switch-on
        # transient, which is not part of a steady absorbed power.
        tail = [p for t, p in series if t >= series[-1][0] - period]
        mean = sum(tail) / len(tail)
        print(f'cycle-averaged absorbed power over the last period '
              f'({len(tail)} frames): {mean:.6e} W/m')
    else:
        print(f'frames span {span:.3e} s, less than one RMF period '
              f'({period:.3e} s): the values above are INSTANTANEOUS and mostly')
        print('reactive. Run to at least one period before quoting an absorbed power.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
