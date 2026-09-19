#!/usr/bin/env python3
"""What the plasma is doing at the antenna, frame by frame, where the run dies.

`run_growth.py` reports max|value| per component and `rmf_diagnostics.py` bins
the column out to RAD_PLASMA = 0.030 - which excludes the antenna at 0.036 and
the wall at 0.050, i.e. exactly where the formation run's peaks sit. This tool
reads the whole disc and answers the questions the instability plan turns on:

  * Is a density HOLE forming at the winding, how deep, and is the electron
    pressure staying flat (Ohmic heating balanced by expansion) or running away?
  * Does the electron fluid reach the state the scheme cannot carry - a nodal
    pressure going non-positive (silently floored, never logged), or
    |u_e| + sound speed crossing the light speed the timestep is pinned to?
  * Is the charge separation real physics or a discrete Gauss-law violation?
    It reports the per-cell residual div E - rho_c/eps0 against div E itself.
  * Radial profiles of n, T_e, the drift velocities, |B_perp| and E_r ALL THE
    WAY TO THE WALL, with the count of nodes in each bin so an empty bin reads
    as NaN rather than zero.

Everything is measured; nothing is guessed. It refuses without the deck (it
needs the masses, the charge, the light speed for eps0, and the geometry), and a
radial bin with no nodes is NaN, not 0.

Usage:
    python3 scripts/winding_anatomy.py <deck.pin> <run>_*.vtu
    python3 scripts/winding_anatomy.py --nbins 100 <deck.pin> frames...
"""

import argparse
import math
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'scripts'))
sys.path.insert(0, os.path.join(REPO, 'test'))

try:
    import numpy as np
except ImportError:                                  # pragma: no cover
    np = None

import rmf_diagnostics as rd
import vtu

RHO_E, MOM_EX, MOM_EY, MOM_EZ, EN_E = 0, 1, 2, 3, 4
RHO_I, MOM_IX, MOM_IY, MOM_IZ, EN_I = 5, 6, 7, 8, 9
E_X, E_Y, E_Z, B_X, B_Y, B_Z, PHI, PSI = 10, 11, 12, 13, 14, 15, 16, 17


def frame_index(path):
    base = os.path.basename(path)
    stem = base[:-4] if base.endswith('.vtu') else base
    if '_' not in stem or not stem.rsplit('_', 1)[1].isdigit():
        raise SystemExit(f'{base}: frame files must end in _<number>.vtu')
    return int(stem.rsplit('_', 1)[1])


def nodal_state(data):
    """Flat per-cell-node coordinates and 18 components (length 3*ncells).

    These are the values as the solver STORES them - not averaged onto vertices
    - because a floor engages on a nodal value, so counting non-positive nodal
    pressures needs the raw nodal state, not a vertex mean that could hide it.
    """
    cells = np.asarray(data['connectivity']).reshape(-1, 3)
    points = np.asarray(data['Position'])[:, :2]
    xs, ys, comp = [], [], {c: [] for c in range(18)}
    for node in range(3):
        vid = cells[:, node]
        xs.append(points[vid, 0])
        ys.append(points[vid, 1])
        for c in range(18):
            comp[c].append(np.asarray(data[f'solutiondg.{node * 18 + c}']))
    x = np.concatenate(xs)
    y = np.concatenate(ys)
    return x, y, {c: np.concatenate(v) for c, v in comp.items()}


def primitives(comp, params):
    """n_e, n_i, |u_e|, u_e_z, p_e, p_i, T_e[eV] from a nodal/vertex state."""
    me, mi, q, g = params['ME'], params['MI'], params['Q'], params.get('GAMMA', 5.0 / 3.0)
    rho_e = comp[RHO_E]
    rho_i = comp[RHO_I]
    with np.errstate(invalid='ignore', divide='ignore'):
        uex, uey, uez = comp[MOM_EX] / rho_e, comp[MOM_EY] / rho_e, comp[MOM_EZ] / rho_e
        uix, uiy, uiz = comp[MOM_IX] / rho_i, comp[MOM_IY] / rho_i, comp[MOM_IZ] / rho_i
        ke_e = 0.5 * rho_e * (uex ** 2 + uey ** 2 + uez ** 2)
        ke_i = 0.5 * rho_i * (uix ** 2 + uiy ** 2 + uiz ** 2)
        p_e = (g - 1.0) * (comp[EN_E] - ke_e)
        p_i = (g - 1.0) * (comp[EN_I] - ke_i)
        n_e = rho_e / me
        n_i = rho_i / mi
        speed_e = np.sqrt(uex ** 2 + uey ** 2 + uez ** 2)
        a_e = np.sqrt(np.maximum(g * p_e, 0.0) / rho_e)
        t_e = p_e / (n_e * q)
    return dict(n_e=n_e, n_i=n_i, p_e=p_e, p_i=p_i, u_e=speed_e, u_ez=uez,
                a_e=a_e, lam_e=speed_e + a_e, t_e=t_e)


def cell_divE(data):
    """Per-cell div E (constant on a P1 triangle) and the cell-mean charge/eps0.

    On a linear triangle E is affine, so div E is a single number per cell got
    from the three nodal E_x, E_y and the shape-function gradients. The cell-mean
    charge density is the mean of the three nodal e(n_i - n_e). Their difference
    is the discrete Gauss residual the cleaning potential phi is fed.
    """
    cells = np.asarray(data['connectivity']).reshape(-1, 3)
    pts = np.asarray(data['Position'])[:, :2]
    ex = np.vstack([np.asarray(data[f'solutiondg.{n * 18 + E_X}']) for n in range(3)]).T
    ey = np.vstack([np.asarray(data[f'solutiondg.{n * 18 + E_Y}']) for n in range(3)]).T
    v = pts[cells]                                   # (ncell, 3, 2)
    x0, y0 = v[:, 0, 0], v[:, 0, 1]
    x1, y1 = v[:, 1, 0], v[:, 1, 1]
    x2, y2 = v[:, 2, 0], v[:, 2, 1]
    twoA = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
    with np.errstate(invalid='ignore', divide='ignore'):
        # b_i = dN_i/dx, c_i = dN_i/dy for the linear triangle
        b0, b1, b2 = (y1 - y2) / twoA, (y2 - y0) / twoA, (y0 - y1) / twoA
        c0, c1, c2 = (x2 - x1) / twoA, (x0 - x2) / twoA, (x1 - x0) / twoA
        dExdx = ex[:, 0] * b0 + ex[:, 1] * b1 + ex[:, 2] * b2
        dEydy = ey[:, 0] * c0 + ey[:, 1] * c1 + ey[:, 2] * c2
    divE = dExdx + dEydy
    ok = np.abs(twoA) > 0
    return divE, ok


def band_extremes(r, prim, params):
    """Winding-band scalar summaries for one frame."""
    rc, w = params['COIL_R'], params['COIL_W']
    a0 = params['a_e0']
    sel = (r >= rc - w) & (r <= rc + w)
    if not np.any(sel):
        return {}
    ne = prim['n_e'][sel]
    ni = prim['n_i'][sel]
    n0 = params['n_dens']
    imin = int(np.nanargmin(ne))
    return {
        'n_e_min/n0': float(np.nanmin(ne) / n0),
        'n_e_mean/n0': float(np.nanmean(ne) / n0),
        '(ni-ne)/n0@min': float((ni[imin] - ne[imin]) / n0),
        'Te_max': float(np.nanmax(prim['t_e'][sel])),
        'Te_mean': float(np.nanmean(prim['t_e'][sel])),
        'max|ue|/a_e0': float(np.nanmax(prim['u_e'][sel]) / a0),
        'max|uez|/a_e0': float(np.nanmax(np.abs(prim['u_ez'][sel])) / a0),
        'max(lam_e)/c0': float(np.nanmax(prim['lam_e'][sel]) / params['LIGHT']),
    }


def main(argv=None):
    if np is None:
        raise SystemExit('winding_anatomy needs numpy')
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('deck')
    ap.add_argument('frames', nargs='+')
    ap.add_argument('--nbins', type=int, default=100)
    ap.add_argument('--profile-frame', type=int, default=-1,
                    help='frame index to print a full radial profile for (default: last)')
    args = ap.parse_args(argv)

    params = rd.deck_parameters(args.deck)
    for need in ('ME', 'MI', 'Q', 'MU0', 'LIGHT', 'COIL_R', 'COIL_W', 'WALL_RADIUS',
                 'n_dens', 'Te', 'TEND', 'OUT'):
        if need not in params:
            raise SystemExit(f'{args.deck}: deck is missing {need}, needed for winding_anatomy')
    eps0 = 1.0 / (params['MU0'] * params['LIGHT'] ** 2)
    g = params.get('GAMMA', 5.0 / 3.0)
    # electron sound speed of the initial uniform column, the natural yardstick
    params['a_e0'] = math.sqrt(g * params['n_dens'] * params['Q'] * params['Te']
                               / (params['n_dens'] * params['ME']))
    dt_frame = params['TEND'] / params['OUT']
    floors = {
        'p': 0.01 * params['n_dens'] * params['Q'] * params['Te'],   # MIN_PRES
        'rho_e': params['ME'] * 0.01 * params['n_dens'],
        'rho_i': params['MI'] * 0.01 * params['n_dens'],
    }

    frames = sorted(args.frames, key=frame_index)
    print(f'# winding_anatomy: {len(frames)} frames, frame dt={dt_frame:.4e} s, '
          f'a_e0={params["a_e0"]:.4e} m/s, c0={params["LIGHT"]:.3e}, '
          f'eps0={eps0:.4e}, floors p={floors["p"]:.3g}')
    print('# t  ' + '  '.join(['n_e_min/n0', '(ni-ne)/n0@min', 'Te_max', 'Te_mean',
                               'max|uez|/a_e0', 'max(lam_e)/c0', 'nodes[p_e<=0]',
                               'nodes[rho_e<floor]', 'gauss_resid/divE', 'c0|phi|/|Eperp|']))
    profframe = frames[-1] if args.profile_frame < 0 else None
    for path in frames:
        idx = frame_index(path)
        t = idx * dt_frame
        data = vtu.read(path)
        x, y, comp = nodal_state(data)
        r, _, _ = rd.polar(x, y)
        prim = primitives(comp, params)
        ext = band_extremes(r, prim, params)
        # floor engagement over the WHOLE disc (nodal)
        npe = int(np.count_nonzero(prim['p_e'] <= 0.0))
        npi = int(np.count_nonzero(prim['p_i'] <= 0.0))
        nre = int(np.count_nonzero(comp[RHO_E] < floors['rho_e']))
        # Gauss residual per cell
        divE, ok = cell_divE(data)
        cells = np.asarray(data['connectivity']).reshape(-1, 3)
        rho_c_node = params['Q'] * (comp[RHO_I] / params['MI'] - comp[RHO_E] / params['ME'])
        ncell = cells.shape[0]
        rho_c_cell = np.mean(rho_c_node.reshape(3, ncell), axis=0)
        resid = divE - rho_c_cell / eps0
        if np.any(ok):
            dnorm = float(np.sqrt(np.nanmean(divE[ok] ** 2)))
            gratio = (float(np.sqrt(np.nanmean(resid[ok] ** 2))) / dnorm
                      if dnorm > 0 else float('nan'))
        else:
            gratio = float('nan')
        # phi vs E_perp (whole disc, raw nodal)
        eperp = np.hypot(comp[E_X], comp[E_Y])
        cphi = params['LIGHT'] * np.nanmax(np.abs(comp[PHI]))
        rphi = cphi / np.nanmax(eperp) if np.nanmax(eperp) > 0 else float('nan')
        if ext:
            print(f'{t:.4e}  {ext["n_e_min/n0"]:.4f}  {ext["(ni-ne)/n0@min"]:+.4f}  '
                  f'{ext["Te_max"]:.2f}  {ext["Te_mean"]:.2f}  {ext["max|uez|/a_e0"]:.3f}  '
                  f'{ext["max(lam_e)/c0"]:.3f}  {npe:5d}  {nre:5d}  {gratio:.3f}  {rphi:.3f}')

    # one full radial profile
    data = vtu.read(profframe if profframe else frames[-1])
    x, y, comp = nodal_state(data)
    r, ct, st = rd.polar(x, y)
    prim = primitives(comp, params)
    edges = np.linspace(0.0, params['WALL_RADIUS'], args.nbins + 1)
    uer = comp[MOM_EX] / comp[RHO_E] * ct + comp[MOM_EY] / comp[RHO_E] * st
    bperp = np.hypot(comp[B_X], comp[B_Y])
    er = comp[E_X] * ct + comp[E_Y] * st
    print(f'# radial profile, frame {frame_index(profframe if profframe else frames[-1])}: '
          'r  n_e/n0  Te[eV]  u_er  |Bperp|  E_r  count')
    for name, val in (('n_e', prim['n_e'] / params['n_dens']), ('Te', prim['t_e']),
                      ('u_er', uer), ('Bperp', bperp), ('E_r', er)):
        c, mean, _, cnt, _ = rd.radial_profile(r, val, edges)
        if name == 'n_e':
            profiles = {'r': c, 'count': cnt}
        profiles[name] = mean
    for i in range(len(profiles['r'])):
        if profiles['count'][i] == 0:
            continue
        print(f'{profiles["r"][i]:.4f}  {profiles["n_e"][i]:.4f}  {profiles["Te"][i]:8.2f}  '
              f'{profiles["u_er"][i]:+.3e}  {profiles["Bperp"][i]:.3e}  {profiles["E_r"][i]:+.3e}  '
              f'{profiles["count"][i]:5d}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
