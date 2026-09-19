#!/usr/bin/env python3
"""Azimuthal structure of a run's fields on the mesh's node rings near the antenna.

`run_growth.py` says WHICH component is growing and `where_peak.py` says in WHICH
REGION its peak sits. Neither says whether the growth is a smooth, low-m,
drive-locked pattern (the antenna's own m=1/m=2 response, or a coherent physical
mode of a few wavelengths) or a grid-scale, cell-to-cell structure that the mesh
cannot resolve (an unshielded Debye sheath, a species-asymmetric-dissipation
artefact). On the formation deck those are DIFFERENT failures with the same
run_growth table, and this is the tool that tells them apart.

For each output frame and each node ring near the winding it reports, per
component (charge density e(n_i - n_e), phi, E_r, E_theta, n_e, n_i, u_e_theta):

  * |A_0|            the axisymmetric (m=0) amplitude
  * |A_1| .. |A_8|   the low azimuthal harmonics
  * grid-band rms    the r.m.s. of what is left after removing |m| <= 4 -
                     i.e. the cell-to-cell structure the low harmonics miss
  * the residual phase arg(A_m) - m*omega*t of the dominant low harmonic, so a
    pattern rotating rigidly with the drive (constant residual) is told apart
    from a wave that rotates at its own speed (drifting residual)

and, across the whole mesh, the (r, theta) of each component's peak cell per
frame, so a peak pinned to a fixed mesh azimuth (a stitching seam) is told apart
from one that co-rotates with the antenna.

It refuses rather than guesses: it needs the deck (for omega, the species masses
and charge, and the winding radius), it needs at least one ring of >= MIN_RING
nodes within the winding band, and it fits a growth rate only over frames past
RAMP_CLEAR * RISE with at least MIN_FIT of them, returning NaN otherwise.

Usage:
    python3 scripts/ring_spectrum.py <deck.pin> <run>_*.vtu
    python3 scripts/ring_spectrum.py --bands winding,column,wall <deck.pin> frames...
    python3 scripts/ring_spectrum.py --mmax 60 --json out.json <deck.pin> frames...
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

# Component indices in the 18-component two-fluid + PH-Maxwell state.
RHO_E, MOM_EX, MOM_EY, MOM_EZ, EN_E = 0, 1, 2, 3, 4
RHO_I, MOM_IX, MOM_IY, MOM_IZ, EN_I = 5, 6, 7, 8, 9
E_X, E_Y, E_Z, B_X, B_Y, B_Z, PHI, PSI = 10, 11, 12, 13, 14, 15, 16, 17

MIN_RING = 50        # a ring with fewer nodes is not resolved enough to Fourier
RAMP_CLEAR = 3.0     # fit growth only past this many antenna rise times
MIN_FIT = 5          # ... and only with at least this many surviving frames
RING_TOL = 0.35      # nodes within RING_TOL * median-spacing share a ring
LOW_M = 4            # harmonics kept as "smooth"; the rest is the grid band


def vertex_values(data, components):
    """Average the P1 cell-node state onto shared mesh vertices.

    vtu stores the solution per cell node (solutiondg.(node*18 + component)); a
    physical vertex is shared by several cells and each holds its own nodal
    value there. The mesh-continuous field at a vertex is their mean, which is
    what an azimuthal transform on a node ring needs. Returns (x, y, {c: val})
    with x, y the vertex coordinates and each val an array over vertices.
    """
    points = np.asarray(data['Position'])[:, :2]
    cells = np.asarray(data['connectivity']).reshape(-1, 3)
    nvert = len(points)
    count = np.zeros(nvert)
    for node in range(3):
        np.add.at(count, cells[:, node], 1.0)
    if np.any(count == 0):
        # A vertex in the file that no cell references cannot carry a value.
        # Keep it (so indices line up with Position) but mark it empty.
        pass
    out = {}
    for c in components:
        acc = np.zeros(nvert)
        for node in range(3):
            np.add.at(acc, cells[:, node], np.asarray(data[f'solutiondg.{node * 18 + c}']))
        with np.errstate(invalid='ignore', divide='ignore'):
            out[c] = np.where(count > 0, acc / np.where(count > 0, count, 1.0), np.nan)
    return points[:, 0], points[:, 1], out, count


def derived_fields(x, y, cons, params):
    """Physical fields at each vertex from the averaged conserved state.

    Primitives are formed AFTER averaging the conserved variables, which is the
    order the solver's own reconstruction uses; forming them per cell node and
    then averaging would divide small nearly-cancelling momenta by small
    densities first and average the noise.
    """
    me, mi, q = params['ME'], params['MI'], params['Q']
    r, ct, st = rd.polar(x, y)
    rho_e = cons[RHO_E]
    rho_i = cons[RHO_I]
    n_e = rho_e / me
    n_i = rho_i / mi
    with np.errstate(invalid='ignore', divide='ignore'):
        u_ex = cons[MOM_EX] / rho_e
        u_ey = cons[MOM_EY] / rho_e
    e_x, e_y = cons[E_X], cons[E_Y]
    fields = {
        'rho_c': q * (n_i - n_e),
        'phi': cons[PHI],
        'E_r': e_x * ct + e_y * st,
        'E_theta': -e_x * st + e_y * ct,
        'n_e': n_e,
        'n_i': n_i,
        'u_e_theta': -u_ex * st + u_ey * ct,
    }
    return r, fields


def find_rings(r, rlo, rhi):
    """Group vertex radii in [rlo, rhi] into concentric node rings.

    A generated disc mesh puts its nodes on concentric rings; radii within a
    small fraction of the local ring spacing belong to one ring. Returns a list
    of (r_mean, node_index_array), coarsest exclusion applied by the caller.
    """
    sel = np.where((r >= rlo) & (r <= rhi))[0]
    if sel.size == 0:
        return []
    order = sel[np.argsort(r[sel])]
    rs = r[order]
    gaps = np.diff(rs)
    spacing = np.median(gaps[gaps > 0]) if np.any(gaps > 0) else 0.0
    rings = []
    start = 0
    for i in range(1, len(rs)):
        if spacing > 0 and (rs[i] - rs[i - 1]) > RING_TOL * spacing:
            rings.append(order[start:i])
            start = i
    rings.append(order[start:])
    return [(float(np.mean(r[idx])), idx) for idx in rings]


def harmonics(theta, values, mmax):
    """A_m = mean_j values_j exp(-i m theta_j) for m = 0..mmax (nonuniform)."""
    good = np.isfinite(values)
    th = theta[good]
    v = values[good]
    if th.size == 0:
        return np.full(mmax + 1, np.nan, dtype=complex), 0
    ms = np.arange(mmax + 1)
    phase = np.exp(-1j * np.outer(ms, th))
    return (phase @ v) / th.size, th.size


def ring_report(theta, values, mmax):
    """Per-ring spectral summary for one component on one ring."""
    A, n = harmonics(theta, values, mmax)
    if n == 0:
        return None
    # Reconstruct the smooth (|m| <= LOW_M) part at each node and subtract it;
    # the r.m.s. of the remainder is the grid-scale content the harmonics miss.
    ms = np.arange(mmax + 1)
    low = ms <= LOW_M
    good = np.isfinite(values)
    th = theta[good]
    recon = np.real(A[low][0]) + 2.0 * np.real(
        np.exp(1j * np.outer(th, ms[low][1:])) @ A[low][1:])
    resid = values[good] - recon
    grid_rms = float(np.sqrt(np.mean(resid ** 2)))
    return {
        'n': n,
        'A': A,
        'm0': float(abs(A[0])),
        'low': [float(2.0 * abs(A[m])) for m in range(1, min(9, mmax + 1))],
        'grid_rms': grid_rms,
    }


def peak_locus(data, comp):
    """(r, theta_deg, |value|) of the cell whose peak-over-nodes |value| is largest."""
    cells = np.asarray(data['connectivity']).reshape(-1, 3)
    points = np.asarray(data['Position'])[:, :2]
    stack = np.vstack([np.abs(np.asarray(data[f'solutiondg.{node * 18 + comp}']))
                       for node in range(3)])
    cellmax = np.nanmax(stack, axis=0)
    c = int(np.nanargmax(cellmax))
    cen = points[cells[c]].mean(axis=0)
    return (float(math.hypot(*cen)), float(math.degrees(math.atan2(cen[1], cen[0]))),
            float(cellmax[c]))


def _r2(y, yhat):
    y = np.asarray(y, dtype=float)
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else float('nan')


def growth(ts, ys):
    """Fit y(t) over finite positive points and say whether it is exponential.

    Returns (exp_rate, n, shape) where exp_rate is the slope of ln y (an
    exponential rate, 1/s) and shape is 'exp', 'lin' or 'flat' from comparing an
    exponential fit's R^2 with a straight line's. A linearly growing quantity
    fitted as an exponential returns a plausible-looking rate; naming the shape
    stops that rate being mistaken for an instability growth rate.
    """
    ts = np.asarray(ts, dtype=float)
    ys = np.asarray(ys, dtype=float)
    ok = np.isfinite(ys) & (ys > 0)
    k = int(np.count_nonzero(ok))
    if k < MIN_FIT:
        return float('nan'), k, 'n/a'
    t, y = ts[ok], ys[ok]
    ln = np.polyfit(t, np.log(y), 1)
    exp_rate = float(ln[0])
    r2_exp = _r2(np.log(y), np.polyval(ln, t))
    lin = np.polyfit(t, y, 1)
    r2_lin = _r2(y, np.polyval(lin, t))
    if abs(y[-1] / y[0]) < 1.5:
        shape = 'flat'
    elif r2_lin > r2_exp + 0.02:
        shape = 'lin'
    else:
        shape = 'exp'
    return exp_rate, k, shape


def band_bounds(name, params):
    """(rlo, rhi) for a named radial band, from the deck's geometry."""
    rc, w = params['COIL_R'], params['COIL_W']
    wall = params['WALL_RADIUS']
    if name == 'winding':
        return rc - w, rc + w
    if name == 'column':
        return 0.5 * params['RAD_PLASMA'], params['RAD_PLASMA']
    if name == 'wall':
        return 0.95 * wall, wall
    raise SystemExit(f'unknown band "{name}"; choose from winding, column, wall')


def frame_index(path):
    base = os.path.basename(path)
    stem = base[:-4] if base.endswith('.vtu') else base
    if '_' not in stem or not stem.rsplit('_', 1)[1].isdigit():
        raise SystemExit(f'{base}: frame files must end in _<number>.vtu')
    return int(stem.rsplit('_', 1)[1])


def main(argv=None):
    if np is None:
        raise SystemExit('ring_spectrum needs numpy')
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('deck')
    ap.add_argument('frames', nargs='+')
    ap.add_argument('--bands', default='winding',
                    help='comma-separated: winding, column, wall (default winding)')
    ap.add_argument('--mmax', type=int, default=60)
    ap.add_argument('--components', default='rho_c,phi,E_r,E_theta,n_e,n_i,u_e_theta')
    args = ap.parse_args(argv)

    params = rd.deck_parameters(args.deck)
    for need in ('omega', 'ME', 'MI', 'Q', 'COIL_R', 'COIL_W', 'WALL_RADIUS', 'RAD_PLASMA'):
        if need not in params:
            raise SystemExit(f'{args.deck}: deck is missing {need}, needed for ring_spectrum')
    omega = 2.0 * math.pi * params['omega']          # deck omega is a frequency in Hz
    comps = args.components.split(',')
    bands = args.bands.split(',')

    frames = sorted(args.frames, key=frame_index)
    if len(frames) < 2:
        raise SystemExit('need at least two frames to see structure grow')
    tend = params.get('TEND')
    nout = params.get('OUT')
    if tend is None or nout is None:
        raise SystemExit(f'{args.deck}: deck is missing TEND/OUT, needed to time the frames')
    dt_frame = tend / nout

    series = {}                                       # (band, ring_r, comp) -> [(t, rep)]
    peak_series = {name: [] for name in ('e-rho', 'phi', 'E_x', 'E_y')}
    ring_hint = None
    for path in frames:
        idx = frame_index(path)
        t = idx * dt_frame
        data = vtu.read(path)
        x, y, cons, _ = vertex_values(data, list(range(18)))
        r, fields = derived_fields(x, y, cons, params)
        _, ct, st = rd.polar(x, y)
        theta = np.arctan2(y, x)
        for band in bands:
            rlo, rhi = band_bounds(band, params)
            rings = find_rings(r, rlo, rhi)
            rings = [ring for ring in rings if ring[1].size >= MIN_RING]
            if band == 'winding' and not rings and ring_hint is None:
                ring_hint = (rlo, rhi)
            for rmean, idxs in rings:
                for comp in comps:
                    rep = ring_report(theta[idxs], fields[comp][idxs], args.mmax)
                    if rep is not None:
                        series.setdefault((band, round(rmean, 5), comp), []).append((t, rep))
        # peak loci are of the raw stored components, named honestly, so a
        # fixed-azimuth peak (a mesh seam) can be read straight off.
        cmap = {'e-rho': RHO_E, 'phi': PHI, 'E_x': E_X, 'E_y': E_Y}
        for name in peak_series:
            peak_series[name].append((t, peak_locus(data, cmap[name])))

    if not series:
        hint = ''
        if ring_hint is not None:
            hint = (f'\nNo node ring in the winding band [{ring_hint[0]:.4f}, '
                    f'{ring_hint[1]:.4f}] had >= {MIN_RING} nodes. On a coarse mesh '
                    'the winding is too sparsely sampled to Fourier; run at k=1 or '
                    'pass --bands column.')
        raise SystemExit('ring_spectrum found no usable ring.' + hint)

    ramp_clear = RAMP_CLEAR * params.get('RISE', 0.0)
    print(f'# ring_spectrum: {len(frames)} frames, omega={omega:.4e} rad/s, '
          f'frame dt={dt_frame:.4e} s, ramp_clear={ramp_clear:.3e} s')
    for (band, rmean, comp), pts in sorted(series.items()):
        ts = [t for t, _ in pts]
        n = pts[0][1]['n']
        m0 = [rep['m0'] for _, rep in pts]
        grid = [rep['grid_rms'] for _, rep in pts]
        past = [i for i, t in enumerate(ts) if t >= ramp_clear]
        s_m0, k0, sh0 = growth([ts[i] for i in past], [m0[i] for i in past])
        s_gr, kg, shg = growth([ts[i] for i in past], [grid[i] for i in past])
        # dominant low harmonic 3<=m<=30 in the last frame, and its phase drift
        last = pts[-1][1]['A']
        band_hi = min(30, args.mmax)
        if band_hi >= 3:
            dom = 3 + int(np.nanargmax([abs(last[m]) for m in range(3, band_hi + 1)]))
            drift = []
            for t, rep in pts:
                drift.append(np.angle(rep['A'][dom]) - dom * omega * t)
            drift = np.unwrap(drift)
            dphi = float((drift[-1] - drift[0]) / (ts[-1] - ts[0])) if ts[-1] > ts[0] else float('nan')
        else:
            dom, dphi = -1, float('nan')
        print(f'{band:7s} r={rmean:.5f} N={n:4d} {comp:9s} '
              f'm0={s_m0:+.2e}/s[{sh0}] grid={s_gr:+.2e}/s[{shg}] '
              f'm0:{m0[0]:.2e}->{m0[-1]:.2e} grid:{grid[0]:.2e}->{grid[-1]:.2e} '
              f'dom_m={dom} dphi/omega={dphi / omega if omega else float("nan"):+.2f}')

    print('# peak loci (r, theta_deg) per frame:')
    for comp, pts in peak_series.items():
        locs = ' '.join(f'({r:.4f},{th:+06.1f})' for _, (r, th, _) in pts)
        print(f'{comp:9s}: {locs}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
