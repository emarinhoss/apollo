#!/usr/bin/env python3
"""RMF current-drive diagnostics: the quantities Phase 3 compares with the literature.

docs/rmf-frc-model-assessment.md, Phase 3, asks for a set of comparisons against
the rotating-magnetic-field literature - the penetration threshold, the
penetrated-limit driven field, the skin-depth-limited current layer, and, before
any of those, which way the drive turns. Each of them is a number extracted from
a run. This module extracts them.

It is split in two on purpose. The functions in the first half take plain arrays
and parameters and return numbers; they know nothing about files, decks, or
Apollo. Everything they compute has a closed form for some field, which is how
test/test_rmf_diagnostics.py checks them - in under a second, without a solver.
The second half reads a deck and a set of .vtu frames and calls them. A
diagnostic that has never been checked against a known answer is not evidence,
and these are the numbers a validation phase would rest on.

    python3 scripts/rmf_diagnostics.py <deck.pin> <case_0.vtu> [more .vtu ...]

SIGN CONVENTIONS, stated because Phase 3 item 0 is entirely a question of sign.

  * (x, y, z) is right-handed and +z is out of the r-theta plane. The 2-D runs
    are a slice through an infinite cylinder whose axis is +z.
  * theta increases counter-clockwise from +x; e_theta = (-sin theta, cos theta)
    points counter-clockwise.
  * A counter-clockwise current loop (J_theta > 0) produces B along +z inside
    itself. So a driven current REVERSES a +z bias field only if J_theta < 0.
  * J = sum over species of q_s n_s u_s, so for electrons J_theta < 0 means
    u_e,theta > 0: electrons dragged COUNTER-CLOCKWISE.
  * In the synchronous limit electrons are dragged in the sense the RMF turns.
    Field reversal against a +z bias therefore needs a COUNTER-CLOCKWISE RMF.

That chain is what makes the rotation sense of a drive a physical claim rather
than a convention, and it is asserted in test/cxx/test_rmf_rotation_sense.cc.

WHY THE CURRENT COMES FROM THE FLUID MOMENTA, NOT FROM curl B. Both are
available. J = sum_s (q_s/m_s)(rho u)_s is local, exact, and needs no
derivatives - it is the same expression WxCurrentSrc feeds back into Ampere's
law (src/hyperapps/multifluid/wxcurrentsrc.h), so it is the current the run
actually used. Taking curl of a P1 discontinuous field would differentiate
across element faces where the solution is not continuous.

READ CYCLE AVERAGES, NOT INSTANTANEOUS VALUES. Every quantity here except the
profiles carries an oscillation at omega (and at 2*omega, since the ponderomotive
and J x B terms are quadratic in the drive) on top of the secular part that the
literature's expressions describe. scripts/rmf_antenna_power.py learned this the
expensive way for the antenna power. The helpers below therefore report a cycle
average when the frames span at least one RMF period and say so plainly when
they do not.
"""

import argparse
import math
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'test'))

try:
    import numpy as np
except ImportError:  # pragma: no cover
    np = None

# numpy 2 renamed trapz to trapezoid and removed the old name. CI installs
# numpy from pip, so it gets 2.x; the repository's own environment may not.
_trapezoid = None
if np is not None:
    _trapezoid = getattr(np, 'trapezoid', None) or getattr(np, 'trapz', None)


# Component indices in the 18-element two-fluid state. Confirmed against the
# InpRange/OutRange blocks of examples/unstructuredDG/multifluid/rmf_frc/frc2d.pin.
RHO_E, MOM_EX, MOM_EY, MOM_EZ, EN_E = 0, 1, 2, 3, 4
RHO_I, MOM_IX, MOM_IY, MOM_IZ, EN_I = 5, 6, 7, 8, 9
E_X, E_Y, E_Z = 10, 11, 12
B_X, B_Y, B_Z = 13, 14, 15
PHI, PSI = 16, 17

MU0 = 4.0e-7 * math.pi


# ---------------------------------------------------------------------------
# The numeric core. Arrays in, numbers out; no files, no deck, no Apollo.
# ---------------------------------------------------------------------------

def polar(x, y):
    """(r, cos theta, sin theta) for nodal coordinates.

    Returned as the trig rather than the angle because every caller wants the
    unit vectors, and going through atan2 and back costs accuracy near the axis.
    On the axis itself r is zero and the direction is undefined; cos and sin are
    set to zero there so that any projection through them is zero rather than a
    NaN. Callers must mask r = 0 out rather than relying on that.
    """
    r = np.hypot(x, y)
    safe = np.where(r > 0.0, r, 1.0)
    return r, x / safe, y / safe


def current_density(comps, m_e, m_i, q):
    """(J_x, J_y, J_z) from the fluid momenta.

    J = sum over species of (q_s/m_s)(rho u)_s, with the electron charge -q and
    the ion charge +q, matching the two `currents` source blocks in the deck.
    """
    jx = (-q / m_e) * comps[MOM_EX] + (q / m_i) * comps[MOM_IX]
    jy = (-q / m_e) * comps[MOM_EY] + (q / m_i) * comps[MOM_IY]
    jz = (-q / m_e) * comps[MOM_EZ] + (q / m_i) * comps[MOM_IZ]
    return jx, jy, jz


def azimuthal(vx, vy, cos_t, sin_t):
    """The e_theta component of a planar vector field: counter-clockwise positive."""
    return -vx * sin_t + vy * cos_t


def radial(vx, vy, cos_t, sin_t):
    """The e_r component of a planar vector field: outward positive."""
    return vx * cos_t + vy * sin_t


def radial_profile(r, values, edges):
    """Bin `values` by radius. Returns (centres, mean, std, count).

    Empty bins come back as NaN in mean and std and 0 in count, rather than
    being dropped, so that profiles from different frames stay index-aligned and
    can be averaged over a cycle.
    """
    edges = np.asarray(edges, dtype=float)
    centres = 0.5 * (edges[:-1] + edges[1:])
    idx = np.digitize(r, edges) - 1
    n = len(centres)
    mean = np.full(n, np.nan)
    std = np.full(n, np.nan)
    count = np.zeros(n, dtype=int)
    for b in range(n):
        sel = values[idx == b]
        count[b] = sel.size
        if sel.size:
            mean[b] = float(np.mean(sel))
            std[b] = float(np.std(sel))
    return centres, mean, std, count


def axial_field_on_axis(r, bz, r_axis):
    """Mean B_z inside r < r_axis.

    The meshes here have no node exactly on the axis, and one node would be a
    noisy estimate in any case, so this averages a small disc. r_axis has to be
    small against the scale on which the driven field varies and large enough to
    hold enough nodes; the caller is told the count so it can judge.
    """
    sel = r < r_axis
    if not np.any(sel):
        return float('nan'), 0
    return float(np.mean(bz[sel])), int(np.count_nonzero(sel))


def rotation_parameter(r, u_theta, omega, r_min, r_max):
    """zeta = u_e,theta / (omega r), the electron rotation parameter.

    zeta = 1 is synchronous rotation - electrons carried around with the RMF -
    and zeta = 0 is no drive. It is the quantity the penetrated-limit expression
    for the driven field is written in.

    Evaluated only for r_min < r < r_max: the denominator vanishes on the axis,
    and the outermost cells sit in the boundary condition's own layer.
    """
    sel = (r > r_min) & (r < r_max)
    if not np.any(sel):
        return float('nan'), 0
    return float(np.mean(u_theta[sel] / (omega * r[sel]))), int(np.count_nonzero(sel))


def penetrated_limit_bz(r_centres, zeta_profile, n_e, q, omega, a):
    """B_z on axis for a measured zeta(r), in the fully-driven (penetrated) limit.

    Each shell of thickness dr carries a surface current J_theta dr, and inside
    an infinite cylinder every shell outside the field point contributes
    mu0 J_theta dr uniformly, so

        B_z(0) = mu0 * integral from 0 to a of J_theta(r) dr
               = -mu0 n e omega * integral from 0 to a of zeta(r) r dr

    using J_theta = -e n u_theta = -e n omega r zeta. With zeta constant this is
    the -mu0 n e omega a^2 / 2 * zeta quoted in the assessment; the integral form
    is used here so that a measured, non-uniform zeta(r) can be fed in.

    The sign is the point: zeta > 0 (electrons counter-clockwise) gives
    B_z < 0, which OPPOSES a +z bias. See the module docstring.

    NaN bins - radii the run put no nodes in - are skipped, not treated as zero.
    """
    r_centres = np.asarray(r_centres, dtype=float)
    zeta_profile = np.asarray(zeta_profile, dtype=float)
    good = np.isfinite(zeta_profile) & (r_centres <= a)
    if np.count_nonzero(good) < 2:
        return float('nan')
    integral = _trapezoid(zeta_profile[good] * r_centres[good], r_centres[good])
    return float(-MU0 * n_e * q * omega * integral)


def penetration_fraction(r, b_perp, a, inner=0.3, outer=0.9):
    """|B_perp| averaged over an inner disc, divided by its value near the edge.

    The RMF literature's "penetrated" state has the transverse field reaching
    the axis essentially undiminished; the screened state has it confined to a
    layer of order delta at the edge. This ratio is 1 for the first and near 0
    for the second. It is a shape measure, deliberately normalised by the field
    the run itself has at the edge rather than by the deck's B_omega, so that it
    says whether the field got in rather than how big the drive was.

    IT CAN EXCEED 1, AND THAT IS NOT A PENETRATED STATE. Before the run settles,
    the transverse field arrives as a wave and rings the domain, so the interior
    can transiently hold more field than the edge: the shipped deck reads 1.20,
    0.63 and 1.24 on successive frames over its first 0.15 of an RMF period.
    Only a value cycle-averaged over a whole period means what this docstring
    says; anything else is a snapshot of a transient. See cycle_average, and the
    warning main() prints when the frames do not span a period.
    """
    inner_sel = r < inner * a
    outer_sel = (r > outer * a) & (r <= a)
    if not np.any(inner_sel) or not np.any(outer_sel):
        return float('nan')
    edge = float(np.mean(b_perp[outer_sel]))
    if edge == 0.0:
        return float('nan')
    return float(np.mean(b_perp[inner_sel]) / edge)


def current_layer_thickness(r_centres, j_theta, a, floor_fraction=1e-3,
                            min_r2=0.9, max_fraction=0.5):
    """e-folding length of |J_theta| inward from the edge, to compare with delta.

    In the skin-depth-limited state the driven current sits in a layer at the
    plasma edge whose thickness is the resistive skin depth
    delta = sqrt(2 eta / (mu0 omega)). This fits log|J_theta| against (a - r)
    over the outer half of the profile and returns the decay length.

    IT RETURNS NaN RATHER THAN A NUMBER WHENEVER THE PROFILE IS NOT A SKIN
    LAYER, and the three ways that is detected are all load-bearing. Run
    unguarded against the first 1.6% of an RMF period of the shipped deck, an
    earlier version of this function reported layers of 40 mm and 59 mm in a
    30 mm column - fits to a transient that was not decaying at all. A number
    like that reads as a measurement and would have been compared against
    delta.

      * a non-negative slope: the current grows inward, so there is no layer.
        This is what a PENETRATED run looks like and is not a failure.
      * a decay length longer than max_fraction of the column: a layer that
        does not fit inside the plasma is not a layer. The default of half the
        radius is generous - the deck's own delta is a/24.
      * a poor fit: R^2 below min_r2 means the profile is not an exponential,
        so its "decay length" describes nothing.
    """
    r_centres = np.asarray(r_centres, dtype=float)
    j = np.abs(np.asarray(j_theta, dtype=float))
    good = np.isfinite(j) & (r_centres <= a) & (r_centres > 0.5 * a)
    if np.count_nonzero(good) < 3:
        return float('nan')
    peak = np.nanmax(j[good])
    if not np.isfinite(peak) or peak <= 0.0:
        return float('nan')
    # Keep the part of the profile that is above the noise floor; taking a log
    # of a value that has decayed into round-off fits the round-off.
    good &= j > floor_fraction * peak
    if np.count_nonzero(good) < 3:
        return float('nan')
    depth = a - r_centres[good]
    logj = np.log(j[good])
    slope, intercept = np.polyfit(depth, logj, 1)
    if slope >= 0.0:
        return float('nan')            # grows inward: not a skin layer

    length = -1.0 / slope
    if length > max_fraction * a:
        return float('nan')            # wider than the column: not a layer

    resid = logj - (slope * depth + intercept)
    spread = float(np.sum((logj - np.mean(logj)) ** 2))
    if spread <= 0.0:
        return float('nan')
    r2 = 1.0 - float(np.sum(resid ** 2)) / spread
    if r2 < min_r2:
        return float('nan')            # not an exponential: the length means nothing
    return float(length)


def cycle_average(times, values, period):
    """Mean of `values` over the last whole RMF period covered by `times`.

    Returns (average, spans_a_period). When the samples do not cover a period
    the plain mean is returned with the flag False, and the caller must say so:
    an average over a fraction of a cycle is dominated by whichever part of the
    cycle it happened to catch, which for these diagnostics is the difference
    between a driven current and its reactive swing.
    """
    times = np.asarray(times, dtype=float)
    values = np.asarray(values, dtype=float)
    if times.size == 0:
        return float('nan'), False
    span = float(times[-1] - times[0])
    if span < period:
        return float(np.mean(values)), False
    sel = times >= times[-1] - period
    return float(np.mean(values[sel])), True


def skin_depth(eta, omega):
    """delta = sqrt(2 eta / (mu0 omega)), the classical resistive skin depth."""
    return math.sqrt(2.0 * eta / (MU0 * omega))


# ---------------------------------------------------------------------------
# The shell: deck parsing and frame reading.
# ---------------------------------------------------------------------------

# Everything a diagnostic or a cost estimate needs out of a deck. Named here so
# that a deck missing one of them fails at parse time with a list, rather than
# with a KeyError three functions later.
REQUIRED_DECK_KEYS = ('MU0', 'Q', 'ME', 'MI', 'omega', 'RAD_PLASMA',
                      'Te', 'LIGHT', 'TEND', 'ETA', 'n_dens', 'OUT')


def deck_parameters_from_text(text, what='deck'):
    """Evaluate the Python macro block of a .pin and return its namespace.

    Read from the deck rather than passed on the command line so that a
    diagnostic cannot silently describe a run with another run's numbers. The
    macro section of a .pin is Python, so the assignments are evaluated in an
    empty namespace with math available - the same thing scripts/wxinpparse.py
    does before handing the file to the solver.

    Lines that do not evaluate are skipped rather than fatal: the macro block
    also holds imports and expressions that depend on things only the solver
    defines, and none of those are wanted here.
    """
    macro = text.split('<warpx>')[0]     # everything before it is the macro block
    ns = {'math': math, 'PI': math.pi}
    for line in macro.splitlines():
        line = line.split('#')[0].strip()
        if not re.match(r'^[A-Za-z_][A-Za-z0-9_]*\s*=', line):
            continue
        try:
            exec(line, ns)          # noqa: S102 - the deck is the input, by design
        except Exception:
            continue
    missing = [k for k in REQUIRED_DECK_KEYS if k not in ns]
    if missing:
        raise SystemExit(f'{what}: deck is missing {", ".join(missing)}; '
                         f'this is not an rmf_frc deck')
    ns.setdefault('P_ORDER', 1)
    return ns


def deck_parameters(path):
    """deck_parameters_from_text for a file on disk."""
    with open(path) as fh:
        return deck_parameters_from_text(fh.read(), what=path)


def frame_time(index, params):
    """The time of output frame `index`.

    The solver writes one frame before the time loop and one per output
    interval, so frame i of a complete run is at TEND*i/OUT
    (src/solvers/apsolver.cc:186, 216).

    Taken from the deck's OUT, NOT from how many frames happen to be present.
    Dividing by the number of frames found is right only for a run that
    finished, and an unfinished run is the normal case here: the solver has no
    checkpoint/restart, so anything interrupted leaves a partial set of frames
    that are still perfectly good to analyse. Timing them by their own count
    would stretch them to fill TEND and silently misreport every time, the
    cycle-average window included.
    """
    tend, out = params.get('TEND'), params.get('OUT')
    if tend is None or not out:
        return float('nan')
    return tend * index / out


def analyse(frames, params, nbins=40):
    """Per-frame diagnostics for a list of (index, {name: array}) frames."""
    import deckrun

    q, m_e, m_i = params['Q'], params['ME'], params['MI']
    omega = 2.0 * math.pi * params['omega']
    a = params['RAD_PLASMA']
    edges = np.linspace(0.0, a, nbins + 1)

    wanted = (RHO_E, MOM_EX, MOM_EY, MOM_EZ,
              RHO_I, MOM_IX, MOM_IY, MOM_IZ, B_X, B_Y, B_Z)

    rows = []
    for index, data in frames:
        x, y, comps = deckrun.nodal(data, wanted)
        r, cos_t, sin_t = polar(x, y)

        # The TOTAL current, both species. The ions carry little of it at these
        # frequencies - omega is 8.7 ion gyrofrequencies - but "little" is a
        # result, not an assumption to build into the instrument.
        jx, jy, _jz = current_density(comps, m_e, m_i, q)
        j_theta = azimuthal(jx, jy, cos_t, sin_t)

        rho_e = comps[RHO_E]
        safe_rho = np.where(rho_e > 0.0, rho_e, np.nan)
        u_theta = azimuthal(comps[MOM_EX] / safe_rho, comps[MOM_EY] / safe_rho,
                            cos_t, sin_t)

        b_perp = np.hypot(comps[B_X], comps[B_Y])
        bz_axis, n_axis = axial_field_on_axis(r, comps[B_Z], 0.1 * a)
        zeta, n_zeta = rotation_parameter(r, u_theta, omega, 0.1 * a, 0.9 * a)
        centres, j_prof, _s, _c = radial_profile(r, j_theta, edges)
        _c2, zeta_prof, _s2, _c3 = radial_profile(
            r, np.where(r > 0, u_theta / (omega * np.where(r > 0, r, 1.0)), np.nan),
            edges)

        rows.append(dict(
            index=index,
            bz_axis=bz_axis, n_axis=n_axis,
            zeta=zeta, n_zeta=n_zeta,
            penetration=penetration_fraction(r, b_perp, a),
            layer=current_layer_thickness(centres, j_prof, a),
            j_theta_mean=float(np.nanmean(j_prof)),
            bz_predicted=penetrated_limit_bz(
                centres, zeta_prof, params.get('n_dens', float('nan')), q, omega, a),
        ))
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('deck', help='the .pin the run came from')
    parser.add_argument('frames', nargs='+', help='.vtu output frames')
    args = parser.parse_args(argv)

    if np is None:
        raise SystemExit('numpy is required')
    import vtu

    params = deck_parameters(args.deck)
    period = 1.0 / params['omega']
    names = sorted(args.frames, key=lambda p: int(os.path.basename(p)
                                                  .rsplit('_', 1)[1].split('.')[0]))
    frames = [(int(os.path.basename(p).rsplit('_', 1)[1].split('.')[0]), vtu.read(p))
              for p in names]

    rows = analyse(frames, params)
    times = [frame_time(r['index'], params) for r in rows]

    delta = skin_depth(params['ETA'], 2.0 * math.pi * params['omega'])
    spanned = (times[-1] - times[0]) / period
    print(f'RMF period {period:.4e} s; frames span {spanned:.3f} periods')
    if spanned < 1.0:
        # Said here as well as beside the averages below, because the per-frame
        # table is what gets read and copied out of.
        print('*** THESE FRAMES DO NOT SPAN AN RMF PERIOD. Every column below is a\n'
              '*** snapshot of a transient, not a driven state: the penetration ratio\n'
              '*** can exceed 1 while the field rings, and zeta has not had a period\n'
              '*** in which to be driven. Do not compare any of it with the\n'
              '*** literature, which describes cycle-averaged steady states.')
    print(f'resistive skin depth delta = {delta * 1e3:.3f} mm, '
          f'lambda = a/delta = {params["RAD_PLASMA"] / delta:.1f}')
    print()
    print(f'{"t [s]":<12}{"B_z axis [T]":<15}{"zeta":<12}'
          f'{"penetration":<14}{"layer [mm]":<12}')
    for t, row in zip(times, rows):
        layer = row['layer']
        # Not "penetrated": the fit is refused for three different reasons (see
        # current_layer_thickness) and only one of them is penetration. The
        # penetration column beside it is what says which.
        layer_s = f'{layer * 1e3:.3f}' if math.isfinite(layer) else 'no layer'
        print(f'{t:<12.4e}{row["bz_axis"]:<15.6e}{row["zeta"]:<12.4f}'
              f'{row["penetration"]:<14.4f}{layer_s:<12}')

    print()
    for key, label in (('bz_axis', 'B_z on axis [T]'), ('zeta', 'zeta')):
        avg, whole = cycle_average(times, [r[key] for r in rows], period)
        note = 'cycle average' if whole else ('MEAN OVER A PARTIAL CYCLE - '
                                              'not a cycle average, do not '
                                              'compare it with one')
        print(f'{label}: {avg:.6e}   ({note})')
    return 0


if __name__ == '__main__':
    sys.exit(main())
