#!/usr/bin/env python3
"""Parameter scans for the RMF example, and what they cost before you start one.

Phase 3 of docs/rmf-frc-model-assessment.md asks for scans: item 1 sweeps the
RMF amplitude B_omega across the penetration threshold, item 3 looks at the
low-amplitude end. Each point is a full run of the two-fluid deck. This script
runs such a scan, applies scripts/rmf_diagnostics.py to each point, and collects
the result.

ITS FIRST JOB IS TO REFUSE TO SURPRISE YOU. `--dry-run` prints what the scan
would cost and stops. That is not a courtesy: one RMF period of the shipped deck
is about 46,000 steps and takes hours, a scan is several runs of several
periods, and the solver has no checkpoint/restart (docs/known-issues.md section
6), so a run that is interrupted is a run that is lost. The cost model below is
measured, not guessed, and is checked against real timings in
test/test_rmf_scan.py.

    python3 scripts/rmf_scan.py deck.pin --param Bomega \\
        --values 10e-4,20e-4,30e-4,50e-4 --periods 5 --dry-run

Drop --dry-run to run it. Points already completed in the output directory are
skipped, so an interrupted scan can be restarted and will only redo what is
missing. That is the only form of restart available here.

THE COST MODEL

Apollo's timestep is, from src/subsolvers/wxnodaldg2dmethod.cc:225,

    dt = (2/3) * cfl * dtscale * rMin / maxSpeed

with dtscale the smallest inscribed-circle radius over the mesh's triangles
(src/lib/wxnodaldggeometry2d.cc:274), rMin the node spacing on the reference
element, and maxSpeed the largest characteristic speed in the state.

The part worth knowing is what maxSpeed turns out to be. For the shipped
hydrogen deck it is NOT dominated by the reduced speed of light, though it looks
as though it should be: LIGHT = 3.0e6 m/s while the electron sound speed
sqrt(gamma k T_e / m_e) at 30 eV is 2.966e6 m/s, 1.1% below it. Measured, over a
scan of LIGHT from 3.0e6 down to 1.0e5, dt rises by 1.16% and then does not move
again - it saturates on the electrons. So the usual lever for making a reduced-c
two-fluid run affordable is worth about one per cent here, and the cost of Phase
3 is set by the electron fluid, not by the field solver. See the assessment.

The lever that does work is the pair (T_e, c) together: every dimensionless
parameter the RMF penetration literature is written in - gamma, lambda,
omega/omega_ci, omega/omega_ce, B_omega/B_bias, the Debye length in cells, and
omega_pe*dt - is invariant under T_e -> T_e/s^2 with c -> c/s, while dt grows by
s. Only beta moves, by 1/s^2. `--speedup` applies that scaling and says what it
did to beta, because a scan that quietly changed the plasma it was scanning
would be worthless.
"""

import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'scripts'))
sys.path.insert(0, os.path.join(REPO, 'test'))

MU0 = 4.0e-7 * math.pi
GAMMA = 5.0 / 3.0

# Measured on the 4-core container this was developed on: 731 steps of the
# shipped deck on optimizedCircle2.msh (7792 triangles) in 339 s of wall clock,
# with the machine otherwise idle. Per triangle so it can be scaled to another
# mesh.
#
# TREAT THE STEP COUNT AS EXACT AND THE HOURS AS AN ORDER-OF-MAGNITUDE. The
# timestep model below reproduces the solver's dt to six significant figures, so
# the step count is exact; this constant is not. The same binary on the same
# mesh measured 0.464 s/step idle, 0.569 under moderate load and 0.720 under
# heavy load on the same afternoon - a spread of 79%. Use --sec-per-step to
# substitute a number measured on the machine that will actually run the scan.
#
# There is no separate fixed-cost term because there is nothing to put in it:
# timed directly, a 3-step run of this deck takes 2.6 s and a 41-step run 30 s,
# which fits a setup cost of 0.4 s. Anything that looks like a large fixed cost
# in a two-point fit is load variation between the two points.
SECONDS_PER_STEP_PER_TRIANGLE = 339.0 / 731.0 / 7792.0


def mesh_dtscale(path):
    """Smallest inscribed-circle radius, min over triangles of area/(perimeter/2).

    The same quantity src/lib/wxnodaldggeometry2d.cc:274 computes, recomputed
    here rather than read out of a run so that a scan can be costed before any
    run exists.
    """
    with open(path) as fh:
        text = fh.read()
    n_nodes = int(re.search(r'\$Nodes\n(\d+)', text).group(1))
    nodes = {}
    for line in text.split('$Nodes\n')[1].split('\n')[1:n_nodes + 1]:
        p = line.split()
        nodes[int(p[0])] = (float(p[1]), float(p[2]))
    n_elem = int(re.search(r'\$Elements\n(\d+)', text).group(1))
    best, count = float('inf'), 0
    for line in text.split('$Elements\n')[1].split('\n')[1:n_elem + 1]:
        p = line.split()
        if len(p) < 2 or p[1] != '2':
            continue
        a, b, c = (nodes[int(i)] for i in p[-3:])
        area = abs(0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])))
        per = sum(math.hypot(u[0] - v[0], u[1] - v[1])
                  for u, v in ((a, b), (b, c), (c, a)))
        best = min(best, area / (per / 2.0))
        count += 1
    return best, count


def max_wave_speed(params):
    """The larger of the reduced speed of light and the electron sound speed.

    Those are the two fast waves of the 18-component system AT REST. The ion
    sound speed is smaller by sqrt(m_i/m_e) = 43 and the Alfven speed by 229, so
    neither can set the step.

    THIS IS AN ESTIMATE FOR A COLD START AND IT DRIFTS LOW AS THE RUN PROCEEDS.
    Apollo recomputes dt every step (wxnodaldg2dmethod.cc:553) from the Euler
    flux's Lax-Friedrichs speed, which is |u| + sqrt(gamma p/rho) - the flow
    speed ADDS to the sound speed rather than competing with it
    (wxeulereqn.cc:1110). So a driven run gets slower as it spins up: at
    synchronous electron rotation, u_theta = omega*a = 1.50e5 m/s at the plasma
    edge, so the electron characteristic reaches 3.116e6 against the deck's
    LIGHT of 3.0e6 and dt falls 3.7%. Cost estimates from this function are
    therefore LOWER bounds: by about that much for a fully penetrated run, and
    by nothing at all for a run that stays screened.

    That same arithmetic is why the shipped deck has a physics problem rather
    than only a cost one: 3.12e6 is above its LIGHT of 3.0e6. See the assessment,
    section 3.10.
    """
    c_se = math.sqrt(GAMMA * params['Te'] * params['Q'] / params['ME'])
    return max(params['LIGHT'], c_se), c_se


def timestep(params, dtscale, p_order):
    """dt as the solver computes it. Verified against runs in test/test_rmf_scan.py."""
    if p_order != 1:
        raise SystemExit(
            f'P_ORDER = {p_order}: the reference-element node spacing rMin was only '
            f'measured for P_ORDER = 1 (where it is 2), and this model would be a '
            f'guess at any other order. Measure it before trusting a cost estimate.')
    r_min = 2.0
    cfl = 1.0 / (2.0 * (p_order + 1) - 1)
    v_max, _c_se = max_wave_speed(params)
    return (2.0 / 3.0) * cfl * dtscale * r_min / v_max


def estimate(params, mesh_path, seconds_per_step=None):
    """(dt, steps, wall-clock seconds, notes) for one run of this deck."""
    dtscale, n_tri = mesh_dtscale(mesh_path)
    dt = timestep(params, dtscale, int(params.get('P_ORDER', 1)))
    steps = math.ceil(params['TEND'] / dt)
    per_step = (seconds_per_step if seconds_per_step is not None
                else SECONDS_PER_STEP_PER_TRIANGLE * n_tri)
    v_max, c_se = max_wave_speed(params)
    return dt, steps, steps * per_step, {
        'triangles': n_tri, 'dtscale': dtscale, 'max_speed': v_max,
        'electron_sound_speed': c_se,
        'limited_by': 'light' if params['LIGHT'] >= c_se else 'electron sound speed',
        'seconds_per_step': per_step,
    }


def deck_number(value):
    """Format a number so the deck parser reads it as a REAL, not an integer.

    Apollo's parser types a literal with no decimal point or exponent as an
    integer, and get<REAL> on an integer throws std::bad_cast, which the solver
    reports as "unexpected error on MPI rank 0: std::bad_cast" with no mention
    of which key. The test/cxx tests carry a comment warning about this; the
    first programmatically generated deck here walked straight into it anyway,
    writing "c0 = 1000000" for LIGHT = 1.0e6 and aborting the run at setup. A
    scan generates every one of its decks, so it would have hit this at every
    point.
    """
    text = repr(float(value))
    if '.' not in text and 'e' not in text and 'n' not in text:  # nan/inf
        text += '.0'
    return text


def apply_speedup(text, s):
    """Rewrite a deck for the (T_e -> T_e/s^2, c -> c/s) scaling.

    Returns (new_text, note). Refuses rather than silently doing nothing if
    either assignment is not where it is expected, since a scan that thought it
    had scaled and had not would compare runs at different parameters.
    """
    new, n_te = re.subn(r'^(Te\s*=\s*)([0-9.eE+-]+)',
                        lambda m: f'{m.group(1)}{deck_number(float(m.group(2)) / (s * s))}',
                        text, count=1, flags=re.M)
    new, n_c = re.subn(r'^(LIGHT\s*=\s*)([0-9.eE+-]+)',
                       lambda m: f'{m.group(1)}{deck_number(float(m.group(2)) / s)}',
                       new, count=1, flags=re.M)
    if n_te != 1 or n_c != 1:
        raise SystemExit('--speedup: could not find both "Te =" and "LIGHT =" at the '
                         'start of a line in this deck; refusing to half-apply it')
    return new, (f'Te/{s * s:g} and LIGHT/{s:g}: dt x{s:g}, beta /{s * s:g}. '
                 f'gamma, lambda, omega/omega_ci, B_omega/B_bias, the Debye length '
                 f'in cells and omega_pe*dt are all unchanged.')


def set_param(text, name, value):
    new, n = re.subn(rf'^({name}\s*=\s*)([0-9.eE+-]+)',
                     lambda m: f'{m.group(1)}{deck_number(value)}', text, count=1,
                     flags=re.M)
    if n != 1:
        raise SystemExit(f'could not find "{name} =" at the start of a line in the deck')
    return new


def run_point(deck_text, deck_name, mesh_path, outdir, binary):
    """One scan point in its own directory. Returns the directory."""
    os.makedirs(outdir, exist_ok=True)
    pin = os.path.join(outdir, deck_name)
    with open(pin, 'w') as fh:
        fh.write(deck_text)
    shutil.copy(mesh_path, outdir)
    env = dict(os.environ, PYTHONPATH=os.path.join(REPO, 'scripts'))
    subprocess.run([sys.executable, os.path.join(REPO, 'scripts', 'wxinpparse.py'),
                    '-i', deck_name], cwd=outdir, env=env, check=True,
                   stdout=subprocess.DEVNULL)
    inp = deck_name.replace('.pin', '.inp')
    with open(os.path.join(outdir, 'solver.log'), 'w') as log:
        subprocess.run([binary, '-i', inp], cwd=outdir, stdout=log,
                       stderr=subprocess.STDOUT, check=True)
    return outdir


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('deck')
    ap.add_argument('--param', default='Bomega',
                    help='deck macro to scan (default Bomega)')
    ap.add_argument('--values', help='comma-separated values for --param')
    ap.add_argument('--periods', type=float, default=None,
                    help='run length in RMF periods; overrides the deck TEND')
    ap.add_argument('--speedup', type=float, default=None,
                    help='apply the (Te/s^2, c/s) scaling for an s-fold larger dt')
    ap.add_argument('--outdir', default='rmf_scan')
    ap.add_argument('--mesh', default=None, help='defaults to the deck Gridname')
    ap.add_argument('--binary', default=os.path.join(REPO, 'src', 'build-opt', 'apollo'))
    ap.add_argument('--sec-per-step', type=float, default=None,
                    help='override the measured cost constant')
    ap.add_argument('--dry-run', action='store_true',
                    help='print the cost of the scan and stop')
    args = ap.parse_args(argv)

    import rmf_diagnostics as rd

    with open(args.deck) as fh:
        base_text = fh.read()
    note = None
    if args.speedup:
        base_text, note = apply_speedup(base_text, args.speedup)

    mesh = args.mesh
    if mesh is None:
        m = re.search(r"Gridname\s*=\s*'([^']+)'", base_text)
        if not m:
            raise SystemExit('no Gridname in the deck; pass --mesh')
        mesh = os.path.join(os.path.dirname(os.path.abspath(args.deck)), m.group(1))

    values = ([float(v) for v in args.values.split(',')] if args.values
              else [None])

    print(f'deck   {args.deck}')
    print(f'mesh   {mesh}')
    if note:
        print(f'scaled {note}')
    print()

    plan, total = [], 0.0
    for v in values:
        text = base_text if v is None else set_param(base_text, args.param, v)
        params = rd.deck_parameters_from_text(text)
        if args.periods is not None:
            text = set_param(text, 'TEND', args.periods / params['omega'])
            params = rd.deck_parameters_from_text(text)
        dt, steps, secs, info = estimate(params, mesh, args.sec_per_step)
        total += secs
        plan.append((v, text, params, dt, steps, secs, info))

    first = plan[0][6]
    print(f'{first["triangles"]} triangles, smallest inscribed radius '
          f'{first["dtscale"] * 1e3:.4f} mm')
    print(f'timestep limited by the {first["limited_by"]}: '
          f'max wave speed {first["max_speed"]:.4e} m/s '
          f'(electron sound speed {first["electron_sound_speed"]:.4e})')
    print()
    print(f'{args.param:<14}{"TEND [s]":<13}{"periods":<10}{"dt [s]":<13}'
          f'{"steps":<11}{"wall clock":<12}')
    for v, _t, params, dt, steps, secs, _i in plan:
        periods = params['TEND'] * params['omega']
        print(f'{("(deck)" if v is None else f"{v:.4g}"):<14}'
              f'{params["TEND"]:<13.4e}{periods:<10.2f}{dt:<13.4e}'
              f'{steps:<11d}{secs / 3600.0:<12.2f}h')
    print(f'\ntotal {total / 3600.0:.2f} hours ({total / 86400.0:.2f} days) '
          f'for {len(plan)} point(s)')

    if args.dry_run:
        print('\n--dry-run: nothing was run.')
        return 0

    if total > 6 * 3600:
        print(f'\nRefusing to start automatically: this scan is estimated at '
              f'{total / 3600.0:.1f} hours and the solver cannot be restarted if it '
              f'is interrupted (docs/known-issues.md section 6). Run the points '
              f'individually, or shorten it with --periods, or make it affordable '
              f'with --speedup.')
        return 2

    os.makedirs(args.outdir, exist_ok=True)
    results = []
    for v, text, params, _dt, _steps, secs, _i in plan:
        name = 'deck' if v is None else f'{args.param}_{v:.6g}'
        point = os.path.join(args.outdir, name)
        done = os.path.join(point, 'diagnostics.json')
        if os.path.exists(done):
            print(f'{name}: already done, skipping')
            results.append(json.load(open(done)))
            continue
        print(f'{name}: running, estimated {secs / 60.0:.1f} min ...', flush=True)
        started = time.time()
        run_point(text, os.path.basename(args.deck), mesh, point, args.binary)
        elapsed = time.time() - started

        import vtu
        names = sorted((f for f in os.listdir(point) if f.endswith('.vtu')),
                       key=lambda p: int(p.rsplit('_', 1)[1].split('.')[0]))
        frames = [(int(n.rsplit('_', 1)[1].split('.')[0]),
                   vtu.read(os.path.join(point, n))) for n in names]
        rows = rd.analyse(frames, params)
        out = {'value': v, 'param': args.param, 'elapsed_s': elapsed,
               'estimated_s': secs, 'rows': rows,
               'TEND': params['TEND'], 'periods': params['TEND'] * params['omega']}
        with open(done, 'w') as fh:
            json.dump(out, fh, indent=2)
        results.append(out)
        print(f'{name}: done in {elapsed / 60.0:.1f} min '
              f'(estimated {secs / 60.0:.1f})')

    summary = os.path.join(args.outdir, 'scan.json')
    with open(summary, 'w') as fh:
        json.dump(results, fh, indent=2)
    print(f'\nwrote {summary}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
