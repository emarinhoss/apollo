#!/usr/bin/env python3
"""Verification tests: does Apollo solve the equations it claims to solve?

The isentropic vortex is the standard verification case for a compressible Euler
solver. It is an exact, smooth, time-dependent solution: the vortex translates at
the free-stream velocity without changing shape, so the difference between the
computed solution and the analytic one at t=1 is entirely discretization error.

That makes two things checkable that no smoke test can reach:

  * accuracy - the error is small, and it is small for *every* numerical flux,
    not just the default. A flux that computes the wrong momentum component
    still runs; it just gives a wrong answer, and only a comparison against a
    known solution notices. Both defects fixed in this branch (the Roe flux
    rotation writing nflux[1] twice, and abs() truncating an HLL wave speed)
    are of exactly that kind.

  * convergence order - refining the mesh reduces the error at the rate the
    method promises. A scheme can be accurate on one mesh by luck; the slope is
    what says the spatial operator is right.

These need a built solver, so they skip when there is not one. Build with
`cd src && scons build-opt` (or set APOLLO_BIN).
"""

import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(REPO, 'scripts')
VORTEX = os.path.join(REPO, 'examples', 'unstructuredDG', 'euler', 'isentropicVortex')

try:
    import numpy as np
except ImportError:
    np = None

try:
    import vtu
except ImportError:
    vtu = None


def apollo_binary():
    # Absolute, always: each case runs in its own scratch directory, so a
    # relative APOLLO_BIN would be resolved against that instead of the repo.
    if os.environ.get('APOLLO_BIN'):
        return os.path.abspath(os.environ['APOLLO_BIN'])
    for variant in ('build-opt', 'build-debug'):
        path = os.path.join(REPO, 'src', variant, 'apollo')
        if os.path.isfile(path) and os.access(path, os.X_OK):
            return path
    return None


# Parameters of the vortex, matching examples/.../isentropicVortex.pin.
XO, YO, BETA, GAMMA, UO, VO = 5.0, 0.0, 5.0, 1.4, 1.0, 0.0


def exact_solution(x, y, t):
    """(rho, rho*u, rho*v, E) of the analytic vortex at time t."""
    pi = math.pi
    xc, yc = XO + UO * t, YO + VO * t
    r2 = (x - xc) ** 2 + (y - yc) ** 2
    envelope = np.exp(1.0 - r2)
    u = UO - BETA * envelope * (y - yc) / (2 * pi)
    v = VO + BETA * envelope * (x - xc) / (2 * pi)
    base = 1.0 - (GAMMA - 1.0) * BETA ** 2 * np.exp(2.0 * (1.0 - r2)) / (16 * GAMMA * pi ** 2)
    rho = base ** (1.0 / (GAMMA - 1.0))
    p = rho ** GAMMA
    energy = p / (GAMMA - 1.0) + 0.5 * rho * (u * u + v * v)
    return rho, rho * u, rho * v, energy


def build_deck(mesh, flux=None, tend=None):
    """Take the shipped vortex deck and vary only what the test needs.

    Reproducing the deck from scratch would test a configuration nobody runs -
    and get it subtly wrong, since the boundary condition and solver-sequence
    blocks carry keys that are easy to omit. Editing the real one keeps the test
    honest about what users actually run.
    """
    with open(os.path.join(VORTEX, 'isentropicVortex.pin')) as fh:
        deck = fh.read()

    deck = re.sub(r"Gridname\s*=\s*'[^']*'", "Gridname = '%s'" % mesh, deck, count=1)
    if tend is not None:
        deck = re.sub(r'Time\s*=\s*\[[^\]]*\]', 'Time = [0.0, %r]' % tend, deck, count=1)
    if flux:
        deck = deck.replace('Kind = eulerEqn',
                            'Kind = eulerEqn\n\t\tNumerical_Flux = %s' % flux, 1)
    # Quieter logs make failures easier to read, and the run marginally faster.
    deck = re.sub(r'Verbosity\s*=\s*\w+', 'Verbosity = warning', deck, count=1)
    return deck


class VortexRunner:
    """Runs one vortex case and reports its nodal L2 error."""

    def __init__(self, workdir):
        self.workdir = workdir

    def run(self, mesh, flux=None, tend=1.0):
        case = os.path.join(self.workdir, f'{os.path.splitext(mesh)[0]}_{flux or "default"}')
        os.makedirs(case, exist_ok=True)
        shutil.copy(os.path.join(VORTEX, mesh), case)

        pin = os.path.join(case, 'vortex.pin')
        with open(pin, 'w') as fh:
            fh.write(build_deck(mesh, flux=flux, tend=tend))

        env = dict(os.environ, PYTHONPATH=SCRIPTS)
        pre = subprocess.run(
            [sys.executable, os.path.join(SCRIPTS, 'wxinpparse.py'), '-i', 'vortex.pin'],
            cwd=case, env=env, capture_output=True, text=True, timeout=300)
        if pre.returncode != 0:
            raise RuntimeError(f'preprocessing failed: {pre.stdout}{pre.stderr}')

        run = subprocess.run([apollo_binary(), '-i', 'vortex.inp'],
                             cwd=case, capture_output=True, text=True, timeout=1800)
        log = run.stdout + run.stderr
        if run.returncode != 0:
            raise RuntimeError(f'solver exited {run.returncode}:\n{log[-2000:]}')
        for bad in ('NaN', 'Negative pressure', 'Negative density'):
            if bad in log:
                raise RuntimeError(f'solver reported "{bad}":\n{log[-2000:]}')

        outputs = sorted(f for f in os.listdir(case) if f.endswith('.vtu'))
        if not outputs:
            raise RuntimeError(f'no .vtu output:\n{log[-2000:]}')
        return self._l2_error(os.path.join(case, outputs[-1]), tend)

    @staticmethod
    def _l2_error(path, t):
        """RMS error over every nodal degree of freedom, per conserved variable.

        For P1 on triangles the three nodes are the cell vertices, and the 15
        cell-data arrays are laid out node-major: solutiondg.(node*5 + component).
        """
        data = vtu.read(path)
        points = data['Position']
        cells = data['connectivity'].reshape(-1, 3)

        errors = {}
        for component, name in ((0, 'rho'), (1, 'rhou'), (2, 'rhov'), (4, 'E')):
            computed, analytic = [], []
            for node in range(3):
                vid = cells[:, node]
                x, y = points[vid, 0], points[vid, 1]
                computed.append(data[f'solutiondg.{node * 5 + component}'])
                analytic.append(exact_solution(x, y, t)[{0: 0, 1: 1, 2: 2, 4: 3}[component]])
            diff = np.concatenate(computed) - np.concatenate(analytic)
            errors[name] = float(np.sqrt(np.mean(diff ** 2)))
        return errors


@unittest.skipIf(np is None, 'numpy is required')
@unittest.skipIf(vtu is None, 'test/vtu.py not importable')
@unittest.skipIf(apollo_binary() is None,
                 'no Apollo binary; build with "cd src && scons build-opt" '
                 'or set APOLLO_BIN')
class TestVortexAccuracy(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix='apollo-vortex-')
        cls.runner = VortexRunner(cls.tmp)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_every_flux_solves_the_vortex(self):
        """Each numerical flux must reproduce the analytic solution.

        A flux that writes the wrong momentum component runs happily and returns
        a wrong answer; the tolerance here is loose enough to admit any
        reasonable scheme on a 1110-cell mesh and tight enough to catch that.

        Wave is not in this list: applyWavePropagationFluxes passes its states
        to the Riemann solver in an order inconsistent with the flux evaluation
        above it and still aborts on this case. Add it here when that is fixed.
        """
        tolerance = 5.0e-2
        for flux in (None, 'LF', 'HLL', 'Roe'):
            with self.subTest(flux=flux or 'default'):
                errors = self.runner.run('grid2.msh', flux=flux)
                for name, err in errors.items():
                    self.assertLess(
                        err, tolerance,
                        f'{flux or "default"} flux: L2 error in {name} is {err:.3e}, '
                        f'above {tolerance:.0e}')
                print(f'\n    {flux or "default":8s} ' +
                      '  '.join(f'{k}={v:.3e}' for k, v in errors.items()))

    def test_refining_the_mesh_reduces_the_error(self):
        """The spatial operator must actually converge.

        290 -> 1110 -> 4454 triangles is close enough to uniform refinement for
        the observed order to be meaningful. P1 DG is formally 2nd order; the
        bar here is 1.5, low enough not to be brittle about mesh quality and
        boundary treatment, high enough that a broken operator cannot pass.
        """
        meshes = [('grid.msh', 290), ('grid2.msh', 1110), ('grid3.msh', 4454)]
        results = []
        for mesh, ncells in meshes:
            errors = self.runner.run(mesh, tend=0.5)
            results.append((ncells, errors['rho']))
            print(f'\n    {mesh:12s} {ncells:5d} cells   L2(rho) = {errors["rho"]:.4e}')

        for (n0, e0), (n1, e1) in zip(results, results[1:]):
            self.assertLess(e1, e0, f'error grew from {n0} to {n1} cells')
            # h scales as 1/sqrt(ncells) in 2D.
            order = math.log(e0 / e1) / math.log(math.sqrt(n1 / n0))
            print(f'    {n0} -> {n1} cells: observed order {order:.2f}')
            self.assertGreater(
                order, 1.5,
                f'observed convergence order {order:.2f} between {n0} and {n1} '
                f'cells is below 1.5; the spatial operator is not converging '
                f'at the rate P1 DG should')


if __name__ == '__main__':
    unittest.main()
