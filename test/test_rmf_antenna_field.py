#!/usr/bin/env python3
"""Verification: does the RMF antenna produce the field theory says it does?

The rotating-magnetic-field drive in examples/.../rmf_frc/antenna is an antenna
current, not a prescribed field, which is the point of it - the plasma can
screen it and load it. The price is that nothing about the field is guaranteed
by construction any more, so it has to be checked.

With the plasma removed the check is exact. A z-directed surface current
k cos(theta - psi) on a shell of radius s, inside a perfectly conducting
cylinder at r = b, produces a uniform transverse field mu0 k (1 - s^2/b^2)/2
everywhere inside itself. Superposing the shells of the winding, the field
inside the antenna is

    |B| = B_rmf (1 - exp(-t/rise))     uniform, rotating at omega

and vacuum.pin is that configuration: same antenna, same wall, no plasma and no
source that couples the fluids to the field. So three things are checkable
against closed form rather than against a previous run:

  * the magnitude, which says setup() converted B_rmf into the right current,
    including the conducting wall's image-current factor;
  * the uniformity, which says the winding really is m = 1 - a monopole or a
    quadrupole component would show up as radial or azimuthal structure;
  * the rotation rate, which says the drive rotates rather than oscillating.

The tolerances are set by the O((omega b / c)^2) = 1% quasi-static correction
(the deck runs c at c/100, so retardation across the domain is not entirely
negligible) and by the coarse mesh, not by anything adjustable.

Needs a built solver, so it skips when there is not one. Build with
`cd src && scons build-opt` (or set APOLLO_BIN).
"""

import math
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(REPO, 'scripts')
ANTENNA = os.path.join(REPO, 'examples', 'unstructuredDG', 'multifluid',
                       'rmf_frc', 'antenna')

try:
    import numpy as np
except ImportError:
    np = None

try:
    import vtu
except ImportError:
    vtu = None


def apollo_binary():
    if os.environ.get('APOLLO_BIN'):
        return os.path.abspath(os.environ['APOLLO_BIN'])
    for variant in ('build-opt', 'build-debug'):
        path = os.path.join(REPO, 'src', variant, 'apollo')
        if os.path.isfile(path) and os.access(path, os.X_OK):
            return path
    return None


# Must match vacuum.pin. Duplicated rather than parsed: a test that reads its
# expectations out of the file under test cannot fail when that file is wrong.
B_RMF, FREQ, RISE = 50.0e-4, 7.95e5, 1.0e-7
COIL_R, COIL_W = 0.036, 0.005
TEND, NOUT = 3.0e-7, 12
OMEGA = 2.0 * math.pi * FREQ

# The eighteen-component two-fluid state, laid out node-major in the .vtu cell
# data: solutiondg.(node*18 + component). 13 and 14 are B_x and B_y.
NEQ, BX, BY = 18, 13, 14


def run_vacuum(workdir):
    """Run vacuum.pin and return [(t, x, y, Bx, By)] for every output frame."""
    case = os.path.join(workdir, 'vacuum')
    os.makedirs(case, exist_ok=True)
    for name in ('vacuum.pin', 'vacuumDisc.msh'):
        shutil.copy(os.path.join(ANTENNA, name), case)

    env = dict(os.environ, PYTHONPATH=SCRIPTS)
    pre = subprocess.run(
        [sys.executable, os.path.join(SCRIPTS, 'wxinpparse.py'), '-i', 'vacuum.pin'],
        cwd=case, env=env, capture_output=True, text=True, timeout=300)
    if pre.returncode != 0:
        raise RuntimeError(f'preprocessing failed: {pre.stdout}{pre.stderr}')

    run = subprocess.run([apollo_binary(), '-i', 'vacuum.inp'],
                         cwd=case, capture_output=True, text=True, timeout=3600)
    log = run.stdout + run.stderr
    if run.returncode != 0:
        raise RuntimeError(f'solver exited {run.returncode}:\n{log[-2000:]}')
    for bad in ('NaN', 'Negative pressure', 'Negative density'):
        if bad in log:
            raise RuntimeError(f'solver reported "{bad}":\n{log[-2000:]}')

    # Sorted by frame INDEX, not by name: a plain sort gives _0, _1, _10, _11,
    # _12, _2, ... which silently shuffles time order, and everything below
    # that takes "the later frames" or measures a sweep between consecutive
    # frames then reads a scrambled sequence.
    def index_of(name):
        return int(name.rsplit('_', 1)[1].split('.')[0])

    frames = []
    for name in sorted((f for f in os.listdir(case) if f.endswith('.vtu')),
                       key=index_of):
        index = index_of(name)
        data = vtu.read(os.path.join(case, name))
        points, cells = data['Position'], data['connectivity'].reshape(-1, 3)
        xs, ys, bxs, bys = [], [], [], []
        for node in range(3):
            vid = cells[:, node]
            xs.append(points[vid, 0])
            ys.append(points[vid, 1])
            bxs.append(data[f'solutiondg.{node * NEQ + BX}'])
            bys.append(data[f'solutiondg.{node * NEQ + BY}'])
        frames.append((TEND * index / NOUT,
                       np.concatenate(xs), np.concatenate(ys),
                       np.concatenate(bxs), np.concatenate(bys)))
    if len(frames) < 4:
        raise RuntimeError(f'expected output frames, got {len(frames)}')
    return frames


@unittest.skipIf(np is None, 'numpy is required')
@unittest.skipIf(vtu is None, 'test/vtu.py not importable')
@unittest.skipIf(apollo_binary() is None,
                 'no Apollo binary; build with "cd src && scons build-opt" '
                 'or set APOLLO_BIN')
class TestRMFAntennaField(unittest.TestCase):
    """One run, three questions asked of it."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix='apollo-rmf-antenna-')
        cls.frames = run_vacuum(cls.tmp)
        # Inside the winding's inner edge, where the closed form applies, and
        # away from the axis where the polar mesh's smallest cells are.
        r_in = COIL_R - 0.5 * COIL_W
        cls.interior = [
            (t, bx, by, (np.hypot(x, y) > 0.004) & (np.hypot(x, y) < r_in))
            for (t, x, y, bx, by) in cls.frames]

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    # Frames before this are still ringing; see settled().
    SETTLE = 1.5

    @classmethod
    def settled(cls):
        """Frames at t >= 1.5*RISE.

        The closed form is the steady rotating solution. Switching the antenna
        on also rings the conducting cylinder at its own modes - about 44 ns
        here - and nothing damps that, there being no plasma. Measured on the
        shipped deck, the error in |B| against the closed form runs:

            t/RISE   0.25   0.5    0.75   1.0    1.25   1.5    1.75   2.0
            |B| err  9.8%   6.5%   6.2%   4.9%   1.6%   0.2%   0.9%   1.9%
            spread   3.8%   9.6%   6.0%   1.8%   3.1%   0.7%   3.1%   1.8%

        so everything from 1.5 rise times on sits well inside the tolerances
        below, and everything before it does not. This is a property of the
        physics rather than a frame count, so it is written as one - shortening
        or lengthening the deck moves the window with it.
        """
        settled = [f for f in cls.interior if f[0] >= cls.SETTLE * RISE]
        assert len(settled) >= 4, (
            f'only {len(settled)} frames at t >= {cls.SETTLE}*RISE; the deck '
            f'needs a longer TEND or more outputs for this test to mean anything')
        return settled

    def test_magnitude_matches_the_closed_form(self):
        """|B| inside the winding must be B_rmf times the switch-on envelope."""
        checked = 0
        for t, bx, by, mask in self.settled():
            expected = B_RMF * (1.0 - math.exp(-t / RISE))
            got = float(np.mean(np.hypot(bx[mask], by[mask])))
            self.assertLess(
                abs(got - expected) / expected, 0.05,
                f't = {t:.3e}: |B| = {got:.6e}, closed form {expected:.6e}')
            checked += 1
        self.assertGreater(checked, 3, 'too few frames were actually checked')

    def test_field_is_uniform(self):
        """The winding is m = 1, so the field it encloses has no structure."""
        for t, bx, by, mask in self.settled():
            mag = np.hypot(bx[mask], by[mask])
            spread = float(np.std(mag) / np.mean(mag))
            self.assertLess(
                spread, 0.06,
                f't = {t:.3e}: |B| varies by {spread:.3f} of its mean across '
                f'the interior; the closed-form field is uniform')

    def test_field_rotates_at_the_drive_frequency(self):
        """Rotating, not oscillating, and at omega rather than some fraction."""
        late = self.settled()
        angles, times = [], []
        for t, bx, by, mask in late:
            angles.append(math.atan2(float(np.mean(by[mask])),
                                     float(np.mean(bx[mask]))))
            times.append(t)

        swept = 0.0
        for a, b in zip(angles, angles[1:]):
            d = b - a
            while d > math.pi:
                d -= 2 * math.pi
            while d < -math.pi:
                d += 2 * math.pi
            swept += d
        expected = OMEGA * (times[-1] - times[0])
        self.assertLess(
            abs(abs(swept) - expected) / expected, 0.10,
            f'field swept {abs(swept):.4f} rad over {times[-1] - times[0]:.3e} s; '
            f'omega*dt = {expected:.4f} rad')
        # An oscillating field would sweep nothing on average. Expressed
        # relative to the expected sweep so that shortening the deck cannot
        # quietly turn this into a tautology, or into a guaranteed failure.
        self.assertGreater(abs(swept), 0.5 * expected,
                           f'field swept only {abs(swept):.4f} rad against an '
                           f'expected {expected:.4f}; this is not a rotating drive')


if __name__ == '__main__':
    unittest.main()
