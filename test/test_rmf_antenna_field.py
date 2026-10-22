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
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ANTENNA = os.path.join(REPO, 'examples', 'unstructuredDG', 'multifluid',
                       'rmf_frc', 'antenna')

try:
    import numpy as np
except ImportError:
    np = None

try:
    import deckrun
except ImportError:
    deckrun = None

# Must match vacuum.pin. Duplicated rather than parsed: a test that reads its
# expectations out of the file under test cannot fail when that file is wrong.
B_RMF, FREQ, RISE = 50.0e-4, 7.95e5, 1.0e-7
COIL_R, COIL_W = 0.036, 0.005
TEND, NOUT = 3.0e-7, 12
OMEGA = 2.0 * math.pi * FREQ

# The eighteen-component two-fluid state; 13 and 14 are B_x and B_y.
BX, BY = 13, 14


@unittest.skipIf(np is None, 'numpy is required')
@unittest.skipIf(deckrun is None or deckrun.vtu is None, 'test helpers not importable')
@unittest.skipIf(deckrun is not None and deckrun.apollo_binary() is None,
                 'no Apollo binary; build with "cd src && scons build-opt" '
                 'or set APOLLO_BIN')
class TestRMFAntennaField(unittest.TestCase):
    """One run, three questions asked of it."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix='apollo-rmf-antenna-')
        case = deckrun.run(ANTENNA, 'vacuum.pin', cls.tmp,
                           aux_files=('vacuumDisc.msh',))
        # Inside the winding's inner edge, where the closed form applies, and
        # away from the axis where the polar mesh's smallest cells are.
        r_in = COIL_R - 0.5 * COIL_W
        cls.interior = []
        for index, data in deckrun.frames(case):
            x, y, f = deckrun.nodal(data, (BX, BY))
            r = np.hypot(x, y)
            cls.interior.append((TEND * index / NOUT, f[BX], f[BY],
                                 (r > 0.004) & (r < r_in)))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    # Frames before this are still ringing; see settled().
    SETTLE = 1.5

    @classmethod
    def settled(cls):
        """Frames at t >= 1.5*RISE.

        The closed form is the steady rotating solution. Switching the antenna
        on also rings the conducting cylinder at its own modes - the lowest one
        the m = 1 drive can excite is j_11 = 3.8317, a 27 ns period; see the
        note in vacuum.pin - and nothing damps that, there being no plasma.
        Measured on the
        shipped deck, the error in |B| against the closed form runs:

            t/RISE  0.25  0.5   0.75  1.0   1.25  1.5   1.75  2.0   2.5   3.0
            |B| err  9.8   6.5   6.2   4.9   1.6   0.2   0.9   1.9   0.4   0.3  %
            spread   3.8   9.6   6.0   1.8   3.1   0.7   3.1   1.8   3.0   0.9  %

        so everything from 1.5 rise times on sits well inside the tolerances
        below, and everything before it does not. Over that window the field
        also turns 0.7491 rad against omega*dt = 0.7493, a rotation-rate error
        of 0.03%. This is a property of the physics rather than a frame count,
        so it is written as one - shortening or lengthening the deck moves the
        window with it.
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
        # SIGNED, not |swept|. Which way the field turns is the physically
        # decisive property: electrons are dragged in the sense of the rotation,
        # so J_theta = -e n u_theta takes the opposite sign to it, and that
        # decides whether the driven axial field opposes the bias (formation) or
        # reinforces it. The antenna's cos(theta - omega t - phase) turns
        # counter-clockwise, so the sweep is positive; twoFluidSimplifiedRMFBC
        # turns the other way. Comparing magnitudes would let a flip through.
        expected = OMEGA * (times[-1] - times[0])
        self.assertLess(
            abs(swept - expected) / expected, 0.10,
            f'field swept {swept:+.4f} rad over {times[-1] - times[0]:.3e} s; '
            f'omega*dt = {expected:+.4f} rad (sign included: the antenna turns '
            f'counter-clockwise)')
        # An oscillating field would sweep nothing on average. Expressed
        # relative to the expected sweep so that shortening the deck cannot
        # quietly turn this into a tautology, or into a guaranteed failure.
        self.assertGreater(swept, 0.5 * expected,
                           f'field swept only {swept:+.4f} rad against an '
                           f'expected {expected:+.4f}; this is not a rotating '
                           f'drive, or it is turning the wrong way')


if __name__ == '__main__':
    unittest.main()
