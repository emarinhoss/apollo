#!/usr/bin/env python3
"""Verification: does the edge-driven RMF boundary condition apply its field?

examples/unstructuredDG/multifluid/rmf_frc/frc2d.pin prescribes a rotating
transverse field at the plasma edge with twoFluidSimplifiedRMFBC. Until now
nothing checked, at solver level, that the field inside the domain is the field
the deck asked for - and two defects in that class went unnoticed until someone
read it: the phase was applied to one Cartesian component only, which turns a
rotating drive into a linearly polarised one, and the induced E_z carried the
right radial amplitude with the azimuthal dependence dropped, satisfying
neither component of Faraday's law.

Both fail this test - but only because it was designed against them after the
claim that it caught them was made and found to be false. Reintroducing each
defect into the boundary condition and running the shipped deck left all three
original checks passing. Two things fixed that, and neither is incidental
(the numbers below are from rebuilding each defect and re-running, which is the
step the original claim skipped):

  * vacuum.pin sets PHASE = PI/2. At PHASE = 0 the defective and correct phase
    expressions are the same expression - the difference is adding +0.0 - so no
    run of that deck could distinguish them even in principle. At PI/2 the
    defect fails the rotation check (0.0000 rad swept against -0.4496 expected:
    a linearly polarised field does not turn) and the magnitude check (38.5%
    against 5%).
  * test_induced_ez_varies_around_the_azimuth checks the angular STRUCTURE of
    E_z. The E_z defect moves |B| by only 0.24%, against a 5% tolerance and a
    0.25% quasi-static error, so no magnitude tolerance can separate it; on
    angular structure it gives 0.013 against a floor of 0.30, where the correct
    field gives 1.115.

test/cxx/test_rmf_boundary.cc also catches both, at unit level and in far less
time; this test is what checks that the field actually reaches the interior of a
run.

With the plasma removed and every fluid-field coupling source gone, the answer
is known: the interior carries the same uniform transverse field the boundary
does, of magnitude B_rmf(1 - exp(-t/rise)), rotating at omega. That is the
quasi-static limit rather than an exact solution - a spatially uniform B_perp
has curl B = 0, so Ampere would need dE/dt = 0 while the induced E_z rotates -
and it is good to O((omega a/c)^2), about 0.25% at this deck's numbers. The
tolerances below are that plus discretisation, not fitted.

The rotation check is SIGNED. This boundary condition writes
B = (-B_t sin(omega t + phase), -B_t cos(omega t + phase)), whose direction
angle decreases: it turns clockwise, the opposite way from the antenna source in
rmf_frc/antenna. Which way an RMF turns decides whether the driven azimuthal
current opposes the bias field or reinforces it, so it is pinned here.

There is room to spare against all four. Measured over the seven frames of the
settled window:

    magnitude    0.19% worst    against 5%
    uniformity   0.17% worst    against 6%
    rotation     0.25%          against 10%
    E_z ratio    1.115 lowest   against a floor of 0.30 (the defect gives 0.013)

The frames before the window are what the window exists for: at t = 15 ns the
magnitude error is 17% and the spread 11%, which is the cavity ringing the
switch-on excites.

Runs in about 4.6 minutes. It is in CI on every push, which is why the deck is
no longer and no finer than the margins above allow; see vacuum.pin.

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
RMF = os.path.join(REPO, 'examples', 'unstructuredDG', 'multifluid', 'rmf_frc')

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
B_RMF, FREQ, RISE = 50.0e-4, 7.95e5, 6.0e-8
RAD = 0.030
TEND, NOUT = 1.8e-7, 12
OMEGA = 2.0 * math.pi * FREQ

EZ, BX, BY = 12, 13, 14   # E_z, B_x, B_y in the 18-component state
SETTLE = 1.5             # rise times; see settled()


@unittest.skipIf(np is None, 'numpy is required')
@unittest.skipIf(deckrun is None or deckrun.vtu is None, 'test helpers not importable')
@unittest.skipIf(deckrun is not None and deckrun.apollo_binary() is None,
                 'no Apollo binary; build with "cd src && scons build-opt" '
                 'or set APOLLO_BIN')
class TestRMFBoundaryField(unittest.TestCase):
    """One plasma-free run, three questions asked of it."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix='apollo-rmf-bc-')
        case = deckrun.run(RMF, 'vacuum.pin', cls.tmp, aux_files=('vacuumDisc.msh',))
        cls.samples = []
        for index, data in deckrun.frames(case):
            x, y, f = deckrun.nodal(data, (EZ, BX, BY))
            r = np.hypot(x, y)
            # Away from the axis, where the polar mesh's smallest cells are,
            # and away from the boundary layer itself.
            mask = (r > 0.004) & (r < 0.8 * RAD)
            cls.samples.append((TEND * index / NOUT, f[BX], f[BY], mask,
                                f[EZ], x, y, r))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    @classmethod
    def settled(cls):
        """Frames at t >= SETTLE*RISE.

        Switching the boundary field on rings the domain at its own modes - the
        lowest an m = 1 drive can excite in a circle of radius a is
        j_11 c/a = 3.83e8 rad/s, a 16 ns period - and with no plasma nothing
        damps it. RISE is 3.7 of those periods, the same ratio the antenna deck
        was validated at, and this window starts well after the ramp, which is
        what puts the ringing below the tolerances.
        """
        out = [s for s in cls.samples if s[0] >= SETTLE * RISE]
        assert len(out) >= 4, (
            f'only {len(out)} frames at t >= {SETTLE}*RISE; vacuum.pin needs a '
            f'longer TEND or more outputs for this test to mean anything')
        return out

    def test_magnitude_matches_the_applied_field(self):
        """|B_perp| inside must be B_rmf times the switch-on envelope."""
        for t, bx, by, mask, _ez, _x, _y, _r in self.settled():
            expected = B_RMF * (1.0 - math.exp(-t / RISE))
            got = float(np.mean(np.hypot(bx[mask], by[mask])))
            self.assertLess(
                abs(got - expected) / expected, 0.05,
                f't = {t:.3e}: |B| = {got:.6e}, applied field {expected:.6e}')

    def test_field_is_uniform(self):
        """The boundary applies a spatially uniform field, so the interior is."""
        for t, bx, by, mask, _ez, _x, _y, _r in self.settled():
            mag = np.hypot(bx[mask], by[mask])
            spread = float(np.std(mag) / np.mean(mag))
            self.assertLess(
                spread, 0.06,
                f't = {t:.3e}: |B| varies by {spread:.3f} of its mean across the '
                f'interior; the applied field is uniform')

    def test_field_rotates_clockwise_at_the_drive_frequency(self):
        """Rotating, at omega, and in the sense this boundary condition writes.

        Signed: twoFluidSimplifiedRMFBC turns clockwise, so the sweep is
        negative. Comparing magnitudes would let a handedness flip through, and
        handedness is what decides whether the drive reverses the bias field or
        reinforces it.
        """
        late = self.settled()
        times = [s[0] for s in late]
        angles = [math.atan2(float(np.mean(s[2][s[3]])), float(np.mean(s[1][s[3]])))
                  for s in late]

        swept = 0.0
        for a, b in zip(angles, angles[1:]):
            d = b - a
            while d > math.pi:
                d -= 2 * math.pi
            while d < -math.pi:
                d += 2 * math.pi
            swept += d

        expected = -OMEGA * (times[-1] - times[0])
        self.assertLess(
            abs(swept - expected) / abs(expected), 0.10,
            f'field swept {swept:+.4f} rad over {times[-1] - times[0]:.3e} s; '
            f'-omega*dt = {expected:+.4f} rad (sign included: this boundary '
            f'condition turns clockwise)')
        self.assertLess(swept, 0.5 * expected,
                        f'field swept only {swept:+.4f} rad against an expected '
                        f'{expected:+.4f}; this is not a rotating drive, or it '
                        f'is turning the wrong way')

    def test_induced_ez_varies_around_the_azimuth(self):
        """E_z must be m = 1, not a function of radius alone.

        This is the check that guards the second Phase 0 defect. The induced
        axial field is E_z = x dB_y/dt - y dB_x/dt, which on a circle of fixed
        radius goes as cos(theta - something): it changes sign twice around.
        The defect carried the right radial amplitude with the azimuthal
        dependence dropped, giving an E_z that is CONSTANT on every such circle.

        It is deliberately not a tolerance on a magnitude. Reintroducing that
        defect and running this deck moves |B| by 0.24% against a 5% tolerance -
        no defensible tolerance separates that from the O((omega a/c)^2) = 0.25%
        quasi-static error the deck itself cites. The angular structure does
        separate them, with nothing to tune: one varies and the other does not.
        """
        for t, _bx, _by, _mask, ez, x, y, r in self.settled():
            band = (r > 0.4 * RAD) & (r < 0.6 * RAD)
            self.assertGreater(int(np.count_nonzero(band)), 50,
                               'too few nodes in the sampling annulus')
            vals = ez[band]
            # Spread around the circle against the typical magnitude. For
            # E_z ~ cos(theta - phi) sampled uniformly in theta, std is
            # 1/sqrt(2) of the amplitude and mean|.| is 2/pi of it, so the
            # ratio is pi/(2 sqrt(2)) = 1.111; the measured value is 1.115,
            # the excess being the r-dependence across the band. For an E_z
            # that depends only on r it is that r-spread alone - the defective
            # form gives 0.013. The threshold sits between them, 3.7x below one
            # and 23x above the other, so it is not a fitted number.
            ratio = float(np.std(vals) / np.mean(np.abs(vals)))
            self.assertGreater(
                ratio, 0.3,
                f't = {t:.3e}: E_z varies by only {ratio:.3f} of its own '
                f'magnitude around the annulus; the induced field is m = 1 and '
                f'should give about 1.11. An E_z that is a function of radius '
                f'alone would give nearly zero here.')


if __name__ == '__main__':
    unittest.main()
