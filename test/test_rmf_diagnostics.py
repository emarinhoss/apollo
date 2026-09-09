#!/usr/bin/env python3
"""Verification: do the Phase 3 RMF diagnostics measure what they claim to?

scripts/rmf_diagnostics.py extracts the numbers Phase 3 of
docs/rmf-frc-model-assessment.md compares against the RMF current-drive
literature - the driven axial field, the electron rotation parameter zeta, the
penetration of the transverse field, and the thickness of the current layer.
Those numbers are the evidence the whole phase would rest on, so the instrument
is checked here before any of them is believed.

Every check below feeds an ANALYTIC field whose answer is known in closed form
and asserts the diagnostic returns it. Nothing here runs the solver: the point
is to separate "the diagnostic is wrong" from "the physics is surprising", and
the only way to do that is to hand it a case where the physics is not in doubt.
It runs in under a second, so it is in CI on every push.

The most important check is test_penetrated_limit_sign_opposes_a_plus_z_bias.
Phase 3 item 0 - which way does the drive turn, and does the driven field oppose
the bias - is entirely a question of sign, and a sign convention that is only
written in a comment is one refactor from being wrong. Here it is executable.
"""

import math
import os
import sys
import unittest

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'scripts'))

try:
    import numpy as np
except ImportError:
    np = None

try:
    import rmf_diagnostics as rd
except ImportError:
    rd = None


MU0 = 4.0e-7 * math.pi
Q = 1.6e-19
ME = 1.67e-27 / 1836.0
MI = 1.67e-27
A = 0.030                       # plasma radius, as in the deck
OMEGA = 2.0 * math.pi * 7.95e5


def disc(n_r=60, n_t=180, a=A):
    """A structured polar cloud of nodes, avoiding r = 0 exactly.

    Deliberately structured rather than a real mesh: a diagnostic that only
    works on one mesh is not a diagnostic, and a bin-average has a clean
    expected value on a cloud whose radial distribution is known.
    """
    r = np.linspace(a / (2 * n_r), a, n_r)
    t = np.linspace(0.0, 2 * math.pi, n_t, endpoint=False)
    rr, tt = np.meshgrid(r, t, indexing='ij')
    return (rr.ravel() * np.cos(tt.ravel()),
            rr.ravel() * np.sin(tt.ravel()))


@unittest.skipIf(np is None, 'numpy is required')
@unittest.skipIf(rd is None, 'scripts/rmf_diagnostics.py did not import')
class TestGeometry(unittest.TestCase):

    def test_polar_returns_the_unit_vectors(self):
        x = np.array([1.0, 0.0, -2.0, 0.0])
        y = np.array([0.0, 3.0, 0.0, -1.0])
        r, c, s = rd.polar(x, y)
        np.testing.assert_allclose(r, [1.0, 3.0, 2.0, 1.0])
        np.testing.assert_allclose(c, [1.0, 0.0, -1.0, 0.0], atol=1e-15)
        np.testing.assert_allclose(s, [0.0, 1.0, 0.0, -1.0], atol=1e-15)

    def test_polar_does_not_divide_by_zero_on_the_axis(self):
        """A node exactly on the axis must not produce NaN or an exception."""
        r, c, s = rd.polar(np.array([0.0]), np.array([0.0]))
        self.assertEqual(float(r[0]), 0.0)
        self.assertTrue(math.isfinite(float(c[0])) and math.isfinite(float(s[0])))

    def test_azimuthal_and_radial_decompose_a_known_field(self):
        """A rigid rotation is purely azimuthal; a source flow purely radial."""
        x, y = disc(n_r=8, n_t=32)
        r, c, s = rd.polar(x, y)
        # u = omega_0 * z_hat cross r  ->  purely azimuthal, magnitude omega_0 r
        w0 = 3.0
        ux, uy = -w0 * y, w0 * x
        np.testing.assert_allclose(rd.azimuthal(ux, uy, c, s), w0 * r, rtol=1e-12)
        np.testing.assert_allclose(rd.radial(ux, uy, c, s), 0.0, atol=1e-12)
        # u = k * r_hat  ->  purely radial
        k = 2.0
        np.testing.assert_allclose(rd.radial(k * x, k * y, c, s), k * r, rtol=1e-12)
        np.testing.assert_allclose(rd.azimuthal(k * x, k * y, c, s), 0.0, atol=1e-12)


@unittest.skipIf(np is None, 'numpy is required')
@unittest.skipIf(rd is None, 'scripts/rmf_diagnostics.py did not import')
class TestCurrent(unittest.TestCase):

    def test_current_from_momenta_matches_q_n_u(self):
        """J = sum_s q_s n_s u_s, with the electron term negative."""
        n = 1e20
        ue, ui = 4.0e5, 1.0e3
        comps = {rd.MOM_EX: np.array([n * ME * ue]), rd.MOM_EY: np.array([0.0]),
                 rd.MOM_EZ: np.array([0.0]),
                 rd.MOM_IX: np.array([n * MI * ui]), rd.MOM_IY: np.array([0.0]),
                 rd.MOM_IZ: np.array([0.0])}
        jx, jy, jz = rd.current_density(comps, ME, MI, Q)
        self.assertAlmostEqual(float(jx[0]), Q * n * (ui - ue), delta=1e-9 * Q * n * ue)
        self.assertEqual(float(jy[0]), 0.0)
        self.assertEqual(float(jz[0]), 0.0)

    def test_electrons_moving_counter_clockwise_give_a_clockwise_current(self):
        """The sign that decides field reversal. See the module's conventions."""
        x, y = disc(n_r=6, n_t=24)
        r, c, s = rd.polar(x, y)
        n = 1e20
        # Electrons rigidly rotating counter-clockwise: u_theta = +omega r.
        uex, uey = -OMEGA * y, OMEGA * x
        comps = {rd.MOM_EX: n * ME * uex, rd.MOM_EY: n * ME * uey,
                 rd.MOM_EZ: np.zeros_like(x),
                 rd.MOM_IX: np.zeros_like(x), rd.MOM_IY: np.zeros_like(x),
                 rd.MOM_IZ: np.zeros_like(x)}
        jx, jy, _ = rd.current_density(comps, ME, MI, Q)
        j_theta = rd.azimuthal(jx, jy, c, s)
        u_theta = rd.azimuthal(uex, uey, c, s)
        self.assertTrue(np.all(u_theta > 0), 'electrons should be counter-clockwise')
        self.assertTrue(np.all(j_theta < 0),
                        'counter-clockwise electrons must carry a CLOCKWISE current; '
                        'if this fails, the charge sign or e_theta is wrong')
        np.testing.assert_allclose(j_theta, -Q * n * OMEGA * r, rtol=1e-12)


@unittest.skipIf(np is None, 'numpy is required')
@unittest.skipIf(rd is None, 'scripts/rmf_diagnostics.py did not import')
class TestProfiles(unittest.TestCase):

    def test_radial_profile_recovers_a_linear_function(self):
        x, y = disc(n_r=50, n_t=120)
        r, _c, _s = rd.polar(x, y)
        edges = np.linspace(0.0, A, 21)
        centres, mean, _std, count = rd.radial_profile(r, 3.0 * r + 1.0, edges)
        self.assertTrue(np.all(count > 0), 'every bin should hold nodes')
        # Bin means of a linear function sit within half a bin width of the centre.
        np.testing.assert_allclose(mean, 3.0 * centres + 1.0,
                                   atol=3.0 * (edges[1] - edges[0]))

    def test_empty_bins_come_back_as_nan_not_zero(self):
        """A radius the run has no nodes at must not read as a measured zero."""
        r = np.array([0.001, 0.002])
        centres, mean, _s, count = rd.radial_profile(r, np.array([1.0, 1.0]),
                                                     np.linspace(0.0, A, 11))
        self.assertTrue(np.isnan(mean[5]), 'empty bin should be NaN')
        self.assertEqual(int(count[5]), 0)
        self.assertEqual(len(centres), 10)

    def test_axial_field_on_axis_recovers_a_uniform_field(self):
        x, y = disc()
        r, _c, _s = rd.polar(x, y)
        bz = np.full_like(x, -0.0123)
        value, n = rd.axial_field_on_axis(r, bz, 0.1 * A)
        self.assertAlmostEqual(value, -0.0123, places=12)
        self.assertGreater(n, 0)

    def test_rotation_parameter_recovers_rigid_rotation(self):
        """u_theta = zeta * omega * r must give exactly that zeta."""
        x, y = disc()
        r, _c, _s = rd.polar(x, y)
        for zeta0 in (1.0, 0.35, -0.6):
            u_theta = zeta0 * OMEGA * r
            got, n = rd.rotation_parameter(r, u_theta, OMEGA, 0.1 * A, 0.9 * A)
            self.assertGreater(n, 0)
            self.assertAlmostEqual(got, zeta0, places=10,
                                   msg=f'zeta = {zeta0} was not recovered')


@unittest.skipIf(np is None, 'numpy is required')
@unittest.skipIf(rd is None, 'scripts/rmf_diagnostics.py did not import')
class TestPenetratedLimit(unittest.TestCase):
    """The closed form the literature's penetrated-limit comparison uses."""

    def test_constant_zeta_gives_the_closed_form(self):
        """B_z(0) = -mu0 n e omega zeta a^2 / 2 for uniform zeta."""
        n = 1e20
        centres = np.linspace(0.0, A, 400)
        for zeta0 in (1.0, 0.4):
            got = rd.penetrated_limit_bz(centres, np.full_like(centres, zeta0),
                                         n, Q, OMEGA, A)
            want = -MU0 * n * Q * OMEGA * zeta0 * A * A / 2.0
            self.assertAlmostEqual(got / want, 1.0, places=4,
                                   msg=f'zeta = {zeta0}: got {got:.6e}, want {want:.6e}')

    def test_magnitude_matches_the_assessment_table(self):
        """The assessment quotes mu0 n e omega a^2 / 2 = 450 G for this deck."""
        got = rd.penetrated_limit_bz(np.linspace(0.0, A, 400),
                                     np.ones(400), 1e20, Q, OMEGA, A)
        self.assertAlmostEqual(abs(got) * 1e4, 452.0, delta=5.0,
                               msg=f'|B_z| = {abs(got) * 1e4:.1f} G, assessment says ~450 G')

    def test_penetrated_limit_sign_opposes_a_plus_z_bias(self):
        """THE item 0 assertion, in code rather than in a comment.

        Electrons dragged counter-clockwise (zeta > 0) carry a clockwise current,
        whose axial field is -z, which OPPOSES a +z bias - that is field
        reversal. Reverse the rotation and the driven field reinforces the bias
        instead. Whether a given deck's drive turns counter-clockwise is a
        separate question, asserted for the drives themselves in
        test/cxx/test_rmf_rotation_sense.cc; this is the link between the two.
        """
        centres = np.linspace(0.0, A, 200)
        b_ccw = rd.penetrated_limit_bz(centres, np.full(200, +0.8), 1e20, Q, OMEGA, A)
        b_cw = rd.penetrated_limit_bz(centres, np.full(200, -0.8), 1e20, Q, OMEGA, A)
        self.assertLess(b_ccw, 0.0,
                        'counter-clockwise electrons must drive B_z < 0, which is '
                        'what reverses a +z bias')
        self.assertGreater(b_cw, 0.0,
                           'clockwise electrons must drive B_z > 0, reinforcing a '
                           '+z bias rather than reversing it')
        self.assertAlmostEqual(b_ccw, -b_cw, places=12)

    def test_missing_bins_are_skipped_not_counted_as_zero(self):
        """NaN bins must not drag the integral down as if zeta were zero there."""
        centres = np.linspace(0.0, A, 100)
        full = np.ones(100)
        holed = full.copy()
        holed[40:45] = np.nan
        a = rd.penetrated_limit_bz(centres, full, 1e20, Q, OMEGA, A)
        b = rd.penetrated_limit_bz(centres, holed, 1e20, Q, OMEGA, A)
        self.assertAlmostEqual(a / b, 1.0, places=3,
                               msg='a gap in the profile changed the answer by more '
                                   'than interpolation across it would')


@unittest.skipIf(np is None, 'numpy is required')
@unittest.skipIf(rd is None, 'scripts/rmf_diagnostics.py did not import')
class TestPenetrationAndLayer(unittest.TestCase):

    def test_uniform_field_is_fully_penetrated(self):
        x, y = disc()
        r, _c, _s = rd.polar(x, y)
        self.assertAlmostEqual(
            rd.penetration_fraction(r, np.full_like(x, 5.0e-3), A), 1.0, places=10)

    def test_a_screened_field_reads_near_zero(self):
        """B_perp confined to a layer of thickness delta at the edge."""
        x, y = disc(n_r=200, n_t=64)
        r, _c, _s = rd.polar(x, y)
        delta = rd.skin_depth(0.5e-5, OMEGA)
        b = 5.0e-3 * np.exp(-(A - r) / delta)
        frac = rd.penetration_fraction(r, b, A)
        self.assertLess(frac, 1e-3,
                        f'a field screened to a {delta * 1e3:.2f} mm layer in a '
                        f'{A * 1e3:.0f} mm column read as {frac:.4f} penetrated')

    def test_layer_thickness_recovers_the_skin_depth(self):
        """J_theta ~ exp(-(a-r)/delta) must give back delta."""
        centres = np.linspace(0.0, A, 300)
        for delta in (1.0e-3, 2.5e-3, 5.0e-4):
            j = 1.0e6 * np.exp(-(A - centres) / delta)
            got = rd.current_layer_thickness(centres, j, A)
            self.assertAlmostEqual(got / delta, 1.0, places=3,
                                   msg=f'delta = {delta:.1e} came back as {got:.3e}')

    def test_a_penetrated_profile_reports_no_layer(self):
        """A flat or inward-growing J_theta has no skin layer; say so, do not fit one.

        Returning a fitted length here would be worse than returning nothing: it
        would be compared against delta and would look like a measurement.
        """
        centres = np.linspace(0.0, A, 100)
        self.assertTrue(math.isnan(
            rd.current_layer_thickness(centres, np.full(100, 1.0e6), A)))
        self.assertTrue(math.isnan(
            rd.current_layer_thickness(centres, 1.0e6 * (2.0 - centres / A), A)))

    def test_a_layer_wider_than_the_column_is_refused(self):
        """The defect real output exposed: 40 mm and 59 mm "layers" in a 30 mm column.

        A profile can decay inward slowly enough that the fitted length exceeds
        the plasma radius. That is not a skin layer, and reporting the number
        invites comparing it against delta - which for this deck is a/24.
        """
        centres = np.linspace(0.0, A, 100)
        gentle = 1.0e6 * np.exp(-(A - centres) / (5.0 * A))
        self.assertTrue(math.isnan(rd.current_layer_thickness(centres, gentle, A)),
                        'a decay length of 5a was reported as a layer')
        # Just inside the cut it is still accepted, so the guard is a bound and
        # not a blanket refusal.
        ok = 1.0e6 * np.exp(-(A - centres) / (0.3 * A))
        self.assertAlmostEqual(rd.current_layer_thickness(centres, ok, A) / (0.3 * A),
                               1.0, places=3)

    def test_a_profile_that_is_not_exponential_is_refused(self):
        """A "decay length" fitted to something that does not decay exponentially
        describes nothing, however neatly the least-squares comes out."""
        centres = np.linspace(0.0, A, 60)
        rng = np.random.default_rng(20260909)
        noise = 1.0e6 * (1.0 + 0.8 * rng.standard_normal(60)) ** 2 + 1.0
        self.assertTrue(math.isnan(rd.current_layer_thickness(centres, noise, A)))

    def test_skin_depth_closed_form(self):
        """delta = sqrt(2 eta / (mu0 omega)); the assessment quotes 1.26 mm."""
        self.assertAlmostEqual(rd.skin_depth(0.5e-5, OMEGA) * 1e3, 1.262, places=3)


@unittest.skipIf(np is None, 'numpy is required')
@unittest.skipIf(rd is None, 'scripts/rmf_diagnostics.py did not import')
class TestCycleAverage(unittest.TestCase):

    def test_a_full_cycle_of_a_sinusoid_averages_to_zero(self):
        period = 1.0 / 7.95e5
        t = np.linspace(0.0, 3.0 * period, 601)
        avg, whole = rd.cycle_average(t, np.sin(2 * math.pi * t / period), period)
        self.assertTrue(whole)
        self.assertAlmostEqual(avg, 0.0, places=6)

    def test_a_dc_offset_survives_the_average(self):
        """The secular part is what the literature's expressions describe."""
        period = 1.0 / 7.95e5
        t = np.linspace(0.0, 2.0 * period, 401)
        avg, whole = rd.cycle_average(
            t, 0.25 + np.sin(2 * math.pi * t / period), period)
        self.assertTrue(whole)
        self.assertAlmostEqual(avg, 0.25, places=6)

    def test_a_partial_cycle_is_flagged(self):
        """Half a cycle of a sinusoid averages to about 2/pi, not to 0.

        That gap is the whole reason for the flag: the same quantity reads 0.64
        or 0.00 depending only on which part of the cycle the frames caught. The
        expected value is 2/pi for the continuous average; an evenly spaced
        discrete mean including both endpoints is a per-cent below it, so the
        tolerance is set to that rather than to floating point.
        """
        period = 1.0 / 7.95e5
        t = np.linspace(0.0, 0.5 * period, 101)
        avg, whole = rd.cycle_average(t, np.sin(2 * math.pi * t / period), period)
        self.assertFalse(whole, 'half a period must not be reported as a cycle average')
        self.assertAlmostEqual(avg, 2.0 / math.pi, delta=0.02)


@unittest.skipIf(rd is None, 'scripts/rmf_diagnostics.py did not import')
class TestFrameTiming(unittest.TestCase):
    """Frame times come from the deck's OUT, not from how many frames exist.

    The solver has no checkpoint/restart, so a partial set of frames is the
    normal thing to be handed. Timing them by their own count would stretch them
    to fill TEND: the last frame of an interrupted run would be reported at the
    end time it never reached, and every cycle average taken from it would be
    over the wrong window.
    """

    PARAMS = {'TEND': 2.0e-6, 'OUT': 10}

    def test_a_complete_run_is_evenly_spaced_to_TEND(self):
        times = [rd.frame_time(i, self.PARAMS) for i in range(11)]
        self.assertAlmostEqual(times[0], 0.0)
        self.assertAlmostEqual(times[-1], 2.0e-6)
        for a, b in zip(times, times[1:]):
            self.assertAlmostEqual(b - a, 2.0e-7, places=15)

    def test_a_partial_run_keeps_the_same_times(self):
        """The killer case: four frames of a ten-frame run are still frames 0-3."""
        full = [rd.frame_time(i, self.PARAMS) for i in range(11)]
        partial = [rd.frame_time(i, self.PARAMS) for i in range(4)]
        self.assertEqual(partial, full[:4])
        self.assertAlmostEqual(partial[-1], 6.0e-7,
                               msg='frame 3 of a TEND = 2e-6, OUT = 10 run is at '
                                   '6e-7 s whether or not the run finished')

    def test_a_deck_without_OUT_gives_nan_rather_than_a_wrong_time(self):
        self.assertTrue(math.isnan(rd.frame_time(3, {'TEND': 1.0})))
        self.assertTrue(math.isnan(rd.frame_time(3, {'TEND': 1.0, 'OUT': 0})))


if __name__ == '__main__':
    unittest.main()
