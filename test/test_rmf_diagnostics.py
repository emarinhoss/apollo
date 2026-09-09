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
        centres, mean, _std, count, outside = rd.radial_profile(r, 3.0 * r + 1.0, edges)
        self.assertEqual(outside, 0)
        self.assertTrue(np.all(count > 0), 'every bin should hold nodes')
        # Bin means of a linear function sit within half a bin width of the centre.
        np.testing.assert_allclose(mean, 3.0 * centres + 1.0,
                                   atol=3.0 * (edges[1] - edges[0]))

    def test_points_on_the_outer_boundary_are_kept(self):
        """np.digitize puts r >= edges[-1] past the last bin.

        On optimizedCircle2.msh the plasma boundary is at exactly RAD_PLASMA and
        the nodal radii computed from it spread over the last few ulps either
        side, so 396 of 608 boundary nodes fell out of every bin. That ring is
        where the RMF enters and where the skin layer lives.
        """
        edges = np.linspace(0.0, A, 11)
        r = np.array([0.0, 0.5 * A, A * (1 - 1e-16), A, A * (1 + 1e-16)])
        _c, _m, _s, count, outside = rd.radial_profile(r, np.ones(5), edges)
        self.assertEqual(outside, 0, 'a boundary node was dropped')
        self.assertEqual(int(count.sum()), 5)
        self.assertEqual(int(count[-1]), 3, 'the three boundary points belong in '
                                            'the last bin')

    def test_points_genuinely_outside_are_dropped_and_counted(self):
        """Kept separate from the rounding case: 2a is not a boundary node."""
        edges = np.linspace(0.0, A, 11)
        r = np.array([0.5 * A, 2.0 * A, -1.0])
        _c, _m, _s, count, outside = rd.radial_profile(r, np.ones(3), edges)
        self.assertEqual(outside, 2)
        self.assertEqual(int(count.sum()), 1)

    def test_empty_bins_come_back_as_nan_not_zero(self):
        """A radius the run has no nodes at must not read as a measured zero."""
        r = np.array([0.001, 0.002])
        centres, mean, _s, count, _out = rd.radial_profile(
            r, np.array([1.0, 1.0]), np.linspace(0.0, A, 11))
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

    def test_axial_field_on_axis_averages_the_axis_and_not_the_domain(self):
        """A uniform field cannot show which region was averaged.

        Fed a constant B_z, every possible mask returns that constant, so the
        test above passes with the selection pointed anywhere - at the whole
        domain, or at everything OUTSIDE the axis disc. That matters because
        this function produces the headline number of Phase 3 item 0: on a
        field-reversed column B_z on axis is the opposite sign to B_z at the
        edge, so averaging the wrong region does not give a slightly wrong
        answer, it gives the wrong sign.
        """
        x, y = disc(n_r=120, n_t=120)
        r, _c, _s = rd.polar(x, y)
        # Reversed core, positive edge: a caricature of a formed FRC.
        bz = np.where(r < 0.4 * A, -1.0e-2, +6.0e-3)
        value, n = rd.axial_field_on_axis(r, bz, 0.1 * A)
        self.assertGreater(n, 0)
        self.assertAlmostEqual(value, -1.0e-2, places=12,
                               msg='the axis average picked up the edge field')
        # And it must not simply be returning the global mean.
        self.assertNotAlmostEqual(value, float(np.mean(bz)), places=6)

    def test_rotation_parameter_reads_the_annulus_it_says_it_does(self):
        """Rigid rotation gives the same zeta on any window, so it proves nothing.

        u_theta = zeta*omega*r returns zeta over any selection whatsoever, so
        the test above cannot see the r-window. A zeta that varies with radius
        can: here zeta is 1 in the sampled annulus and 0 outside it, so a
        function reading the wrong region returns a diluted number.
        """
        x, y = disc(n_r=200, n_t=64)
        r, _c, _s = rd.polar(x, y)
        band = (r > 0.1 * A) & (r < 0.9 * A)
        u_theta = np.where(band, 1.0 * OMEGA * r, 0.0)
        got, n = rd.rotation_parameter(r, u_theta, OMEGA, 0.1 * A, 0.9 * A)
        self.assertGreater(n, 0)
        self.assertAlmostEqual(got, 1.0, places=12,
                               msg='zeta was averaged over more than the stated '
                                   'annulus, so the zeros outside it diluted it')

    def test_rotation_parameter_divides_by_the_local_radius(self):
        """zeta = u_theta/(omega r) pointwise, not u_theta/(omega <r>)."""
        x, y = disc(n_r=200, n_t=64)
        r, _c, _s = rd.polar(x, y)
        # Solid-body u_theta = k r gives zeta = k/omega everywhere; a CONSTANT
        # u_theta does not, and dividing by a mean radius would not notice.
        u_theta = np.full_like(r, 1.0e5)
        got, _n = rd.rotation_parameter(r, u_theta, OMEGA, 0.1 * A, 0.9 * A)
        sel = (r > 0.1 * A) & (r < 0.9 * A)
        want = float(np.mean(1.0e5 / (OMEGA * r[sel])))
        self.assertAlmostEqual(got, want, places=12)
        self.assertNotAlmostEqual(got, 1.0e5 / (OMEGA * float(np.mean(r[sel]))),
                                  places=6)

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

    def test_bin_centres_integrate_over_the_whole_column(self):
        """The grid analyse() actually uses, not the one the tests used to use.

        radial_profile returns bin CENTRES, which stop half a bin short of
        r = a - and the integrand is zeta*r, so that outer sliver carries more
        weight than anywhere else. Integrating between the first and last centre
        read 2.5% low for uniform zeta and 24% low for a screened profile, and
        every test passed a full-range linspace and so never saw it.
        """
        n = 1e20
        for nbins in (20, 40, 400):
            edges = np.linspace(0.0, A, nbins + 1)
            centres = 0.5 * (edges[:-1] + edges[1:])
            got = rd.penetrated_limit_bz(centres, np.ones(nbins), n, Q, OMEGA, A)
            want = -MU0 * n * Q * OMEGA * A * A / 2.0
            self.assertAlmostEqual(got / want, 1.0, places=6,
                                   msg=f'{nbins} bin centres: {got / want:.6f} of the '
                                       f'closed form')

    def test_a_screened_profile_is_integrated_to_within_a_per_cent(self):
        """The case that was 24% low: all the current sits in the outer sliver."""
        n, delta = 1e20, 1.262e-3
        edges = np.linspace(0.0, A, 41)
        centres = 0.5 * (edges[:-1] + edges[1:])
        got = rd.penetrated_limit_bz(centres, np.exp(-(A - centres) / delta),
                                     n, Q, OMEGA, A)
        fine = np.linspace(0.0, A, 20001)
        want = -MU0 * n * Q * OMEGA * float(
            np.trapezoid(np.exp(-(A - fine) / delta) * fine, fine))
        self.assertAlmostEqual(got / want, 1.0, places=1,
                               msg=f'{got / want:.4f} of a converged integral')

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

    # A power law that decays inward over the fit window without the extreme
    # dynamic range that would trip the noise floor first. It is the dangerous
    # case: the full-window fit gives 1.48 mm, within 20% of this deck's actual
    # skin depth of 1.26 mm, so the number it would report is not obviously
    # wrong to anyone reading it.
    PLAUSIBLE_POWER_LAW = staticmethod(
        lambda centres: 1.0 / (A - centres + 1.0e-3) ** 3)

    def test_a_power_law_is_refused_although_it_fits_beautifully(self):
        """The case that defeated the R^2 guard this replaced.

        A power law is smooth and monotone in log space, so R^2 came out 1.00000
        for one and the function reported a layer. This one would report 1.48 mm
        against a real delta of 1.26 mm - not a number anyone would question.
        The two-half consistency test refuses it, because a power law's decay
        length varies through the window (0.91 mm in the outer half, 2.51 mm in
        the inner) and an exponential's does not.
        """
        centres = np.linspace(0.0, A, 100)
        self.assertTrue(
            math.isnan(rd.current_layer_thickness(
                centres, self.PLAUSIBLE_POWER_LAW(centres), A)),
            'a power law was accepted as an exponential layer')

    def test_an_exponential_on_a_pedestal_is_refused(self):
        """Two scales in one profile has no single decay length."""
        centres = np.linspace(0.0, A, 100)
        j = 1.0e6 * (np.exp(-(A - centres) / 1.262e-3) + 0.1)
        self.assertTrue(math.isnan(rd.current_layer_thickness(centres, j, A)))

    def test_the_consistency_guard_is_what_rejects_the_power_law(self):
        """Guard against the guard being dead code.

        The R^2 check this replaced was exercised by no test at all: setting its
        threshold to -1e9 left the whole suite passing, which is how it survived
        being useless. Loosening this one has to change the answer, or the same
        thing has happened again.
        """
        centres = np.linspace(0.0, A, 100)
        power = self.PLAUSIBLE_POWER_LAW(centres)
        self.assertTrue(math.isnan(rd.current_layer_thickness(centres, power, A)))
        loosened = rd.current_layer_thickness(centres, power, A, max_ratio=1e12)
        self.assertFalse(math.isnan(loosened),
                         'disabling the consistency guard changed nothing, so it is '
                         'not the guard that rejects a power law')
        self.assertAlmostEqual(loosened * 1e3, 1.48, delta=0.05,
                               msg='the number the old guard would have reported')

    def test_a_profile_too_spiky_to_check_is_refused(self):
        """A near-singular profile leaves too few samples above the noise floor.

        |J| ~ (a-r)^-2 with no offset spans eight decades over the fit window, so
        the floor filter leaves two points - enough to fit a line through, not
        enough to establish that it is an exponential. Refused for that reason
        rather than by the consistency test, and that distinction is why this is
        a separate case from the one above.
        """
        centres = np.linspace(0.0, A, 100)
        spike = 1.0e-2 / np.maximum(A - centres, 1e-5) ** 2
        self.assertTrue(math.isnan(rd.current_layer_thickness(centres, spike, A)))
        self.assertTrue(math.isnan(
            rd.current_layer_thickness(centres, spike, A, max_ratio=1e12)),
            'this one must be refused even with the consistency guard disabled')

    def test_the_fit_uses_the_outer_half_only(self):
        """The layer lives at the plasma edge; the interior is a different problem.

        Fed a profile that is exponential outside 0.5a and flat inside it - a
        skin layer sitting on a penetrated core - restricting the fit to the
        outer half recovers delta. Widening it to the whole column mixes two
        scales and the consistency guard refuses, which is the right answer for
        a fit but the wrong answer for the question. Without this, removing the
        r > 0.5a restriction changed nothing any test could see.
        """
        centres = np.linspace(0.0, A, 200)
        # 3 mm, not the deck's 1.26: at 1.26 mm the flat core sits at 7e-6 of
        # the peak and the noise floor removes it anyway, so the window makes no
        # difference and the test could not see it either.
        delta = 3.0e-3
        j = 1.0e6 * np.where(centres > 0.5 * A,
                             np.exp(-(A - centres) / delta),
                             math.exp(-(0.5 * A) / delta))
        got = rd.current_layer_thickness(centres, j, A)
        self.assertAlmostEqual(got / delta, 1.0, places=2,
                               msg=f'the outer-half fit gave {got * 1e3:.4f} mm '
                                   f'for a {delta * 1e3:.3f} mm layer')

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

    def test_a_cosine_averages_to_zero_too(self):
        """The endpoint double-count, which a sine cannot reveal.

        A closed window over one period holds both endpoints, which are the same
        phase. For a sine they sit on the zero crossings and contribute nothing,
        which is why the two original tests passed while the bug was live; for a
        cosine they sit on the extremum and leak amplitude/(N+1) - measured at
        +0.004975 = 1/201 on this grid.
        """
        period = 1.0 / 7.95e5
        t = np.linspace(0.0, 3.0 * period, 601)
        for name, fn in (('cos', np.cos), ('sin', np.sin)):
            avg, whole = rd.cycle_average(t, fn(2 * math.pi * t / period), period)
            self.assertTrue(whole)
            self.assertAlmostEqual(avg, 0.0, places=9,
                                   msg=f'{name} over whole periods should average to 0, '
                                       f'got {avg:.6e}')

    def test_a_phase_shifted_wave_averages_to_zero_at_any_phase(self):
        """Nothing about the answer may depend on where the window happens to start."""
        period = 1.0 / 7.95e5
        t = np.linspace(0.0, 2.0 * period, 401)
        for phase in (0.0, 0.3, 1.0, 2.5, math.pi):
            avg, whole = rd.cycle_average(
                t, np.sin(2 * math.pi * t / period + phase), period)
            self.assertTrue(whole)
            self.assertAlmostEqual(avg, 0.0, places=9,
                                   msg=f'phase {phase}: got {avg:.6e}')

    def test_too_few_samples_per_period_is_not_a_cycle_average(self):
        """Averaging a period sampled twice is aliasing, not averaging.

        A scan that rewrites TEND for more periods but leaves the deck's OUT at
        10 writes 10/periods frames per period. At five periods that is two.
        """
        period = 1.0 / 7.95e5
        t = np.linspace(0.0, 5.0 * period, 11)          # 2 frames per period
        _avg, whole = rd.cycle_average(t, np.cos(2 * math.pi * t / period), period)
        self.assertFalse(whole, 'two samples per period was called a cycle average')
        t_ok = np.linspace(0.0, 5.0 * period, 201)      # 40 per period
        _avg2, whole2 = rd.cycle_average(t_ok, np.cos(2 * math.pi * t_ok / period), period)
        self.assertTrue(whole2)

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
