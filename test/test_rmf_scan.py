#!/usr/bin/env python3
"""Verification: does the Phase 3 cost model predict the timestep the solver picks?

scripts/rmf_scan.py refuses to start a scan it estimates at more than a few
hours, and prints what a scan would cost before running it. Both depend on a
model of Apollo's timestep, and a cost model that is wrong in the cheap
direction is worse than none: it would talk someone into launching a run that
cannot finish, on a solver with no checkpoint/restart.

So the model is checked here against timesteps the solver actually chose. The
reference values below were measured by running
examples/unstructuredDG/multifluid/rmf_frc/frc2d.pin on optimizedCircle2.msh
with src/build-opt/apollo and reading the "Simulation dt" it printed:

    Te [eV]   LIGHT [m/s]   dt printed [s]
    30        3.0e6         2.73658e-11
    30        1.0e6         2.76826e-11
    30        3.0e5         2.76826e-11
    30        1.0e5         2.76826e-11
    10        1.732e6       4.74004e-11
    3.3333    1.0e6         8.20975e-11

The interesting row is the third: dropping the speed of light by a further
factor of ten changes nothing at all. dt has saturated on the ELECTRON SOUND
SPEED, sqrt(gamma k T_e/m_e) = 2.9657e6 m/s at 30 eV, which is 1.1% below the
deck's 3.0e6 m/s "speed of light". Reducing c is the standard way to make a
two-fluid run affordable and here it is worth 1.16% in total. That is the single
fact that sets what Phase 3 can cost, so it is pinned by
test_reducing_the_speed_of_light_saturates below rather than left in a comment.

Runs in milliseconds; no solver needed.
"""

import math
import os
import sys
import unittest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'scripts'))

try:
    import rmf_scan
except ImportError:
    rmf_scan = None

MESH = os.path.join(REPO, 'examples', 'unstructuredDG', 'multifluid', 'rmf_frc',
                    'optimizedCircle2.msh')
DECK = os.path.join(REPO, 'examples', 'unstructuredDG', 'multifluid', 'rmf_frc',
                    'frc2d.pin')

# (Te [eV], LIGHT [m/s], dt the solver printed [s])
MEASURED = [
    (30.0, 3.0e6, 2.73658e-11),
    (30.0, 1.0e6, 2.76826e-11),
    (30.0, 3.0e5, 2.76826e-11),
    (30.0, 1.0e5, 2.76826e-11),
    (10.0, 1.732e6, 4.74004e-11),
    (3.3333, 1.0e6, 8.20975e-11),
]


@unittest.skipIf(rmf_scan is None, 'scripts/rmf_scan.py did not import')
@unittest.skipIf(not os.path.exists(MESH), 'the rmf_frc mesh is not present')
class TestCostModel(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.dtscale, cls.n_tri = rmf_scan.mesh_dtscale(MESH)

    def params(self, te, light):
        return {'Te': te, 'LIGHT': light, 'Q': 1.6e-19, 'ME': 1.67e-27 / 1836.0}

    def test_mesh_geometry_matches_the_shipped_mesh(self):
        """dtscale is min over triangles of area/(perimeter/2), as the solver does."""
        self.assertEqual(self.n_tri, 7792)
        self.assertAlmostEqual(self.dtscale * 1e3, 0.18472, places=4)

    def test_predicts_every_measured_timestep(self):
        for te, light, dt_measured in MEASURED:
            got = rmf_scan.timestep(self.params(te, light), self.dtscale, 1)
            self.assertAlmostEqual(
                got / dt_measured, 1.0, places=5,
                msg=f'Te = {te} eV, LIGHT = {light:.3e}: model {got:.6e}, '
                    f'solver printed {dt_measured:.6e}')

    def test_reducing_the_speed_of_light_saturates(self):
        """The finding that decides what Phase 3 costs.

        Below the electron sound speed, lowering c buys nothing: the timestep
        stops moving. Measured over 3.0e6 -> 1.0e5, a factor of 30, dt changes
        by 1.16% and then not at all.
        """
        p = self.params(30.0, 3.0e6)
        c_se = math.sqrt(rmf_scan.GAMMA * p['Te'] * p['Q'] / p['ME'])
        self.assertAlmostEqual(c_se, 2.9657e6, delta=1e3)
        self.assertLess(c_se, 3.0e6,
                        'at 30 eV the electron sound speed should sit just below '
                        "the deck's reduced speed of light")

        dt_full = rmf_scan.timestep(p, self.dtscale, 1)
        dts = [rmf_scan.timestep(self.params(30.0, c), self.dtscale, 1)
               for c in (1.0e6, 3.0e5, 1.0e5, 1.0e3)]
        for dt in dts:
            self.assertAlmostEqual(dt, dts[0], places=15,
                                   msg='below the electron sound speed dt must not '
                                       'depend on the speed of light at all')
        self.assertLess(dts[0] / dt_full, 1.02,
                        f'a 3000x cut in c bought {100 * (dts[0] / dt_full - 1):.2f}%; '
                        f'if this is ever large, the electron sound speed stopped '
                        f'being the constraint and the assessment needs revisiting')

    def test_the_limiting_speed_is_named_correctly(self):
        _v, c_se = rmf_scan.max_wave_speed(self.params(30.0, 3.0e6))
        v_hot, _ = rmf_scan.max_wave_speed(self.params(30.0, 3.0e6))
        v_cold, _ = rmf_scan.max_wave_speed(self.params(30.0, 1.0e6))
        self.assertAlmostEqual(v_hot, 3.0e6)
        self.assertAlmostEqual(v_cold, c_se)

    def test_higher_orders_are_refused_rather_than_guessed(self):
        """rMin was measured at P_ORDER = 1 only; a guess would misprice a scan."""
        with self.assertRaises(SystemExit):
            rmf_scan.timestep(self.params(30.0, 3.0e6), self.dtscale, 2)


@unittest.skipIf(rmf_scan is None, 'scripts/rmf_scan.py did not import')
@unittest.skipIf(not os.path.exists(DECK), 'the rmf_frc deck is not present')
class TestScaling(unittest.TestCase):
    """The (Te -> Te/s^2, c -> c/s) rewrite --speedup applies."""

    def setUp(self):
        with open(DECK) as fh:
            self.text = fh.read()
        self.dtscale, _ = rmf_scan.mesh_dtscale(MESH)

    def test_speedup_multiplies_the_timestep_by_s(self):
        import rmf_diagnostics as rd
        base = rd.deck_parameters_from_text(self.text)
        dt0 = rmf_scan.timestep(base, self.dtscale, 1)
        for s in (2.0, 3.0, 5.0):
            new, _note = rmf_scan.apply_speedup(self.text, s)
            p = rd.deck_parameters_from_text(new)
            self.assertAlmostEqual(
                rmf_scan.timestep(p, self.dtscale, 1) / dt0, s, places=6,
                msg=f'--speedup {s} should give exactly {s}x the timestep')

    def test_speedup_preserves_the_rmf_parameters(self):
        """gamma, lambda and the Debye length in cells must not move.

        Those are what the penetration literature's expressions are written in;
        a speedup that changed them would be scanning a different plasma.
        """
        import rmf_diagnostics as rd
        base = rd.deck_parameters_from_text(self.text)
        new, _n = rmf_scan.apply_speedup(self.text, 3.0)
        scaled = rd.deck_parameters_from_text(new)

        def invariants(p):
            eps0 = 1.0 / (p['MU0'] * p['LIGHT'] ** 2)
            nu_ei = p['n_dens'] * p['Q'] ** 2 * p['ETA'] / p['ME']
            w = 2.0 * math.pi * p['omega']
            return dict(
                # omega_ce in the ROTATING field, which is the literature's
                # convention (see the assessment, section 3.2). Only invariance
                # is being checked here, so the choice does not change the
                # verdict - but it should not read as if it did.
                gamma=p['Q'] * 50e-4 / p['ME'] / nu_ei,
                lam=p['RAD_PLASMA'] / math.sqrt(2 * p['ETA'] / (p['MU0'] * w)),
                debye=math.sqrt(eps0 * p['Te'] * p['Q'] / (p['n_dens'] * p['Q'] ** 2)),
                skin=p['LIGHT'] / math.sqrt(
                    p['n_dens'] * p['Q'] ** 2 / (eps0 * p['ME'])),
            )

        a, b = invariants(base), invariants(scaled)
        for key in a:
            self.assertAlmostEqual(b[key] / a[key], 1.0, places=9,
                                   msg=f'{key} moved under --speedup')

    def test_speedup_reduces_beta_and_that_is_the_price(self):
        """beta is the one thing the scaling does move: by 1/s^2. Say so, do not hide it."""
        import rmf_diagnostics as rd
        base = rd.deck_parameters_from_text(self.text)
        new, note = rmf_scan.apply_speedup(self.text, 3.0)
        scaled = rd.deck_parameters_from_text(new)
        beta = lambda p: 2 * p['MU0'] * p['n_dens'] * p['Te'] * p['Q'] / (60e-4) ** 2
        self.assertAlmostEqual(beta(scaled) / beta(base), 1.0 / 9.0, places=9)
        self.assertIn('beta', note, 'the note must tell the user beta changed')

    def test_every_written_number_is_a_float_literal(self):
        """Apollo's parser types "1000000" as an integer and get<REAL> then throws.

        The failure is std::bad_cast at setup, naming no key. This is not
        hypothetical: the first deck this harness generated set LIGHT = 1.0e6,
        wrote "LIGHT = 1000000", and the run aborted before its first step with
        c0 = 1000000 in the generated .inp.
        """
        for v in (1.0e6, 3.0e5, 30.0, 2.0, 1e-3, 0.0, 5e-4, 1e20, 1.5):
            text = rmf_scan.deck_number(v)
            self.assertTrue('.' in text or 'e' in text,
                            f'{v!r} was written as {text!r}, which the deck parser '
                            f'would read as an integer')
            self.assertEqual(float(text), float(v))

    def test_speedup_writes_float_literals(self):
        """The whole deck after a rewrite, not just the helper."""
        new, _n = rmf_scan.apply_speedup(self.text, 3.0)
        for name in ('Te', 'LIGHT'):
            line = [l for l in new.splitlines() if l.startswith(name + ' ')][0]
            value = line.split('=', 1)[1].split('#')[0].strip()
            self.assertTrue('.' in value or 'e' in value,
                            f'{line.strip()!r} would be parsed as an integer')

    def test_set_param_writes_float_literals(self):
        for value in (1.0e6, 30.0, 2.0):
            new = rmf_scan.set_param(self.text, 'Bomega', value)
            line = [l for l in new.splitlines() if l.startswith('Bomega')][0]
            got = line.split('=', 1)[1].split('#')[0].strip()
            self.assertTrue('.' in got or 'e' in got,
                            f'{line.strip()!r} would be parsed as an integer')

    def test_a_deck_without_the_expected_assignments_is_refused(self):
        """Half-applying the scaling would silently compare different plasmas."""
        with self.assertRaises(SystemExit):
            rmf_scan.apply_speedup('Te = 30.\n# no LIGHT here\n', 3.0)


if __name__ == '__main__':
    unittest.main()
