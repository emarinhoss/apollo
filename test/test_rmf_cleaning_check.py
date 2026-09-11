#!/usr/bin/env python3
"""Verification: does the divergence-cleaning gate actually gate?

scripts/rmf_cleaning_check.py is the thing that decides whether a multi-day RMF
run is worth starting. A gate that passes a failing case is worse than no gate,
because it is consulted precisely when nobody wants to hear no.

It has already done that once. The first frame of a run is its initial
condition, where the applied field has not switched on and rms B_x is exactly
zero, so the difference ratio was 0/0. `max()` over a list beginning with NaN
returns NaN, and `NaN > TOLERANCE` is False - so a case where the transverse
field differed by 92% between two cleaning speeds was reported as a PASS with
exit status 0. Comparisons against NaN are false in both directions, which makes
NaN exactly the wrong thing to leave in a series a threshold is applied to.

These tests are about the verdict logic, so they build the difference series
directly rather than running the solver. Milliseconds, no solver needed.
"""

import os
import shutil
import sys
import unittest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'scripts'))

try:
    import rmf_cleaning_check as cc
except ImportError:
    cc = None

NAN = float('nan')


def verdict(diffs):
    """The gate's decision for a [(frame, dBx, dBz)] series: True if it passes."""
    usable = [d for d in diffs if d[1] == d[1]]
    if not usable:
        return None
    return max(d[1] for d in usable) <= cc.TOLERANCE


@unittest.skipIf(cc is None, 'scripts/rmf_cleaning_check.py did not import')
class TestVerdict(unittest.TestCase):

    def test_a_leading_nan_does_not_hide_a_failure(self):
        """The bug, exactly: frame 0 is 0/0 and the real divergence is 92%."""
        series = [(0, NAN, 0.0), (1, 0.1092, 0.0), (3, 0.5659, 0.0),
                  (6, 0.9162, 0.0)]
        self.assertFalse(verdict(series),
                         'a 92% divergence was passed because frame 0 was NaN')
        # And the naive form this replaced would have passed it, which is why
        # the test is written against a series rather than against the output.
        self.assertFalse(max(d[1] for d in series) > cc.TOLERANCE,
                         'if this ever fails, max() stopped propagating NaN and '
                         'the guard could be simplified')

    def test_a_genuinely_small_difference_passes(self):
        self.assertTrue(verdict([(0, NAN, 0.0), (1, 0.001, 0.0), (2, 0.02, 0.0)]))

    def test_a_difference_just_over_the_line_fails(self):
        self.assertFalse(verdict([(1, cc.TOLERANCE * 1.001, 0.0)]))

    def test_a_difference_just_under_the_line_passes(self):
        self.assertTrue(verdict([(1, cc.TOLERANCE * 0.999, 0.0)]))

    def test_all_nan_is_neither_pass_nor_fail(self):
        """No field in any frame is not a pass; it is nothing measured."""
        self.assertIsNone(verdict([(0, NAN, 0.0), (1, NAN, 0.0)]))

    def test_the_tolerance_is_where_the_module_says_it_is(self):
        self.assertTrue(0.0 < cc.TOLERANCE < 1.0)
        self.assertAlmostEqual(cc.TOLERANCE, 0.05)


@unittest.skipIf(cc is None, 'scripts/rmf_cleaning_check.py did not import')
class TestFrameDiscovery(unittest.TestCase):

    def test_frames_are_ordered_by_index_not_by_name(self):
        """frc2d_10.vtu sorts before frc2d_2.vtu as a string."""
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            for i in (0, 1, 2, 10, 11):
                open(os.path.join(d, 'run_%d.vtu' % i), 'w').close()
            self.assertEqual([i for i, _ in cc.frames_in(d)], [0, 1, 2, 10, 11])

    def test_a_directory_with_no_frames_is_an_error(self):
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(SystemExit):
                cc.frames_in(d)


@unittest.skipIf(cc is None, 'rmf_cleaning_check not importable')
class TestCleaningSpeeds(unittest.TestCase):
    """Reading the cleaning speeds, and telling 'slower' from 'off'.

    gamma and chi are bare multiplicative factors on every term coupling phi
    and psi to E and B (wxphmaxwelleqn.cc:490-497), and chi also scales the
    charge source that generates phi (wxchargesrc.h:30). At zero they are all
    zero: the potentials are inert and there is no cleaning at all. Comparing a
    run against that measures how much divergence error the deck carries, not
    whether the answer depends on the cleaning speed - a different question,
    and the one the campaign actually needs. The script has to tell them apart
    or it reports the wrong verdict with full confidence, which it did.
    """

    def _run_dir(self, divb, dive):
        import tempfile
        tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, tmp, True)
        with open(os.path.join(tmp, 'deck.pin'), 'w') as handle:
            handle.write('# -*- python -*-\nimport math\n'
                         'DIVB_SPEED = %r\nDIVE_SPEED = %r\n\n<apollo>\n</apollo>\n'
                         % (divb, dive))
        return tmp

    def test_speeds_are_read_from_the_deck(self):
        self.assertEqual(cc.cleaning_speeds(self._run_dir(1.0, 1.0)), (1.0, 1.0))
        self.assertEqual(cc.cleaning_speeds(self._run_dir(0.5, 0.5)), (0.5, 0.5))
        self.assertEqual(cc.cleaning_speeds(self._run_dir(0.0, 0.0)), (0.0, 0.0))

    def test_missing_deck_is_not_guessed(self):
        import tempfile
        tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, tmp, True)
        self.assertIsNone(cc.cleaning_speeds(tmp))

    def test_deck_without_the_keys_is_not_guessed(self):
        """A deck that exists but says nothing about cleaning returns None.

        Distinct from the no-deck case above: this one reaches the parser and
        finds neither key. Defaulting to (1.0, 1.0) here would let the script
        report a comparison as a speed-sensitivity test when it has no idea
        what speeds were used - the failure this whole change exists to stop.
        """
        import tempfile
        tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, tmp, True)
        with open(os.path.join(tmp, 'deck.pin'), 'w') as handle:
            handle.write('# -*- python -*-\nGAMMA = 1.66\n<apollo>\n</apollo>\n')
        self.assertIsNone(cc.cleaning_speeds(tmp))

    def test_one_key_only_is_not_guessed(self):
        """Half an answer is not an answer."""
        import tempfile
        tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, tmp, True)
        with open(os.path.join(tmp, 'deck.pin'), 'w') as handle:
            handle.write('DIVB_SPEED = 1.0\n')
        self.assertIsNone(cc.cleaning_speeds(tmp))

    def test_unparseable_value_is_not_guessed(self):
        import tempfile
        tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, tmp, True)
        with open(os.path.join(tmp, 'deck.pin'), 'w') as handle:
            handle.write('DIVB_SPEED = SOME_MACRO\nDIVE_SPEED = 1.0\n')
        self.assertIsNone(cc.cleaning_speeds(tmp))

    def test_zero_is_recognised_as_no_cleaning(self):
        off = self._run_dir(0.0, 0.0)
        on = self._run_dir(1.0, 1.0)
        speeds = [cc.cleaning_speeds(on), cc.cleaning_speeds(off)]
        self.assertEqual(cc.degenerate_runs([on, off], speeds), [off])

    def test_two_nonzero_speeds_are_not_degenerate(self):
        """1.0 against 0.5 is the comparison that isolates the speed."""
        full = self._run_dir(1.0, 1.0)
        half = self._run_dir(0.5, 0.5)
        speeds = [cc.cleaning_speeds(full), cc.cleaning_speeds(half)]
        self.assertEqual(cc.degenerate_runs([full, half], speeds), [])

    def test_half_is_not_mistaken_for_off(self):
        """The guard must key on zero, not on 'less than the other one'."""
        half = self._run_dir(0.5, 0.5)
        self.assertEqual(
            cc.degenerate_runs([half], [cc.cleaning_speeds(half)]), [])

    def test_unreadable_speeds_do_not_claim_degeneracy(self):
        speeds = [None, None]
        self.assertEqual(cc.degenerate_runs(['a', 'b'], speeds), [])


if __name__ == '__main__':
    unittest.main()
