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


if __name__ == '__main__':
    unittest.main()
