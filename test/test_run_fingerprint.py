#!/usr/bin/env python3
"""Tests for scripts/run_fingerprint.py.

The tool exists to answer one question - "did these two builds produce the same
answer?" - and the only thing a caller reads is the exit status. So these tests
are about the exit status, not the prose.

They are written as mutations: take a fingerprint that compares clean against
itself, break exactly one thing, and require that the tool notices. That is the
convention this repository settled on after `current_layer_thickness` shipped an
R-squared check that accepted a power law with R2 = 1.00000 and that no test
touched.

The mutations here are not hypothetical. Every one of the four "cannot compare"
cases below - a renamed array, a length mismatch, a changed mesh, non-finite
values - passed with exit 0 in the tool's first draft: it printed the refusal
and then returned success because the arrays that *did* line up agreed. A CI
step would have read that as agreement.
"""

import copy
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
TOOL = os.path.join(ROOT, 'scripts', 'run_fingerprint.py')

AGREE = 0        # the two runs agree to round-off
DIFFER = 1       # they are comparable, and they disagree
REFUSED = 2      # they cannot be compared at all


def _array(n=4, lo=-1.0, hi=1.0, total=0.25, nonfinite=0):
    return {'n': n, 'n_nonfinite': nonfinite, 'min': lo, 'max': hi,
            'sorted_sum': total, 'abs_sorted_sum': abs(total) + 1.0}


def _fingerprint():
    """Two frames, a mesh array and two solution arrays each."""
    frame = {
        'file': 'run_0.vtu',
        'arrays': {
            'Position': _array(12, -0.05, 0.05, 0.0),
            'solutiondg.0': _array(4, 9.0e-11, 9.1e-11, 3.6e-10),
            'solutiondg.12': _array(4, -67.5, 65.1, -2.4),
        },
    }
    second = copy.deepcopy(frame)
    second['file'] = 'run_1.vtu'
    return {'source': '/nowhere', 'banner': 'Apollo test, PETSc 3.19.6',
            'frames': [frame, second]}


class TestRunFingerprint(unittest.TestCase):

    def _compare(self, left, right):
        with tempfile.TemporaryDirectory() as tmp:
            lp = os.path.join(tmp, 'left.json')
            rp = os.path.join(tmp, 'right.json')
            with open(lp, 'w') as h:
                json.dump(left, h)
            with open(rp, 'w') as h:
                json.dump(right, h)
            done = subprocess.run(
                [sys.executable, TOOL, '--compare', lp, rp],
                capture_output=True, text=True)
        return done.returncode, done.stdout + done.stderr

    def _mutated(self, mutate):
        right = _fingerprint()
        mutate(right)
        return self._compare(_fingerprint(), right)

    # ---- the two comparable outcomes ------------------------------------

    def test_identical_runs_agree(self):
        """The baseline. If this fails nothing below means anything."""
        code, out = self._compare(_fingerprint(), _fingerprint())
        self.assertEqual(code, AGREE, out)

    def test_a_real_difference_is_reported(self):
        """A 1e-6 relative shift is far above round-off and must fail."""
        def mutate(fp):
            fp['frames'][1]['arrays']['solutiondg.12']['sorted_sum'] *= 1 + 1e-6
        code, out = self._mutated(mutate)
        self.assertEqual(code, DIFFER, out)

    def test_round_off_is_not_reported(self):
        """The other direction: a guard that fires on everything is useless."""
        def mutate(fp):
            fp['frames'][1]['arrays']['solutiondg.12']['sorted_sum'] *= 1 + 1e-13
        code, out = self._mutated(mutate)
        self.assertEqual(code, AGREE, out)

    # ---- the four that must REFUSE, not report a small number -----------

    def test_frame_count_mismatch_refuses(self):
        def mutate(fp):
            fp['frames'] = fp['frames'][:1]
        code, out = self._mutated(mutate)
        self.assertEqual(code, REFUSED, out)
        self.assertIn('frame count differs', out)

    def test_renamed_array_refuses(self):
        """Every other array still agrees; that must not read as agreement."""
        def mutate(fp):
            arrays = fp['frames'][0]['arrays']
            arrays['solutiondg.99'] = arrays.pop('solutiondg.12')
        code, out = self._mutated(mutate)
        self.assertEqual(code, REFUSED, out)
        self.assertIn('array names differ', out)

    def test_length_mismatch_refuses(self):
        def mutate(fp):
            fp['frames'][0]['arrays']['solutiondg.12']['n'] = 999
        code, out = self._mutated(mutate)
        self.assertEqual(code, REFUSED, out)
        self.assertIn('length differs', out)

    def test_changed_mesh_refuses(self):
        """Comparing fields across two different meshes is meaningless."""
        def mutate(fp):
            fp['frames'][0]['arrays']['Position']['max'] *= 1.01
        code, out = self._mutated(mutate)
        self.assertEqual(code, REFUSED, out)
        self.assertIn('MESH DIFFERS', out)

    def test_non_finite_values_refuse(self):
        """A run that produced NaN has not produced a comparable answer."""
        def mutate(fp):
            fp['frames'][1]['arrays']['solutiondg.0']['n_nonfinite'] = 7
        code, out = self._mutated(mutate)
        self.assertEqual(code, REFUSED, out)
        self.assertIn('non-finite', out)

    # ---- the reported difference must be the largest one ----------------

    def test_the_worst_array_is_the_one_named(self):
        def mutate(fp):
            fp['frames'][0]['arrays']['solutiondg.0']['sorted_sum'] *= 1 + 1e-8
            fp['frames'][1]['arrays']['solutiondg.12']['sorted_sum'] *= 1 + 1e-3
        code, out = self._mutated(mutate)
        self.assertEqual(code, DIFFER, out)
        worst = [l for l in out.split('\n') if l.startswith('worst over')]
        self.assertTrue(worst, out)
        self.assertIn('solutiondg.12', worst[0])

    def test_mesh_arrays_do_not_mask_a_field_difference(self):
        """Mesh arrays are listed apart so a field change is not buried."""
        def mutate(fp):
            fp['frames'][1]['arrays']['solutiondg.12']['sorted_sum'] *= 1 + 1e-4
        code, out = self._mutated(mutate)
        self.assertEqual(code, DIFFER, out)
        self.assertNotIn('MESH DIFFERS', out)

    # ---- a different partitioning is not a different mesh ---------------
    #
    # This is the distinction the tool got wrong twice. PETSc writes one
    # <Piece> per rank and a boundary vertex appears in every piece touching
    # it, so Position gets LONGER on more ranks (measured on the RMF antenna
    # deck: 18303 points on 1 rank, 18588 on 2, identical mesh) and
    # connectivity indexes into that renumbered list. Cell data is not
    # duplicated - 11989 either way - so the solution arrays remain
    # comparable, which is the whole basis of known-issues 15.

    def test_more_points_is_a_partitioning_note_not_a_refusal(self):
        def mutate(fp):
            for frame in fp['frames']:
                frame['arrays']['Position']['n'] = 18588
        code, out = self._mutated(mutate)
        self.assertEqual(code, AGREE, out)
        self.assertIn('PARTITIONED DIFFERENTLY', out)
        self.assertNotIn('MESH DIFFERS', out)

    def test_repartitioning_does_not_hide_a_field_difference(self):
        """The reason the note must not become a blanket excuse."""
        def mutate(fp):
            for frame in fp['frames']:
                frame['arrays']['Position']['n'] = 18588
            fp['frames'][1]['arrays']['solutiondg.12']['sorted_sum'] *= 1.9
        code, out = self._mutated(mutate)
        self.assertEqual(code, DIFFER, out)

    def test_rank_array_alone_is_not_a_mesh_difference(self):
        """'Rank' differs by construction across rank counts."""
        def mutate(fp):
            for frame in fp['frames']:
                frame['arrays']['Rank'] = _array(4, 0.0, 1.0, 2.0)
        left = _fingerprint()
        for frame in left['frames']:
            frame['arrays']['Rank'] = _array(4, 0.0, 0.0, 0.0)
        right = _fingerprint()
        mutate(right)
        code, out = self._compare(left, right)
        self.assertEqual(code, AGREE, out)
        self.assertNotIn('MESH DIFFERS', out)

    # ---- the banner is what tells you which build produced which side ----

    def test_both_banners_are_printed(self):
        left = _fingerprint()
        right = _fingerprint()
        right['banner'] = 'Apollo test, PETSc 3.25.5'
        code, out = self._compare(left, right)
        self.assertEqual(code, AGREE, out)
        self.assertIn('3.19.6', out)
        self.assertIn('3.25.5', out)


if __name__ == '__main__':
    unittest.main(verbosity=2)
