#!/usr/bin/env python3
"""Tests for scripts/run_growth.py.

The tool has one job: when a run dies, name the component that was growing.
So the tests are about the RANKING, and specifically about the baseline it is
measured from.

THE BUG THESE PIN. The first version ranked by last/first with frame 0 as the
denominator. Frame 0 is the initial condition, where the driven fields and both
cleaning potentials are exactly zero, so `first` was 0, the ratio was scored as
0, and every field that grew from nothing sorted LAST.

That is exactly backwards, and it was not theoretical: on a run that died at
1.73 of 20 RMF periods after ~12 hours with no checkpoint, the tool reported
e-rho at 1.5x as the fastest grower and buried phi (298x), E_y (186x) and
E_x (131x) at the bottom. The three it buried were the instability.
"""

import os
import struct
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
TOOL = os.path.join(ROOT, 'scripts', 'run_growth.py')

# Component index -> name, as the solver flattens 18 components x 3 nodes.
PHI, E_X, E_Y, E_Z, PSI, E_RHO = 16, 10, 11, 12, 17, 0


def _vtu(path, peaks, cells=4):
    """A minimal .vtu whose solutiondg.N arrays have the given peak values.

    Written with the reader's own conventions (little-endian, UInt64 header,
    appended raw) so this exercises the real path rather than a stub.
    """
    arrays = {}
    for comp, peak in peaks.items():
        # Node 0 of that component carries the peak; the others are quieter.
        arrays['solutiondg.%d' % comp] = [peak] + [peak * 0.1] * (cells - 1)
    arrays.setdefault('solutiondg.%d' % E_RHO, [1.0] * cells)

    names = sorted(arrays)
    blocks, offset, decls = [], 0, []
    for name in names:
        data = struct.pack('<%dd' % cells, *arrays[name])
        decls.append('<DataArray type="Float64" Name="%s" format="appended" '
                     'offset="%d"/>' % (name, offset))
        blocks.append(struct.pack('<Q', len(data)) + data)
        offset += 8 + len(data)

    payload = b''.join(blocks)
    head = (
        '<?xml version="1.0"?>\n'
        '<VTKFile type="UnstructuredGrid" byte_order="LittleEndian" '
        'header_type="UInt64">\n'
        '<UnstructuredGrid>\n'
        '<Piece NumberOfPoints="%d" NumberOfCells="%d">\n'
        '<CellData>\n%s\n</CellData>\n'
        '</Piece>\n</UnstructuredGrid>\n'
        '<AppendedData encoding="raw">\n_' % (cells, cells, '\n'.join(decls))
    ).encode()
    with open(path, 'wb') as h:
        h.write(head)
        h.write(payload)
        h.write(b'\n</AppendedData>\n</VTKFile>\n')


def _run(frames):
    done = subprocess.run([sys.executable, TOOL] + frames,
                          capture_output=True, text=True)
    return done.returncode, done.stdout + done.stderr


def _ranking(out):
    """The component names in the order the tool ranked them."""
    lines = out.split('\n')
    start = next(i for i, l in enumerate(lines) if l.startswith('growth, each'))
    names = []
    for line in lines[start + 1:]:
        if not line.startswith('  ') or line.strip().startswith('never'):
            break
        names.append(line.split()[0])
    return names


class TestRunGrowth(unittest.TestCase):

    def _write(self, tmp, series):
        """series: list of {component: peak} per frame."""
        paths = []
        for i, peaks in enumerate(series):
            p = os.path.join(tmp, 'run_%d.vtu' % i)
            _vtu(p, peaks)
            paths.append(p)
        return paths

    def test_a_field_growing_from_zero_ranks_first(self):
        """The bug. phi starts at 0, so frame 0 cannot be the baseline."""
        with tempfile.TemporaryDirectory() as tmp:
            frames = self._write(tmp, [
                {PHI: 0.0,     E_RHO: 1.0},     # initial condition
                {PHI: 1.0e-5,  E_RHO: 1.0},
                {PHI: 1.0e-3,  E_RHO: 1.1},
                {PHI: 3.0e-3,  E_RHO: 1.5},     # 300x from frame 1
            ])
            code, out = _run(frames)
        self.assertEqual(code, 0, out)
        self.assertEqual(_ranking(out)[0], 'phi', out)

    def test_the_slow_grower_does_not_win(self):
        """e-rho at 1.5x must not outrank phi at 300x, which is what happened."""
        with tempfile.TemporaryDirectory() as tmp:
            frames = self._write(tmp, [
                {PHI: 0.0,     E_RHO: 1.0},
                {PHI: 1.0e-5,  E_RHO: 1.0},
                {PHI: 3.0e-3,  E_RHO: 1.5},
            ])
            code, out = _run(frames)
        rank = _ranking(out)
        self.assertLess(rank.index('phi'), rank.index('e-rho'), out)

    def test_a_field_that_never_leaves_zero_is_reported_not_ranked(self):
        """Dividing by its baseline is meaningless; say so instead."""
        with tempfile.TemporaryDirectory() as tmp:
            frames = self._write(tmp, [
                {PHI: 0.0, E_Z: 0.0, E_RHO: 1.0},
                {PHI: 1.0e-5, E_Z: 0.0, E_RHO: 1.0},
                {PHI: 1.0e-3, E_Z: 0.0, E_RHO: 1.0},
            ])
            code, out = _run(frames)
        self.assertEqual(code, 0, out)
        self.assertIn('never non-zero', out)
        self.assertIn('E_z', out.split('never non-zero')[1])

    def test_a_flat_field_does_not_outrank_a_growing_one(self):
        with tempfile.TemporaryDirectory() as tmp:
            frames = self._write(tmp, [
                {PHI: 0.0,    PSI: 40.0, E_RHO: 1.0},
                {PHI: 1.0e-5, PSI: 40.0, E_RHO: 1.0},
                {PHI: 3.0e-3, PSI: 45.0, E_RHO: 1.0},
            ])
            code, out = _run(frames)
        rank = _ranking(out)
        self.assertEqual(rank[0], 'phi', out)
        self.assertLess(rank.index('phi'), rank.index('psi'), out)

    def test_an_argument_with_no_frame_number_is_refused_legibly(self):
        """A traceback is a poor answer when the run you were watching died.

        The frames are ordered by the N in `<run>_N.vtu`, and the first version
        read it with an unguarded rsplit: every argument that did not carry one
        - `--help` included - died with an IndexError naming neither the
        argument nor the reason.
        """
        with tempfile.TemporaryDirectory() as tmp:
            odd = os.path.join(tmp, 'formation.vtu')     # no _N
            _vtu(odd, {PHI: 1.0})
            code, out = _run([odd])
        self.assertNotEqual(code, 0, out)
        self.assertNotIn('Traceback', out)
        self.assertIn('formation.vtu', out)
        self.assertIn('cannot tell which frame', out)

    def test_help_prints_usage_rather_than_failing(self):
        code, out = _run(['--help'])
        self.assertEqual(code, 0, out)
        self.assertIn('usage: run_growth.py', out)
        self.assertNotIn('Traceback', out)

    def test_frames_are_ordered_by_number_not_by_name(self):
        """Frame 10 sorts between 1 and 2 by name, which measures the growth
        across the wrong pair of frames."""
        with tempfile.TemporaryDirectory() as tmp:
            series = [{PHI: 0.0}] + [{PHI: 1.0e-5 * (1.7 ** i)} for i in range(11)]
            frames = self._write(tmp, series)
            self.assertTrue(any(f.endswith('_10.vtu') for f in frames), frames)
            code, out = _run(list(reversed(frames)))     # argv order must not matter
        self.assertEqual(code, 0, out)
        rows = [l for l in out.split('\n') if l[:1].isdigit()]
        phis = [float(r.split()[2 + PHI]) for r in rows]
        self.assertEqual(phis, sorted(phis), out)

    def test_non_finite_values_are_counted(self):
        """A frame holding NaN must say so rather than be summarised."""
        with tempfile.TemporaryDirectory() as tmp:
            frames = self._write(tmp, [
                {PHI: 1.0e-5, E_RHO: 1.0},
                {PHI: float('nan'), E_RHO: 1.0},
            ])
            code, out = _run(frames)
        self.assertEqual(code, 0, out)
        # The nonfin column is the second field of each frame row.
        rows = [l for l in out.split('\n') if l[:1].isdigit()]
        self.assertTrue(rows, out)
        self.assertGreater(int(rows[-1].split()[1]), 0, out)


if __name__ == '__main__':
    unittest.main(verbosity=2)
