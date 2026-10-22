#!/usr/bin/env python3
"""Tests for scripts/where_peak.py.

The tool exists to separate two failures that produce the SAME growth table:
growth at the antenna winding (the drive loading up) from growth at the
conducting wall (an unresolved electron sheath, which kills the run). So the
tests are about the LABEL it puts on a radius, and about the two ways it could
put a confident label on nothing:

  - inventing region boundaries when no deck was given;
  - ordering frames by name, so that frame 10 lands between 1 and 2 and the
    "last frames" it prints are not the last frames.
"""

import os
import struct
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
TOOL = os.path.join(ROOT, 'scripts', 'where_peak.py')

WALL, COIL_R, COIL_W, RAD_PLASMA = 0.050, 0.036, 0.005, 0.030

# deck_parameters insists on a real rmf_frc deck, so give it one.
DECK = '''# -*- python -*-
import math
PI = math.pi
MU0 = PI*4.0e-7
Q = 1.6e-19
MP = 1.67e-27
MI = MP
ME = MP/1836
omega = 7.95e5
Te = 30.
LIGHT = 3.0e6
ETA = 0.5e-5
n_dens = 1.e20
TEND = 1.0e-6
OUT = 4
RAD_PLASMA = %r
COIL_R = %r
COIL_W = %r
WALL_RADIUS = %r
<warpx>
''' % (RAD_PLASMA, COIL_R, COIL_W, WALL)


def _vtu(path, cells, peaks, peak_node=0):
    """A .vtu whose cells sit at the given radii, with a chosen peak cell.

    cells : list of radii, one per cell (each cell is a tiny triangle there)
    peaks : {component: (cell index, value)} - that cell gets the value, the
            rest get something smaller.
    peak_node : which of the three DG nodes carries the peak.
    """
    import math as _m
    pts, conn = [], []
    for i, r in enumerate(cells):
        th = 0.3 * i
        cx, cy = r * _m.cos(th), r * _m.sin(th)
        d = 1e-5
        base = len(pts)
        pts += [(cx - d, cy - d, 0.0), (cx + d, cy - d, 0.0), (cx, cy + d, 0.0)]
        conn += [base, base + 1, base + 2]

    arrays = {}
    n = len(cells)
    for comp, (cell, value) in peaks.items():
        for node in range(3):
            vals = [value * 1e-3] * n
            if node == peak_node:
                vals[cell] = value
            arrays['solutiondg.%d' % (node * 18 + comp)] = ('f8', vals)
    arrays['Position'] = ('f8', [c for p in pts for c in p])
    arrays['connectivity'] = ('i8', conn)
    arrays['offsets'] = ('i8', [3 * (i + 1) for i in range(n)])
    arrays['types'] = ('u1', [5] * n)

    order = sorted(arrays)
    blocks, offset, decls = [], 0, {}
    for name in order:
        kind, vals = arrays[name]
        fmt = {'f8': '<%dd', 'i8': '<%dq', 'u1': '<%dB'}[kind]
        data = struct.pack(fmt % len(vals), *vals)
        vtype = {'f8': 'Float64', 'i8': 'Int64', 'u1': 'UInt8'}[kind]
        ncomp = 3 if name == 'Position' else 1
        decls[name] = ('<DataArray type="%s" Name="%s" NumberOfComponents="%d" '
                       'format="appended" offset="%d"/>' % (vtype, name, ncomp, offset))
        blocks.append(struct.pack('<Q', len(data)) + data)
        offset += 8 + len(data)

    sols = '\n'.join(decls[k] for k in order if k.startswith('solutiondg.'))
    head = (
        '<?xml version="1.0"?>\n'
        '<VTKFile type="UnstructuredGrid" byte_order="LittleEndian" '
        'header_type="UInt64">\n<UnstructuredGrid>\n'
        '<Piece NumberOfPoints="%d" NumberOfCells="%d">\n'
        '<Points>\n%s\n</Points>\n'
        '<Cells>\n%s\n%s\n%s\n</Cells>\n'
        '<CellData>\n%s\n</CellData>\n'
        '</Piece>\n</UnstructuredGrid>\n<AppendedData encoding="raw">\n_'
        % (len(pts), n, decls['Position'], decls['connectivity'],
           decls['offsets'], decls['types'], sols)).encode()
    with open(path, 'wb') as h:
        h.write(head)
        h.write(b''.join(blocks))
        h.write(b'\n</AppendedData>\n</VTKFile>\n')


def _run(args):
    done = subprocess.run([sys.executable, TOOL] + args,
                          capture_output=True, text=True)
    return done.returncode, done.stdout + done.stderr


class TestWherePeak(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix='apollo-wherepeak-')
        self.deck = os.path.join(self.tmp, 'formation.pin')
        with open(self.deck, 'w') as h:
            h.write(DECK)
        # One cell in each region the tool knows how to name.
        self.cells = [0.010, 0.0355, 0.0440, 0.0490]   # column, winding, gap, WALL
        self.REGION = {0: 'column', 1: 'winding', 2: 'gap', 3: 'WALL'}

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _frames(self, peak_cell, n=3, comp=10, peak_node=0):
        paths = []
        for i in range(n):
            p = os.path.join(self.tmp, 'run_%d.vtu' % i)
            _vtu(p, self.cells, {comp: (peak_cell, 1.0e3)}, peak_node=peak_node)
            paths.append(p)
        return paths

    # -- the labels, which are the whole point -----------------------------

    def test_each_region_is_named_from_the_deck(self):
        for cell, want in self.REGION.items():
            frames = self._frames(cell)
            code, out = _run([self.deck] + frames + ['--components', 'E_x'])
            self.assertEqual(code, 0, out)
            self.assertIn(want, out, 'cell at r=%.4f should be %r:\n%s'
                          % (self.cells[cell], want, out))

    def test_the_wall_ring_is_counted(self):
        frames = self._frames(3)                       # the wall cell
        code, out = _run([self.deck] + frames + ['--components', 'E_x'])
        self.assertEqual(code, 0, out)
        self.assertIn('3 of   3', out.replace('  3 of', ' 3 of'), out)

    def test_an_interior_peak_is_not_counted_as_wall(self):
        frames = self._frames(1)                       # the winding cell
        code, out = _run([self.deck] + frames + ['--components', 'E_x'])
        self.assertEqual(code, 0, out)
        self.assertIn('winding', out)
        self.assertNotIn('WALL', out.split('wall ring')[0], out)

    # -- the two ways to be confidently wrong ------------------------------

    def test_without_a_deck_no_region_is_named(self):
        """The boundaries are deck values; another deck's geometry differs."""
        frames = self._frames(3)
        code, out = _run(frames + ['--components', 'E_x'])
        self.assertEqual(code, 0, out)
        self.assertIn('regions cannot be named', out)
        for word in ('WALL', 'winding', 'column', 'gap'):
            self.assertNotIn(' %s ' % word, out, 'invented a label: %s' % out)
        self.assertIn('r=0.0490', out)                 # the radius is still a fact

    def test_frames_are_ordered_by_number_not_by_name(self):
        paths = []
        for i in list(range(12)):
            p = os.path.join(self.tmp, 'run_%d.vtu' % i)
            # the peak moves to the wall only at frame 10 and 11
            _vtu(p, self.cells, {10: (3 if i >= 10 else 1, 1.0e3 * (i + 1))})
            paths.append(p)
        code, out = _run([self.deck] + list(reversed(paths))
                         + ['--components', 'E_x', '--frames', '2'])
        self.assertEqual(code, 0, out)
        body = out.split('wall ring')[0]
        self.assertIn('    10 ', body, out)
        self.assertIn('    11 ', body, out)
        self.assertNotIn('     1 ', body, 'frame 1 printed as one of the last:\n%s' % out)

    # -- refusals ----------------------------------------------------------

    def test_a_frame_with_no_number_is_refused_legibly(self):
        p = os.path.join(self.tmp, 'formation.vtu')     # no _N
        _vtu(p, self.cells, {10: (3, 1.0e3)})
        code, out = _run([self.deck, p])
        self.assertNotEqual(code, 0, out)
        self.assertNotIn('Traceback', out)
        self.assertIn('cannot tell which frame', out)

    def test_no_frames_is_refused(self):
        code, out = _run([self.deck])
        self.assertNotEqual(code, 0, out)
        self.assertIn('no .vtu frames', out)

    def test_help_prints_usage(self):
        code, out = _run(['--help'])
        self.assertEqual(code, 0, out)
        self.assertIn('usage: where_peak.py', out)
        self.assertNotIn('Traceback', out)

    def test_an_unknown_component_is_refused(self):
        frames = self._frames(3)
        code, out = _run([self.deck] + frames + ['--components', 'E_q'])
        self.assertNotEqual(code, 0, out)
        self.assertIn('unknown component', out)

    # -- the DG storage layout --------------------------------------------

    def test_a_peak_on_node_2_is_found(self):
        """solutiondg.N is node*18+component: a reader that looks only at the
        first 18 arrays misses any peak living on node 1 or 2."""
        frames = self._frames(3, peak_node=2)
        code, out = _run([self.deck] + frames + ['--components', 'E_x'])
        self.assertEqual(code, 0, out)
        self.assertIn('1.000e+03', out, out)
        self.assertIn('WALL', out)


if __name__ == '__main__':
    unittest.main(verbosity=2)
