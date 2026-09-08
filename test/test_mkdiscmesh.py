#!/usr/bin/env python3
"""The disc mesh generator, scripts/mkdiscmesh.py.

The repository ships .msh files with no .geo alongside them, so before this
generator existed there was no way to change the domain of a case - which is
what Phase 1 of docs/rmf-frc-model-assessment.md needed, to put the RMF antenna
somewhere other than the plasma edge.

Two kinds of thing are checked here. One is that the output is a mesh:
counter-clockwise triangles, no inversions, an area that converges to pi R^2, a
closed boundary. The other is that it is a mesh Apollo's reader will accept -
src/lib/wxpreadgmshgrid.h is strict in ways the format is not, and a file that
violates any of them is rejected at load time with a message that does not say
which rule was broken:

  * the header line must be exactly "2.2 0 8";
  * node numbers must be 1..N in order and element numbers 1..M in order, each
    checked against the loop index;
  * facets must all precede cells, because the reader infers the topological
    dimension by scanning for the highest one and then indexes cells as
    (element index - number of facets).

The generator's own report already counts inverted triangles and prints the
area; these tests exist so that a regression fails the build rather than a run.
"""

import math
import os
import subprocess
import sys
import tempfile
import unittest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GENERATOR = os.path.join(REPO, 'scripts', 'mkdiscmesh.py')


def generate(tmp, name, **kwargs):
    """Run the generator and return the path it wrote."""
    path = os.path.join(tmp, name)
    argv = [sys.executable, GENERATOR, '-o', path]
    for key, value in kwargs.items():
        flag = '--' + key.replace('_', '-')
        if value is True:
            argv.append(flag)
        else:
            argv += [flag, repr(value)]
    done = subprocess.run(argv, capture_output=True, text=True, timeout=300)
    if done.returncode != 0:
        raise AssertionError(f'generator failed:\n{done.stdout}{done.stderr}')
    return path


class Mesh:
    """The parts of a gmsh 2.2 ASCII file these tests care about."""

    def __init__(self, path):
        with open(path) as fh:
            self.lines = fh.read().split('\n')
        self.header = self.lines[1]

        i = self.lines.index('$Nodes')
        count = int(self.lines[i + 1])
        self.node_ids, self.nodes = [], []
        for k in range(i + 2, i + 2 + count):
            parts = self.lines[k].split()
            self.node_ids.append(int(parts[0]))
            self.nodes.append((float(parts[1]), float(parts[2]), float(parts[3])))

        j = self.lines.index('$Elements')
        count = int(self.lines[j + 1])
        self.elem_ids, self.kinds, self.tags, self.conn = [], [], [], []
        for k in range(j + 2, j + 2 + count):
            parts = self.lines[k].split()
            self.elem_ids.append(int(parts[0]))
            self.kinds.append(int(parts[1]))
            ntags = int(parts[2])
            self.tags.append([int(t) for t in parts[3:3 + ntags]])
            self.conn.append([int(v) for v in parts[3 + ntags:]])

    @property
    def triangles(self):
        return [c for k, c in zip(self.kinds, self.conn) if k == 2]

    @property
    def edges(self):
        return [c for k, c in zip(self.kinds, self.conn) if k == 1]

    def signed_area(self, tri):
        (ax, ay, _), (bx, by, _), (cx, cy, _) = (self.nodes[v - 1] for v in tri)
        return 0.5 * ((bx - ax) * (cy - ay) - (cx - ax) * (by - ay))

    def radii(self):
        return [math.hypot(x, y) for (x, y, _) in self.nodes]


class TestMeshIsValid(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix='apollo-mkdiscmesh-')
        cls.path = generate(cls.tmp, 'graded.msh', outer=0.05, ring=0.03,
                            h_inner=1.5e-3, h_outer=3.0e-3)
        cls.mesh = Mesh(cls.path)

    @classmethod
    def tearDownClass(cls):
        import shutil
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_no_inverted_triangles(self):
        """Every triangle counter-clockwise, as the shipped meshes are."""
        areas = [self.mesh.signed_area(t) for t in self.mesh.triangles]
        bad = [a for a in areas if a <= 0.0]
        self.assertEqual(bad, [], f'{len(bad)} triangles are inverted or degenerate')

    def test_area_matches_the_disc(self):
        """The triangles tile the disc, up to the polygonal boundary."""
        area = sum(self.mesh.signed_area(t) for t in self.mesh.triangles)
        exact = math.pi * 0.05 ** 2
        self.assertLess(abs(area - exact) / exact, 1e-3,
                        f'area {area:.8e} against pi R^2 = {exact:.8e}')

    def test_area_converges_under_refinement(self):
        """Halving the spacing must cut the area deficit, not leave it."""
        deficits = []
        for h in (4.0e-3, 2.0e-3, 1.0e-3):
            mesh = Mesh(generate(self.tmp, f'conv_{h}.msh', outer=0.05, h_inner=h))
            area = sum(mesh.signed_area(t) for t in mesh.triangles)
            deficits.append((math.pi * 0.05 ** 2 - area) / (math.pi * 0.05 ** 2))
        self.assertTrue(
            deficits[0] > deficits[1] > deficits[2] > 0,
            f'area deficit did not decrease under refinement: {deficits}')

    def test_boundary_is_a_closed_loop_on_the_outer_radius(self):
        """Every boundary node on r = outer, and in exactly two edges."""
        edges = self.mesh.edges
        self.assertTrue(edges, 'no boundary edges were written')
        seen = {}
        for a, b in edges:
            for v in (a, b):
                seen[v] = seen.get(v, 0) + 1
        self.assertEqual(sorted(set(seen.values())), [2],
                         'the boundary is not a simple closed loop')
        for v in seen:
            x, y, _ = self.mesh.nodes[v - 1]
            self.assertAlmostEqual(math.hypot(x, y), 0.05, places=12)

    def test_conforming_ring(self):
        """A node ring lands exactly on --ring, so a jump there is on a face."""
        self.assertTrue(any(abs(r - 0.03) < 1e-12 for r in self.mesh.radii()),
                        'no node ring landed on r = 0.03')

    def test_grading_actually_coarsens(self):
        """Cells outside the ring must be bigger, in both directions.

        An earlier version forced the per-ring node count to be non-decreasing
        outwards, which left the outer rings carrying the inner spacing
        azimuthally: radially graded, azimuthally not.
        """
        inner, outer = [], []
        for tri in self.mesh.triangles:
            pts = [self.mesh.nodes[v - 1] for v in tri]
            rc = sum(math.hypot(p[0], p[1]) for p in pts) / 3.0
            lengths = [math.hypot(pts[i][0] - pts[(i + 1) % 3][0],
                                  pts[i][1] - pts[(i + 1) % 3][1]) for i in range(3)]
            (inner if rc <= 0.03 else outer).extend(lengths)
        mean_in = sum(inner) / len(inner)
        mean_out = sum(outer) / len(outer)
        self.assertGreater(mean_out, 1.5 * mean_in,
                           f'mean edge {mean_out:.3e} outside against '
                           f'{mean_in:.3e} inside; h-outer was twice h-inner')

    def test_triangle_quality(self):
        """Nothing close to degenerate: the time step is set by the worst cell."""
        worst = 1.0
        for tri in self.mesh.triangles:
            pts = [self.mesh.nodes[v - 1] for v in tri]
            lengths = [math.hypot(pts[i][0] - pts[(i + 1) % 3][0],
                                  pts[i][1] - pts[(i + 1) % 3][1]) for i in range(3)]
            area = abs(self.mesh.signed_area(tri))
            worst = min(worst, 4.0 * math.sqrt(3.0) * area / sum(l * l for l in lengths))
        self.assertGreater(worst, 0.4, f'worst triangle quality is {worst:.4f}')


class TestReaderConstraints(unittest.TestCase):
    """Requirements src/lib/wxpreadgmshgrid.h imposes beyond the format."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix='apollo-mkdiscmesh-fmt-')
        cls.mesh = Mesh(generate(cls.tmp, 'fmt.msh', outer=0.04, ring=0.02,
                                 h_inner=2.0e-3, h_outer=4.0e-3))

    @classmethod
    def tearDownClass(cls):
        import shutil
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_header(self):
        """The reader fscanfs "2.2 %d %d" and rejects any other file type."""
        self.assertEqual(self.mesh.header, '2.2 0 8')

    def test_node_numbering_is_one_based_and_dense(self):
        self.assertEqual(self.mesh.node_ids, list(range(1, len(self.mesh.node_ids) + 1)))

    def test_element_numbering_is_one_based_and_dense(self):
        self.assertEqual(self.mesh.elem_ids, list(range(1, len(self.mesh.elem_ids) + 1)))

    def test_facets_precede_cells(self):
        """Lines first, then triangles: the reader indexes cells off the count."""
        kinds = self.mesh.kinds
        first_triangle = kinds.index(2)
        self.assertNotIn(1, kinds[first_triangle:],
                         'a boundary line appears after the first triangle')

    def test_every_element_carries_two_tags(self):
        for tags in self.mesh.tags:
            self.assertEqual(len(tags), 2)

    def test_boundary_tag_selects_the_first_boundary_condition(self):
        """Physical tag 1 becomes Face Sets value 1, and
        WxNodalDG2dMethod::applyBc looks up boundaryConditions[tag-1]."""
        line_tags = {t[0] for k, t in zip(self.mesh.kinds, self.mesh.tags) if k == 1}
        self.assertEqual(line_tags, {1})

    def test_connectivity_is_in_range(self):
        n = len(self.mesh.node_ids)
        for element in self.mesh.conn:
            for v in element:
                self.assertTrue(1 <= v <= n, f'vertex {v} outside 1..{n}')

    def test_region_tag_marks_the_two_regions(self):
        """--region-tag is for other tools; Apollo's reader drops cell tags."""
        mesh = Mesh(generate(self.tmp, 'tagged.msh', outer=0.04, ring=0.02,
                             h_inner=2.0e-3, h_outer=4.0e-3, region_tag=True))
        tags = {t[0] for k, t in zip(mesh.kinds, mesh.tags) if k == 2}
        self.assertEqual(tags, {1, 2})


if __name__ == '__main__':
    unittest.main()
