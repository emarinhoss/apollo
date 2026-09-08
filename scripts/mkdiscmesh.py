#!/usr/bin/env python3
"""Generate a graded triangular disc mesh in gmsh 2.2 ASCII format.

Some examples ship a .geo beside their .msh (euler/backwardFacingStep does); the
multifluid rmf_frc discs do not, so there was no way to reproduce or vary one.
This is that generator.

It exists for the RMF-FRC problem, which needs a disc whose interior is meshed
at the plasma resolution out to the plasma radius and coarsely beyond it, with
a node ring lying exactly on the plasma edge so that a discontinuous initial
condition there is mesh-conforming:

    python3 scripts/mkdiscmesh.py --outer 0.05 --ring 0.03 \
            --h-inner 9.3e-4 --h-outer 3.0e-3 -o disc.msh

The construction is polar: concentric rings of nodes at radii chosen by
integrating a spacing function h(r), each ring carrying about 2*pi*r/h(r)
nodes, with consecutive rings stitched by an advancing front in angle. That
gives near-equilateral triangles, a spacing that follows h(r), and - unlike a
Delaunay triangulation of scattered points - exact control over where nodes
land, which is what makes the conforming ring possible.

The output is what gmsh itself writes for a 2.2 ASCII mesh, and deliberately
nothing looser:

  - the header line is exactly "2.2 0 8" (ASCII, sizeof(double));
  - node numbers run 1..N in order, element numbers 1..M in order;
  - facets (2-node lines, type 1) all precede cells (3-node triangles, type 2);
  - every element carries two tags, (physical, geometrical).

Apollo loads meshes with PETSc's DMPlexCreateGmsh (src/solvers/apsolver.cc),
which is more permissive than that - it reads 4.1 and binary too. There is also
a hand-written reader in src/lib/wxpreadgmshgrid.h that requires exactly the
above, but nothing calls it; it is one of the unreachable files in src/lib that
docs/known-issues.md records. Writing to the stricter shape costs nothing and
keeps the meshes loadable by either.

Physical tags on boundary lines reach the DMPlex as the "Face Sets" label, which
wxNodalDGgeometry2D reads and WxNodalDG2dMethod::applyBc turns into
boundaryConditions[tag-1] - so tag 1 selects the first entry of the deck's
boundaryConditions list, tag 2 the second. (Confirmed by running: the antenna
example's single tag-1 boundary picks up its one BC.) Cell tags do not survive
into anything Apollo reads, so a two-region mesh cannot be told apart by cell
tag - regions have to be distinguished geometrically, by radius, in the initial
condition and in the source terms. The --region-tag option still writes a tag
per region, for other tools.
"""

import argparse
import math
import sys


def ring_radii(r_outer, h_inner, h_outer, r_ring=None, grade_len=None):
    """Radii of the node rings, from the first ring out to r_outer.

    The spacing function is h_inner inside r_ring and relaxes geometrically to
    h_outer over `grade_len` beyond it. Radii are produced by stepping
    r -> r + h(r), then rescaled so that a ring lands exactly on r_ring and
    exactly on r_outer: the rescaling is a monotone piecewise-linear map, so it
    perturbs the spacing by at most one step's worth and preserves ordering.
    """
    if r_ring is None:
        r_ring = r_outer
    if grade_len is None:
        grade_len = max(h_outer, 0.25 * (r_outer - r_ring))

    def h(r):
        if r <= r_ring or grade_len <= 0.0:
            return h_inner
        # Geometric relaxation towards h_outer with an e-folding of grade_len:
        # smooth, always between the two spacings, and monotone in r.
        w = 1.0 - math.exp(-(r - r_ring) / grade_len)
        return h_inner + (h_outer - h_inner) * w

    radii = []
    r = 0.0
    # Guard against a spacing so small the loop would never terminate.
    if h_inner <= 0.0 or h_outer <= 0.0:
        raise ValueError('spacings must be positive')
    while True:
        r += h(r)
        if r >= r_outer - 0.5 * h(r):
            break
        radii.append(r)
        if len(radii) > 100000:
            raise ValueError('spacing too fine for the requested radius')
    radii.append(r_outer)

    if r_ring < r_outer:
        # Snap the nearest ring onto r_ring and stretch each side to match, so
        # the plasma edge is a mesh line rather than falling mid-element.
        k = min(range(len(radii)), key=lambda i: abs(radii[i] - r_ring))
        old = radii[k]
        if old > 0.0 and k + 1 < len(radii):
            inner = [rr * (r_ring / old) for rr in radii[:k + 1]]
            span_old = r_outer - old
            span_new = r_outer - r_ring
            outer = [r_outer - (r_outer - rr) * (span_new / span_old)
                     for rr in radii[k + 1:]]
            radii = inner + outer
    return radii


def ring_counts(radii, h_inner, h_outer, r_ring, grade_len=None):
    """Number of nodes on each ring: about 2*pi*r/h(r), so that the azimuthal
    spacing tracks the radial one.

    Counts are allowed to DECREASE outwards. An earlier version forced them to
    be non-decreasing, on the theory that the stitching needs it; it does not -
    the advancing front below always steps whichever ring is behind in angle,
    so either ordering works - and the constraint quietly defeated coarsening,
    since a coarse outer ring inherited the fine inner ring's node count and
    ended up with the plasma's azimuthal resolution out to the wall.

    The same smooth h(r) as ring_radii is used, so consecutive rings differ by
    a little rather than by a step, which is what keeps the transition
    triangles well shaped.
    """
    if r_ring is None:
        r_ring = radii[-1]
    if grade_len is None:
        grade_len = max(h_outer, 0.25 * (radii[-1] - r_ring))

    def h(r):
        if r <= r_ring or grade_len <= 0.0:
            return h_inner
        w = 1.0 - math.exp(-(r - r_ring) / grade_len)
        return h_inner + (h_outer - h_inner) * w

    return [max(6, int(round(2.0 * math.pi * r / h(r)))) for r in radii]


def build(radii, counts):
    """Nodes and triangles for the rings, plus the outer boundary edges.

    Returns (nodes, triangles, boundary_edges) with 1-based node ids, every
    triangle counter-clockwise (the shipped meshes are, and
    wxNodalDGgeometry2D::elementArea takes the signed area).
    """
    nodes = [(0.0, 0.0)]          # the centre, node 1
    ring_ids = []
    for r, n in zip(radii, counts):
        ids = []
        for k in range(n):
            th = 2.0 * math.pi * k / n
            nodes.append((r * math.cos(th), r * math.sin(th)))
            ids.append(len(nodes))
        ring_ids.append(ids)

    tris = []
    # Centre fan: the innermost ring closes onto node 1.
    inner = ring_ids[0]
    for k in range(len(inner)):
        tris.append((1, inner[k], inner[(k + 1) % len(inner)]))

    # Stitch consecutive rings by advancing whichever front is behind in angle.
    for a_ids, b_ids in zip(ring_ids, ring_ids[1:]):
        na, nb = len(a_ids), len(b_ids)
        ia = ib = 0
        while ia < na or ib < nb:
            # Angle of the next candidate node on each ring, as a fraction of a
            # full turn; +1.0 stands for "this front is finished".
            fa = (ia + 1) / na if ia < na else 2.0
            fb = (ib + 1) / nb if ib < nb else 2.0
            if fa <= fb:
                # Advance the inner ring. Vertex order (a_i, b_j, a_i+1), not
                # the inner pair first: going out to the outer ring before
                # closing back is what makes the triangle counter-clockwise.
                tris.append((a_ids[ia % na], b_ids[ib % nb], a_ids[(ia + 1) % na]))
                ia += 1
            else:
                # Advance the outer ring: (a_i, b_j, b_j+1), likewise CCW.
                tris.append((a_ids[ia % na], b_ids[ib % nb], b_ids[(ib + 1) % nb]))
                ib += 1

    outer = ring_ids[-1]
    edges = [(outer[k], outer[(k + 1) % len(outer)]) for k in range(len(outer))]
    return nodes, tris, edges


def signed_area(nodes, tri):
    (ax, ay), (bx, by), (cx, cy) = (nodes[i - 1] for i in tri)
    return 0.5 * ((bx - ax) * (cy - ay) - (cx - ax) * (by - ay))


def write_msh(path, nodes, tris, edges, r_ring, region_tag):
    """Write gmsh 2.2 ASCII. Lines first, then triangles - the reader requires it."""
    with open(path, 'w') as fh:
        fh.write('$MeshFormat\n2.2 0 8\n$EndMeshFormat\n')
        fh.write('$Nodes\n%d\n' % len(nodes))
        for i, (x, y) in enumerate(nodes, start=1):
            fh.write('%d %.17g %.17g 0\n' % (i, x, y))
        fh.write('$EndNodes\n')
        fh.write('$Elements\n%d\n' % (len(edges) + len(tris)))
        eid = 0
        for a, b in edges:
            eid += 1
            # Physical tag 1 -> boundaryConditions[0] in the deck.
            fh.write('%d 1 2 1 1 %d %d\n' % (eid, a, b))
        for tri in tris:
            eid += 1
            tag = 1
            if region_tag:
                r = max(math.hypot(*nodes[i - 1]) for i in tri)
                tag = 1 if r <= r_ring + 1e-12 else 2
            fh.write('%d 2 2 %d %d %d %d %d\n' % (eid, tag, tag, tri[0], tri[1], tri[2]))
        fh.write('$EndElements\n')


def report(nodes, tris, radii, r_ring, r_outer):
    """Quality numbers, printed so a bad mesh is obvious before it is used."""
    areas = [signed_area(nodes, t) for t in tris]
    inverted = sum(1 for a in areas if a <= 0.0)
    total = sum(areas)
    quals, edges_in, edges_out = [], [], []
    for t, a in zip(tris, areas):
        p = [nodes[i - 1] for i in t]
        ls = [math.hypot(p[i][0] - p[(i + 1) % 3][0], p[i][1] - p[(i + 1) % 3][1])
              for i in range(3)]
        # Normalised shape quality: 1 for an equilateral triangle, 0 degenerate.
        quals.append(4.0 * math.sqrt(3.0) * abs(a) / sum(l * l for l in ls))
        rc = sum(math.hypot(*q) for q in p) / 3.0
        (edges_in if rc <= r_ring else edges_out).extend(ls)

    def stats(v):
        return (min(v), sum(v) / len(v), max(v)) if v else (0.0, 0.0, 0.0)

    print('nodes            %d' % len(nodes))
    print('triangles        %d  (inverted: %d)' % (len(tris), inverted))
    print('rings            %d' % len(radii))
    print('area             %.8e   exact pi*R^2 = %.8e   (%.4f%% low)'
          % (total, math.pi * r_outer ** 2,
             100.0 * (math.pi * r_outer ** 2 - total) / (math.pi * r_outer ** 2)))
    print('edge len  r<=%.4g   min %.3e  mean %.3e  max %.3e' % ((r_ring,) + stats(edges_in)))
    print('edge len  r> %.4g   min %.3e  mean %.3e  max %.3e' % ((r_ring,) + stats(edges_out)))
    print('quality (1=equilateral)  min %.4f  mean %.4f' % (min(quals), sum(quals) / len(quals)))
    return inverted


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-o', '--output', required=True, help='output .msh path')
    ap.add_argument('--outer', type=float, required=True, help='domain radius b [m]')
    ap.add_argument('--ring', type=float, default=None,
                    help='radius that must carry a node ring, and inside which '
                         '--h-inner applies [m]; default: --outer')
    ap.add_argument('--h-inner', type=float, required=True,
                    help='target edge length inside --ring [m]')
    ap.add_argument('--h-outer', type=float, default=None,
                    help='target edge length outside --ring [m]; default: --h-inner')
    ap.add_argument('--grade-len', type=float, default=None,
                    help='e-folding length over which the spacing relaxes from '
                         'h-inner to h-outer [m]')
    ap.add_argument('--region-tag', action='store_true',
                    help='tag triangles 1 inside --ring and 2 outside (the mesh '
                         'reader ignores cell tags; other tools may not)')
    args = ap.parse_args(argv)

    r_ring = args.outer if args.ring is None else args.ring
    h_outer = args.h_inner if args.h_outer is None else args.h_outer
    if not 0.0 < r_ring <= args.outer:
        ap.error('--ring must be positive and no larger than --outer')

    radii = ring_radii(args.outer, args.h_inner, h_outer, r_ring, args.grade_len)
    counts = ring_counts(radii, args.h_inner, h_outer, r_ring, args.grade_len)
    nodes, tris, edges = build(radii, counts)

    on_ring = [r for r in radii if abs(r - r_ring) < 1e-12]
    if r_ring < args.outer and not on_ring:
        print('warning: no node ring landed on r = %g' % r_ring, file=sys.stderr)

    write_msh(args.output, nodes, tris, edges, r_ring, args.region_tag)
    print('wrote %s' % args.output)
    inverted = report(nodes, tris, radii, r_ring, args.outer)
    return 1 if inverted else 0


if __name__ == '__main__':
    sys.exit(main())
