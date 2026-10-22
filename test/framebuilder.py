"""Write a synthetic P1 .vtu for the ring/anatomy diagnostics' tests.

Shared by test_ring_spectrum.py and test_winding_anatomy.py so the appended-raw
packing lives in one place. A frame is a list of triangles; each triangle is
three (x, y, {component: value}) nodes. Components not named default to zero.
"""

import os
import struct

NEQ = 18


def write_vtu(path, triangles):
    """triangles: list of [ (x,y,{comp:val}), (x,y,{...}), (x,y,{...}) ]."""
    pts, conn = [], []
    node_vals = {c: [] for c in range(NEQ)}
    for tri in triangles:
        base = len(pts)
        for (x, y, vals) in tri:
            pts.append((x, y, 0.0))
        conn += [base, base + 1, base + 2]
        for node in range(3):
            _, _, vals = tri[node]
            for c in range(NEQ):
                node_vals[c].append(float(vals.get(c, 0.0)))

    n = len(triangles)
    arrays = {}
    per = {c: [[0.0] * n for _ in range(3)] for c in range(NEQ)}
    for cell, tri in enumerate(triangles):
        for node in range(3):
            _, _, vals = tri[node]
            for c in range(NEQ):
                per[c][node][cell] = float(vals.get(c, 0.0))
    for c in range(NEQ):
        for node in range(3):
            arrays['solutiondg.%d' % (node * NEQ + c)] = ('f8', per[c][node])

    arrays['Position'] = ('f8', [v for p in pts for v in p])
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


def ring_triangles(radius, ncells, comp_fn):
    """ncells thin triangles ON a ring at `radius`, comp_fn(theta,i) -> {comp:val}.

    All three nodes sit at exactly `radius` (angles theta, theta+/-delta), so the
    3*ncells vertices form ONE concentric ring at that radius - which is what the
    diagnostic's ring finder expects. The same value dict is put on all three
    nodes, so vertex averaging returns comp_fn(theta) around the ring.
    """
    import math
    tris = []
    delta = 1.0e-3
    for i in range(ncells):
        th = 2.0 * math.pi * i / ncells
        vals = comp_fn(th, i)
        node = lambda a: (radius * math.cos(a), radius * math.sin(a), vals)
        tris.append([node(th), node(th + delta), node(th - delta)])
    return tris


def deck_text(tend=1.0e-6, out=4, path=None):
    txt = '''# -*- python -*-
import math
PI = math.pi
MU0 = PI*4.0e-7
Q = 1.6e-19
MP = 1.67e-27
MI = MP
ME = MP/1836
GAMMA = 1.66666666667
omega = 7.95e5
Te = 30.
LIGHT = 3.0e6
ETA = 0.5e-5
n_dens = 1.e20
TEND = %r
OUT = %r
RAD_PLASMA = 0.030
COIL_R = 0.036
COIL_W = 0.005
WALL_RADIUS = 0.050
<warpx>
''' % (tend, out)
    if path:
        with open(path, 'w') as h:
            h.write(txt)
    return txt
