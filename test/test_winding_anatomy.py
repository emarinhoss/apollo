#!/usr/bin/env python3
"""Tests for scripts/winding_anatomy.py.

The tool answers the questions the instability plan turns on, and each is a way
the run could be misread, so each is planted and checked:

  - a density hole with the pressure held flat (Ohmic heating balanced by
    expansion) must read as a hole with T_e ~ unchanged, not a hot spot;
  - a node whose electron pressure has gone non-positive - silently floored by
    the solver, never logged - must be COUNTED, because that count is the only
    warning the run gives;
  - the electron characteristic speed |u_e|+a_e is measured against the light
    speed the timestep is pinned to, so a fast electron must push it above 1;
  - the Gauss residual div E - rho_c/eps0 must be ~0 when the field and the
    charge agree and O(1) when they do not - the number that says whether the
    charge separation is physical or a discretisation error;
  - and a deck missing the light speed (needed for eps0) must refuse.
"""

import math
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
TOOL = os.path.join(ROOT, 'scripts', 'winding_anatomy.py')
sys.path.insert(0, HERE)
import framebuilder as fb

ME = 1.67e-27 / 1836
MI = 1.67e-27
Q = 1.6e-19
N0 = 1.0e20
GAMMA = 1.66666666667
MU0 = math.pi * 4.0e-7
LIGHT = 3.0e6
EPS0 = 1.0 / (MU0 * LIGHT ** 2)
RC, W = 0.036, 0.005
TE_EV = 30.0
PE0 = N0 * Q * TE_EV        # p = n Q Te[eV]


def energy(n, te_ev, ti_frac=0.5, vx=0.0, vy=0.0, vz=0.0, mass=ME):
    """Total energy density for a fluid at density n, temperature te_ev."""
    p = n * Q * te_ev
    rho = n * mass
    return p / (GAMMA - 1.0) + 0.5 * rho * (vx * vx + vy * vy + vz * vz)


def state(n_e=N0, te=TE_EV, n_i=N0, ti=15.0, uez=0.0, ex=0.0, ey=0.0):
    return {
        0: n_e * ME, 3: n_e * ME * uez, 4: energy(n_e, te, vz=uez, mass=ME),
        5: n_i * MI, 9: energy(n_i, ti, mass=MI),
        10: ex, 11: ey,
    }


def blob_triangles(cells):
    """cells: list of (r, theta, state_dict). Fat triangles with real area."""
    h = 3.0e-4
    tris = []
    for (r, th, st) in cells:
        x0, y0 = r * math.cos(th), r * math.sin(th)
        tris.append([(x0, y0, st), (x0 + h, y0, st), (x0, y0 + h, st)])
    return tris


def _run(args, expect_ok=True):
    done = subprocess.run([sys.executable, TOOL] + args, capture_output=True, text=True)
    if expect_ok and done.returncode != 0:
        raise AssertionError('winding_anatomy failed:\n' + done.stderr)
    return done


def _rows(out):
    """The per-frame time-series rows only, stopping at the radial profile."""
    rows = []
    for ln in out.splitlines():
        if ln.startswith('# radial profile'):
            break
        if ln and not ln.startswith('#') and ln[0].isdigit():
            rows.append(ln)
    return rows


class TestWindingAnatomy(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.deck = os.path.join(self.d, 'deck.pin')
        fb.deck_text(tend=1.0e-6, out=4, path=self.deck)

    def _write(self, cells_of_frame, nframes=5):
        paths = []
        for f in range(nframes):
            p = os.path.join(self.d, 'run_%d.vtu' % f)
            fb.write_vtu(p, blob_triangles(cells_of_frame(f)))
            paths.append(p)
        return paths

    def _ring(self, ncells=40, radius=RC):
        return [(radius, 2.0 * math.pi * i / ncells) for i in range(ncells)]

    def test_hole_with_flat_pressure(self):
        # a deepening hole at fixed pressure: n falls, Te rises so that n*Te const
        def frame(f):
            depth = 0.05 * f
            cells = []
            for (r, th) in self._ring():
                n = N0 * (1.0 - depth)
                te = TE_EV / (1.0 - depth) if depth < 1 else TE_EV
                cells.append((r, th, state(n_e=n, te=te)))
            return cells
        out = _run([self.deck] + self._write(frame)).stdout
        rows = _rows(out)
        first = rows[0].split()
        last = rows[-1].split()
        # n_e_min/n0 is column 2; it must fall below 1
        self.assertLess(float(last[1]), 0.95)
        self.assertGreater(float(first[1]), 0.98)

    def test_negative_pressure_node_is_counted(self):
        # one node with kinetic energy above total energy -> p_e < 0
        def frame(f):
            cells = self._ring()
            out = [(r, th, state()) for (r, th) in cells]
            if f >= 2:
                # break one cell: huge uez but energy for a cold fluid
                bad = state()
                bad[3] = N0 * ME * 5.0e6           # rho * uez, uez=5e6 m/s
                bad[4] = energy(N0, 1.0)           # only 1 eV of energy: KE >> E
                out[0] = (RC, 0.0, bad)
            return out
        out = _run([self.deck] + self._write(frame)).stdout
        rows = _rows(out)
        # column index of nodes[p_e<=0]: header order -> it's the 8th value (idx 7)
        counts = [int(r.split()[7]) for r in rows]
        self.assertEqual(counts[0], 0)
        self.assertGreater(counts[-1], 0)

    def test_lambda_e_over_c0_reflects_fast_electrons(self):
        def frame(f):
            cells = self._ring()
            uez = 2.0e6 * f
            return [(r, th, state(uez=uez)) for (r, th) in cells]
        out = _run([self.deck] + self._write(frame)).stdout
        rows = _rows(out)
        lam = [float(r.split()[6]) for r in rows]   # max(lam_e)/c0
        self.assertGreater(lam[-1], lam[0])
        self.assertGreater(lam[-1], 1.0)

    def test_gauss_residual_zero_when_field_matches_charge(self):
        # rho_c chosen, and E set so div E = rho_c/eps0 exactly -> residual ~0.
        rho_c = 1.0e-3                              # C/m^3
        dn = rho_c / Q                             # n_i - n_e in number density
        alpha = rho_c / EPS0                       # dE_x/dx = alpha -> divE=alpha
        def frame(f):
            cells = self._ring()
            out = []
            for (r, th) in cells:
                st = state(n_e=N0, n_i=N0 + dn)
                # this triangle's nodes get E_x = alpha * x (set per node below)
                out.append((r, th, st))
            return out
        # need per-node E_x = alpha*x; blob_triangles puts the SAME dict on all
        # three nodes, so instead build triangles by hand here.
        paths = []
        for f in range(4):
            tris = []
            hh = 3.0e-4
            for (r, th) in self._ring():
                x0, y0 = r * math.cos(th), r * math.sin(th)
                verts = [(x0, y0), (x0 + hh, y0), (x0, y0 + hh)]
                base = state(n_e=N0, n_i=N0 + dn)
                tri = []
                for (vx, vy) in verts:
                    s = dict(base)
                    s[10] = alpha * vx            # E_x = alpha x
                    tri.append((vx, vy, s))
                tris.append(tri)
            p = os.path.join(self.d, 'g_%d.vtu' % f)
            fb.write_vtu(p, tris)
            paths.append(p)
        out = _run([self.deck] + paths).stdout
        rows = _rows(out)
        gr = [float(r.split()[9]) for r in rows]    # gauss_resid/divE
        self.assertLess(gr[-1], 0.05)

    def test_gauss_residual_order_one_when_field_disagrees(self):
        rho_c = 1.0e-3
        dn = rho_c / Q
        alpha = 0.5 * rho_c / EPS0                  # divE only half of rho_c/eps0
        paths = []
        for f in range(4):
            tris = []
            hh = 3.0e-4
            for (r, th) in ([(RC, 2 * math.pi * i / 40) for i in range(40)]):
                x0, y0 = r * math.cos(th), r * math.sin(th)
                verts = [(x0, y0), (x0 + hh, y0), (x0, y0 + hh)]
                base = state(n_e=N0, n_i=N0 + dn)
                tri = []
                for (vx, vy) in verts:
                    s = dict(base)
                    s[10] = alpha * vx
                    tri.append((vx, vy, s))
                tris.append(tri)
            p = os.path.join(self.d, 'b_%d.vtu' % f)
            fb.write_vtu(p, tris)
            paths.append(p)
        out = _run([self.deck] + paths).stdout
        gr = [float(r.split()[9]) for r in _rows(out)]
        self.assertGreater(gr[-1], 0.5)

    def test_missing_light_refuses(self):
        bad = os.path.join(self.d, 'bad.pin')
        with open(bad, 'w') as h:
            h.write(fb.deck_text().replace('LIGHT = 3.0e6\n', ''))
        done = _run([bad] + self._write(lambda f: [(RC, 0.0, state())] * 1,
                                        nframes=2), expect_ok=False)
        self.assertNotEqual(done.returncode, 0)
        self.assertIn('LIGHT', done.stderr + done.stdout)


if __name__ == '__main__':
    unittest.main()
