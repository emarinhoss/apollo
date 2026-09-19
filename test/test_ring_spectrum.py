#!/usr/bin/env python3
"""Tests for scripts/ring_spectrum.py.

The tool has to make one call reliably: is a growing field a smooth, low-m,
drive-locked pattern, or grid-scale, cell-to-cell structure the mesh cannot
resolve? Those are the physical and the numerical instability and the tool is
what tells them apart, so the tests plant each and check the label:

  - an m=0 field that grows LINEARLY must be reported [lin], not as an
    exponential rate (a linear fit dressed as an exponential is exactly the
    "2.2e6/s" that a run's max|value| produces and that must not be mistaken
    for a growth rate);
  - a node-to-node grid field that grows EXPONENTIALLY must be [exp];
  - a planted m=7 pattern must come back as the dominant harmonic;
  - a peak pinned to a fixed azimuth must be reported at that azimuth;
  - and the two ways it could label nothing: too few nodes on a ring to Fourier,
    and a missing deck key, must both refuse rather than guess.
"""

import math
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
TOOL = os.path.join(ROOT, 'scripts', 'ring_spectrum.py')
sys.path.insert(0, HERE)
import framebuilder as fb

ME = 1.67e-27 / 1836
MI = 1.67e-27
N0 = 1.0e20
RC = 0.036


def uniform(extra=None):
    base = {0: N0 * ME, 5: N0 * MI}
    if extra:
        base.update(extra)
    return base


def _run(args, expect_ok=True):
    done = subprocess.run([sys.executable, TOOL] + args,
                          capture_output=True, text=True)
    if expect_ok and done.returncode != 0:
        raise AssertionError('ring_spectrum failed:\n' + done.stderr)
    return done


class TestRingSpectrum(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.deck = os.path.join(self.d, 'deck.pin')
        fb.deck_text(tend=1.0e-6, out=8, path=self.deck)

    def _frames(self, comp_fn_of_frame, nframes=9, ncells=120, radius=RC):
        paths = []
        for f in range(nframes):
            p = os.path.join(self.d, 'run_%d.vtu' % f)
            fb.write_vtu(p, fb.ring_triangles(radius, ncells, comp_fn_of_frame(f)))
            paths.append(p)
        return paths

    def _line(self, out, comp):
        for ln in out.splitlines():
            if comp in ln and ln.startswith('winding'):
                return ln
        raise AssertionError('no winding line for %s in:\n%s' % (comp, out))

    def test_linear_m0_is_labelled_linear(self):
        # phi uniform in theta, growing linearly frame to frame past the ramp
        def fn(f):
            return lambda th, i: uniform({16: 1e-6 * (1 + f)})
        out = _run(['--components', 'phi', self.deck] + self._frames(fn)).stdout
        self.assertIn('[lin]', self._line(out, 'phi '))

    def test_exponential_grid_band_is_labelled_exp(self):
        # phi with a node-to-node alternating (grid) part growing exponentially
        def fn(f):
            amp = 1e-8 * (2.0 ** f)
            return lambda th, i: uniform({16: amp * (-1) ** i})
        out = _run(['--components', 'phi', self.deck] + self._frames(fn)).stdout
        line = self._line(out, 'phi ')
        # the grid field must be seen growing and called exponential
        self.assertIn('grid=+', line)
        self.assertIn('[exp]', line.split('grid=')[1])

    def test_planted_azimuthal_mode_is_recovered(self):
        # phi = cos(m theta) with m=7 (a scalar: no radial projection to mix it
        # into m-1/m+1) -> dominant harmonic 7.
        m = 7
        def fn(f):
            return lambda th, i: uniform({16: 1.0e-6 * math.cos(m * th)})
        out = _run(['--components', 'phi', self.deck] + self._frames(fn)).stdout
        line = self._line(out, 'phi ')
        dom = int(line.split('dom_m=')[1].split()[0])
        self.assertEqual(dom, m)

    def test_peak_azimuth_is_fixed(self):
        # e-rho peak planted in one cell at a fixed angle across frames
        theta_peak = 2.0 * math.pi * 40 / 120     # cell 40 of 120
        def fn(f):
            def g(th, i):
                v = uniform()
                if i == 40:
                    v[0] = N0 * ME * (1.0 + 0.5 * (f + 1))
                return v
            return g
        out = _run([self.deck] + self._frames(fn)).stdout
        row = [ln for ln in out.splitlines() if ln.startswith('e-rho')][0]
        # every non-first frame should sit near the planted azimuth
        want = math.degrees(theta_peak)
        angles = [float(tok.split(',')[1].rstrip(')'))
                  for tok in row.split(':', 1)[1].split()]
        self.assertTrue(all(abs(a - want) < 5.0 for a in angles[1:]),
                        'peak azimuths %s not pinned near %.1f' % (angles[1:], want))

    def test_too_few_ring_nodes_refuses(self):
        def fn(f):
            return lambda th, i: uniform({16: 1e-6 * (1 + f)})
        # 10 cells -> ~30 nodes on the ring, below MIN_RING
        done = _run(['--components', 'phi', self.deck]
                    + self._frames(fn, ncells=10), expect_ok=False)
        self.assertNotEqual(done.returncode, 0)
        self.assertIn('ring', (done.stderr + done.stdout).lower())

    def test_missing_deck_key_refuses(self):
        bad = os.path.join(self.d, 'bad.pin')
        text = fb.deck_text().replace('omega = 7.95e5\n', '')
        with open(bad, 'w') as h:
            h.write(text)
        def fn(f):
            return lambda th, i: uniform({16: 1e-6})
        done = _run(['--components', 'phi', bad] + self._frames(fn), expect_ok=False)
        self.assertNotEqual(done.returncode, 0)
        self.assertIn('omega', done.stderr + done.stdout)


if __name__ == '__main__':
    unittest.main()
