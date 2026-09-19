#!/usr/bin/env python3
"""Does the .pin -> .inp preprocessor actually substitute the deck's Python?

A deck's preamble is Python: it sets P_ORDER, MEQN, densities and so on, and the
body refers to those names in arithmetic. scripts/wxinputparser.py evaluates the
preamble and substitutes the results. When that substitution silently does not
happen, the expressions reach the solver verbatim and it aborts with

    Error: unexpected symbol '/' on line 23 of the input file

naming neither the key nor the reason. That is what happened on Python 3.13:
the parser did `exec(co)` inside a method and then recovered each name with
`eval(name)`, which works only because CPython <= 3.12 let exec() writes into
locals() persist on the frame. PEP 667 made locals() an independent snapshot,
every eval() raised NameError, a bare `except: pass` swallowed them all, and the
substitution dictionary came back empty. Nothing failed loudly; the deck just
stopped being preprocessed.

The interesting property of that bug is that it is invisible on the interpreter
you happen to run. So this checks two things: that substitution works here, and
that every python3.x on this machine produces the same file.
"""

import glob
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SCRIPTS = os.path.join(ROOT, 'scripts')
PARSER = os.path.join(SCRIPTS, 'wxinpparse.py')
DECK = os.path.join(ROOT, 'examples', 'unstructuredDG', 'multifluid', 'rmf_frc',
                    'phase3', '00-divergence-gate', 'edge-cleaning-on.pin')

# A deck in miniature: a preamble that defines names, and a body that does
# arithmetic on them. Kept here so the test still means something if the phase3
# decks move.
MINIMAL_PIN = """# -*- python -*-
P_ORDER = 1
MEQN = 18
SCALE = 4.0

<apollo>
  Simulation = tiny
  <hyperbolic>
    polynomialOrder = P_ORDER
    cfl = 1.0/((2.*(P_ORDER+1)-1))
    slots = MEQN*(P_ORDER+1)*(P_ORDER+2)//2
    half = SCALE/2.
  </hyperbolic>
</apollo>
"""


def interpreters():
    """Every python3.x on this machine, plus the one running the tests."""
    found = {os.path.realpath(sys.executable)}
    for path in glob.glob('/usr/bin/python3.*') + glob.glob('/usr/local/bin/python3.*'):
        if path.endswith('-config'):
            continue
        if os.access(path, os.X_OK):
            found.add(os.path.realpath(path))
    return sorted(found)


def preprocess(pin_text, python=None, name='deck'):
    """Run the .pin through the preprocessor and return the generated .inp."""
    python = python or sys.executable
    tmp = tempfile.mkdtemp()
    try:
        pin = os.path.join(tmp, name + '.pin')
        with open(pin, 'w') as handle:
            handle.write(pin_text)
        env = dict(os.environ, PYTHONPATH=SCRIPTS)
        proc = subprocess.run([python, PARSER, '-i', name + '.pin'],
                              cwd=tmp, env=env, capture_output=True, text=True)
        out = os.path.join(tmp, name + '.inp')
        if not os.path.exists(out):
            raise AssertionError('no .inp produced by %s:\n%s\n%s'
                                 % (python, proc.stdout, proc.stderr))
        with open(out) as handle:
            return handle.read()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def unsubstituted(inp_text):
    """Value lines that still carry a macro name or bare arithmetic.

    Only the left of a '#' and outside quotes matters: exprList and progn hold
    expression strings that are meant to reach the solver unevaluated.
    """
    bad = []
    for number, line in enumerate(inp_text.splitlines(), start=1):
        if '=' not in line or '"' in line or "'" in line:
            continue
        value = line.split('=', 1)[1].strip()
        if not value or value.startswith('['):
            continue
        # A value that parses as a number is substituted, whatever characters it
        # contains: 1e+20 holds a '+' and is perfectly fine.
        try:
            float(value)
            continue
        except ValueError:
            pass
        if any(op in value for op in ('/', '*', '+')) or 'P_ORDER' in value:
            bad.append((number, line.strip()))
    return bad


class TestDeckPreprocessing(unittest.TestCase):

    def test_minimal_deck_is_substituted(self):
        out = preprocess(MINIMAL_PIN)
        self.assertIn('polynomialOrder = 1', out)
        self.assertIn('cfl = 0.3333333333333333', out)
        self.assertIn('slots = 54', out)
        self.assertIn('half = 2.0', out)
        self.assertEqual(unsubstituted(out), [],
                         'arithmetic survived preprocessing: %s' % unsubstituted(out))

    def test_the_check_can_fail(self):
        """The detector must flag an unsubstituted deck, or it proves nothing."""
        broken = 'a = 1.0/((2.*(P_ORDER+1)-1))\n'
        self.assertTrue(unsubstituted(broken),
                        'unsubstituted() did not flag an obviously broken line')

    @unittest.skipUnless(os.path.exists(DECK), 'phase3 decks not present')
    def test_real_deck_is_substituted(self):
        """The deck that actually failed, at the line that actually failed."""
        with open(DECK) as handle:
            out = preprocess(handle.read(), name='edge-cleaning-on')
        self.assertIn('cfl = 0.3333333333333333', out)
        self.assertIn('polynomialOrder = 1', out)
        self.assertEqual(unsubstituted(out), [])

    def test_every_interpreter_agrees(self):
        """The bug was invisible on 3.11 and fatal on 3.13.

        Comparing the interpreters present on this machine catches a
        version-dependent regression even when the default one is unaffected.
        """
        pythons = interpreters()
        if len(pythons) < 2:
            self.skipTest('only one python3 available: %s' % pythons)
        outputs = {}
        for python in pythons:
            outputs[python] = preprocess(MINIMAL_PIN, python=python)
        reference = outputs[pythons[0]]
        for python, out in outputs.items():
            with self.subTest(python=python):
                self.assertEqual(
                    out, reference,
                    '%s produced a different .inp than %s' % (python, pythons[0]))
                self.assertEqual(unsubstituted(out), [])


if __name__ == '__main__':
    unittest.main(verbosity=2)
