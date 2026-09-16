#!/usr/bin/env python3
"""Restart: does a resumed run reproduce the run it claims to continue?

WHY THIS EXISTS. Apollo had no restart, so an interrupted run was a lost run
(docs/known-issues.md section 6 said exactly that). A 20-period formation run
died 1.73 periods in after about 12 hours, and every one of those hours had to
be paid again. The state was on disk in the output frames the whole time; only
a way to read one back was missing.

WHAT MUST HOLD, in order of how badly it hurts to get wrong:

1. A resumed run produces the SAME NUMBERS as running straight through. Not
   close - the same. If it does not, a restart silently changes the answer and
   is worse than no restart at all.
2. Resuming twice must work. The first implementation recorded the RUNNING
   interval count rather than the deck's, so resuming at frame 2 of 6 wrote
   "nout 4"; resuming from that then compared 4 against the deck's 6 and
   refused. One resume worked and a second did not - the wrong half.
3. A checkpoint from a DIFFERENT problem must be refused, not resumed into.

These need a built solver and are skipped without one, so check the skip count.
"""

import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

APOLLO = os.environ.get('APOLLO_BIN', os.path.join(ROOT, 'src', 'build-opt', 'apollo'))
DECK_DIR = os.path.join(ROOT, 'examples', 'unstructuredDG', 'multifluid',
                        'rmf_frc', 'phase3', '00-divergence-gate')
DECK = os.path.join(DECK_DIR, 'antenna-cleaning-on.pin')
MESH = os.path.join(DECK_DIR, 'disc.msh')
PARSER = os.path.join(ROOT, 'scripts', 'wxinpparse.py')

# The gate deck: OUT = 6, so seven frames, about 90 s straight through.
FRAMES = 6


def _have_solver():
    return (os.path.exists(APOLLO) and os.access(APOLLO, os.X_OK)
            and os.path.exists(DECK) and os.path.exists(MESH))


def _stage(work):
    """Copy the deck and mesh into `work` and expand the .pin."""
    shutil.copy(DECK, os.path.join(work, 'gate.pin'))
    shutil.copy(MESH, work)
    env = dict(os.environ, PYTHONPATH=os.path.join(ROOT, 'scripts'))
    subprocess.run([sys.executable, PARSER, '-i', 'gate.pin'],
                   cwd=work, env=env, capture_output=True, check=True)


def _run(work, extra=(), timeout=900):
    return subprocess.run([APOLLO, '-i', 'gate.inp'] + list(extra),
                          cwd=work, capture_output=True, text=True,
                          timeout=timeout)


def _meta(work):
    out = {}
    with open(os.path.join(work, 'gate.checkpoint.meta')) as h:
        for line in h:
            bits = line.split()
            if len(bits) == 2 and not line.startswith('#'):
                out[bits[0]] = bits[1]
    return out


def _write_meta(work, frame, nout, size, time='2.73584905660377360e-09'):
    with open(os.path.join(work, 'gate.checkpoint.meta'), 'w') as h:
        h.write('frame %s\ntime %s\ntend 5.47169811320754719e-09\n'
                'nout %s\nsize %s\n' % (frame, time, nout, size))


@unittest.skipUnless(_have_solver(), 'needs a built solver and the gate deck')
class TestRestartRefusals(unittest.TestCase):
    """The cheap half: a bad checkpoint must be named, not resumed into.

    These need only one short run to produce a checkpoint to corrupt.
    """

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix='apollo-restart-')
        _stage(cls.tmp)
        done = _run(cls.tmp)
        if done.returncode != 0:
            raise unittest.SkipTest('the gate deck did not run: %s'
                                    % done.stdout[-400:])
        cls.size = _meta(cls.tmp)['size']
        # Keep a pristine copy: a SUCCESSFUL resume rewrites the checkpoint,
        # which silently invalidated an earlier version of these tests.
        shutil.copy(os.path.join(cls.tmp, 'gate.checkpoint'),
                    os.path.join(cls.tmp, 'pristine'))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def _refused(self, frame, nout, size):
        shutil.copy(os.path.join(self.tmp, 'pristine'),
                    os.path.join(self.tmp, 'gate.checkpoint'))
        _write_meta(self.tmp, frame, nout, size)
        done = _run(self.tmp, ['-r', 'gate.checkpoint'])
        self.assertNotEqual(done.returncode, 0,
                            'resumed into a bad checkpoint:\n' + done.stdout[-600:])
        return done.stdout + done.stderr

    def test_a_state_of_the_wrong_size_is_refused(self):
        """A different mesh or equation set is not a continuation."""
        out = self._refused(3, FRAMES, 12345)
        self.assertIn('holds a state of', out)

    def test_a_different_output_count_is_refused(self):
        """Frame spacing would differ, so the frames would not line up."""
        out = self._refused(3, 12, self.size)
        self.assertIn('Output_files', out)

    def test_a_frame_past_the_end_is_refused(self):
        out = self._refused(99, FRAMES, self.size)
        self.assertIn('outside', out)

    def test_a_negative_frame_is_refused(self):
        out = self._refused(-1, FRAMES, self.size)
        self.assertIn('outside', out)

    def test_a_finished_run_is_refused(self):
        """frame == nout means there is nothing left to do."""
        out = self._refused(FRAMES, FRAMES, self.size)
        self.assertIn('already finished', out)

    def test_a_missing_sidecar_is_refused(self):
        shutil.copy(os.path.join(self.tmp, 'pristine'),
                    os.path.join(self.tmp, 'gate.checkpoint'))
        meta = os.path.join(self.tmp, 'gate.checkpoint.meta')
        _write_meta(self.tmp, 3, FRAMES, self.size)
        os.rename(meta, meta + '.hidden')
        try:
            done = _run(self.tmp, ['-r', 'gate.checkpoint'])
        finally:
            os.rename(meta + '.hidden', meta)
        self.assertNotEqual(done.returncode, 0, done.stdout[-600:])
        self.assertIn('no checkpoint metadata', done.stdout + done.stderr)

    def test_the_checkpoint_records_the_decks_interval_count(self):
        """Not the running one - that is what made restart non-idempotent."""
        self.assertEqual(_meta(self.tmp)['nout'], str(FRAMES))


@unittest.skipUnless(_have_solver(), 'needs a built solver and the gate deck')
class TestRestartReproduces(unittest.TestCase):
    """The expensive half: does resuming give the same answer? ~3 min."""

    def _fingerprint(self, work):
        sys.path.insert(0, os.path.join(ROOT, 'scripts'))
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            'rf', os.path.join(ROOT, 'scripts', 'run_fingerprint.py'))
        rf = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(rf)
        return rf, rf.fingerprint(work)

    def test_resuming_twice_reproduces_a_straight_run_exactly(self):
        straight = tempfile.mkdtemp(prefix='apollo-straight-')
        resumed = tempfile.mkdtemp(prefix='apollo-resumed-')
        try:
            _stage(straight)
            self.assertEqual(_run(straight).returncode, 0)

            # Run, stop at a frame, resume - twice, because the first
            # implementation worked once and refused the second time.
            _stage(resumed)
            for stop_after in (2, 4):
                proc = subprocess.Popen(
                    [APOLLO, '-i', 'gate.inp'] +
                    (['-r', 'gate.checkpoint'] if stop_after > 2 else []),
                    cwd=resumed, stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT, text=True)
                # Stop once the frame we want is on disk.
                import time as _t
                deadline = _t.time() + 600
                while _t.time() < deadline:
                    try:
                        if int(_meta(resumed)['frame']) >= stop_after:
                            break
                    except (OSError, KeyError, ValueError):
                        pass
                    if proc.poll() is not None:
                        break
                    _t.sleep(1)
                if proc.poll() is None:
                    proc.terminate()
                    proc.wait(timeout=60)
            # Finish it.
            done = _run(resumed, ['-r', 'gate.checkpoint'])
            self.assertEqual(done.returncode, 0, done.stdout[-600:])

            rf, left = self._fingerprint(straight)
            right = rf.fingerprint(resumed)
            worst, rows, refusals, notes = rf.compare(left, right, 1e-12)
            self.assertFalse(refusals, refusals)
            fields = [r for r in rows if r[1] not in rf.MESH_ARRAYS]
            self.assertTrue(fields, 'nothing was compared')
            self.assertEqual(
                fields[0][2], 0.0,
                'a resumed run must reproduce a straight run exactly; worst '
                'was %.3e in %s' % (fields[0][2], fields[0][1]))
        finally:
            shutil.rmtree(straight, ignore_errors=True)
            shutil.rmtree(resumed, ignore_errors=True)


if __name__ == '__main__':
    unittest.main(verbosity=2)
