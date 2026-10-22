#!/usr/bin/env python3
"""run_one.sh: what it does when there is already a run in the results directory.

WHY THIS EXISTS. Before the solver could checkpoint, re-running `run_one.sh`
into a populated directory cost nothing - it recomputed frames that were cheap
to recompute. Now the directory can hold a checkpoint worth hours, and a fresh
run overwrites frame 0 and the checkpoint with it. The script's answer is to
REFUSE unless told which way to go, and a refusal nobody tests is decoration.

So these check the four ways in:

    no checkpoint                  -> run, no -r
    unfinished, APOLLO_RESUME unset-> refuse, exit 3, checkpoint UNTOUCHED
    unfinished, APOLLO_RESUME=1    -> run with -r <run>.checkpoint
    unfinished, APOLLO_RESUME=0    -> delete it and run, no -r

and the post-mortem, which distinguishes a run the wall clock killed (resume it)
from a run the solver stopped on its own (resuming lands back in the same
place). None of this needs a solver: APOLLO_BIN is a stub that records its
argv. The physics of restart is test_restart.py's job, not this file's.
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

from test_run_growth import _vtu, PHI                       # noqa: E402

PHASE3 = os.path.join(ROOT, 'examples', 'unstructuredDG', 'multifluid',
                      'rmf_frc', 'phase3')
RUN_ONE = os.path.join(PHASE3, 'run_one.sh')
DECK = '00-divergence-gate/antenna-cleaning-on.pin'
NAME = 'antenna-cleaning-on'

# A stand-in for the solver. It writes down how it was called, optionally lays
# down frames and a checkpoint, and exits with whatever it was told to.
STUB = r'''#!/bin/bash
printf '%s\n' "$@" > "$ARGV_FILE"
if [[ -n "${STUB_FRAMES:-}" ]]; then
    cp "$STUB_FRAMES"/*.vtu . 2>/dev/null || true
fi
if [[ -n "${STUB_META:-}" ]]; then
    printf 'frame %s\ntime 1.0e-09\ntend 5.0e-09\nnout 6\nsize 658800\n' \
        "$STUB_META" > "$STUB_NAME.checkpoint.meta"
    : > "$STUB_NAME.checkpoint"
fi
exit "${STUB_EXIT:-0}"
'''


def _have_deck():
    return os.path.exists(os.path.join(PHASE3, DECK)) and os.path.exists(RUN_ONE)


@unittest.skipUnless(_have_deck(), 'needs the phase3 decks')
class TestRunOneResume(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix='apollo-run-one-')
        self.results = os.path.join(self.tmp, 'results')
        os.makedirs(self.results)
        self.stub = os.path.join(self.tmp, 'stub-apollo')
        with open(self.stub, 'w') as h:
            h.write(STUB)
        os.chmod(self.stub, 0o755)
        self.argv = os.path.join(self.tmp, 'argv')

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    # -- helpers ------------------------------------------------------------

    def _checkpoint(self, frame, nout=6, body=b'STATE'):
        """Put a checkpoint of our own in the results directory.

        With the deck and mesh beside it, because that is what an interrupted
        run leaves behind and what the resume path reads: it deliberately does
        NOT re-copy them, so that a deck edited in the source folder since
        cannot be resumed into.
        """
        src = os.path.join(PHASE3, os.path.dirname(DECK))
        shutil.copy(os.path.join(PHASE3, DECK), self.results)
        shutil.copy(os.path.join(src, 'disc.msh'), self.results)
        with open(os.path.join(self.results, NAME + '.checkpoint'), 'wb') as h:
            h.write(body)
        with open(os.path.join(self.results, NAME + '.checkpoint.meta'), 'w') as h:
            h.write('frame %d\ntime 1.82e-09\ntend 5.47e-09\nnout %d\nsize 658800\n'
                    % (frame, nout))

    def _frames(self, n):
        d = os.path.join(self.tmp, 'frames')
        os.makedirs(d, exist_ok=True)
        for i in range(n):
            _vtu(os.path.join(d, '%s_%d.vtu' % (NAME, i)), {PHI: 1.0e-5 * (10 ** i)})
        return d

    def _run(self, resume=None, exit_code=0, meta=None, frames=None):
        env = dict(os.environ, APOLLO_BIN=self.stub, ARGV_FILE=self.argv,
                   STUB_EXIT=str(exit_code), STUB_NAME=NAME)
        if resume is not None:
            env['APOLLO_RESUME'] = resume
        else:
            env.pop('APOLLO_RESUME', None)
        if meta is not None:
            env['STUB_META'] = str(meta)
        if frames is not None:
            env['STUB_FRAMES'] = frames
        done = subprocess.run([RUN_ONE, DECK, self.results], env=env,
                              capture_output=True, text=True, timeout=300)
        return done

    def _argv(self):
        with open(self.argv) as h:
            return [l.rstrip('\n') for l in h]

    def _ckpt_bytes(self):
        with open(os.path.join(self.results, NAME + '.checkpoint'), 'rb') as h:
            return h.read()

    # -- the four ways in ---------------------------------------------------

    def test_an_empty_directory_runs_from_the_start(self):
        self._run(meta=6)
        self.assertNotIn('-r', self._argv())

    def test_an_unfinished_run_is_not_overwritten_without_being_told(self):
        """The guard. Exit 3, and the checkpoint still byte-for-byte there."""
        self._checkpoint(2, body=b'IRREPLACEABLE')
        done = self._run()
        self.assertEqual(done.returncode, 3, done.stdout + done.stderr)
        self.assertFalse(os.path.exists(self.argv), 'the solver was started anyway')
        self.assertEqual(self._ckpt_bytes(), b'IRREPLACEABLE')
        out = done.stdout + done.stderr
        self.assertIn('UNFINISHED', out)
        self.assertIn('APOLLO_RESUME=1', out)
        self.assertIn('APOLLO_RESUME=0', out)

    def test_the_refusal_names_the_results_directory_it_was_given(self):
        """Not the default one - following the advice has to reach this run."""
        self._checkpoint(2)
        done = self._run()
        for line in (done.stdout + done.stderr).splitlines():
            if 'APOLLO_RESUME=1' in line:
                self.assertIn(self.results, line, line)
                break
        else:
            self.fail('no APOLLO_RESUME=1 line in the refusal')

    def test_resume_passes_the_checkpoint_to_the_solver(self):
        self._checkpoint(2)
        self._run(resume='1', meta=6)
        argv = self._argv()
        self.assertIn('-r', argv)
        self.assertEqual(argv[argv.index('-r') + 1], NAME + '.checkpoint')

    def test_resume_zero_discards_the_checkpoint_and_starts_over(self):
        self._checkpoint(2)
        stale = os.path.join(self.results, '%s_5.vtu' % NAME)
        _vtu(stale, {PHI: 1.0})
        self._run(resume='0', meta=6)
        self.assertNotIn('-r', self._argv())
        self.assertFalse(os.path.exists(stale),
                         'a stale frame from the discarded run survived')

    def test_a_finished_run_is_not_silently_redone(self):
        self._checkpoint(6)
        done = self._run()
        self.assertEqual(done.returncode, 3, done.stdout + done.stderr)
        self.assertIn('COMPLETE', done.stdout + done.stderr)
        self.assertFalse(os.path.exists(self.argv), 'it ran the whole thing again')

    # -- the post-mortem ----------------------------------------------------

    def test_a_signal_is_reported_as_an_interruption(self):
        """143 is SIGTERM: the wall clock. Resuming is exactly right, and the
        message must not tell them it ran into an instability."""
        done = self._run(exit_code=143, meta=2, frames=self._frames(2))
        out = done.stdout + done.stderr
        self.assertEqual(done.returncode, 143, out)
        self.assertIn('INTERRUPTED', out)
        self.assertIn('APOLLO_RESUME=1', out)
        self.assertNotIn('stopped on its own', out)

    def test_a_nonzero_exit_is_reported_as_the_solver_stopping(self):
        """exit(1) is what the NaN guard does; a resume lands back in it."""
        done = self._run(exit_code=1, meta=2, frames=self._frames(2))
        out = done.stdout + done.stderr
        self.assertEqual(done.returncode, 1, out)
        self.assertIn('DIED', out)
        self.assertIn('stopped on its own', out)

    def test_a_death_runs_the_growth_table_over_the_frames(self):
        """The whole point: the run that most needs analysing used to get none,
        because `set -e` ended the script on the solver's exit status."""
        done = self._run(exit_code=1, meta=2, frames=self._frames(3))
        out = done.stdout + done.stderr
        self.assertIn('what was growing', out)
        self.assertIn('growth, each from its own first non-zero frame', out)
        self.assertTrue(os.path.exists(os.path.join(self.results, 'growth.txt')), out)

    def test_a_death_with_no_frames_says_so_rather_than_failing(self):
        done = self._run(exit_code=1)
        out = done.stdout + done.stderr
        self.assertEqual(done.returncode, 1, out)
        self.assertIn('no checkpoint was written', out)


if __name__ == '__main__':
    unittest.main(verbosity=2)
