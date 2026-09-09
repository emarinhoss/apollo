#!/usr/bin/env python3
"""Does src/lib/petsc_compat.h substitute at the right PETSc versions?

Apollo is built against whatever PETSc the machine has, and the compatibility
header decides per symbol whether to use the real declaration or its own
stand-in. That decision is a version comparison, and a version comparison that
names the wrong release is invisible until someone builds somewhere else:

    Ubuntu 24.04 ships PETSc 3.19  -> built fine
    Ubuntu 22.04 ships PETSc 3.15  -> did not compile at all

because two shims claimed their symbols were removed in 3.18 when they were
removed in 3.13 and 3.14. CI only ran the newer release, so nothing noticed.

This test compiles the real header against a stub <petsc.h> that declares what
each release actually declares (see petsc_compat_stub/petsc.h for where that
table was read from), once per release across the supported range. It needs a
C++ compiler and nothing else - no PETSc, no MPI, no built solver.
"""

import os
import shutil
import subprocess
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
STUB = os.path.join(HERE, 'petsc_compat_stub')
PROBE = os.path.join(STUB, 'probe.cc')
COMPAT = os.path.join(ROOT, 'src', 'lib')

# 3.11 is the floor because the DMPlexGetLabelValue macro is guarded there and
# every PETSc anyone can still install is above it; the ceiling is open-ended
# on purpose, so a release newer than any that existed when this was written
# still has to pass.
SUPPORTED_MINORS = list(range(11, 24))


def compiler():
    for candidate in (os.environ.get('CXX'), 'g++', 'c++', 'clang++'):
        if candidate and shutil.which(candidate):
            return candidate
    return None


def compile_at(minor, extra=()):
    """Syntax-check the probe against the stub emulating PETSc 3.<minor>."""
    cxx = compiler()
    with tempfile.TemporaryDirectory() as tmp:
        cmd = [
            cxx, '-std=c++14', '-fsyntax-only', '-Wall',
            '-I', STUB, '-I', COMPAT,
            '-DSTUB_PETSC_MINOR=%d' % minor,
            *extra,
            PROBE,
        ]
        proc = subprocess.run(cmd, cwd=tmp, capture_output=True, text=True)
    return proc.returncode, proc.stderr


@unittest.skipIf(compiler() is None, 'no C++ compiler')
class TestPetscCompat(unittest.TestCase):

    def test_compiles_against_every_supported_release(self):
        """The header must work on every PETSc from 3.11 up.

        A failure here names the release, which is the whole point: the bug
        this replaces reported itself as 300 lines of template errors on one
        distribution and nothing at all on another.
        """
        for minor in SUPPORTED_MINORS:
            with self.subTest(petsc='3.%d' % minor):
                code, err = compile_at(minor)
                self.assertEqual(
                    code, 0,
                    'petsc_compat.h does not compile against PETSc 3.%d:\n%s'
                    % (minor, err))

    def test_shim_does_not_fire_while_petsc_still_declares_the_symbol(self):
        """The other direction, which is the one a too-early guard breaks.

        Defining a static inline over a declaration still in scope is an error,
        not a silent override, so the releases that predate each removal are
        what pin the lower edge of the guards: 3.12 still has
        DMPlexGetHybridBounds and 3.13 still has the adjacency pair.
        """
        for minor in (11, 12, 13):
            with self.subTest(petsc='3.%d' % minor):
                code, err = compile_at(minor)
                self.assertEqual(code, 0, err)

    def test_the_versions_that_actually_regressed(self):
        """3.13 through 3.17 - the window the old 3.18 guards left broken.

        Ubuntu 22.04's 3.15 and Debian bullseye's 3.14 both live here, so this
        is not a hypothetical range.
        """
        for minor in (13, 14, 15, 16, 17):
            with self.subTest(petsc='3.%d' % minor):
                code, err = compile_at(minor)
                self.assertEqual(
                    code, 0,
                    'PETSc 3.%d is the case that shipped broken:\n%s'
                    % (minor, err))

    def test_the_probe_can_fail(self):
        """Guard against a probe that would pass no matter what.

        If the stub and the header ever stopped being compiled together the
        tests above would pass vacuously, so check that the same command does
        fail when the header is deliberately denied its shims.
        """
        code, err = compile_at(15, extra=('-DPETSC_COMPAT_H',))
        self.assertNotEqual(
            code, 0,
            'the probe compiled with petsc_compat.h suppressed, so it is not '
            'actually testing the header')
        self.assertIn('DMPlexGetHybridBounds', err)


if __name__ == '__main__':
    unittest.main(verbosity=2)
