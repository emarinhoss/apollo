#!/usr/bin/env python3
"""Does the build link against the conda environment it is running inside?

SConstruct builds with a scrubbed environment - only PATH, LM_LICENSE_FILE and
LD_LIBRARY_PATH are passed through - so conda's own

    LDFLAGS="-Wl,-rpath,$CONDA_PREFIX/lib -L$CONDA_PREFIX/lib"

never reaches the linker. The build then compiles against the environment's
headers and links against whatever the system has. On a cluster that produced a
binary needing libgfortran.so.3, a GCC 4 era library, out of a conda GCC 15
toolchain whose runtime is libgfortran.so.5 - so it linked and would not start.

The -L and the rpath do different jobs and both are needed: the first decides
which library is chosen at link time, the second whether it can be found at run
time. A batch job that runs the solver may never activate the environment, and
the SONAME baked into a binary cannot be corrected afterwards - it has to be
relinked.
"""

import os
import sys
import tempfile
import unittest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'src'))

from buildconf.conda_paths import conda_link_paths


class TestCondaLinkPaths(unittest.TestCase):

    def setUp(self):
        self.prefix = tempfile.mkdtemp()
        os.mkdir(os.path.join(self.prefix, 'lib'))
        self.lib = os.path.join(self.prefix, 'lib')

    def tearDown(self):
        import shutil
        shutil.rmtree(self.prefix, ignore_errors=True)

    def test_inert_outside_conda(self):
        """A distribution build must be completely unaffected."""
        self.assertEqual(conda_link_paths({}), ([], []))

    def test_active_environment_supplies_both(self):
        libpath, linkflags = conda_link_paths({'CONDA_PREFIX': self.prefix})
        self.assertEqual(libpath, [self.lib])
        self.assertEqual(linkflags, ['-Wl,-rpath,' + self.lib])

    def test_rpath_is_recorded_not_merely_searched(self):
        """-rpath, not -rpath-link: it has to survive into the binary.

        A SLURM job that runs the solver may not have activated the
        environment, and LD_LIBRARY_PATH is not inherited from the build.
        """
        _libpath, linkflags = conda_link_paths({'CONDA_PREFIX': self.prefix})
        self.assertTrue(linkflags[0].startswith('-Wl,-rpath,'), linkflags)
        self.assertNotIn('rpath-link', linkflags[0])

    def test_conda_build_prefixes(self):
        """conda-build sets PREFIX rather than CONDA_PREFIX, and never activates."""
        self.assertEqual(conda_link_paths({'PREFIX': self.prefix})[0], [self.lib])
        self.assertEqual(conda_link_paths({'BUILD_PREFIX': self.prefix})[0],
                         [self.lib])

    def test_conda_prefix_wins_over_build_prefix(self):
        other = tempfile.mkdtemp()
        os.mkdir(os.path.join(other, 'lib'))
        try:
            libpath, _ = conda_link_paths({'CONDA_PREFIX': self.prefix,
                                           'PREFIX': other})
            self.assertEqual(libpath, [self.lib])
        finally:
            import shutil
            shutil.rmtree(other, ignore_errors=True)

    def test_prefix_without_a_lib_directory_is_ignored(self):
        """Refuse rather than emit -L at a path that does not exist."""
        empty = tempfile.mkdtemp()
        try:
            self.assertEqual(conda_link_paths({'CONDA_PREFIX': empty}), ([], []))
        finally:
            import shutil
            shutil.rmtree(empty, ignore_errors=True)

    def test_nonexistent_prefix_is_ignored(self):
        self.assertEqual(conda_link_paths({'CONDA_PREFIX': '/nope/not/here'}),
                         ([], []))


if __name__ == '__main__':
    unittest.main(verbosity=2)
