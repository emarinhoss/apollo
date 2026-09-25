#!/usr/bin/env python3
"""Does the build look for Eigen where this repository tells people to put it?

`environment.yml` in this repository tells people to build inside a conda
environment and installs Eigen into it. SConstruct's candidate list covered
Debian, Fedora, macOS and an explicit override — and not conda. So the
documented path failed with

    hyperapps/multifluid/MathX.h:20:10: fatal error: Eigen/Dense: No such file

on a cluster, in the middle of several hundred lines of unrelated warnings.

The reason it survived is that SConstruct is not importable, so nothing about
it can be checked except by running a build on a machine that reproduces the
case. The decision now lives in buildconf/eigen_paths.py as a pure function of
the environment, which is what these tests exercise.
"""

import os
import sys
import unittest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'src'))

from buildconf.eigen_paths import eigen_candidates


class TestEigenCandidates(unittest.TestCase):

    def test_conda_prefix_is_searched(self):
        """The case that failed: a conda environment, no system Eigen."""
        got = eigen_candidates({'CONDA_PREFIX': '/home/u/miniforge3/envs/apollo'},
                               'linux')
        self.assertIn('/home/u/miniforge3/envs/apollo/include/eigen3', got)

    def test_conda_is_preferred_over_system(self):
        """Inside an environment, its own Eigen is the one that matches its compiler."""
        got = eigen_candidates({'CONDA_PREFIX': '/env'}, 'linux')
        self.assertLess(got.index('/env/include/eigen3'),
                        got.index('/usr/include/eigen3'))

    def test_explicit_override_wins(self):
        got = eigen_candidates({'EIGEN3_INCLUDE_DIR': '/opt/eigen',
                                'CONDA_PREFIX': '/env'}, 'linux')
        self.assertEqual(got[0], '/opt/eigen')

    def test_conda_build_prefixes(self):
        """conda-build sets PREFIX and BUILD_PREFIX rather than CONDA_PREFIX."""
        got = eigen_candidates({'PREFIX': '/p', 'BUILD_PREFIX': '/bp'}, 'linux')
        self.assertIn('/p/include/eigen3', got)
        self.assertIn('/bp/include/eigen3', got)

    def test_system_paths_still_present(self):
        """The distribution packages must keep working - that is how CI builds."""
        got = eigen_candidates({}, 'linux')
        self.assertEqual(got, ['/usr/include/eigen3', '/usr/local/include/eigen3'])

    def test_homebrew_only_on_macos(self):
        mac = eigen_candidates({}, 'darwin')
        self.assertIn('/opt/homebrew/include/eigen3', mac)
        self.assertNotIn('/opt/homebrew/include/eigen3',
                         eigen_candidates({}, 'linux'))

    def test_no_duplicates(self):
        """/usr/local appears twice on macOS if nothing dedupes."""
        got = eigen_candidates({'CONDA_PREFIX': '/usr'}, 'darwin')
        self.assertEqual(len(got), len(set(got)), got)

    def test_paths_point_at_the_parent_of_Eigen(self):
        """Every candidate must end in eigen3, not include.

        `#include <Eigen/Dense>` resolves against the directory CONTAINING
        Eigen, and Eigen ships one level down from include. Handing the build
        $CONDA_PREFIX/include - which is already on conda's CPATH, and is why
        PETSc and Boost are found without help - would not fix anything.
        """
        for path in eigen_candidates({'CONDA_PREFIX': '/env', 'PREFIX': '/p'},
                                     'darwin'):
            self.assertTrue(path.endswith('eigen3'), path)

    def test_real_environment_does_not_crash(self):
        self.assertTrue(eigen_candidates())


if __name__ == '__main__':
    unittest.main(verbosity=2)
