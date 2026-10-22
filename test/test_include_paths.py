#!/usr/bin/env python3
"""Does the build avoid handing the compiler its own default search paths?

`-I/usr/include` is redundant when it is right and fatal when it is wrong: -I
puts a directory AHEAD of the toolchain's own sysroot. A conda-forge GCC 15
makes _Float128, _Float32 and _Float64 built-in types and ships a matching
glibc; an older system glibc declares them as typedefs. Both are fine alone.
Together, with the system copy promoted, every file that reaches <cmath> dies:

    /usr/include/bits/floatn.h:86:20: error: redeclaration of C++ built-in
                                      type '_Float128' [-fpermissive]

Apollo did that to itself: this_host_config.py sets mpi_base = "/usr" and
config_mpi.py appends $mpi_base/include to CPPPATH. On a distribution build the
path is merely redundant, which is why it went unnoticed; on any toolchain
carrying its own libc headers - conda, spack, a cross-compiler - it is fatal.

The same argument applies to -L one stage later, where the consequence is worse:
linking the distribution's library against conda's headers is an ABI mismatch
that surfaces at run time, not build time.
"""

import os
import sys
import unittest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'src'))

from buildconf.include_paths import (DEFAULT_SYSTEM_INCLUDES,
                                     DEFAULT_SYSTEM_LIBS,
                                     drop_default_system_includes,
                                     drop_default_system_libs)


class TestIncludePaths(unittest.TestCase):

    def test_usr_include_is_dropped(self):
        """The exact case that broke the conda build."""
        kept, dropped = drop_default_system_includes(
            ['/usr/include', '/usr/include/hdf5/openmpi', '/opt/petsc/include'])
        self.assertEqual(dropped, ['/usr/include'])
        self.assertEqual(kept, ['/usr/include/hdf5/openmpi', '/opt/petsc/include'])

    def test_subdirectories_are_kept(self):
        """/usr/include/hdf5/openmpi is a real location, not a default one.

        Dropping it would break the HDF5 include that config_hdf5.py finds, so
        this is the check that stops the fix from becoming a worse bug.
        """
        kept, dropped = drop_default_system_includes(
            ['/usr/include/hdf5/openmpi', '/usr/include/eigen3'])
        self.assertEqual(dropped, [])
        self.assertEqual(len(kept), 2)

    def test_trailing_slash_is_recognised(self):
        _kept, dropped = drop_default_system_includes(['/usr/include/'])
        self.assertEqual(dropped, ['/usr/include/'])

    def test_order_is_preserved(self):
        paths = ['/a', '/usr/include', '/b', '/usr/local/include', '/c']
        kept, _dropped = drop_default_system_includes(paths)
        self.assertEqual(kept, ['/a', '/b', '/c'])

    def test_conda_prefix_is_never_dropped(self):
        """The environment's own include directory is the one that must survive."""
        kept, dropped = drop_default_system_includes(
            ['/home/u/miniforge3/envs/apollo/include', '/usr/include'])
        self.assertEqual(kept, ['/home/u/miniforge3/envs/apollo/include'])
        self.assertEqual(dropped, ['/usr/include'])

    def test_lib_defaults_are_dropped(self):
        kept, dropped = drop_default_system_libs(
            ['/usr/lib', '/usr/lib/petsc/lib', '/env/lib'])
        self.assertEqual(dropped, ['/usr/lib'])
        self.assertEqual(kept, ['/usr/lib/petsc/lib', '/env/lib'])

    def test_usr_local_lib_is_kept(self):
        """Homebrew and hand-built dependencies live there legitimately."""
        kept, dropped = drop_default_system_libs(['/usr/local/lib'])
        self.assertEqual(dropped, [])
        self.assertEqual(kept, ['/usr/local/lib'])

    def test_nothing_to_drop_leaves_everything(self):
        paths = ['/env/include', '/opt/x/include']
        kept, dropped = drop_default_system_includes(paths)
        self.assertEqual(kept, paths)
        self.assertEqual(dropped, [])

    def test_the_two_tables_are_disjoint_in_intent(self):
        """An include default must not be treated as a library default."""
        for path in DEFAULT_SYSTEM_INCLUDES:
            self.assertNotIn(path, DEFAULT_SYSTEM_LIBS)


if __name__ == '__main__':
    unittest.main(verbosity=2)
