##
# Find HDF5.
#
# Apollo does not include <hdf5.h> itself, but PETSc's petscviewerhdf5.h does,
# and solvers/apsolver.h pulls that in transitively. So any PETSc built with
# HDF5 support - which includes the distribution packages - makes HDF5's headers
# a hard requirement for compiling Apollo, even though nothing here calls the
# HDF5 API directly. This module locates them.
#
# The layout is awkward on Debian and Ubuntu: the MPI build of HDF5 installs its
# header to /usr/include/hdf5/openmpi rather than /usr/include, so the compiler
# does not find it without help.
##

Import('warpMConstructionEnv')
import os
import subprocess
import sys


def _candidateIncludeDirs(env):
    """Directories that might hold hdf5.h, most specific first."""
    dirs = []

    # 1. Explicit configuration always wins.
    if env['hdf5_base']:
        dirs.append(os.path.join(env['hdf5_base'], 'include'))

    # 2. The environment variable the HDF5 tools themselves use.
    for var in ('HDF5_DIR', 'HDF5_ROOT'):
        if os.environ.get(var):
            dirs.append(os.path.join(os.environ[var], 'include'))

    # 3. Ask the parallel HDF5 compiler wrapper, if there is one. `h5pcc -show`
    #    prints the full compile line, from which the -I flags can be read.
    for wrapper in ('h5pcc', 'h5cc'):
        try:
            out = subprocess.run([wrapper, '-show'], capture_output=True,
                                 text=True, timeout=10)
        except (OSError, subprocess.SubprocessError):
            continue
        if out.returncode == 0:
            for tok in out.stdout.split():
                if tok.startswith('-I'):
                    dirs.append(tok[2:])

    # 4. pkg-config, under the several names distributions use.
    for pkg in ('hdf5-openmpi', 'hdf5-mpich', 'hdf5-serial', 'hdf5'):
        try:
            out = subprocess.run(['pkg-config', '--cflags-only-I', pkg],
                                 capture_output=True, text=True, timeout=10)
        except (OSError, subprocess.SubprocessError):
            break
        if out.returncode == 0:
            for tok in out.stdout.split():
                if tok.startswith('-I'):
                    dirs.append(tok[2:])

    # 5. Well-known locations. MPI builds come first: Apollo is an MPI program,
    #    and mixing a serial hdf5.h with a parallel libhdf5 breaks at link time.
    dirs += [
        '/usr/include/hdf5/openmpi',
        '/usr/include/hdf5/mpich',
        '/usr/include/hdf5/serial',
        '/usr/lib64/openmpi/include',
    ]
    if sys.platform == 'darwin':
        for prefix in ('/opt/homebrew', '/usr/local'):
            dirs.append(os.path.join(prefix, 'opt', 'hdf5-mpi', 'include'))
            dirs.append(os.path.join(prefix, 'opt', 'hdf5', 'include'))

    # Preserve order, drop duplicates and anything that is not a directory.
    seen = set()
    out = []
    for d in dirs:
        if d and d not in seen and os.path.isdir(d):
            seen.add(d)
            out.append(d)
    return out


conf = Configure(warpMConstructionEnv)

if conf.CheckCXXHeader('hdf5.h'):
    print("HDF5: hdf5.h already on the default include path")
else:
    found = None
    for d in _candidateIncludeDirs(warpMConstructionEnv):
        if not os.path.isfile(os.path.join(d, 'hdf5.h')):
            continue
        conf.env.AppendUnique(CPPPATH=[d])
        if conf.CheckCXXHeader('hdf5.h'):
            found = d
            break
        # Not usable after all; do not leave a bad -I behind.
        conf.env['CPPPATH'] = [p for p in conf.env['CPPPATH'] if p != d]

    if found:
        print("HDF5: found hdf5.h in", found)
        libDir = os.path.join(os.path.dirname(os.path.dirname(found)), 'lib')
        if os.path.isdir(libDir):
            conf.env.AppendUnique(LIBPATH=[libDir])
    else:
        # Not fatal: a PETSc built without HDF5 support never includes it.
        print("")
        print("WARNING: hdf5.h was not found.")
        print("  Apollo itself does not need it, but PETSc's petscviewerhdf5.h")
        print("  includes it, so the build will fail with")
        print("    fatal error: hdf5.h: No such file or directory")
        print("  if your PETSc was built with HDF5 support (the distribution")
        print("  packages are). Install the parallel development package:")
        print("    Debian/Ubuntu  sudo apt-get install libhdf5-openmpi-dev")
        print("    RHEL/Fedora    sudo dnf install hdf5-openmpi-devel")
        print("    macOS          brew install hdf5-mpi")
        print("  or point the build at an existing install:")
        print("    scons build-opt hdf5_base=/path/to/hdf5")
        print("")

warpMConstructionEnv = conf.Finish()
