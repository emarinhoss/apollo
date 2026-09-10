##
# Add Blas libraries to build warpMConstructionEnvironment
##

Import('warpMConstructionEnv')
import os
import platform
import sys

# Add specified library paths to search paths for libs and headers
# Use system-detected paths instead of hardcoded absolute paths
if warpMConstructionEnv['blas_base'] == '':
    # Let SCons use default system library search paths
    # Don't set to absolute path - let the linker find it

    # On macOS, check for Homebrew OpenBLAS
    if sys.platform == 'darwin':
        homebrew_prefix = '/opt/homebrew'  # Apple Silicon
        if not os.path.isdir(homebrew_prefix):
            homebrew_prefix = '/usr/local'  # Intel Mac

        openblas_include = os.path.join(homebrew_prefix, 'opt', 'openblas', 'include')
        openblas_lib = os.path.join(homebrew_prefix, 'opt', 'openblas', 'lib')

        if os.path.isdir(openblas_include):
            warpMConstructionEnv.AppendUnique(CPPPATH=[openblas_include])
            print(f"Added OpenBLAS include path: {openblas_include}")
        if os.path.isdir(openblas_lib):
            warpMConstructionEnv.AppendUnique(LIBPATH=[openblas_lib])
            print(f"Added OpenBLAS library path: {openblas_lib}")
elif warpMConstructionEnv['blas_base'] != '':
    # User specified a custom BLAS location
    warpMConstructionEnv.AppendUnique(LIBPATH=(warpMConstructionEnv['blas_base']))

    # Add include directory
    blas_include = os.path.join(warpMConstructionEnv['blas_base'], 'include')
    if os.path.isdir(blas_include):
        warpMConstructionEnv.AppendUnique(CPPPATH=[blas_include])

    # Add architecture-specific subdirectory if it exists (Ubuntu/Debian convention)
    arch_subdir = platform.machine() + '-linux-gnu'
    arch_path = os.path.join(warpMConstructionEnv['blas_base'], arch_subdir)
    if os.path.isdir(arch_path):
        warpMConstructionEnv.AppendUnique(LIBPATH=(arch_path))

# Configure libraries using some autoconf like functionality of scons
conf = Configure(warpMConstructionEnv)
if warpMConstructionEnv['blas_base']:
    print("Blas_base set as", warpMConstructionEnv['blas_base'])
else:
    print("Using system default paths for BLAS")

# Try to find BLAS library (try multiple common names)
blas_found = False

# Try OpenBLAS first (most common optimized BLAS)
if conf.CheckLib('openblas'):
    print("Found OpenBLAS library")
    blas_found = True
# Try standard BLAS + CBLAS
elif conf.CheckLib('cblas') or conf.CheckLib('blas'):
    print("Found standard BLAS library")
    if conf.CheckLib('cblas'):
        print("Found CBLAS library")
    blas_found = True
# Try combined BlasLapack (PETSc-provided)
elif conf.CheckLib('BlasLapack'):
    print("Found BlasLapack library")
    blas_found = True

if not blas_found:
    print("WARNING: No BLAS library found. Code will use slower fallback implementation.")
    print("For best performance, install OpenBLAS: sudo apt-get install libopenblas-dev")
    # Don't exit - allow compilation without BLAS, will use fallback

# lib/wxcubature2d.cc includes <cblas.h> under #ifdef USE_BLAS. Having a BLAS
# *library* to link is not enough: the header has to be present too, or the
# build fails at compile time on a machine that only ships the runtime package.
if blas_found and not conf.CheckCHeader('cblas.h'):
    print("WARNING: a BLAS library was found but <cblas.h> was not; building "
          "without BLAS. Install the development package (Debian/Ubuntu: "
          "libopenblas-dev, RHEL/Fedora: openblas-devel).")
    blas_found = False

# Check for gfortran, but do NOT add -lgfortran to the link.
#
# A shared BLAS records its own dependency: libopenblas.so.0 carries
# NEEDED libgfortran.so.5, so the linker pulls the right one transitively and
# the explicit flag is redundant. It is also actively harmful, because -lgfortran
# resolves through the unversioned libgfortran.so development symlink, and which
# one that is depends on which lib directory the linker reaches first. In a
# conda environment the symlink lives with the compiler rather than in
# $CONDA_PREFIX/lib, so -lgfortran fell through to an old system copy and
# produced a binary needing BOTH:
#
#     libgfortran.so.5 => .../envs/apollo/lib/libgfortran.so.5   (transitive)
#     libgfortran.so.3 => not found                              (this flag)
#
# autoadd=0 keeps the probe - the note below is still worth printing - without
# putting the library on the link line. If a static BLAS ever does need it, the
# link fails loudly with undefined _gfortran_* symbols, which is a better
# failure than silently binding the wrong runtime.
if not conf.CheckLib('gfortran', autoadd=0):
    print("Note: gfortran not found. Only matters for a static BLAS; a shared "
          "one records its own dependency.")

warpMConstructionEnv = conf.Finish() # replace the environment with the one modified by Configure's auto-conf actions

# Record the result so SConstruct defines USE_BLAS only when the code can
# actually compile against it, and defines it identically for every build
# variant. Set after Finish() so it lands on the environment SConstruct holds.
warpMConstructionEnv['APOLLO_HAVE_BLAS'] = blas_found
