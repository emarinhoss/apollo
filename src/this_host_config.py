##
# Host-specific defaults for the Apollo build.
#
# Every name set here is an SCons build variable (see the vars.AddVariables call
# in SConstruct) and can be overridden on the command line, which always wins:
#
#     scons build-opt petsc_base=/opt/petsc mpi_base=/usr arch=x86-64-v3
#
# Leave a variable unset (or empty) to let SCons autodetect it. Detection order
# for PETSc is: this file -> $PETSC_DIR (with $PETSC_ARCH if set) -> the
# compiler's default search paths. Set values here only when autodetection
# cannot find your installation; do not commit a path that exists on just one
# machine, because it silently overrides everyone else's $PETSC_DIR.
##

import os

########################################
#  Common settings
########################################

# Base directory of the MPI installation. '/usr' suits a distro package
# (openmpi-bin / libopenmpi-dev on Debian and Ubuntu, openmpi-devel on RHEL).
# The build uses the mpicc/mpicxx wrappers from $mpi_base/bin.
mpi_base = "/usr"

# Base directory of the PETSc installation. Empty means "use $PETSC_DIR, or the
# system paths if that is unset". Uncomment and edit only if neither works:
# petsc_base = "/opt/petsc"
petsc_base = ""

# Base directory of BLAS and Boost. Empty means "use the system paths", which is
# right for distro packages and for Homebrew on macOS.
# blas_base = ""
# boost_base = ""

########################################
#  Optimized-build tuning
########################################

# -march=/-mtune= value for build-opt. 'native' tunes for this machine and the
# resulting binary may SIGILL on an older node; on a heterogeneous cluster set
# the oldest supported ISA (for example 'x86-64-v2') or 'none'.
# arch = "native"

# -ffast-math is off by default because it deletes Apollo's NaN guards; see the
# comment on the 'fastmath' variable in SConstruct before turning it on.
# fastmath = "no"
