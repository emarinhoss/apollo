#!/usr/bin/env bash
#
# Build Apollo in a Google Colab session (or any Ubuntu box with sudo).
#
# In a Colab cell:
#
#     !git clone https://github.com/emarinhoss/apollo.git
#     !bash apollo/examples/unstructuredDG/multifluid/rmf_frc/phase3/colab_setup.sh
#
# Takes a few minutes, most of it apt. The build is 85 objects and took 42 s on
# four cores; Colab's free tier gives two, so expect roughly ninety seconds.
#
# WHAT COLAB CAN AND CANNOT DO HERE
#
# It can run 00-divergence-gate, which is the gate on everything else and takes
# a few minutes. It can run short probes, measure your own seconds-per-step for
# the cost model, and do all of the analysis.
#
# It cannot run the campaign. The shortest long run in these folders is 35 hours
# on one core, against a Colab session that is capped around 12 and disconnects
# sooner if the tab is idle. There is no checkpoint/restart, so a disconnect is
# not a pause, it is a loss. Use a cluster for 01, 02 and 03 — slurm_array.sh is
# a starting point.

set -euo pipefail

SUDO=""
if [[ "$(id -u)" -ne 0 ]]; then SUDO="sudo"; fi

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
APOLLO="${APOLLO:-$(cd "$HERE/../../../../.." && pwd)}"

echo "=== system packages"
$SUDO apt-get update -qq
# The same set README.md lists. PETSc comes from the distribution, so no source
# build - but note the version differs: 22.04 gives 3.15, 24.04 gives 3.19. Both
# work; the version is reported below because it decides which compatibility
# shims apply.
$SUDO apt-get install -y -qq --no-install-recommends \
    build-essential \
    libopenmpi-dev openmpi-bin \
    libpetsc-real-dev \
    libhdf5-openmpi-dev \
    libboost-dev \
    libeigen3-dev \
    libgsl-dev \
    libopenblas-dev

# Which PETSc this actually got. Apollo builds against 3.11 and up, but the
# version varies by distribution - Ubuntu 22.04 (what Colab runs) ships 3.15,
# 24.04 ships 3.19 - and the compatibility shims in src/lib/petsc_compat.h
# switch on it. Printing it costs nothing and turns "300 lines of template
# errors" into "you have version X", which is the difference between a
# five-minute problem and an afternoon.
export PETSC_DIR="${PETSC_DIR:-/usr/lib/petsc}"
PETSC_VERSION="$(sed -n 's/#define PETSC_VERSION_\(MAJOR\|MINOR\|SUBMINOR\) *//p' \
    "$PETSC_DIR/include/petscversion.h" 2>/dev/null | paste -sd. - || true)"
echo "=== petsc ${PETSC_VERSION:-UNKNOWN} in $PETSC_DIR"

PETSC_MINOR="$(printf '%s' "$PETSC_VERSION" | cut -d. -f2)"
if [[ -z "$PETSC_MINOR" ]]; then
    echo "WARNING: could not read a PETSc version from $PETSC_DIR." >&2
    echo "         If the build fails below, that is the first thing to check." >&2
elif (( PETSC_MINOR < 11 )); then
    echo "PETSc $PETSC_VERSION is older than 3.11, which is the oldest release" >&2
    echo "the compatibility header is tested against (test/test_petsc_compat.py)." >&2
    echo "The build below will probably fail. Use a newer PETSc, or the conda" >&2
    echo "environment in environment.yml, which pins 3.19." >&2
    exit 1
fi

echo "=== python packages"
python3 -m pip install -q -r "$APOLLO/requirements.txt"

echo "=== build"
cd "$APOLLO/src"
# -j$(nproc) rather than a fixed number: Colab gives two cores, a cluster login
# node may give many, and scons will happily use whatever it is told.
scons build-opt -j"$(nproc)"

echo
echo "=== check"
"$APOLLO/src/build-opt/apollo" --help >/dev/null 2>&1 || true
if [[ -x "$APOLLO/src/build-opt/apollo" ]]; then
    echo "solver:  $APOLLO/src/build-opt/apollo"
else
    echo "BUILD FAILED: no binary at $APOLLO/src/build-opt/apollo" >&2
    exit 1
fi

# The solver-free tests. They take under a second and they check the analysis
# path this session will actually use, so a failure here is worth seeing before
# any run rather than after one.
python3 "$APOLLO/test/test_rmf_diagnostics.py" 2>&1 | tail -3
python3 "$APOLLO/test/test_rmf_scan.py" 2>&1 | tail -3

cat <<EOF

=== next

Measure this machine, so the cost estimates mean something here:

    cd $HERE
    ./run_one.sh 00-divergence-gate/edge-cleaning-on.pin

Then the gate, which is the whole reason to start with 00:

    ./run_one.sh 00-divergence-gate/edge-cleaning-off.pin
    python3 $APOLLO/scripts/rmf_cleaning_check.py \\
        results/00-divergence-gate/edge-cleaning-on \\
        results/00-divergence-gate/edge-cleaning-off

Read $HERE/README.md before starting anything longer.
EOF
