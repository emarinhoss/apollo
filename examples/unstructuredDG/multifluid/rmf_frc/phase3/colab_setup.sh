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

echo "=== system packages"
$SUDO apt-get update -qq
# The same set README.md lists, which is enough on Ubuntu 22.04 and 24.04. PETSc
# comes from the distribution: no source build.
$SUDO apt-get install -y -qq --no-install-recommends \
    build-essential \
    libopenmpi-dev openmpi-bin \
    libpetsc-real-dev \
    libhdf5-openmpi-dev \
    libboost-dev \
    libeigen3-dev \
    libgsl-dev \
    libopenblas-dev

echo "=== python packages"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
APOLLO="${APOLLO:-$(cd "$HERE/../../../../.." && pwd)}"
python3 -m pip install -q -r "$APOLLO/requirements.txt"

echo "=== build"
export PETSC_DIR="${PETSC_DIR:-/usr/lib/petsc}"
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
