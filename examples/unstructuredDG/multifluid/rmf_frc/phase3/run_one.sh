#!/usr/bin/env bash
#
# Run one Phase 3 deck in its own directory, then analyse it.
#
#     ./run_one.sh 00-divergence-gate/edge-cleaning-on.pin [results_dir]
#
# The directory matters: every run writes frames named after the deck's
# `Simulation` key, not after the file, so two runs sharing a directory
# overwrite each other's output without a word. This copies the deck and its
# mesh somewhere of their own first.
#
# ONE RANK, DELIBERATELY. Apollo writes about 1/N of the cells to each .vtu when
# run on N ranks, and which cells is up to the partitioner, so the analysis of a
# multi-rank run describes a fraction of the domain. Set APOLLO_RANKS to
# override if you are only after timing; the diagnostics will warn you.
#
# Environment:
#   APOLLO        checkout root      (default: five levels up from this script)
#   APOLLO_BIN    solver binary      (default: $APOLLO/src/build-opt/apollo)
#   APOLLO_RANKS  MPI ranks          (default: 1 — read the paragraph above)

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
APOLLO="${APOLLO:-$(cd "$HERE/../../../../.." && pwd)}"
APOLLO_BIN="${APOLLO_BIN:-$APOLLO/src/build-opt/apollo}"
RANKS="${APOLLO_RANKS:-1}"

if [[ $# -lt 1 ]]; then
    sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'
    exit 2
fi

DECK="$1"
[[ -f "$DECK" ]] || DECK="$HERE/$1"
if [[ ! -f "$DECK" ]]; then
    echo "no such deck: $1" >&2
    exit 2
fi
DECK="$(cd "$(dirname "$DECK")" && pwd)/$(basename "$DECK")"
NAME="$(basename "$DECK" .pin)"
SRC_DIR="$(dirname "$DECK")"

RESULTS="${2:-$HERE/results/$(basename "$SRC_DIR")/$NAME}"
mkdir -p "$RESULTS"

if [[ ! -x "$APOLLO_BIN" ]]; then
    echo "no solver at $APOLLO_BIN" >&2
    echo "build it first:  cd $APOLLO/src && scons build-opt" >&2
    exit 2
fi

# The mesh the deck names, not every mesh in the folder: 00-divergence-gate
# holds two, and copying both would leave it ambiguous which was used.
MESH="$(sed -n "s/.*Gridname *= *'\\([^']*\\)'.*/\\1/p" "$DECK" | head -1)"
if [[ -z "$MESH" || ! -f "$SRC_DIR/$MESH" ]]; then
    echo "cannot find the mesh '$MESH' named by $DECK" >&2
    exit 2
fi

cp "$DECK" "$SRC_DIR/$MESH" "$RESULTS/"
cd "$RESULTS"

echo "=== $NAME"
echo "    deck    $DECK"
echo "    mesh    $MESH"
echo "    results $RESULTS"
echo "    ranks   $RANKS"

PYTHONPATH="$APOLLO/scripts" python3 "$APOLLO/scripts/wxinpparse.py" -i "$NAME.pin"

started=$(date +%s)
if [[ "$RANKS" -gt 1 ]]; then
    echo "    WARNING: $RANKS ranks. Each .vtu will hold roughly 1/$RANKS of the"
    echo "             cells; the diagnostics below will say so. Timing only."
    # shellcheck disable=SC2086
    mpirun ${APOLLO_MPI_FLAGS:---oversubscribe} -np "$RANKS" \
        "$APOLLO_BIN" -i "$NAME.inp" 2>&1 | tee solver.log
else
    "$APOLLO_BIN" -i "$NAME.inp" 2>&1 | tee solver.log
fi
elapsed=$(( $(date +%s) - started ))

steps=$(grep -c 'Simulation dt' solver.log || true)
frames=$(ls ./*.vtu 2>/dev/null | wc -l)
echo "    done in ${elapsed}s: $steps steps, $frames frames"

echo
echo "=== diagnostics"
python3 "$APOLLO/scripts/rmf_diagnostics.py" "$NAME.pin" ./*.vtu | tee diagnostics.txt
