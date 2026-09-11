#!/usr/bin/env bash
#
# SLURM job array for one Phase 3 folder: one task per deck, one core per task.
#
#     sbatch --array=0-8 slurm_array.sh 02-threshold-scan
#     sbatch --array=0-1 slurm_array.sh 01-c-sensitivity
#     sbatch --array=0-2 slurm_array.sh 00-divergence-gate
#
# The array index picks the deck, so set the range to match the folder: nine
# decks means --array=0-8. The script tells you the count if you get it wrong.
#
# ONE CORE PER TASK, AND THAT IS STILL THE POINT. The reason used to be that a
# frame held only 1/N of the cells on N ranks; that was a reader bug and it is
# fixed. The reason now is that the solver's answer depends on the rank count
# (known-issues.md 15), which the reader bug had been hiding. A scan is
# independent runs, so the throughput belongs at the array level regardless.
# See README.md in this folder.
#
# EDIT THESE BEFORE SUBMITTING. Partition names, account strings and module
# names are site-specific and no default can be right.

#SBATCH --job-name=apollo-rmf
#SBATCH --output=slurm-%x-%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=48:00:00
##SBATCH --partition=REPLACE_ME
##SBATCH --account=REPLACE_ME

set -euo pipefail

# --- site setup: pick ONE of these -----------------------------------------
#
# a) modules, if your site provides PETSc built with MPI and real scalars
# module load gcc openmpi petsc hdf5 gsl eigen boost python
#
# b) conda, from the environment.yml at the repository root
# source "$(conda info --base)/etc/profile.d/conda.sh"
# conda activate apollo
#
# Nothing is loaded by default, because loading the wrong PETSc is worse than
# loading none: the build will succeed and the run will not.
# ---------------------------------------------------------------------------

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
APOLLO="${APOLLO:-$(cd "$HERE/../../../../.." && pwd)}"
export APOLLO
export APOLLO_BIN="${APOLLO_BIN:-$APOLLO/src/build-opt/apollo}"
export APOLLO_RANKS=1

FOLDER="${1:-02-threshold-scan}"
[[ -d "$HERE/$FOLDER" ]] || { echo "no such folder: $HERE/$FOLDER" >&2; exit 2; }

# Sorted so the array index means the same thing on every task and every
# resubmission. Do not replace this with a glob in task order.
mapfile -t DECKS < <(find "$HERE/$FOLDER" -maxdepth 1 -name '*.pin' | sort)
N=${#DECKS[@]}

IDX="${SLURM_ARRAY_TASK_ID:-0}"
if (( IDX >= N )); then
    echo "array index $IDX but $FOLDER has only $N decks (use --array=0-$((N-1)))" >&2
    exit 2
fi

DECK="${DECKS[$IDX]}"
echo "task $IDX of $N: $(basename "$DECK")"
echo "host $(hostname), started $(date -Is)"

if [[ ! -x "$APOLLO_BIN" ]]; then
    echo "no solver at $APOLLO_BIN — build it before submitting:" >&2
    echo "    cd $APOLLO/src && scons build-opt" >&2
    exit 2
fi

# There is no checkpoint/restart, so a task that runs out of wall clock is a
# task whose work is gone. Say so where it will be read: in the log, at the top.
echo "NOTE: no checkpoint/restart. If this task hits its --time limit the run is"
echo "      lost, not resumable. Estimate first:"
echo "      python3 $APOLLO/scripts/rmf_scan.py $DECK --dry-run"

exec "$HERE/run_one.sh" "$DECK" "$HERE/results/$FOLDER/$(basename "$DECK" .pin)"
