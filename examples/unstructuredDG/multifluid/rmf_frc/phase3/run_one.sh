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
# ONE RANK, DELIBERATELY. The output of a multi-rank run is complete - the old
# claim that a frame held 1/N of the cells was a bug in test/vtu.py, which kept
# only the last of PETSc's per-rank <Piece> elements, and that is fixed. What
# fixing it revealed is worse: the solver's answer depends on the rank count
# (known-issues.md 15). One rank and two give different fields on this deck and
# on the Maxwell pulse alike. APOLLO_RANKS still works and is fine for timing,
# but do not analyse what it produces.
#
# Environment:
#   APOLLO        checkout root      (default: five levels up from this script)
#   APOLLO_BIN    solver binary      (default: $APOLLO/src/build-opt/apollo)
#   APOLLO_RANKS  MPI ranks          (default: 1 — read the paragraph above)
#   APOLLO_RESUME 1 = continue from a checkpoint in the results directory if
#                 there is one, 0 = start over and delete it. Unset, a run that
#                 would overwrite an unfinished checkpoint is REFUSED.
#
# RESUMING. The solver writes <run>.checkpoint at every output frame and takes
# `-r <file>` to continue from it (docs/known-issues.md 6). Before that existed
# this script had nothing to lose by re-running into a populated directory; now
# it does, because a fresh run overwrites frame 0 and the checkpoint with it.
# That is why the default is to refuse rather than to pick one for you: the
# thing at risk is hours of compute, and both answers are one flag away.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
APOLLO="${APOLLO:-$(cd "$HERE/../../../../.." && pwd)}"
APOLLO_BIN="${APOLLO_BIN:-$APOLLO/src/build-opt/apollo}"
RANKS="${APOLLO_RANKS:-1}"

if [[ $# -lt 1 ]]; then
    # The whole comment header, however long it grows. It used to be `sed -n
    # '2,20p'`, a line range, and the header outgrew it: adding the paragraph
    # on resuming pushed the list of environment variables past line 20, so
    # `./run_one.sh` with no arguments stopped printing the one part of the
    # header that is an instruction rather than a caveat.
    awk 'NR > 1 && /^#/ { sub(/^# ?/, ""); print; next } NR > 1 { exit }' "$0"
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

# How to say "this same command again" in the messages below. Without $2 the
# advice would send a run whose results directory was given explicitly back to
# the default one, which is a different directory and usually an empty one.
SELF="$0 $1"
if [[ $# -ge 2 ]]; then SELF="$0 $1 $2"; fi

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

# What is already in the results directory, before anything is copied over it.
CKPT="$RESULTS/$NAME.checkpoint"
META="$CKPT.meta"
CK_FRAME=""; CK_NOUT=""; CK_TIME=""
if [[ -f "$CKPT" && -f "$META" ]]; then
    CK_FRAME="$(awk '$1=="frame"{print $2}' "$META")"
    CK_NOUT="$(awk '$1=="nout"{print $2}' "$META")"
    CK_TIME="$(awk '$1=="time"{print $2}' "$META")"
fi

RESUME="${APOLLO_RESUME:-}"
DO_RESUME=0
if [[ -n "$CK_FRAME" ]]; then
    if [[ "$CK_FRAME" -ge "$CK_NOUT" ]]; then
        # Finished. Resuming is refused by the solver itself; re-running would
        # silently redo however many hours this took.
        if [[ "$RESUME" == "0" ]]; then
            echo "the run here is COMPLETE ($CK_FRAME of $CK_NOUT frames);"
            echo "APOLLO_RESUME=0, so running it again from the start."
        elif [[ "$RESUME" == "1" ]]; then
            echo "=== $NAME: already complete ($CK_FRAME of $CK_NOUT frames), nothing to resume"
            cd "$RESULTS"
            echo
            echo "=== diagnostics"
            python3 "$APOLLO/scripts/rmf_diagnostics.py" "$NAME.pin" ./*.vtu \
                | tee diagnostics.txt
            exit 0
        else
            echo "the run in $RESULTS is COMPLETE: $CK_FRAME of $CK_NOUT frames." >&2
            echo "Starting over would overwrite it. Choose:" >&2
            echo "    APOLLO_RESUME=1 $SELF     # keep it, just re-analyse" >&2
            echo "    APOLLO_RESUME=0 $SELF     # discard it and run again" >&2
            exit 3
        fi
    elif [[ "$RESUME" == "1" ]]; then
        DO_RESUME=1
    elif [[ "$RESUME" == "0" ]]; then
        echo "discarding the checkpoint at frame $CK_FRAME of $CK_NOUT (APOLLO_RESUME=0)"
    else
        echo "an UNFINISHED run is in $RESULTS: frame $CK_FRAME of $CK_NOUT, t = $CK_TIME." >&2
        echo "Starting over would overwrite its checkpoint and lose it. Choose:" >&2
        echo "    APOLLO_RESUME=1 $SELF     # continue it" >&2
        echo "    APOLLO_RESUME=0 $SELF     # discard it and start over" >&2
        exit 3
    fi
fi

if [[ "$DO_RESUME" -eq 1 ]]; then
    # Deliberately NOT re-copying the deck and mesh: the ones in the results
    # directory are what the checkpoint was made from, and a deck edited in the
    # source folder since would resume into a different problem. The solver
    # catches a changed Output_files and a changed state size; it cannot catch
    # a changed TEND or a changed source term.
    cd "$RESULTS"
else
    if [[ -n "$CK_FRAME" ]]; then
        # A shorter re-run would otherwise leave the old run's later frames
        # sitting in the directory, and ./*.vtu below would analyse both.
        # $RESULTS-relative, not ./-relative: this runs BEFORE the cd below,
        # so `./` here is whatever directory the caller happened to be in.
        rm -f "$CKPT" "$META" "$RESULTS/$NAME"_*.vtu
    fi
    cp "$DECK" "$SRC_DIR/$MESH" "$RESULTS/"
    cd "$RESULTS"
    # tee appends below, so that a resumed run keeps the log of the half that
    # came before it. A fresh run has to start the file over, or the step count
    # reported at the end counts a previous run's steps too.
    : > solver.log
fi

echo "=== $NAME"
echo "    deck    $DECK"
echo "    mesh    $MESH"
echo "    results $RESULTS"
echo "    ranks   $RANKS"
if [[ "$DO_RESUME" -eq 1 ]]; then
    echo "    resume  frame $CK_FRAME of $CK_NOUT, t = $CK_TIME"
fi

PYTHONPATH="$APOLLO/scripts" python3 "$APOLLO/scripts/wxinpparse.py" -i "$NAME.pin"

# Expanded below as ${RESUME_ARGS[@]+"${RESUME_ARGS[@]}"}: an EMPTY array
# expanded as "${arr[@]}" is an unbound variable under `set -u` on bash before
# 4.4, which is what a CentOS 7 login node still ships. Every non-resumed run
# would have died on the line that runs the solver.
RESUME_ARGS=()
if [[ "$DO_RESUME" -eq 1 ]]; then
    RESUME_ARGS=(-r "$NAME.checkpoint")
fi

started=$(date +%s)
status=0
if [[ "$RANKS" -gt 1 ]]; then
    # This warning used to say each .vtu would hold roughly 1/$RANKS of the
    # cells. That was wrong - a bug in test/vtu.py, which kept only the last
    # rank's <Piece> - and the header of this file has said so since it was
    # fixed, while these two lines went on printing the retracted claim to
    # anyone who actually ran it. The real reason is worse; see known-issues 15.
    echo "    WARNING: $RANKS ranks. The output will be COMPLETE - but the"
    echo "             solver's answer depends on the rank count, so it will"
    echo "             not match the one-rank answer. Timing only; do not"
    echo "             analyse these numbers. See docs/known-issues.md 15."
    # shellcheck disable=SC2086
    mpirun ${APOLLO_MPI_FLAGS:---oversubscribe} -np "$RANKS" \
        "$APOLLO_BIN" -i "$NAME.inp" ${RESUME_ARGS[@]+"${RESUME_ARGS[@]}"} 2>&1 \
        | tee -a solver.log || status=$?
else
    "$APOLLO_BIN" -i "$NAME.inp" ${RESUME_ARGS[@]+"${RESUME_ARGS[@]}"} 2>&1 \
        | tee -a solver.log || status=$?
fi
elapsed=$(( $(date +%s) - started ))

steps=$(grep -c 'Simulation dt' solver.log || true)
# `find`, not `ls ./*.vtu`: with no frames at all the glob does not match, ls
# exits non-zero, and under `set -e` with pipefail the ASSIGNMENT carries that
# status and ends the script - exit 2, no message, right where the post-mortem
# below was about to explain what happened. A run that dies before its first
# output frame is exactly the run that needs the explanation.
frames=$(find . -maxdepth 1 -name '*.vtu' | wc -l)
echo "    done in ${elapsed}s: $steps steps, $frames frames"

# THE POST-MORTEM IS THE POINT OF THIS BRANCH. `set -e` used to end the script
# on the line above the moment the solver returned non-zero - so the one run
# that most needed analysing, the one that died, was the one that got none.
# Somebody had to go and run the tools by hand afterwards.
if [[ "$status" -ne 0 ]]; then
    echo
    # Two different endings, and the advice is opposite. A status of 128+N is a
    # SIGNAL: the wall clock, a Ctrl-C, a disconnected session - the run was
    # interrupted and resuming is exactly right. Anything below 128 is the
    # solver deciding to stop: the NaN guard at wxnodaldg2dmethod.cc:530 calls
    # exit(1), and a deck it cannot parse dies the same way. Resuming into that
    # arrives back where it left off, which is where the trouble was.
    if [[ "$status" -ge 128 ]]; then
        echo "=== the run was INTERRUPTED (signal $((status - 128))) after ${elapsed}s"
    else
        echo "=== the run DIED (exit $status) after ${elapsed}s"
    fi
    if [[ -f "$NAME.checkpoint.meta" ]]; then
        f="$(awk '$1=="frame"{print $2}' "$NAME.checkpoint.meta")"
        n="$(awk '$1=="nout"{print $2}' "$NAME.checkpoint.meta")"
        t="$(awk '$1=="time"{print $2}' "$NAME.checkpoint.meta")"
        echo "    a checkpoint is on disk: frame $f of $n, t = $t"
        echo "    continue it with:  APOLLO_RESUME=1 $SELF"
        if [[ "$status" -lt 128 ]]; then
            echo "    but the solver stopped on its own, so a resume will run"
            echo "    back into whatever stopped it unless something changes."
        fi
    else
        echo "    no checkpoint was written, so there is nothing to resume from"
    fi
    if [[ "$frames" -gt 1 ]]; then
        # run_growth, not rmf_diagnostics: the physics table bins over
        # r = 0..a and can read healthy to the last frame of a run that went
        # unstable somewhere its bins do not cover. This asks which of the 18
        # components was growing instead.
        echo
        echo "=== what was growing (scripts/run_growth.py)"
        python3 "$APOLLO/scripts/run_growth.py" ./*.vtu | tee growth.txt || true
    fi
    exit "$status"
fi

echo
echo "=== diagnostics"
python3 "$APOLLO/scripts/rmf_diagnostics.py" "$NAME.pin" ./*.vtu | tee diagnostics.txt
