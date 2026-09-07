#!/usr/bin/env bash
#
# Smoke-test Apollo against its own bundled examples.
#
# Each case is preprocessed (.pin -> .inp, which is what the solver actually
# reads) and run to completion in a scratch directory. A case passes if the
# solver exits 0 and prints no NaN, negative-pressure or exception diagnostic.
#
# This is deliberately cheap - the cases below finish in a few seconds each - so
# it can run on every push. It is the check that catches "the solver still
# starts and still integrates", which is exactly what a data race in the element
# loop broke without any compile error.
#
#   test/run_examples.sh [-b <apollo binary>] [-j <mpi ranks>] [case ...]
#
# With no case arguments every case in CASES is run.

set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
APOLLO="${APOLLO_BIN:-$REPO/src/build-opt/apollo}"
RANKS=1

# Cases are "<label>:<path to .pin relative to examples/>".
CASES=(
  "advection:unstructuredDG/advection/advection.pin"
  "euler-isentropic-vortex:unstructuredDG/euler/isentropicVortex/isentropicVortex.pin"
  "maxwell-circular-pulse:unstructuredDG/maxwell/circularPulse/circularPulse.pin"
)

while getopts ":b:j:" opt; do
  case "$opt" in
    b) APOLLO="$OPTARG" ;;
    j) RANKS="$OPTARG" ;;
    *) echo "usage: $0 [-b <apollo binary>] [-j <ranks>] [case ...]" >&2; exit 2 ;;
  esac
done
shift $((OPTIND - 1))

if [[ ! -x "$APOLLO" ]]; then
  echo "error: no Apollo binary at $APOLLO" >&2
  echo "       build one with 'cd src && scons build-opt', or pass -b <path>." >&2
  exit 2
fi

if [[ $# -gt 0 ]]; then
  wanted=("$@")
  selected=()
  for want in "${wanted[@]}"; do
    for case_spec in "${CASES[@]}"; do
      [[ "${case_spec%%:*}" == "$want" ]] && selected+=("$case_spec")
    done
  done
  if [[ ${#selected[@]} -eq 0 ]]; then
    echo "error: no case matched: $*" >&2
    printf 'known cases: %s\n' "$(printf '%s ' "${CASES[@]%%:*}")" >&2
    exit 2
  fi
  CASES=("${selected[@]}")
fi

# Extra mpirun flags. Open MPI refuses to run as root without being told to,
# which is the normal situation inside a container, and CI runners have fewer
# cores than the rank counts worth testing.
MPI_FLAGS=()
if [[ "$RANKS" -gt 1 ]]; then
  read -r -a MPI_FLAGS <<< "${APOLLO_MPI_FLAGS:-}"
  if [[ ${#MPI_FLAGS[@]} -eq 0 ]]; then
    MPI_FLAGS=(--oversubscribe)
    [[ "$(id -u)" -eq 0 ]] && MPI_FLAGS+=(--allow-run-as-root)
  fi
fi

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

# Diagnostics the solver prints when the solution has gone bad. It does not
# always exit non-zero on these, so grep for them as well as checking the status.
BAD_PATTERNS='NaN|NAN|Negative pressure|Negative density|Exception caught|not found'

pass=0
fail=0
failed_cases=()

for case_spec in "${CASES[@]}"; do
  label="${case_spec%%:*}"
  pin="$REPO/examples/${case_spec#*:}"
  dir="$WORK/$label"

  printf '%-28s ' "$label"

  if [[ ! -f "$pin" ]]; then
    echo "SKIP (no such deck: $pin)"
    continue
  fi

  mkdir -p "$dir"
  # The mesh files a deck names are looked up relative to the working directory.
  cp "$(dirname "$pin")"/*.pin "$(dirname "$pin")"/*.msh "$dir/" 2>/dev/null

  pin_name="$(basename "$pin")"
  inp_name="${pin_name%.pin}.inp"

  # .pin decks are Python-macro templates; the solver reads the expanded .inp.
  if ! (cd "$dir" && PYTHONPATH="$REPO/scripts" python3 "$REPO/scripts/wxinpparse.py" \
        -i "$pin_name" >preprocess.log 2>&1); then
    echo "FAIL (preprocessing)"
    sed 's/^/    /' "$dir/preprocess.log" | tail -15
    fail=$((fail + 1)); failed_cases+=("$label"); continue
  fi

  if [[ "$RANKS" -gt 1 ]]; then
    runner=(mpirun "${MPI_FLAGS[@]}" -np "$RANKS" "$APOLLO")
  else
    runner=("$APOLLO")
  fi

  start=$(date +%s)
  (cd "$dir" && "${runner[@]}" -i "$inp_name" >run.log 2>&1)
  status=$?
  elapsed=$(( $(date +%s) - start ))

  if [[ $status -ne 0 ]]; then
    echo "FAIL (exit $status, ${elapsed}s)"
    tail -15 "$dir/run.log" | sed 's/^/    /'
    fail=$((fail + 1)); failed_cases+=("$label"); continue
  fi

  if grep -qE "$BAD_PATTERNS" "$dir/run.log"; then
    echo "FAIL (bad solution, ${elapsed}s)"
    grep -nE "$BAD_PATTERNS" "$dir/run.log" | head -5 | sed 's/^/    /'
    fail=$((fail + 1)); failed_cases+=("$label"); continue
  fi

  echo "ok (${elapsed}s)"
  pass=$((pass + 1))
done

echo
echo "$pass passed, $fail failed"
if [[ $fail -gt 0 ]]; then
  printf 'failed: %s\n' "${failed_cases[*]}"
  exit 1
fi
