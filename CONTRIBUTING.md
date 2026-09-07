# Contributing to Apollo

Apollo is a Discontinuous Galerkin solver used for plasma physics and CFD
research. The thing that makes it different from ordinary software is that a bug
here usually does not crash: it produces a number that is wrong. Most of what
follows is aimed at that.

## Getting set up

`README.md` has the full dependency list. On Ubuntu everything is apt-installable
and no source builds are needed. Then:

```bash
cd src
scons build-opt          # or build-debug
```

If a dependency is missing the build says which one and how to install it. Do not
commit a change to `src/this_host_config.py` that hardcodes a path on your own
machine: that file's values override everyone else's environment, and a stale
`petsc_base` there is why the build used to fail out of the box.

## Before you open a pull request

```bash
python3 -m unittest discover -s test    # ~1 s
test/run_examples.sh                    # ~15 s, needs a built binary
```

CI runs both, plus a syntax check over the SCons scripts and a full build of both
variants. Getting them green locally first is much faster than round-tripping.

## Things worth knowing

**Both build variants must agree numerically.** `build-debug` and `build-opt`
have to take the same code paths, or a discrepancy you reproduce under the
debugger is not the one production produced. This is why `USE_BLAS` is derived
from what the build actually found and applied to both. If you add a
`#ifdef`-gated numerical path, gate both variants the same way.

**Do not turn on `-ffast-math`.** It implies `-ffinite-math-only`, under which the
compiler may assume no value is ever NaN and folds `x != x` to `false`. That is
the idiom this codebase uses to detect a diverged solution, so the flag silently
deletes those guards. The `fastmath=yes` build variable exists for people who
have thought about this; it adds `-fno-finite-math-only` to keep the guards.

**Be careful with threads.** Apollo's parallelism is MPI. The one OpenMP loop
that exists is disabled by default because it calls PETSc per element and PETSc
is not thread-safe unless built `--with-threadsafety`. If you add an OpenMP
region: no PETSc calls inside it, and no scratch buffer that lives on the object
rather than on the stack. A shared `_meqn`-sized member array used as scratch is
exactly the bug that made 4-thread runs produce negative pressure on a smooth
problem.

**Prefer failing loudly.** A NaN that reaches the output file costs more than a
run that stops. If you add a validity check, make sure it survives optimization -
`std::isnan` is clearer than `x != x` and just as fast.

**Error paths run on every rank.** `exit()` from one rank leaves the others
blocked in a collective. If you add an error path in parallel code, report from
the rank that saw the problem, not just from rank 0, and prefer `MPI_Abort` over
`exit()`.

## Python

The scripts under `scripts/` are Python 3. A previous mechanical 2-to-3 pass left
every one of them broken - four would not even parse - so `test/` now pins the
specific failure modes: mangled `print(A) + B`, removed NumPy aliases,
Python-2-only builtins, mixed tab and space indentation, and undefined names.
Run the suite before committing; if it flags something in code you did not touch,
it has found a real bug.

## Style

There is no enforced formatter, so match the file you are editing. The codebase
mixes two prefixes - `wx` from WarpM, which Apollo descends from, and `ap` for
newer code. New files should use `ap`. Indent Python with spaces.

## Commit messages

Say what changed and why it was wrong before. For a numerical change, say how you
convinced yourself the new answer is the right one - which case you ran, and what
you compared against.
