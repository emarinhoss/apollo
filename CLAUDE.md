# Working on Apollo

Apollo is a discontinuous-Galerkin solver for hyperbolic PDEs on unstructured
triangular meshes, built on PETSc's DMPlex. Its largest application is a
two-fluid (five-moment) plasma model coupled to perfectly-hyperbolic Maxwell —
eighteen components: electrons 0–4, ions 5–9, **E** 10–12, **B** 13–15, and the
cleaning potentials φ and ψ at 16 and 17.

## Build and test

```bash
cd src && PETSC_DIR=/usr/lib/petsc scons build-opt   # 85 objects, ~40 s on 4 cores
PETSC_DIR=/usr/lib/petsc make -C ../test/cxx        # C++ unit tests, 86 checks
cd .. && python3 -m unittest discover -s test       # Python tests
```

**Mind the working directory.** Those three lines run in sequence, and only the
first is from `src/`. Get it wrong and neither failure points at the cause:

- `make -C test/cxx` after `cd src` says
  `make: *** test/cxx: No such file or directory.  Stop.` — naming a path that
  *does* exist in the repo, just not relative to `src/`. So the error reads as
  nonsense: you look, `test/cxx` is right there.
- `python3 -m unittest discover -s test` from `src/` is worse. It prints
  `Ran 0 tests in 0.000s` / `OK` and exits 0. A green run that tested nothing.

Both verified by running them.

**`PETSC_DIR` on the scons line is not optional**, and this file used to say it
was: "`PETSC_DIR` defaults to `/usr/lib/petsc`". It does not. The build looks at
`petsc_base`, then `$PETSC_DIR`, then the compiler's own search paths — and
Debian's `libpetsc-real-dev` installs under `/usr/lib/petsc`, which is not one
of them. Without it the build stops at

    Checking for C++ header file petsc.h... no
    ERROR: the PETSc header 'petsc.h' was not found.

which at least names the fix. Both CI build jobs set `PETSC_DIR` explicitly;
only this file assumed otherwise. Inside a conda environment it is
`$CONDA_PREFIX`. There is a `build-debug` variant.
Apollo does **not** use CMake. See `README.md` for dependencies,
`requirements.txt` for the Python side, `environment.yml` for a conda
environment on a cluster.

**238 Python tests, and 4 modules are 99% of the runtime.** Measured on a
4-core container with a `build-opt` binary present:

| | tests | time |
| --- | --- | --- |
| `test_rmf_antenna_field.py` | 3 | 645 s |
| `test_restart.py` | 8 | 260 s |
| `test_rmf_bc_field.py` | 4 | 215 s |
| `test_vortex_accuracy.py` | 2 | 20 s |
| everything else (17 modules) | 221 | **8 s** |
| whole suite | 238 | ~1150 s |

So run the 209 while you work and the whole suite before you push. The 209 need
no solver, and naming them is the only reliable way to select them — a glob is
not, because the fast and slow modules interleave alphabetically:

```bash
cd test && python3 -m unittest \
    test_rmf_diagnostics test_rmf_scan test_rmf_cleaning_check \
    test_mkdiscmesh test_python_tooling test_deck_preprocess \
    test_vtu_reader test_eigen_paths test_include_paths \
    test_conda_paths test_petsc_compat test_run_fingerprint \
    test_run_growth test_run_one test_where_peak \
    test_ring_spectrum test_winding_anatomy
# Ran 221 tests in 8s -- OK
```

`cd test` first. From the repository root the stdlib's own `test` package wins
the import, so `python3 -m unittest test.test_rmf_scan` dies with
`ModuleNotFoundError: No module named 'test.test_rmf_scan'` — pointing at a file
that is plainly sitting there. Use `discover -s test`, or `cd test` and name the
modules as above. `test_petsc_compat` is in the fast set; it needs a C++
compiler but no PETSc.

The 17 that do need a solver skip themselves without one, so **a green run does
not by itself mean they ran** — check the skip count. With a binary present the
suite reports `Ran 238 tests` and no skips.

## Things that will cost you a day if you do not know them

Most are in `docs/known-issues.md`, which is worth reading before trusting any
result. The ones that bite hardest while editing:

- **A deck value written without a decimal point aborts the run.** The parser
  types `3000000` as an integer, `get<REAL>` throws `std::bad_cast`, and the
  message names neither the key nor the file. **The mirror is just as fatal**:
  `Output_files` is read with `get<int>`, so `OUT = 6.0` dies the same way. Use
  `rmf_scan.deck_number()` and `deck_integer()` when generating decks. (§12)
- **Run one rank for anything whose numbers you will use — the solver's answer
  depends on the rank count.** (§15) This entry used to give a different reason:
  that MPI output was incomplete, one frame holding ~1/N of the cells. That was
  wrong. PETSc writes one `<Piece>` per rank and all of them are in the file;
  `test/vtu.py` flattened them so the last piece won, and the counts behind the
  old claim (7792, 3896, 1961) were that piece's size. The reader is fixed —
  a two-rank file reads back 4304 cells against the one-rank 4304 — and fixing
  it made the real problem visible for the first time: one rank and two ranks
  give materially different fields, on the Maxwell pulse as well as the RMF
  decks. A multi-rank run now produces a complete, plausible file that disagrees
  with the serial answer, which is worse than one that was obviously short.
- **The PETSc version decides whether it compiles at all.** `petsc_compat.h`
  substitutes for plex calls upstream removed at 3.13 and 3.14; the guards said
  3.18 for both, so Apollo built on 24.04 (PETSc 3.19) and would not compile on
  22.04 (PETSc 3.15) — including in Colab. Fixed, and pinned by
  `test/test_petsc_compat.py`, which compiles the header against a stub
  `<petsc.h>` for every release from 3.11 up. If you touch a version guard, run
  it: it fails in both directions.
- **The Python version decides whether a deck is preprocessed.** A `.pin`'s
  preamble is Python, and `wxinputparser.py` evaluates it to substitute names
  into the body. It used to `exec(co)` inside a method and recover each name
  with `eval(name)`, which works only because CPython ≤ 3.12 let exec() writes
  into `locals()` persist. PEP 667 ended that in 3.13: substitution silently
  stopped and the solver aborted with `unexpected symbol '/' on line 23`, naming
  neither key nor cause. Fixed by using an explicit namespace, and pinned by
  `test/test_deck_preprocess.py`, which also compares every `python3.x` on the
  machine — because the bug was invisible on the interpreter you happened to run.
  CI runs the Python gate on 3.11 and 3.13.
- **The `.vtu` layout is PETSc's, not Apollo's, so it changes with PETSc.**
  Output goes through `PETSCVIEWERVTK`, and the length prefix before each
  appended block is 4 bytes on PETSc 3.15 (no `header_type` attribute, so the
  VTK default UInt32) and 8 bytes on 3.19 (`header_type="UInt64"`).
  `test/vtu.py` assumed 8 and could not read a 3.15 run at all — it failed with
  numpy's "buffer size must be a multiple of element size", naming neither file
  nor cause, and every diagnostic in `scripts/` goes through it. It now reads
  the attribute and refuses with a diagnosable message when a block does not
  add up. Pinned by `test/test_vtu_reader.py`, which builds both layouts, and by
  a CI step that reads back a `.vtu` the job's own PETSc actually wrote — in
  **both** build jobs. It used to run only on the 22.04 job, so the sole check
  that the reader matches its writer was performed against exactly one PETSc,
  the oldest one; and that runner is being deprecated, which would have taken
  the check with it.
- **Restart exists now; use it on anything long.** A run writes
  `<runName>.checkpoint` at every output frame and `-r <file>` resumes from it,
  reproducing a straight-through run *exactly* — `0.000e+00` across all 378
  solution arrays, asserted by `test/test_restart.py`. It costs one rolling
  file (5.27 MB on the phase3 mesh), and it refuses rather than resuming into a
  different mesh, a different `Output_files`, or a finished run. This entry used
  to read "an interrupted run is a lost run", which cost someone 12 hours of a
  140-hour formation run. Frames written before the feature existed cannot be
  resumed from. For the phase3 folders, `APOLLO_RESUME=1 ./run_one.sh <deck>`
  continues a run and plain `./run_one.sh` **refuses** to overwrite an
  unfinished one; `slurm_array.sh` sets it, so resubmitting an array after a
  wall-clock kill picks every task up where it stopped. (§6)
- **The two-fluid slope limiter produces NaN** after a few tens of steps, so no
  multifluid case can be limited. (§2)
- **On PETSc 3.22 and newer the Euler limiter silently does nothing.** The
  limiters write back into the TS solution vector (`applyLimiter(in,in)`,
  `wxnodaldg2dmethod.cc:250`), which `TSComputeRHSFunction` locks read-only.
  PETSc always forbade this; 3.22 started *enforcing* it in optimized builds,
  where it used to compile to nothing. Apollo discards the error code, so the
  write is refused, a traceback prints once per RHS evaluation, and the run
  continues unlimited. Seven shipped Euler decks are affected; no RMF deck is.
  A compile-only gate cannot catch this — it is a runtime check inside PETSc,
  not a symbol that comes and goes. (§17)
- **The formation run dies at the ANTENNA, not anywhere the mesh is refined.**
  `03-formation/formation.pin` reaches 1.73 of 20 RMF periods. Every peak sits at
  r = 0.0350-0.0398; the winding is at 0.036 and the mesh refines only r < 0.030,
  so λ_D, c/ω_pe and the resistive skin depth are 0.27, 0.36 and 0.84 cells there
  against 0.40, 0.52 and 1.23 in the column. Neither refinement test could detect
  an improvement, and both were flawed the same way - every arm stayed
  under-resolved. Use `scripts/where_peak.py` before assuming where a run failed.
  (§19)
- **You cannot pass PETSc options on Apollo's command line**, and the log will
  not tell you. `getopt_long(argc, argv, "hi:r:o:")` exits 2 on any `-ts_*`,
  `-dm_*` or `-log_view`; use `PETSC_OPTIONS=` or a `.petscrc`. The banner's
  "Time integration done using: ssp" prints before the options are read and says
  `ssp` for every scheme, so a dropped variable produces an identical log. Only
  `-ts_view` shows `Scheme:`. (§21)
- **`Numerical_Flux = Wave` aborts**, and `eigenSystem()` is dead code that is
  also wrong where it can be read. (§3, §4)

## Conventions this codebase has settled on

- **Verify before claiming.** Several findings in `docs/` were once confidently
  wrong, and the corrections are recorded in place rather than quietly edited
  out. If you write "this test would catch X", reintroduce X and run it.
- **Tests are checked by mutation.** A guard nobody's test exercises is
  decoration: `current_layer_thickness` shipped an R² check that accepted a
  power law with R² = 1.00000 and that no test touched. Before adding a guard,
  break it deliberately and confirm the suite notices.
- **Diagnostics refuse rather than guess.** A number that looks plausible is
  worse than no number: a "40 mm skin layer" in a 30 mm column, a cycle average
  taken over a third of a cycle, a field on an axis that is not in the file.
  Return NaN and say why.
- **Numbers in documentation are measured, and say so.** Where something could
  not be established — the coefficient in the penetration threshold, for one —
  the documents say that rather than supplying a remembered value.

## Layout

| path | what |
| --- | --- |
| `src/lib` | geometry, deck parsing, logging, MPI helpers |
| `src/hyper`, `src/hyperapps` | equation sets, sources and boundary conditions |
| `src/solvers`, `src/subsolvers` | the DG scheme and time stepping |
| `examples/unstructuredDG/` | input decks, one directory per physics module |
| `scripts/` | deck preprocessing and post-processing |
| `test/cxx/` | C++ unit tests; `test/` the Python ones |
| `docs/` | the known-issues list and the RMF-FRC model assessment |

`scripts/` also holds several near-duplicate legacy visualisation scripts of
uncertain provenance (§14) — prefer `test/vtu.py` and `test/deckrun.py`, which
are small and current.

## The RMF-FRC work

`docs/rmf-frc-model-assessment.md` assesses `examples/.../multifluid/rmf_frc`
against the rotating-magnetic-field current-drive literature and carries a
phased plan. Phases 0–2 are done; Phase 3, the validation campaign, is
instrumented and priced but blocked on three findings about the shipped decks
(§3.10, known-issues §13 for the edge-driven decks, and known-issues §19 for the
antenna deck, whose formation run dies at 1.73 of its 20 periods). `examples/.../rmf_frc/phase3/` holds ready-to-run
input folders and a README explaining what to run and in what order.

Tools: `scripts/rmf_diagnostics.py` (what a run produced),
`scripts/rmf_scan.py` (parameter scans and what they cost),
`scripts/rmf_cleaning_check.py` (the gate on any long run),
`scripts/rmf_make_phase3_runs.py` (regenerates the run folders),
`scripts/run_fingerprint.py` (reduces a run to a few kB so two builds can be
compared - use it before trusting a long run on a PETSc this repository has not
executed, which today means anything above 3.19),
`scripts/run_growth.py` (per-frame max|value| for all 18 components - what to
reach for when a run DIES, because rmf_diagnostics answers "what physics did
this produce" and its radial bins can read healthy to the last frame),
`scripts/where_peak.py` (WHERE each component peaks - at the conducting wall, at
the antenna winding, or in the column. run_growth says which component is
growing and not where, and on this deck those are different failures with
indistinguishable growth tables: a peak at the winding is the drive loading up,
a peak in a wall cell is an unresolved electron sheath that kills the run).
