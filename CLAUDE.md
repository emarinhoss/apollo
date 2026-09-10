# Working on Apollo

Apollo is a discontinuous-Galerkin solver for hyperbolic PDEs on unstructured
triangular meshes, built on PETSc's DMPlex. Its largest application is a
two-fluid (five-moment) plasma model coupled to perfectly-hyperbolic Maxwell —
eighteen components: electrons 0–4, ions 5–9, **E** 10–12, **B** 13–15, and the
cleaning potentials φ and ψ at 16 and 17.

## Build and test

```bash
cd src && scons build-opt                        # 85 objects, about 40 s on 4 cores
PETSC_DIR=/usr/lib/petsc make -C ../test/cxx     # C++ unit tests, 86 checks
cd .. && python3 -m unittest discover -s test    # Python tests
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

`PETSC_DIR` defaults to `/usr/lib/petsc`; inside a conda environment it is
`$CONDA_PREFIX`. There is a `build-debug` variant.
Apollo does **not** use CMake. See `README.md` for dependencies,
`requirements.txt` for the Python side, `environment.yml` for a conda
environment on a cluster.

**161 Python tests, and 3 of them are 72% of the runtime.** Measured on a
4-core container with a `build-opt` binary present:

| | tests | time |
| --- | --- | --- |
| `test_rmf_antenna_field.py` | 3 | 645 s |
| `test_rmf_bc_field.py` | 4 | 215 s |
| `test_vortex_accuracy.py` | 2 | 20 s |
| everything else (11 modules) | 152 | **3.4 s** |
| whole suite | 161 | 894 s |

So run the 152 while you work and the whole suite before you push. The 152 need
no solver, and naming them is the only reliable way to select them — a glob is
not, because the fast and slow modules interleave alphabetically:

```bash
cd test && python3 -m unittest \
    test_rmf_diagnostics test_rmf_scan test_rmf_cleaning_check \
    test_mkdiscmesh test_python_tooling test_deck_preprocess \
    test_vtu_reader test_eigen_paths test_include_paths \
    test_conda_paths test_petsc_compat
# Ran 152 tests in 2.598s -- OK
```

`cd test` first. From the repository root the stdlib's own `test` package wins
the import, so `python3 -m unittest test.test_rmf_scan` dies with
`ModuleNotFoundError: No module named 'test.test_rmf_scan'` — pointing at a file
that is plainly sitting there. Use `discover -s test`, or `cd test` and name the
modules as above. `test_petsc_compat` is in the fast set; it needs a C++
compiler but no PETSc.

The 9 that do need a solver skip themselves without one, so **a green run does
not by itself mean they ran** — check the skip count. With a binary present the
suite reports `Ran 161 tests` and no skips.

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
  add up. Pinned by `test/test_vtu_reader.py`, which builds both layouts, and
  by a CI step that reads back a `.vtu` the 22.04 job's PETSc actually wrote.
- **There is no checkpoint or restart.** An interrupted run is a lost run. (§6)
- **The two-fluid slope limiter produces NaN** after a few tens of steps, so no
  multifluid case can be limited. (§2)
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
instrumented and priced but blocked on two findings about the shipped decks
(§3.10 and known-issues §13). `examples/.../rmf_frc/phase3/` holds ready-to-run
input folders and a README explaining what to run and in what order.

Tools: `scripts/rmf_diagnostics.py` (what a run produced),
`scripts/rmf_scan.py` (parameter scans and what they cost),
`scripts/rmf_cleaning_check.py` (the gate on any long run),
`scripts/rmf_make_phase3_runs.py` (regenerates the run folders).
