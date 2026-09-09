# Working on Apollo

Apollo is a discontinuous-Galerkin solver for hyperbolic PDEs on unstructured
triangular meshes, built on PETSc's DMPlex. Its largest application is a
two-fluid (five-moment) plasma model coupled to perfectly-hyperbolic Maxwell —
eighteen components: electrons 0–4, ions 5–9, **E** 10–12, **B** 13–15, and the
cleaning potentials φ and ψ at 16 and 17.

## Build and test

```bash
cd src && scons build-opt          # 85 objects, about 40 s on 4 cores
PETSC_DIR=/usr/lib/petsc make -C test/cxx    # C++ unit tests, 86 checks
python3 -m unittest discover -s test         # Python tests
```

`PETSC_DIR` defaults to `/usr/lib/petsc`. There is a `build-debug` variant.
Apollo does **not** use CMake. See `README.md` for dependencies,
`requirements.txt` for the Python side, `environment.yml` for a conda
environment on a cluster.

The fast tests need no solver and run in under a second:
`test/test_rmf_diagnostics.py`, `test/test_rmf_scan.py`,
`test/test_mkdiscmesh.py`, `test/test_python_tooling.py`. So does
`test/test_petsc_compat.py`, which needs a C++ compiler but no PETSc. The ones that do need
a solver skip themselves without one, so **a green run does not by itself mean
they ran** — check the skip count.

## Things that will cost you a day if you do not know them

Most are in `docs/known-issues.md`, which is worth reading before trusting any
result. The ones that bite hardest while editing:

- **A deck value written without a decimal point aborts the run.** The parser
  types `3000000` as an integer, `get<REAL>` throws `std::bad_cast`, and the
  message names neither the key nor the file. **The mirror is just as fatal**:
  `Output_files` is read with `get<int>`, so `OUT = 6.0` dies the same way. Use
  `rmf_scan.deck_number()` and `deck_integer()` when generating decks. (§12)
- **MPI output is incomplete.** One `.vtu` per frame holds roughly 1/N of the
  cells on N ranks, and which cells is up to the partitioner. Diagnostics on a
  multi-rank run describe a fraction of the domain and do not fail. Run one rank
  when the output matters. (§13)
- **The PETSc version decides whether it compiles at all.** `petsc_compat.h`
  substitutes for plex calls upstream removed at 3.13 and 3.14; the guards said
  3.18 for both, so Apollo built on 24.04 (PETSc 3.19) and would not compile on
  22.04 (PETSc 3.15) — including in Colab. Fixed, and pinned by
  `test/test_petsc_compat.py`, which compiles the header against a stub
  `<petsc.h>` for every release from 3.11 up. If you touch a version guard, run
  it: it fails in both directions.
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
