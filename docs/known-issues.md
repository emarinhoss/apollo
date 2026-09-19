# Known issues

Findings from an end-to-end review that are **not** fixed, with enough detail to
act on. Everything here was confirmed by reading the code and, where a `repro`
block is given, by running it — none of it is speculation. Items are ordered by
how much damage they can do to a result.

Fixed issues are not listed; see the git history.

---

## 1. The bundled collisional–radiative example does not run

```
cd examples/unstructuredDG/multifluid/collisionalRadiative
PYTHONPATH=../../../../scripts python3 ../../../../scripts/wxinpparse.py -i collisonalRadiative.pin
../../../../src/build-opt/apollo -i collisonalRadiative.inp
```

reaches the atomic-data stage and then aborts:

```
file:: could not be opened
no elastic collision xs data for Ar+1
double free or corruption (out)
Signal: Aborted (6)
```

Two things are visible in that: a file is opened with an empty name and the
failure is not treated as fatal, and the subsequent path corrupts the heap.
`KRates.cpp` reads its atomic data with unchecked `fopen`/`fscanf` throughout;
that is where to start. (The directory *creation* path in the same file was
rewritten separately — see the git history — but the data *reading* path was
not.)

---

## 2. The two-fluid slope limiter produces NaN energy

`tuAliabadiLimiter` (`src/hyperapps/multifluid/aptualiabadilimiter.cc`) is the
only slope limiter that handles the **two-fluid** system, and it does not work.
(`eulerLimiterHW`, `WxHestavenWarburtonEulerLimiter`, limits the five-component
Euler state and is enabled in several shipped Euler decks; it cannot be applied
to the eighteen-component multifluid state.) It is commented out in the
one deck that references it (`examples/unstructuredDG/multifluid/rmf_frc/frc2d.pin`,
`#Limiter = [twoFluidLimiter]`), which is presumably why this has gone unnoticed.

Three defects that stopped it before it could take a step have been fixed:

- it calls boundary conditions with a null `AreaInts`, and every RMF boundary
  condition dereferenced it;
- `wxNodalDGgeometry2D` sized `_ETETF` by `_Ktotal` while indexing it by local
  cell id — a heap overflow present in every run, not only limited ones;
- the height-0 stratum contains more than the mesh's cells: `createMesh` calls
  `DMPlexConstructGhostCells`, which appends one ghost cell per boundary facet
  (192 of them on the shipped `rmf_frc` mesh, taking 7792 to 7984), and the
  limiter asked those for their normals;
- it handed boundary conditions a two-element coordinate array. The DG scheme
  passes `xc[5]` = (t, x, y, -, dt) and every boundary condition here reads
  `xc[0]` as the time and `xc[1]`, `xc[2]` as the position; the limiter declared
  `XC[2]` and filled it with (x, y). So on that path an RMF boundary condition
  took the node's x-coordinate for the time - an envelope evaluated at t = 0.03 s
  and a phase omega*t of about 1.5e5 radians - and read `xc[2]` past the end of
  the array. Nothing else set the limiter's clock either, so
  `WxNodalDG2dMethod::applyLimiter` now does before calling it.

What remains is in the limiter's own numerics: it runs for a few tens of steps
and then stops with `*** NaN energy in Euler limiter ***`, on the shipped
`rmf_frc` deck as much as on the antenna deck beside it. The coordinate defect
above looked like a promising explanation - a boundary condition fed a nonsense
time returns a nonsense ghost state - but it is not the cause: with it fixed the
run still stops at the same step 28.

```repro
REPO=$(git rev-parse --show-toplevel)
cd $(mktemp -d) && cp $REPO/examples/unstructuredDG/multifluid/rmf_frc/{*.pin,*.msh} .
sed -i 's/^TEND = 2.e-6.*/TEND = 4.e-9/;s/^     #Limiter/     Limiter/' frc2d.pin
PYTHONPATH=$REPO/scripts python3 $REPO/scripts/wxinpparse.py -i frc2d.pin
$REPO/src/build-opt/apollo -i frc2d.inp | tail -3
```

This matters beyond that example: with it broken, no **multifluid** case can be
limited at all, which is part of why the RMF antenna deck cannot carry a plasma
edge and so has no vacuum region around its column — see
`docs/rmf-frc-model-assessment.md`, Phase 1.

---

## 3. `eigenSystem()` is dead code, and wrong where it can be read

Every equation object defines `eigenSystem(d, q, ev, lev, rev)`, filling in
eigenvalues and left and right eigenvectors of the flux Jacobian. Nothing uses
them. The only caller is `WxHyperbolicEqnSet::eigenSystem`
(`src/hyper/wxhyperboliceqnset.cc:428`), and nothing calls that:

```repro
grep -rn "eigenSystem" src --include=*.cc --include=*.h | grep -v build
```

returns only the definitions and that one forwarding call.

That matters because it is the natural building block for anything needing a
characteristic decomposition - which is exactly what
`docs/rmf-frc-model-assessment.md` Phase 2 proposed to build on - and it is
unexercised, so nothing would notice if it were wrong. In
`WxPHMaxwellEqn::eigenSystem` it is: lines 567-568 read

```cpp
  gamma = q[6];
  kappa = q[7];
```

and use those as the divergence-cleaning wave speeds. `q[6]` and `q[7]` are the
cleaning POTENTIALS phi and psi - solution components, varying in space and time
- not the speeds, which are the `_gamma` and `_chi` members set from the deck.
The eigenvalues it returns are therefore `c*phi` and `c*psi`, which are not wave
speeds of anything.

It is also written only for axis-aligned directions (`d` = 0, 1, 2), so on an
unstructured mesh a caller would have to rotate into the face frame around it.

Phase 2 went around it: the Maxwell characteristic decomposition is written out
in closed form in `src/hyperapps/maxwell/apmaxwellcharacteristics.h` and checked
against a Jacobian assembled by calling `flux()` itself
(`test/cxx/test_maxwell_characteristics.cc`). Fixing `eigenSystem` would still be
worth doing before anything else relies on it - a Roe flux for Maxwell, say -
and that test shows what the answer should be.

---

## 4. `Numerical_Flux = Wave` aborts

`WxEulerEqn::applyWavePropagationFluxes` (`src/hyperapps/euler/wxeulereqn.cc`)
forms the jump as `df[m] = qM[m] - qP[m]` and then calls
`this->rp(0, qP, qM, 0, 0, df, ...)`, passing `qP` as the left state and `qM` as
the right — the opposite order to the flux evaluation directly above it. On the
isentropic vortex it aborts with negative pressure even after the rotation bug
in the same function was fixed.

Deciding which convention `rp()` expects needs someone who knows the wave
decomposition. Once fixed, add `'Wave'` to the flux list in
`test/test_vortex_accuracy.py::test_every_flux_solves_the_vortex`, which is
already written to cover it.

---

## 5. Roughly two thirds of `src/lib` is not built

A transitive include closure from the 89 sources the SConscripts actually
compile reaches 222 files. **129 of the 188 files in `src/lib` are in none of
them.** Among the unreachable:

| Group | Files | Note |
| --- | --- | --- |
| NDG++ value types | `Vec_Type.h`, `Mat_COL.h`, `VecObj_Type.h`, `MatObj_Type.h`, `ArrayGen.h`, `ArrayMacros.h`, `Region1D/2D.h`, `MappedRegion1D/2D.h`, `Index.h`, `Registry_Type.h`, … | Cannot link even if included: `umERROR`, `umWARNING`, `gVecData` are declared and never defined |
| HDF5 I/O | `wxhdf5io.cc/.h`, `wxhdf5iotmpl.*`, `wxhdf5traits.*`, `wxiobase.cc` | Why there is no checkpointing (§6) |
| FEM geometry | `wxfemgeometry.*`, `wxfemquadrature.*`, `wxfemshapefuncs.h`, `wx3dfemgeom.*` | |
| Box/grid machinery | `wxbox*.cc/h`, `wxgridbox.*`, `wxgridrange.*`, `wxsplitbox.*` | |
| Dependency graph | `wmdependencygraph.*`, `wmnametree.*`, `wmindexer.*` | |
| Quadrature tables | `dataN02.h` … `dataN16.h`, `data3dN01.h` … `data3dN09.h` | Reachable only through the commented-out dispatch in `wxpnodaldgfunctions.h` |
| Dead vendors | `Blas_ACML.h`, `acml.h` | AMD ACML has been discontinued since 2015 |

`src/lapack_lite/` is likewise excluded — its `SConscript` call is commented out
in `src/SConscript`.

This is not a bug on its own, but it is a large maintenance cost: the majority
of what looks like the core library is inert, and a reader cannot tell which
half matters without doing this analysis. Each group is worth a decision:
revive it, or delete it and let the history hold it. Note that
`src/lib/nudg/`, which was a byte-identical copy of much of this, has already
been removed.

**Reproduce:**

```bash
python3 - <<'PY'
import os, re, subprocess
incdirs = ['lib','lapack_lite','solvers','subsolvers','hyper','hyperapps','xapollo','.']
sources = []
for sc in subprocess.run(['find','.','-name','SConscript'],
                         capture_output=True, text=True).stdout.split():
    d = os.path.dirname(sc)
    for m in re.finditer(r"'([^']+\.(?:cc|cpp|c))'", re.sub(r'#.*', '', open(sc).read())):
        p = os.path.normpath(os.path.join(d, m.group(1)))
        if os.path.isfile(p): sources.append(p)
inc = re.compile(r'^\s*#\s*include\s*[<"]([^>"]+)[>"]', re.M)
seen, queue = set(sources), list(sources)
while queue:
    f = queue.pop()
    for i in inc.findall(open(f, errors='replace').read()):
        for d in [os.path.dirname(f)] + incdirs:
            c = os.path.normpath(os.path.join(d, i))
            if os.path.isfile(c):
                if c not in seen: seen.add(c); queue.append(c)
                break
lib = [f for f in sorted(os.listdir('lib')) if f.endswith(('.h','.cc'))]
un = [f for f in lib if os.path.normpath('lib/'+f) not in seen]
print(f'{len(sources)} sources reach {len(seen)} files; {len(un)}/{len(lib)} of src/lib unreachable')
PY
```

(run from `src/`)

---

## 6. Checkpoint/restart — FIXED, and here is how to use it

**This entry used to say there was no restart and that an interrupted run was a
lost run.** That was true and it cost someone twelve hours: a 20-period
formation run died 1.73 periods in, and the whole thing had to be paid again.

A run now writes a checkpoint at every output frame, and `-r` resumes from it:

```bash
apollo -i formation.inp                      # writes <runName>.checkpoint each frame
apollo -i formation.inp -r formation.checkpoint    # picks up where it stopped
```

**A resumed run reproduces a straight run exactly.** Not closely - the
fingerprint of a twice-interrupted, twice-resumed gate run against the same deck
run straight through is `0.000e+00` across all 378 solution arrays and all seven
frames. `test/test_restart.py` asserts that, and it is the assertion that
matters: a restart that changes the answer is worse than no restart.

**What it costs.** One file, rewritten each frame, holding one solution vector -
5.27 MB on the phase3 mesh, against 5.4 MB for a single `.vtu`. A per-frame
checkpoint would have been 1.2 GB over 240 frames and doubled what the run
already writes, so the file is deliberately rolling. The price of that choice:
a crash *during* the checkpoint write loses the checkpoint. The write is
ordered after the `.vtu` so the frames stay the record of last resort.

**Why PETSc binary and not the `.vtu`.** The frames hold the same numbers, and
reading one back would let you resume a run that predates this feature. But
`VecView`/`VecLoad` round-trip a `Vec` bit-exactly and handle a different rank
count on the way in, whereas parsing our own output means re-deriving the DMPlex
ordering by hand - and getting that subtly wrong yields a plausible wrong answer
rather than an error. Frames written before this feature existed cannot be
resumed from; that is the accepted cost.

**Why the state is reloaded *after* `init()`, not instead of it.** `init()`
creates the solution vector, lays down the initial condition, and sizes `_dt`
from it — and Apollo sizes `dt` once and never revisits it
(`apsolver.cc`, `_dt = fmin(_dt, suggestedDt)` sits in setup, not in the time
loop). A resumed run must integrate with the same `dt` as the run it continues,
so `init()` has to happen. The old commented-out `load()` branch called
*instead* of `init()` and would have left `dt` unsized.

**It refuses rather than resuming into a different problem.** Each of these is
tested, and each would otherwise produce numbers that look like a continuation
and are not:

| the checkpoint | the refusal |
| --- | --- |
| a state of a different length | `holds a state of N values but this problem has M` |
| a different `Output_files` | `was written by a run with Output_files = N` — the frame spacing would differ |
| a frame outside `0..OUT` | `names frame N, which is outside 0..OUT` |
| a frame equal to `OUT` | `that run already finished` |
| no `.meta` sidecar | `no checkpoint metadata beside <path>` |

**Resuming twice works**, which is not free: the first implementation recorded
the *running* interval count rather than the deck's, so resuming at frame 2 of 6
wrote `nout 4`, and resuming from that compared 4 against the deck's 6 and
refused. One resume worked and a second did not. `_noutDeck` is now kept apart
from `_nout` for exactly this reason, and `test_restart.py` interrupts and
resumes twice.

**Through the phase3 runner.** `examples/.../rmf_frc/phase3/run_one.sh` puts
each run in its own directory, so the checkpoint is there too:

```bash
./run_one.sh 03-formation/formation.pin                  # dies at hour 12
APOLLO_RESUME=1 ./run_one.sh 03-formation/formation.pin  # continues it
```

Plain `./run_one.sh` on a directory that already holds an unfinished checkpoint
**refuses**, exit 3, and prints both commands. That refusal is the point: a
fresh run writes frame 0 and the checkpoint before anything else, so re-running
by reflex is what destroys the thing you would have resumed from.
`APOLLO_RESUME=0` discards it deliberately. `slurm_array.sh` sets
`APOLLO_RESUME=1`, so resubmitting an array after a wall-clock kill continues
every task where it stopped. `test/test_run_one.py` pins all of it against a
stub solver, in under a second and with no PETSc.

`run_one.sh` also now survives the solver returning non-zero, which `set -e`
used to turn into a silent exit: it says whether the run was **interrupted**
(status ≥ 128, a signal — resume it) or **stopped by the solver itself**
(the NaN guard at `wxnodaldg2dmethod.cc:530` calls `exit(1)` — a resume lands
back in the same place), and runs `scripts/run_growth.py` over the frames. The
run that most needed analysing was the one that got none.

---

## 7. Error paths bypass MPI

67 `exit()`/`abort()` calls remain in `src/`, against 4 `MPI_Abort`. Calling
`exit()` from one rank leaves its peers blocked in whatever collective they were
in; the job hangs until the scheduler kills it, and any diagnostic is lost.
Several sit inside per-element loops (the NaN guards in `wxnodaldg2dmethod.cc`
and `wxeulereqn.cc`), which is also the worst place to call it from.

The entry point now aborts through MPI correctly; the fix for the rest is to
throw a `WxExcept` and let it reach that handler.

---

## 8. PETSc return codes are discarded

555 PETSc calls, 63 `CHKERRQ`. A failing `DMPlexCreateFromFile`, `VecGetArray`
or `MatAssemblyEnd` is not noticed, and execution continues with an object that
was never created. `ApSolver::createMesh` and `ApSolver::OutputVTK`
(`src/solvers/apsolver.cc`) additionally open with `PetscFunctionBeginUser` and
never call `PetscFunctionReturn`, so they do not appear correctly in PETSc's own
stack traces either.

Converting these to `PetscCall(...)` is mechanical but touches most of the
codebase, so it wants to be its own change.

---

## 9. Deprecated PETSc APIs

`src/subsolvers/wxpetsctimestepping.h` calls `TSSetDuration`,
`TSSetInitialTimeStep` and `TSGetTimeStepNumber`, all deprecated since PETSc 3.8
(2017), and `PETSC_NULL`, deprecated in 3.19. They still compile and account for
151 of the build's warnings. Replacements are `TSSetMaxSteps`/`TSSetMaxTime`,
`TSSetTimeStep`, `TSGetStepNumber` and `PETSC_NULLPTR`.

---

## 10. Build warnings

A clean `scons build-opt` emits 786 warnings (down from 821 at the start of the
review; the uninitialized-variable reports specifically went from 35 to 5).

| Category | Count | Comment |
| --- | --- | --- |
| `-Wsign-compare` | 185 | `unsigned` loop counters against `int` bounds; mostly benign, individually cheap to fix |
| `-Wmisleading-indentation` | 175 | Worth reading each one: indentation that lies about control flow is how a guard silently stops guarding |
| `-Wdeprecated-declarations` | 151 | §9, plus `std::auto_ptr` in the vendored muParser |
| `-Wunused-variable` | 93 | |
| `-Wreorder` | 72 | Member initialiser lists that do not match declaration order — currently harmless, a real bug the moment one member's initialiser reads another |
| `-Woverloaded-virtual` | 47 | A derived `init()` hiding the base's `init(PetscReal, Vec)`. §11 |
| `-Wmaybe-uninitialized` | 5 | |
| `-Wformat-security` | 4 | `sprintf(dst, variable)` in the CR module; a `%` in the data is a crash |

`-Werror` is not realistic yet, but `-Wextra` on new code would be.

---

## 11. `init()` signature mismatches hide the base virtual

`WxObject::init(PetscReal, Vec)` is virtual; `ApFVM2dScheme::init()` and
`ApDomainDecompCheck::init()` declare a no-argument `init()`, which hides rather
than overrides it. The base's do-nothing body runs instead, so those schemes'
initial conditions are never applied. `apfvm2dscheme.cc` is currently commented
out of `src/subsolvers/SConscript`, so this is latent — but it is exactly the
trap that `override` exists to catch, and the codebase uses it nowhere.

---

## 12. A deck value written without a decimal point aborts the run

`WxCryptSetLexer::_scan_number` (`src/lib/wxcryptsetlexer.cc:213-223`) types a
number with no `.` and no exponent as `WX_INT` and stores it through `atoi`;
everything else goes through `atof` as `WX_REAL`. `WxCrypt::get<REAL>`
(`src/lib/wxcrypt.h:73-81`) then does `wx_any_cast<double>` on a stored `int`,
which throws `std::bad_cast`. The solver reports

    Apollo: unexpected error on MPI rank 0: std::bad_cast

and calls `MPI_ABORT`. **The message names neither the key nor the file**, so a
deck with one integer-looking constant in it fails at setup with nothing to go
on. The value need not look unusual: `LIGHT = 1000000` is fatal where
`LIGHT = 1000000.0` and `LIGHT = 1.0e6` are fine.

This bites hardest where decks are generated rather than typed. A parameter scan
writes a deck per point, and the obvious formatters produce integer-looking text
for round values: `'%.10g' % 1.0e6` is `1000000` and `'%.10g' % 3.0e5` is
`300000`. (Bare `'%g'` is not the hazard — it switches to exponential at
`1e+06` — but any format with enough precision to keep small values readable
renders large round ones without a point.) So a scan fails at every point,
identically, with an error that points nowhere. `scripts/rmf_scan.py` carries a
`deck_number()` helper for exactly this reason, and `test/test_rmf_scan.py`
pins it.

```repro
REPO=$(git rev-parse --show-toplevel)
cd $(mktemp -d) && cp $REPO/examples/unstructuredDG/multifluid/rmf_frc/{frc2d.pin,optimizedCircle2.msh} .
sed -i 's/^LIGHT = 3.0e6.*/LIGHT = 3000000/;s/^TEND = 2.e-6.*/TEND = 1.e-10/' frc2d.pin
PYTHONPATH=$REPO/scripts python3 $REPO/scripts/wxinpparse.py -i frc2d.pin
$REPO/src/build-opt/apollo -i frc2d.inp 2>&1 | grep -E "bad_cast|error"
```

**How loudly it fails depends on where the value is used**, which is worth
knowing before hunting one of these. A value read with `get<REAL>` (like `c0`)
throws `std::bad_cast` and names nothing. A value that reaches an *expression*
is caught by the expression parser instead, which does say what is wrong:

    Undefined name 'num_dens' in expression. ... note that a constant written
    without a decimal point is read as an integer and is not exported to
    expressions.

So `n_dens = 100000000000000000000` aborts with that message rather than
producing a wrong density; an earlier version of this section claimed it was
silently wrong, and it is not. The routine does carry its author's note that
`atoi` will overflow on a large literal (`wxcryptsetlexer.cc:136-138`, "FIX
THIS: MAIN ISSUE IS THAT AN OVERFLOW MAY OCCUR"), but every path tried here
fails loudly and no case was found where the overflow is silent. The fix for the
whole family is the same: parse every numeric literal as a double and let
`get<int>` narrow, or name the offending key in the `bad_cast` path as the
expression parser already does.

---

## 13. The edge-driven RMF boundary drives the divergence-cleaning potential without bound

`twoFluidSimplifiedRMFBC` writes the applied transverse field into the ghost
state at every boundary face. The interior field on the other side of that face
is whatever the run has produced, so the ghost state carries a jump in the normal
component of **B** — which is exactly what the perfectly-hyperbolic Maxwell
system's cleaning potential ψ (component 17) exists to transport away. Here it is
being fed faster than it can carry, and it grows without bound.

Measured on the shipped `rmf_frc/frc2d.pin` over 210 steps (`DIVB_SPEED =
DIVE_SPEED = 1.0`, everything else as shipped):

| frame | rms\|ψ\| | rms\|E\| | ψ / E |
| --- | --- | --- | --- |
| 1 | 58.0 | 152.3 | 0.38 |
| 2 | 189.0 | 107.6 | 1.76 |
| 3 | 367.1 | 104.6 | 3.51 |
| 4 | 581.4 | 104.1 | 5.59 |
| 5 | 826.3 | 93.5 | 8.84 |
| 6 | 1097.5 | 90.8 | **12.08** |

ψ is still rising linearly at the end. φ (component 16) stays identically zero
throughout, so this is specifically the magnetic cleaning potential. The
antenna-driven deck fails the other way round — φ and the in-plane **E** grow
while ψ stays flat — and for a different reason; see §19. It is the
boundary condition that produces it: with `Bomega = 0` it is identically zero.

**It changes the answer.** ψ enters the induction equation in the slot **E**
occupies (`wxphmaxwelleqn.cc`: `f[ibx] = gamma*q[17]` against `f[iby] = -q[iez]`),
so a large ψ is not a diagnostic quantity sitting to one side. Two runs differing
in nothing but the cleaning speed diverge fast:

| frame | rms ΔB_x / rms B_x | rms ΔB_z / B_bias |
| --- | --- | --- |
| 1 | 0.109 | 0.00002 |
| 3 | 0.566 | 0.0002 |
| 6 | **0.916** | 0.0009 |

After 210 steps — 0.005 of one RMF period — the transverse field is essentially
a different field depending on a numerical parameter. B_z is far better behaved,
which is why a short run looks healthy: the quantity Phase 3 item 0 watches is
the last one to notice.

**The antenna deck does not have this problem.** Driving the same plasma with an
antenna current instead of an edge field gives ψ/E of 0.03 to 0.09 over the same
210 steps, not 12, because it never writes a discontinuous **B** into a ghost
state. Any Phase 3 item needing more than a few hundred steps should use
`rmf_frc/antenna/frc2d.pin`.

```repro
REPO=$(git rev-parse --show-toplevel)
for DS in 1.0 0.0; do
  d=$(mktemp -d); cd $d
  cp $REPO/examples/unstructuredDG/multifluid/rmf_frc/{frc2d.pin,optimizedCircle2.msh} .
  sed -i -e "s/^DIVB_SPEED = 1.0/DIVB_SPEED = $DS/" -e "s/^DIVE_SPEED = 1.0/DIVE_SPEED = $DS/" \
         -e 's/^TEND = 2.e-6.*/TEND = 5.4732e-9/' -e 's/^OUT  = 10/OUT  = 6/' frc2d.pin
  PYTHONPATH=$REPO/scripts python3 $REPO/scripts/wxinpparse.py -i frc2d.pin
  $REPO/src/build-opt/apollo -i frc2d.inp > log.txt 2>&1
  echo "$DS -> $d"
done
# then compare component 17 and component 13 between the two directories with
# test/vtu.py and test/deckrun.py.
```

The mechanism to attack is the normal-component jump the ghost state creates.
**Characteristic injection is not the remedy, and this was tested rather than
assumed.** Switching the boundary condition to impose only the incoming Riemann
invariants (`c0 = LIGHT` in the `coil1` block, the Phase 2 path in
`src/hyperapps/maxwell/apmaxwellcharacteristics.h`) reproduces the ψ/E column
above to four decimal places, frame for frame. That is the expected answer for
the reason Phase 2 recorded: at χ = γ = 1 — which is what `DIVB_SPEED =
DIVE_SPEED = 1.0` sets — Lax–Friedrichs is exactly upwind, so the ghost's
outgoing part cannot influence the flux and the two ghost states give identical
answers. The over-specification Phase 2 removed is real but is not what is
driving ψ here.

What is left to try, in order of cost: give ψ somewhere to go by damping it
(the usual GLM treatment adds a source −ψ/τ, which this equation set does not
have); make the applied normal component consistent with the interior rather
than imposed, so there is no jump to clean; or use the antenna deck, which
avoids the question entirely.

---

## 14. Duplicated and uncertain post-processing scripts

Several files in `scripts/` are near-duplicates with no indication of which is
canonical: `wxdata.py` / `wxdata_3949.py`, `wxunsdgdata.py` / `wxunsdgdata2.py`,
`SGC.py` / `SGC_MB.py`, and the four-way `wxxdmf*.py` family. All of them now
parse and import, so a maintainer can diff them and delete the losers.

`scripts/SGC.py` also carries genuinely dead code: `Node.getLocalID` referenced
an `owners` attribute that is never assigned anywhere (it is now an explicit
`NotImplementedError`), and `buildLineFaces`/`buildRectangleFaces`/
`buildHexahedronFaces` are unreferenced.

## 15. The solver's answer depends on the number of MPI ranks

The same case run on one rank and on two produces materially different fields.
This is not a partitioning artefact in the reader and it is not floating-point
noise: it is the solution.

**How it was hidden.** Until now this comparison could not be made. PETSc writes
one `<Piece>` per rank into each `.vtu`, and `test/vtu.py` walked `DataArray`
elements into a flat dictionary, so only the last rank's copy of each array
survived. A multi-rank frame therefore *looked* truncated, and the repository
recorded that as a defect in the writer — "one `.vtu` per frame holds roughly
1/N of the cells". It does not: the counts behind that claim (7792, 3896, 1961
cells on 1, 2, 4 ranks) were the size of the last piece. With the reader fixed a
two-rank frame reads back 4304 cells against the one-rank file's 4304, and the
runs can be compared for the first time.

**repro** — the Maxwell circular pulse, chosen because it is a different module
and a different boundary condition from the RMF work, so this is not specific to
`twoFluidSimplifiedRMFBC`:

```bash
APOLLO_EXAMPLES_WORK=/tmp/one test/run_examples.sh -b src/build-opt/apollo maxwell-circular-pulse
APOLLO_EXAMPLES_WORK=/tmp/two test/run_examples.sh -b src/build-opt/apollo -j 2 maxwell-circular-pulse
python3 - <<'PY'
import sys, glob; sys.path.insert(0, 'test')
import vtu, numpy as np
one = sorted(glob.glob('/tmp/one/*/*.vtu')); two = sorted(glob.glob('/tmp/two/*/*.vtu'))
for i in (0, 3, 5):
    a, b = vtu.read(one[i])['solutiondg.0'], vtu.read(two[i])['solutiondg.0']
    print(i, a.size, b.size, np.sort(a.ravel()).sum(), np.sort(b.ravel()).sum())
PY
```

Observed: frame 0 is identical on both — the initial condition is laid down the
same way — and the runs have diverged by frame 3.

```
0  4304 4304   +0.000000e+00   +0.000000e+00
3  4304 4304   +1.291388e+02   -2.405673e+01
5  4304 4304   +1.393984e+00   -2.393447e+01
```

Values are compared as sorted multisets, so cell ordering and partitioning play
no part; the cell counts agree exactly. The difference exceeds the magnitude of
the solution itself, which rules out summation-order noise in a linear problem
integrated for six frames.

The same is visible in the RMF gate deck: on one rank the penetration ratio is
0.0000 at every frame over 210 steps, on two ranks it climbs to 0.71, and B_z on
axis moves from 6.000000e-03 to 6.003557e-03.

**Consequence.** Run one rank for anything whose numbers you intend to use. That
was already the standing advice, but for a reason that turned out to be false;
this is the real one, and it is worse, because a multi-rank run now produces a
complete, plausible, coherent file that simply disagrees with the serial answer.

**Not diagnosed here.** The obvious suspects are the halo exchange for the DG
face terms and whether a partition-interior face is being treated as a domain
boundary — `twoFluidSimplifiedRMFBC` writing the applied field at a rank
interface would look exactly like the RMF observation above — but the Maxwell
case uses a different boundary condition and shows it too, so a shared cause in
the face pairing or the ghost exchange is more likely than either boundary
condition. `wxnodaldggeometry2d.cc:FacePair2d` and the `DMPlexGetHybridBounds`
shim's `-1` returns for face bounds are worth reading first.

## 16. The CFL timestep is set by the last equation in the deck's list, not the fastest

`WxHyperbolicEqnSet::DGnumericalFlux` passes **one** `maxSpeed` pointer to every
equation in the set, in the order the deck lists them:

```c
// src/hyper/wxhyperboliceqnset.cc:186-203
for (i=_eqnSys.begin(); i!=_eqnSys.end(); ++i)
  (*i)->DGnumericalFlux(normals, qM+mloc, qP+mloc, nflux+mloc, maxSpeed);
```

and each equation **assigns** to it rather than accumulating a maximum:

```c
// src/hyperapps/maxwell/wxphmaxwelleqn.cc:523
*maxSpeed = lambda;                       // lambda = dmax(chi*c0, gamma*c0, c0)
// src/hyperapps/euler/wxeulereqn.cc:1118
*maxSpeed = lambda;                       // lambda = |u| + sqrt(gamma p / rho)
```

There is no comparison at either site. Whatever the last equation writes is what
reaches the CFL, `dt = (2/3)*cfl*dtscale*rMin/maxSpeed`
(`wxnodaldg2dmethod.cc:553`). The earlier equations' wave speeds are computed,
used for their own Lax-Friedrichs dissipation, and then discarded.

**On the shipped decks this happens to be nearly right, which is why it has not
bitten.** Every rmf_frc deck lists `Equations = [eulerElc, eulerIon, maxwell]`,
so Maxwell is last and dt comes from `dmax(chi*c0, gamma*c0, c0)`. With
c0 = 3.0e6 m/s against an electron sound speed of 2.9657e6 m/s, the number that
survives is the larger one anyway — by 1.2%.

**Two ways it is wrong.**

*The electron fluid is not represented in the timestep at all.* Section 3.10 of
the RMF assessment finds that in the penetrated state — the regime the whole
Phase 3 campaign exists to measure — the fastest electron characteristic reaches
1.039 c0 = 3.117e6 m/s. That speed never reaches the CFL: dt is computed as if
the fastest wave in the problem were 3.000e6 m/s, so the step is about 3.9% too
large exactly where the physics of interest lives. What was recorded as a
modelling concern ("the reduced speed of light is level with the electron
speeds") is also a numerical one.

*Reordering the deck's `Equations` list silently changes the timestep.* It reads
as a list of things to solve, not as something with a significant order. Writing
`Equations = [maxwell, eulerElc, eulerIon]` — the same physics — leaves the ion
Euler speed last, and at Te = Ti = 30 eV with a hydrogen mass ratio that is
6.92e4 m/s against 3.0e6:

| last equation | lambda that sets dt | dt relative to shipped |
| --- | --- | --- |
| `maxwell` (as shipped) | 3.000e6 m/s | 1.0x |
| `eulerElc` | 2.966e6 m/s | 1.01x |
| `eulerIon` | 6.921e4 m/s | **43.3x** |

A 43x timestep will not survive contact with the light waves it is not
resolving, so this fails loudly rather than quietly — but it fails for a reason
nothing in the deck or its documentation would lead anyone to suspect.

**repro** — read it off the source; no run is needed. `grep -n 'maxSpeed' ` on
the two equation files shows the bare assignments, and
`wxhyperboliceqnset.cc:186-203` shows the shared pointer and the loop order.

**Not fixed here.** The obvious repair is to make the two sites accumulate
(`*maxSpeed = dmax(*maxSpeed, lambda)`) with the caller zeroing it first, but
that changes the timestep of every multifluid run in the repository — including
every result quoted in `docs/rmf-frc-model-assessment.md` — so it is a decision
about revalidation, not a one-line patch. Note also that the current behaviour
is what makes dt exactly state-independent for the RMF decks, which is why the
cleaning-speed gate in `scripts/rmf_cleaning_check.py` can compare two runs
without a timestep confound; a fix would need that gate re-examined.

---

## 17. The slope limiters write into a vector PETSc locks, and from 3.22 that is enforced

`WxNodalDG2dMethod::step` limits **in place**:

```c
// src/subsolvers/wxnodaldg2dmethod.cc:248-251
if(_haveLimiter){
    applyLimiter(in,in);
    isInfinityOrNAN(in, "NAN/INF encountered in limiter vector of DG step-function.\n");
}
```

`in` is the TS solution vector, handed to Apollo by `TSComputeRHSFunction` for
the duration of one right-hand-side evaluation. Both limiters end by writing the
limited state back into it:

```c
// src/hyperapps/euler/wxhesthavenwarburtoneulerlimiter.cc:487-488
// src/hyperapps/multifluid/aptualiabadilimiter.cc:451-452
DMLocalToGlobalBegin(_dm, local_out, INSERT_VALUES, q_limited);
DMLocalToGlobalEnd(_dm, local_out, INSERT_VALUES, q_limited);
```

An RHS function must not modify its input, and PETSc has always said so:
`TSComputeRHSFunction` wraps the callback in `VecLockReadPush(U)` /
`VecLockReadPop(U)`. **What changed is whether the rule is enforced.** Through
PETSc 3.21 the whole lock API sat behind `#if defined(PETSC_USE_DEBUG)`, with an
`#else` branch defining `VecSetErrorIfLocked(x, arg)` as `PETSC_SUCCESS` — so in
an optimized build the check compiled to nothing and the illegal write silently
succeeded. PETSc 3.22 removed the guard ("Make `VecLock` API active in optimized
mode", `doc/changes/322.md`), and from 3.22 the lock is live in every build.

**So this is a pre-existing defect that 3.22 exposes, not one it creates.** Every
Euler limiter result this repository has ever produced was mutating the TS state
vector mid-evaluation.

**It does not abort.** `DMLocalToGlobalBegin/End` are called bare — no `ierr =`,
no `CHKERRQ` — and nothing up the chain through `ApSolver::ComputeRHSforTS`
(`src/solvers/apsolver.cc:461-495`) checks a `PetscErrorCode` either. On 3.22+
the write is refused with `PETSC_ERR_ARG_WRONGSTATE` ("Vector ... was locked for
read-only access in TSComputeRHSFunction()"), PETSc prints a traceback, and the
run **continues with the limiter doing nothing**. A wrong answer with a noisy
log, which is the worst of the three possible outcomes.

**Who is affected.** Only decks that set a `Limiter` key; `_haveLimiter` defaults
to false (`wxnodaldg2dmethod.cc:131`). Seven shipped decks do:

```
examples/unstructuredDG/euler/backwardFacingStep/backwardFacingStep.pin
examples/unstructuredDG/euler/explosionTest/explosion.pin
examples/unstructuredDG/euler/forwardFacingStep/forwardFacingStep.pin
examples/unstructuredDG/euler/forwardFacingStep/HLLflux/forwardFacingStep.pin
examples/unstructuredDG/euler/forwardFacingStep/Roeflux/forwardFacingStep.pin
examples/unstructuredDG/euler/scramjet/scramjet_inlet.pin
examples/unstructuredDG/euler/test/runtestcase.pin
```

No RMF or multifluid deck in `examples/.../rmf_frc/` sets one — §2 explains why
(the two-fluid limiter produces NaN), so the multifluid path is unaffected in
practice. `phase3/03-formation/formation.pin` has no `Limiter` key; the only
occurrences of the word are prose at lines 186-187.

**repro** — no run needed for the source facts:

```bash
grep -n 'applyLimiter(in,in)' src/subsolvers/wxnodaldg2dmethod.cc
grep -rln '^\s*Limiter\s*=' examples/          # the seven decks
grep -cw PETSC_USE_DEBUG $PETSC_DIR/include/petscconf.h  # 0 means optimized
```

The `-w` in the third command is load-bearing: without it the pattern also
matches `#define PETSC_USE_DEBUGGER "gdb"`, which every `petscconf.h` carries,
so a plain `grep -c PETSC_USE_DEBUG` reports 1 on an optimized build and reads
as "this is a debug build". Measured on PETSc 3.19.6: naive 1, `-w` 0.

On a PETSc ≥ 3.22, running any of the seven decks shows the traceback once per
RHS evaluation.

**Not fixed here.** The repair is ownership, not a version guard: limit into a
scratch vector and feed *that* to the rest of `step()`, never writing back into
`in`. That is a small change, but it alters the numbers every limited Euler deck
produces — on 3.21 and earlier they were being limited, and the fix keeps them
limited while changing *when* the limited state is visible within a stage — so it
needs those cases revalidated rather than a one-line patch. It also interacts
with §2: if the two-fluid limiter is ever fixed, this must be fixed first or the
multifluid path inherits the same silent no-op.

**Version note.** Nothing here depends on Apollo's own version guards, so
`petsc_compat.h` cannot help: the behaviour is a runtime check inside PETSc, not
a symbol that appears or disappears. `test/test_petsc_compat.py` will not catch
it, and neither will any compile-only gate.

---

## 18. The DG hot path runs on GSL's reference CBLAS, not OpenBLAS

Every DG volume and surface integral goes through one `cblas_dgemm`:

```c
// src/lib/wxcubature2d.cc:459
cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
            rows, meqn, cols,
            1.0, (double*)A, cols, (double*)x, meqn,
            0.0, (double*)y, meqn);
```

The binary links both `libgslcblas` (pulled in by GSL) and `libopenblas`, and
`libgslcblas` comes first in `DT_NEEDED`, so the dynamic linker binds
`cblas_dgemm` to GSL's **reference** implementation. `libopenblas` is loaded and
never used for this symbol.

**repro**:

```bash
objdump -p src/build-opt/apollo | grep NEEDED | grep -iE 'blas|gsl'
#   NEEDED   libgsl.so.27
#   NEEDED   libgslcblas.so.0     <- wins
#   NEEDED   libopenblas.so.0

LD_DEBUG=bindings LD_DEBUG_OUTPUT=/tmp/b ./src/build-opt/apollo --help >/dev/null 2>&1
grep -h cblas_dgemm /tmp/b.* | head -1
#   binding file ./src/build-opt/apollo [0] to
#   /lib/x86_64-linux-gnu/libgslcblas.so.0 [0]: normal symbol `cblas_dgemm'
```

**It costs nothing measurable, which was not the expectation.** The obvious
worry is speed, and it was measured rather than assumed: the RMF antenna gate
deck (162 steps) takes **88 s either way** - once as linked, once with
`LD_PRELOAD=libopenblas.so.0`, verified by `LD_DEBUG` to have actually
rebound the symbol. Apollo's matrices are per-element and small, so the call
overhead dominates and OpenBLAS's blocking and vectorisation have nothing to
work with. Do not "fix" this expecting a speedup; there isn't one here.

**What it does change is the numbers, and that is the useful part.** Swapping
only the BLAS - same binary, same commit, same PETSc, same machine, one rank,
same deck - moves the answer:

```bash
python3 scripts/run_fingerprint.py results-as-linked results-with-openblas
#   worst over 378 solution arrays: 2.310e-08  (frame 1, solutiondg.52)
#     driven by:      sorted_sum
#     absolute:       8.483e-14
```

That is a calibration worth having. It is what "the same code, one library
different" costs on this deck, so a cross-machine comparison showing a few
times that much needs no further explanation. When a PETSc 3.19.6 build here
was compared against a PETSc 3.25.5 build on a cluster the worst figure was
6.562e-08 - 2.8x this, from a comparison that also changed compiler, `-march`
(SConstruct defaults `arch` to `native`, so two hosts give two ISAs), and the
commit.

**Not fixed here.** Putting `-lopenblas` ahead of `-lgslcblas` would bind the
faster implementation, but it changes the numbers of every result in this
repository for no measured gain, so it is a revalidation decision rather than a
link-order tweak. Note also that `apversion.cc` prints a bare `BLAS` in the
banner and does not record which one, so a `.vtu` produced here cannot be
attributed to an implementation after the fact.

---

## 19. The formation run's instability grows at the antenna, which the mesh does not refine

`examples/.../rmf_frc/phase3/03-formation/formation.pin` dies with
`*** NaN after Mass matrix multiplication RHS ***` at t = 2.17e-06 s — 1.73 of
the 20 RMF periods it asks for, about 12 hours in on one core. `run_growth.py`
over its 21 frames names the components: φ 298×, E_y 186×, E_x 131×, against
E_z 1.06× (the physical RMF drive), ψ 1.49× and e-rho 1.25×. So the growth is in
the in-plane fields and the E-cleaning potential, and the drive itself is flat.

**It is at the antenna.** `scripts/where_peak.py` on that run's own frames:

```
frames whose peak is in the wall ring (r > 0.95 x WALL_RADIUS):
  E_x      0 of  21    E_y   0 of 21    E_z   0 of 21
  phi      3 of  21    psi   0 of 21    e-rho 0 of 21
```

Every peak sits at r = 0.0350–0.0398, through to the frame before death. The
winding itself spans `COIL_R` ± `COIL_W`/2 = 0.0335–0.0385, so the peaks sit on
it and a little outboard of it, drifting outward over the last frames
(E_x: 0.0375, 0.0375, 0.0375, 0.0389, 0.0375, 0.0375). φ reaches the wall only in
the last three, which is the failure spreading rather than its origin.

**And the antenna is outside the refined region.** The shipped mesh is generated
with `--ring 0.030`, so cells are `h-inner` = 9.3e-4 inside the plasma column and
grade to `h-outer` = 1.5e-3 beyond it. `COIL_R` = 0.036 is beyond it. Measured on
the shipped mesh — regenerated by `mkdiscmesh.py` with the command in the deck's
header, which produces a file **byte-identical** to the shipped `disc.msh`, so
these are measurements of the mesh the runs actually used:

| region | mean h [m] | λ_D/h | (c/ω_pe)/h | δ_resistive/h |
| --- | --- | --- | --- | --- |
| column, r < 0.030 | 1.026e-3 | 0.397 | 0.518 | 1.230 |
| winding, r = 0.036 ± 0.0025 | 1.426e-3 | 0.286 | 0.373 | 0.885 |
| **where the peak sits, r = 0.036–0.040** | **1.496e-3** | **0.272** | **0.355** | **0.844** |
| wall ring, r > 0.047 | 1.854e-3 | 0.220 | 0.287 | 0.681 |

with λ_D = 4.072e-4 m, c/ω_pe = 5.317e-4 m and δ = 1.262e-3 m for this deck's
reduced `LIGHT` and its `ETA`.

The deck's own header anticipates the problem and gives one of the two numbers
against the wrong cell size. It says the collisionless skin depth "is BELOW the
0.9-1.5 mm cells", which is fair — it quotes the range. But it then says the
resistive skin depth "is resolved at about 1.4 cells", and 1.4 is
1.262e-3 / 9.3e-4: the **column's** cell size, for a layer that lives at the
winding. At the winding it is 0.885 — under one cell. Its closing advice,
"Refine before quoting a penetration threshold", is right; what it understates is
that the region needing refinement is the one the mesh deliberately coarsens.

**Two refinement tests, and neither could detect an improvement.**

1. *Uniform*, scaling the whole mesh by k = 1, 2, 2.5, 2.75, 3, 3.25, 3.5, 4 —
   a factor 4 in h. Fitted φ growth rates, all under `-ts_ssp_type rk104` so the
   integrator of §20 plays no part, and all fitted past 3× the antenna ramp:
   2.21e6 (k=1, the real run), 1.57e6, 2.37e6, 3.57e6, 5.90e6, 2.22e6, 2.24e6,
   1.06e6 s⁻¹. Power-law fit against h: exponent −0.18 with r = **−0.15** — no
   correlation. The real run sits in the middle of the ladder's spread.
2. *Winding only*, three meshes identical except in the antenna annulus.
   `h-inner` is the same in all three, but `dt` is **not** quite: it is
   7.145e-11 in the first two and 7.655e-11 in the third, 7% larger, because
   moving the conforming ring changed which cell is the smallest. All three are
   run under `rk104`, where §20's amplification is absent at either value, and
   the rates below are per second rather than per step, so the difference does
   not enter:

   | arm | winding h | (c/ω_pe)/h | E_x rate | φ rate | final E_x | final φ |
   | --- | --- | --- | --- | --- | --- | --- |
   | coarse | 3.633e-3 | 0.146 | 1.84e6 | 1.11e6 | 1.01e3 | 3.93e-4 |
   | mid | 2.858e-3 | 0.186 | 1.77e6 | 1.57e6 | 1.62e3 | 3.92e-4 |
   | fine | 2.024e-3 | 0.263 | **2.13e6** | **1.77e6** | **8.02e2** | **1.28e-4** |

   The finest arm has the **highest** fitted growth rate and the **lowest**
   amplitude at the same physical time. Those disagree, and both differences are
   inside the factor-5.6 scatter of test 1.

**So refinement is untested, not disproved, and the reason is a flaw in both
tests**: every arm of both is deeply under-resolved — (c/ω_pe)/h never exceeds
0.52 and λ_D/h never exceeds 0.40. A sweep conducted entirely on one side of a
threshold says nothing about crossing it. (§20 records the same mistake made in
the other direction.)

**What crossing it would cost.** λ_D = 1 cell needs h = 4.07e-4, which is 2.28×
finer than the shipped `h-inner`: 5.2× the cells and 2.3× the steps, so **12×**
a run that already takes 12 hours. The cheap alternative — `--ring 0.040`, which
merely puts the antenna inside the existing refined region — costs +25% cells and
−0.8% on `dt`, about +26% in total, and moves (c/ω_pe)/h only from 0.355 to
~0.57. On the evidence above that is not enough to expect a difference. The
deck's conforming node ring at r = 0.030 is vestigial for this deck, whose
initial condition is uniform, so moving it has no other effect.

### The same deck has a *different* fatal mode at 3× the shipped cell size

Worth separating, because its growth table is indistinguishable and it is what
`where_peak.py` exists to tell apart. On a uniformly 3× mesh both integrators die
at t ≈ 5.48e-07 — 4752 steps with `rks2`, 4747 with `rk104`, five steps apart —
and there every peak is in **one cell at the conducting wall** (r = 0.0485),
whose electrons accelerate 1.2e5 → 2.1e5 → 5.4e5 m/s over three frames while the
cell empties (e-density 89.8 → 58.2 × floor). It is not a mesh defect: that mesh
is the best-conditioned of the ladder (minimum angle 35.9°, boundary ring 45.2°).
It is not a knife-edge: `OUT` = 15, 30 and 60 all die at t = 5.483–5.489e-07. And
it is not the first of a trend: the meshes at 2.75× and 3.25× both ran to
t = 1.4e-06, 2.6× past it, with no NaN. Why 3× specifically is **not explained**.

### What the structure is, measured on the surviving frames

Two diagnostics were written to separate a smooth, drive-locked pattern from
grid-scale structure the mesh cannot resolve, and to read the winding state the
column-only `rmf_diagnostics` bins cannot reach:

```bash
python3 scripts/ring_spectrum.py   03-formation/formation.pin <results>/*.vtu
python3 scripts/winding_anatomy.py 03-formation/formation.pin <results>/*.vtu
```

Run on a k=1 formation frame set (the `rk104`, `cfl` x4 proxy that reaches
0.45 of a period and reproduces the shipped `rks2` run's winding fields to four
digits) and on the k=2 and k=2.75 sets, they show the following. These are
measurements of those runs, not of the 140-hour run itself, whose frames are on
the cluster.

- **The winding density hole is shallow and decelerating, and the pressure is
  flat.** `n_e` at the winding falls to 0.971 of the column by 0.45 periods, in
  steps that shrink frame to frame (0.0074, 0.0098, 0.0052, 0.0033 per 0.11
  period); `T_e` rises 30.0 -> 31.4 eV at the hottest node and the electron
  pressure over the whole disc holds to +-1.3%. There is no local hot spot and
  no 1/n thermal runaway at this resolution: the hole is a bounded, quasi-static
  J x B / Ohmic depression, not the thing that grows.

- **What grows exponentially is grid-scale; the smooth part grows linearly.**
  On the winding node rings the m = 0 amplitude of `phi` and of the charge
  density `e(n_i - n_e)` grows **linearly** in time, while the r.m.s. of the
  cell-to-cell (high-m) remainder grows **exponentially**, at 3-7e6 s^-1 (k=1)
  and 3.9-4.8e6 s^-1 (k=2). Fitting a single exponential to `run_growth.py`'s
  global max|value| returns ~2e6 s^-1 - which is that linear m = 0 growth dressed
  as a rate. **The "2.2e6 s^-1 growth rate" of this run is a fit artefact, not a
  mode rate**; the real exponential is in the unresolved charge structure.

- **The charge the fields carry is increasingly inconsistent with Gauss's law.**
  The per-cell residual `div E - rho_c/eps0`, normalised by `div E`, grows from
  0.28 to 0.64 over the run: the electrostatic field and the charge density it
  should satisfy drift apart at the cell scale, which is exactly what an
  under-resolved Debye layer (lambda_D/h = 0.3 here) does. `phi`, whose only
  source is that residual, records it - it leads `E` but is a symptom, and the
  in-plane `E` growth is independent of `DIVE_SPEED` over the first 0.15 period.

- **The peak sits at fixed mesh azimuths, not on the rotating drive.**
  `where_peak.py` (now printing theta) shows `e-rho` pinned at theta ~ -120 deg
  and `phi` at theta ~ 0/-120 deg for frame after frame while the antenna
  rotates - the ring-count stitch seams of the generated mesh (`mkdiscmesh`
  closes each ring pair near theta = 0 and emits mismatch cells at the winding
  near +-120 deg). The failure nucleates where the mesh is locally worst.

- **The electron-fluid state limit is the endgame, not the cause.** The
  return-current drift `u_ez` stays at 0.12 of the electron sound speed and
  `|u_e| + a_e` at 1.12 c0 through 0.45 periods, with no node reaching a
  non-positive pressure or the density floor. Those limits (which the amplitude
  ladder hits at 1.9-3.8 a_e, and which end the k=3 wall cell) are how the run
  finally NaNs once a cell has hollowed, not why the structure grew.

So on the evidence available the growth is a **numerical charge-separation
structure at the unresolved Debye scale, seeded by a bounded physical hole and
sited by the mesh seams**, not a resolved plasma instability - but that is a
hypothesis with a decisive test that has not been run (below), not a conclusion.

**A correction this measurement forces.** The "clean at 4x `cfl`" figure once
read as a spatial-CFL margin is not one: at `cfl` x4 the default `rks2`
integrator amplifies the electron oscillation (net +5.8e-3 per step, §20) and a
run dies at ~1100 steps from that, while the same mesh under `rk104` runs clean
- so the number measured the integrator, not the flux. The spatial CFL margin is
>= 4x under `rk104` and unmeasured above it.

**What is still open.** Whether resolving `lambda_D` removes the growth is
untested: every mesh tried stayed at `lambda_D/h` <= 0.42. **The tanh-ramp
reproducer is NOT the cheap version of this test, and §22 is why**: that death
is hydrodynamic, survives deleting every electromagnetic source, and is
insensitive to the mesh, so it can neither confirm nor refute what is happening
at the antenna. It was the obvious cheap proxy and it is the wrong one. The
remaining test of this hypothesis is refinement of the antenna deck itself.
The superseded reasoning kept below for the record: the antenna-off tanh-ramp
reproducer (the deck header's own "dies in 30-260 steps with or without the
antenna" case) on a small disc across
`lambda_D/h` = 0.44, 0.68, 1.02, 1.63; if the deaths retreat as the layer is
resolved, the mechanism is the unresolved sheath and the cure is resolution or a
positivity-preserving scheme; if they do not, it is the two-fluid model's own
vacuum-expansion or the Lax-Friedrichs species asymmetry. Neither refinement nor
a physical-mode search (an azimuthal m-spectrum with a real frequency) has been
done at `lambda_D/h` >= 1.

### repro

```bash
# Regenerate the shipped mesh and coarsened versions of it (k scales both h).
python3 scripts/mkdiscmesh.py -o k1.msh --outer 0.05 --ring 0.03 \
        --h-inner 9.3e-4 --h-outer 1.5e-3          # == the shipped disc.msh
python3 scripts/mkdiscmesh.py -o k3.msh --outer 0.05 --ring 0.03 \
        --h-inner 2.79e-3 --h-outer 4.5e-3

# Run the deck on one of them, then ask run_growth WHAT grew and where_peak WHERE.
python3 scripts/run_growth.py <results>/*.vtu
python3 scripts/where_peak.py 03-formation/formation.pin <results>/*.vtu
```

`where_peak.py` needs the deck because `WALL_RADIUS`, `COIL_R`, `COIL_W` and
`RAD_PLASMA` are deck values; without one it prints radii and declines to name
regions.

---

## 20. The default PETSc time integrator amplifies every oscillatory mode, and only resistivity hides it

`wxpetsctimestepping.h:33` sets `TSSetType(solver, TSSSP)` and leaves the scheme
to PETSc's default, which is `TSSSPRKS2` with five stages (`ssp.c`: `TSCreate_SSP`
ends with `TSSSPSetType(ts, TSSSPRKS2)`, which sets `nstages = 5` and
`default_adapt_type = TSADAPTNONE`). Transcribing that stepper's vector
operations gives

    R(z) = ((s-1)/s + z/s)(1 + z/(s-1))^(s-1) + 1/s,   s = 5

and on the imaginary axis, in exact rational arithmetic,

    |R(iy)|² − 1 = y⁴/32 + y⁶/640 + y⁸/20480 + y¹⁰/1638400

Every coefficient is strictly positive, so **|R(iy)| > 1 for every y ≠ 0**: the
imaginary-axis stability interval is exactly zero. That is a proof, not a
numerical scan. `TSSSPRK104` is the opposite — stable out to y = 4.9215.

The two-fluid decks carry an undamped electron oscillation at ω_pe, and the
`currents`/`lorentzForces` source pair is integrated inside the same Runge–Kutta,
so it passes through that stability function. On the shipped formation mesh
ω_pe·dt = 5.642e9 × 3.57232e-11 = 0.2015, giving +2.58e-05 per step.

**Resistivity beats it, by about ten times.** `ETA` = 0.5e-5 damps the same
oscillation at ν/2 = n e² η / (2 mₑ) = 7.04e6 s⁻¹, which is −2.51e-04 per step at
the same `dt`. Net: **−2.26e-04 per step — damped.** Setting amplification equal
to damping gives a threshold at 2.14× the shipped cell size, since coarsening by
k multiplies ω_pe·dt by k (growth ∝ k⁴) but damping only by k.

Measured across that threshold, with everything else held fixed:

| mesh | ω_pe·dt | net rate/step | `rks2` | `rk104` | E_x growth |
| --- | --- | --- | --- | --- | --- |
| 2× shipped | 0.403 | −8.7e-05 (damped) | survived 11112 steps | survived 11112 steps | 9.2× / 9.2× |
| 3× shipped | 0.655 | +2.1e-03 (grows) | NaN at 4752 | NaN at 4747 | 636× / 79× |

At 2× the two integrators are **identical to every digit printed** — E_x 9.2×,
E_y 11.7×, φ 16.6×, final values 1.619e+03 / 1.134e+03 / 3.920e-04 in both.

Read the 3× row carefully: `rk104` dies there too, five steps from `rks2`, but
*not* from this mechanism — it grows E_x 8× less and still dies at the same
physical time, because at 3× the deck has the separate wall-sheath failure of
§19. The integrator claim in that row is the 636× against 79×, not the NaN. The
shipped mesh is finer still. So on the decks as shipped the integrator choice
changes nothing, and `-ts_ssp_type rk104` will not rescue a formation run; it
costs 2.1× per step (ten RHS evaluations against five) for numbers that do not
move.

Above the threshold it is exactly as advertised: `rks2` died at 299, 661, 2221
and 4752 steps on 6×, 5×, 4× and 3× meshes, and the ratio of the first two
matches the ratio of predicted rates to 6.5%. `rk104` ran the 6×, 5× and 4×
cases to t = 1.2e-06 with no NaN.

**This matters to anyone who coarsens the mesh, lowers `ETA`, or raises `LIGHT`.**
All three move the balance the wrong way, and the failure looks like physics.
The ladder above spans 3× to 6× — entirely above the threshold — which is why it
cannot say anything about the shipped mesh except through the threshold
calculation and the 2×/3× crossing.

---

## 21. Apollo's argument parser rejects every PETSc option

`apsimulation.cc:228` parses `argv` with `getopt_long(argc, argv, "hi:r:o:", ...)`
and its `default:` case calls `usage(); exit(2)`. So anything PETSc-shaped dies
before `PetscInitialize` is reached:

```
$ apollo -i formation.inp -ts_ssp_type rk104
apollo: invalid option -- 't'
*** Welcome to Apollo ***
... usage ...                                   # exit 2
```

That is every `-ts_*`, `-dm_*`, `-log_view`, `-options_left` and `-malloc_debug`.
`--` does not help: `getopt` stops there and the remainder becomes a positional,
which hits the "unexpected argument" exit on the next line. `run_one.sh` has no
`"$@"` passthrough either.

Two routes work, both read by `PetscInitialize`:

```bash
PETSC_OPTIONS="-ts_ssp_type rk104 -ts_view" apollo -i formation.inp
printf -- '-ts_ssp_type rk104\n-ts_view\n' > .petscrc   # CURRENT directory
```

**The confirmation trap.** The only integrator line Apollo logs is
`Time integration done using: ` + `TSGetType`, printed in the constructor —
before `TSSetFromOptions` is called at `wxpetsctimestepping.h:81` — and
`TSGetType` returns `"ssp"` for every scheme anyway. So a dropped, misspelled or
unexported `PETSC_OPTIONS` produces a log **identical** to a correct run, and the
run silently uses the default. `-ts_view` is what distinguishes them:

```
TS Object: 1 MPI process
  type: ssp
    Scheme: rk104
```

Refuse to interpret any run whose log does not carry the `Scheme:` line you
expected. The §20 comparison was set up through the command line first and would
have returned a confident false result — both arms `rks2`, both growing, read as
"the integrator is exonerated" — had the option not been rejected loudly.

---

## 22. The plasma–vacuum ramp fails on positivity at the expansion front, not on the Debye length

`03-formation/formation.pin`'s header explains at length why the deck has no
vacuum annulus, and the explanation is wrong. It says:

> Apollo carries explicit charge separation, so the electron dynamics have two
> scales that have to be resolved. [...] A density ramp at the column edge
> excites exactly them: the runs diverge at the ramp, near r = 0.035 [...]
> A hard step fails on the first residual evaluation; a 1.5 mm tanh ramp lasts
> about 30 steps, a 3 mm ramp about 150, a 5 mm ramp at least 244.

The step counts are right. The attribution to λ_D and ω_pe is not, and neither
is "at the ramp". Ten runs of the ramp reproducer, on one core each, say so.

**The ramp death does not care about the Debye length.** Holding the ramp width
at 3.0 mm and refining the mesh across the threshold the whole §19 question
turns on:

| arm | λ_D/h (from the run's own dt) | triangles | dies at |
| --- | --- | --- | --- |
| L1 | 0.45 | 16 809 | 5.3475e-09 |
| L2 | 0.76 | 31 456 | 5.3219e-09 |
| L3 | **1.14** | 60 497 | 5.2913e-09 |

A 2.3× refinement, 3.6× the cells, crossing λ_D/h = 1 — and the death time moves
by **1%**. Against that, halving the ramp width at the *finest* mesh (C1,
λ_D/h = 1.14, W = 1.29 mm) brings the death forward by **42%**, to 3.0613e-09.
The clock is set by a physical length, the ramp width, not by the cell.

**It does not need the electromagnetic sources at all.** Deleting the whole
`Sources` list — no Lorentz force, no currents, no charge sources, no
resistivity, so no charge separation anywhere in the run — still dies, at
3.1544e-09. A failure that survives the removal of every term that could
produce charge separation is not a Debye-shielding failure. (An empty
`Sources = []` does not parse; the key has to be deleted.)

**It does not care about the antenna, and it is not at the wall.** With
`Bomega = 0.0` and `rmfAntenna` dropped it dies at 5.3215e-09 against
5.3228e-09 with the drive running — the header's "with or without the antenna"
was asserted but never run, and it holds to three figures. Moving the
conducting wall from r = 0.060 to r = 0.090 changes the death time by **0%**
(5.3475e-09), so the wall-reflected expansion is not it either; that hypothesis
was tested and refuted rather than assumed.

**Where it actually fails, and why.** The negative-pressure nodes appear at the
same radius and the same *local density* whatever the mesh or the wall:

| arm | wall | fails at r | local n/n₀ there |
| --- | --- | --- | --- |
| G0 | 0.060 | 0.0434–0.0440 | 6.2e-05 – 1.8e-04 |
| L3 | 0.060 | 0.0425–0.0436 | 1.6e-04 – 3.4e-04 |
| C4 | 0.090 | 0.0440 | 1.0e-04 – 1.2e-04 |

That local density is `VAC_FRAC` = 1e-4: the failure is where the expansion
front arrives at the undisturbed background, travelling at 2.63e6 m/s = 0.89 of
the electron sound speed. There the flow is fast and the gas is thin, so
`p = (γ-1)(E - ρv²/2)` is a difference of nearly equal numbers — which is
exactly what the deck's own `constantResistivity` comment warns about, in a
scheme with no limiter and no positivity fix. The ramp itself is *static*: its
θ-mean density profile is unchanged from the first frame to the last.

**The decisive test is an intervention, not a correlation.** Raise the
background from 1e-4 to 1e-2 of the column and the run **survives to 3.0e-08**,
840 steps, no NaN — 5.6× past the baseline death, on the same mesh with the same
ramp. Nothing about λ_D changed; the thing the front expands into did.

**What this costs and what it buys.** The vacuum annulus that Phase 1 wanted is
not blocked by the Debye length, so the "three times the resolution for nine
times the cells" the header offers as the price would have bought nothing. It is
blocked by the absence of a positivity-preserving scheme (§2: the only two-fluid
limiter produces NaN). A deck that needs a tenuous region can have one today by
keeping it above roughly 1% of the column density.

**This does NOT explain §19.** The antenna instability is a different failure and
the two are easy to tell apart with `winding_anatomy.py`: the formation run's
Gauss residual `|div E - ρ/ε₀|/|div E|` grows 0.28 → 0.64 and c₀φ/E⊥ sits at
0.25–0.33, while this ramp case stays at 0.06 → 0.12 and 0.011 → 0.045. The ramp
death is hydrodynamic; the antenna death is electrostatic. §19's standing
hypothesis is untouched by this result — it only loses the ramp deaths as
supporting evidence, which were never evidence for it.

The superseded explanation is still in `antenna/frc2d.pin`'s header (and so in
every deck generated from it) and in the assessment's §3.10 discussion.

### repro

```bash
# the ramp reproducer, antenna off, on a uniform-in-the-ramp mesh
python3 scripts/mkdiscmesh.py -o disc.msh --outer 0.060 --ring 0.035 \
        --h-inner 9.3e-4 --h-outer 1.5e-3
# edit a copy of antenna/frc2d.pin: Bomega = 0.0, drop rmfAntenna from Sources,
# EDGE_W = 3.0e-3, VAC_FRAC = 1.e-4, TEND = 3.0e-8, OUT = 60
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 apollo -i frc2d.inp
python3 scripts/winding_anatomy.py frc2d.pin frc2d_*.vtu   # neg-p nodes, Gauss residual
```

`winding_anatomy.py` reports the negative-pressure nodes **two frames before the
NaN** (486 of them), so this failure announces itself if anything is looking.
