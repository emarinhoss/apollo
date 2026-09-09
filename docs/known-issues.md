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

## 6. There is no checkpoint/restart

`--restart` now fails with a message rather than being silently ignored, but the
feature itself does not exist: the load path in `apolloMain()` is commented out,
and the HDF5 I/O layer that would back it is in the unreachable set above (§5).
For a code that runs multi-day fusion simulations this is the largest missing
capability.

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
throughout, so this is specifically the magnetic cleaning potential. It is the
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
