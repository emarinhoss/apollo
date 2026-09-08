# Known issues

Findings from an end-to-end review that are **not** fixed, with enough detail to
act on. Everything here was confirmed by reading the code and, where a `repro`
block is given, by running it — none of it is speculation. Items are ordered by
how much damage they can do to a result.

Fixed issues are not listed; see the git history.

---

## 1. Multifluid collisional–radiative source terms

### 1.1 The electron–ion thermal equilibration mass ratio looks inverted

`src/hyperapps/multifluid/wxbraginskiifrictionsrc.h:209` and
`src/hyperapps/multifluid/ap2dmomentumxfer.h:193`, identically:

```cpp
REAL Q_delta = 3*_mi/_me*ne*nue*(Pe/ne-Pi/ni);
```

Braginskii's electron–ion energy equilibration is

> Q_Δ = 3 (m_e / m_i) n_e ν_e (T_e − T_i)

i.e. the coefficient is `3*m_e/m_i`. As written the ratio is inverted, which for
hydrogen is a factor of (m_i/m_e)² ≈ 3.4 × 10⁶ — equilibration would be
effectively instantaneous instead of the slowest process in the system.

**Not changed here** because this module has no test coverage and no reference
solution in the repository, so a change could not be validated. Before changing
it, check the coefficient against the resistivity-form siblings in the same
directory (which carry an additional 1/0.51) and confirm the intended
normalisation of `nue`.

### 1.2 The friction energy terms do not conserve energy (error ∝ ion velocity)

All four friction sources (`constantResistivity`, `anisotropicResistivity`,
`braginskiiFriction`, `MomentumXfer2D`) share the same energy terms. With
`q[0..4]` the electron fluid, `q[5..9]` the ion fluid, `u = ue - ui`, and
`R = -alpha*u` the friction on electrons:

```cpp
s[3] = -(Rux*ui+Ruy*vi+Ruz*wi)-Q_delta;   // electron energy
s[7] = Q_delta;                            // ion energy
```

For total-energy variables the sources are `s_e = R·u_e + Q_e` and
`s_i = -R·u_i + Q_i`, with `Q_e + Q_i = -R·(u_e - u_i) = alpha|u|²` the frictional
heat, which Braginskii deposits in the electrons (`Q_e = alpha|u|² - Q_delta`,
`Q_i = Q_delta`). Substituting: `s_e = +R·u_i - Q_delta`, `s_i = -R·u_i + Q_delta`,
summing to zero. So the code has the right structure, the wrong sign on the
electron work term, and no ion work term; total energy acquires a spurious source
`-R·u_i`.

An earlier version of this note said Ohmic heating was absent altogether. That
was wrong for a total-energy formulation: with ions at rest the code and the
correct form agree (both `-Q_delta`), because friction thermalises the electron
drift *within* the electron fluid — the `alpha|u|²` heating is exactly the
kinetic energy the momentum source removes. The error is proportional to `u_i`:
negligible on the 2 μs `rmf_frc` runs (ion gyroperiod 11 μs), and growing as the
ions spin up, which the RMF literature says they do quickly.

**Fix** (proposed in `docs/rmf-frc-model-assessment.md`, Phase 0): in all four
files set `s[3] = +(Rux*ui+Ruy*vi+Ruz*wi) - Q_delta` and
`s[7] = -(Rux*ui+Ruy*vi+Ruz*wi) + Q_delta`, and add a test that a uniform box
with a relative drift and no fields conserves total energy to round-off.
**Not changed here** because the module has no test to protect the change; the
test should land with it.

### 1.3 The bundled collisional–radiative example does not run

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

## 2. `Numerical_Flux = Wave` aborts

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

## 3. Roughly two thirds of `src/lib` is not built

A transitive include closure from the 89 sources the SConscripts actually
compile reaches 222 files. **129 of the 188 files in `src/lib` are in none of
them.** Among the unreachable:

| Group | Files | Note |
| --- | --- | --- |
| NDG++ value types | `Vec_Type.h`, `Mat_COL.h`, `VecObj_Type.h`, `MatObj_Type.h`, `ArrayGen.h`, `ArrayMacros.h`, `Region1D/2D.h`, `MappedRegion1D/2D.h`, `Index.h`, `Registry_Type.h`, … | Cannot link even if included: `umERROR`, `umWARNING`, `gVecData` are declared and never defined |
| HDF5 I/O | `wxhdf5io.cc/.h`, `wxhdf5iotmpl.*`, `wxhdf5traits.*`, `wxiobase.cc` | Why there is no checkpointing (§4) |
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

## 4. There is no checkpoint/restart

`--restart` now fails with a message rather than being silently ignored, but the
feature itself does not exist: the load path in `apolloMain()` is commented out,
and the HDF5 I/O layer that would back it is in the unreachable set above (§3).
For a code that runs multi-day fusion simulations this is the largest missing
capability.

---

## 5. Error paths bypass MPI

67 `exit()`/`abort()` calls remain in `src/`, against 4 `MPI_Abort`. Calling
`exit()` from one rank leaves its peers blocked in whatever collective they were
in; the job hangs until the scheduler kills it, and any diagnostic is lost.
Several sit inside per-element loops (the NaN guards in `wxnodaldg2dmethod.cc`
and `wxeulereqn.cc`), which is also the worst place to call it from.

The entry point now aborts through MPI correctly; the fix for the rest is to
throw a `WxExcept` and let it reach that handler.

---

## 6. PETSc return codes are discarded

555 PETSc calls, 63 `CHKERRQ`. A failing `DMPlexCreateFromFile`, `VecGetArray`
or `MatAssemblyEnd` is not noticed, and execution continues with an object that
was never created. `ApSolver::createMesh` and `ApSolver::OutputVTK`
(`src/solvers/apsolver.cc`) additionally open with `PetscFunctionBeginUser` and
never call `PetscFunctionReturn`, so they do not appear correctly in PETSc's own
stack traces either.

Converting these to `PetscCall(...)` is mechanical but touches most of the
codebase, so it wants to be its own change.

---

## 7. Deprecated PETSc APIs

`src/subsolvers/wxpetsctimestepping.h` calls `TSSetDuration`,
`TSSetInitialTimeStep` and `TSGetTimeStepNumber`, all deprecated since PETSc 3.8
(2017), and `PETSC_NULL`, deprecated in 3.19. They still compile and account for
151 of the build's warnings. Replacements are `TSSetMaxSteps`/`TSSetMaxTime`,
`TSSetTimeStep`, `TSGetStepNumber` and `PETSC_NULLPTR`.

---

## 8. Build warnings

A clean `scons build-opt` emits 786 warnings (down from 821 at the start of the
review; the uninitialized-variable reports specifically went from 35 to 5).

| Category | Count | Comment |
| --- | --- | --- |
| `-Wsign-compare` | 185 | `unsigned` loop counters against `int` bounds; mostly benign, individually cheap to fix |
| `-Wmisleading-indentation` | 175 | Worth reading each one: indentation that lies about control flow is how a guard silently stops guarding |
| `-Wdeprecated-declarations` | 151 | §7, plus `std::auto_ptr` in the vendored muParser |
| `-Wunused-variable` | 93 | |
| `-Wreorder` | 72 | Member initialiser lists that do not match declaration order — currently harmless, a real bug the moment one member's initialiser reads another |
| `-Woverloaded-virtual` | 47 | A derived `init()` hiding the base's `init(PetscReal, Vec)`. §9 |
| `-Wmaybe-uninitialized` | 5 | |
| `-Wformat-security` | 4 | `sprintf(dst, variable)` in the CR module; a `%` in the data is a crash |

`-Werror` is not realistic yet, but `-Wextra` on new code would be.

---

## 9. `init()` signature mismatches hide the base virtual

`WxObject::init(PetscReal, Vec)` is virtual; `ApFVM2dScheme::init()` and
`ApDomainDecompCheck::init()` declare a no-argument `init()`, which hides rather
than overrides it. The base's do-nothing body runs instead, so those schemes'
initial conditions are never applied. `apfvm2dscheme.cc` is currently commented
out of `src/subsolvers/SConscript`, so this is latent — but it is exactly the
trap that `override` exists to catch, and the codebase uses it nowhere.

---

## 10. Duplicated and uncertain post-processing scripts

Several files in `scripts/` are near-duplicates with no indication of which is
canonical: `wxdata.py` / `wxdata_3949.py`, `wxunsdgdata.py` / `wxunsdgdata2.py`,
`SGC.py` / `SGC_MB.py`, and the four-way `wxxdmf*.py` family. All of them now
parse and import, so a maintainer can diff them and delete the losers.

`scripts/SGC.py` also carries genuinely dead code: `Node.getLocalID` referenced
an `owners` attribute that is never assigned anywhere (it is now an explicit
`NotImplementedError`), and `buildLineFaces`/`buildRectangleFaces`/
`buildHexahedronFaces` are unreferenced.
