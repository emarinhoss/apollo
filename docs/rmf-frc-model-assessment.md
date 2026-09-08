# RMF-FRC model assessment

An assessment of `examples/unstructuredDG/multifluid/rmf_frc` — the rotating-magnetic-field
(RMF) field-reversed-configuration (FRC) formation problem — against the published physics,
with a plan for closing the gaps.

**How this was done.** The deck, the boundary conditions, the source terms and the meshes were
read in full and the deck was run. The literature was surveyed by web search; the network
policy of the environment this was written in blocks every publisher and preprint host
(AIP, Cambridge, IOP, arXiv, OSTI, ADS), so statements attributed to a paper below come
from its abstract or from sentences quoted in search results, not from full text. Where a
formula is stated from the standard literature rather than a fetched source, that is said.
The dimensionless numbers are computed from the deck's own values.

Verdict in one paragraph: **the model class is right, the parameters are in the right regime,
and the formation scenario is exactly the one the literature studies — but the way the RMF is
applied at the boundary is not physically consistent with how an RMF interacts with a plasma,
and the friction source terms do not conserve energy.** The first is the item that most limits
what the example can be compared against; the second is small on the two-microsecond runs
shipped here and grows as the ions spin up.

---

## 1. What the literature says the problem is

**Foundational theory.** The RMF current-drive mechanism was proposed by Blevin and Thonemann
(1962) and worked out for a plasma cylinder by Jones and Hugrass, *J. Plasma Phys.* **26**, 441
(1981) — "steady-state solutions of magneto-fluid equations show that, provided the amplitude
and rotation frequency of the field are suitably chosen, the penetration is not limited by the
usual classical skin effect. The enhanced penetration of the rotating field is accompanied by
the generation of a unidirectional azimuthal electron current which is totally absent in a
purely resistive plasma cylinder" — and numerically by Hugrass and Grimm, *J. Plasma Phys.*
**26**, 455 (1981), who established that "classical RMF theory predicts a threshold condition for
magnetic field penetration into the plasma." The mechanism requires the frequency ordering
ω_ci ≪ ω ≪ ω_ce: electrons are magnetised and dragged around by the rotating transverse field,
ions are not.

The penetration parameters are (form standard in this literature; the exact threshold curve
is in Hugrass & Grimm 1981 and, empirically, Milroy 1999):

- γ = ω_ce / ν_ei — how magnetised the electrons are against collisions;
- λ = r_s / δ, with δ = √(2η / μ₀ω) the classical skin depth — how many skin depths the plasma
  is wide;
- penetration when γ is large compared with λ (order unity ratio); the Michigan thruster papers
  state the same dependence on "the electron cyclotron frequency, electron-ion collision
  frequency, thruster radius, and the classical skin depth of the RMF."

Once penetrated, the electrons rotate near-synchronously, J_θ ≈ −n e ω r, and Guo, Hoffman and
Milroy (*Phys. Plasmas* **14**, 112502, 2007) define ζ as "the ratio of average electron rotation
frequency to RMF frequency." Milroy's fixed-ion model (*Phys. Plasmas* **6**, 2771, 1999) gives
"empirical expressions … to characterize the critical RMF magnitude required for full
penetration and the rate of RMF penetration", and finds that "in the presence of strong
anisotropic plasma resistivity, the direction and magnitude of the axial bias field can have a
strong influence on the penetration."

**The formation problem — the deck's scenario.** Milroy's r-θ MHD model (*Phys. Plasmas* **7**,
4135, 2000): "For the formation problem, a RMF is applied to a plasma column with an initially
uniform axial magnetic field and background plasma density, and the RMF-induced current
reverses this bias field, forming a FRC." Experimentally, Guo et al. (*Phys. Plasmas* **9**, 185,
2002): "The RMF creates an FRC by driving an azimuthal current which reverses an initial
positive bias field. The FRC then expands radially, compressing the initial axial bias flux and
raising the plasma density, until a balance is reached between the RMF drive force and the
electron–ion friction." That is precisely what `frc2d.pin` sets up: uniform column, uniform
bias field, RMF switched on, flux conserver outside.

**Screening — the point that matters for the boundary condition.** From the TCS programme
(Hoffman, Guo, Slough et al.): "The uniform transverse RMF in vacuum is shielded by the
conducting plasma, resulting in a mostly azimuthal field near the FRC separatrix with a very
small radial component." The field at the plasma edge is *not* the vacuum RMF; it is the vacuum
RMF plus the field of the currents the plasma drives in response, and that sum is what the
plasma actually sees. Modern simulations therefore put the drive in through the antenna, not
at the plasma edge: Milroy, Kim and Sovinec (NIMROD, *J. Comput. Phys.* **195**, 355, 2004, and
later) added "boundary conditions to capture the effects of a finite length RMF antenna" and
"modifications to the radial boundary conditions [that] capture most of the effects of
multiple discrete coils", noting that "the Hall term is a zeroth order effect."

**Energetics and torque.** Hoffman et al. (*Phys. Plasmas* **13**, 012507, 2006): "A balance
between the RMF applied torque and electron-ion friction will determine the peak plasma
density"; "measurements of total absorbed power and comparisons of applied RMF torque to torque
on the electrons due to electron-ion friction allowed the separation of classical Ohmic and
anomalous heating to be inferred." Power "is absorbed by the plasma due to oscillating axial
currents (which create the azimuthal torque), proportional to B_ω², and due to the azimuthal
FRC currents, proportional to B_e²." Guo et al. 2007: "A large fraction of the RMF power is
absorbed by an anomalous mechanism directly proportional to the square of the RMF magnitude."

**Ion spin-up.** "The RMF torque on the electrons is quickly transferred to the ions, but ion
spin-up is limited in these low density experiments, presumably by ion-neutral friction";
"the ion rotation is determined by a balance between electron-ion friction, the end shorting
effect, and ion drag against neutrals"; "ion spin-up can substantially reduce or cancel the RMF
current drive effect" (TCS/TCSU papers; ions reached ~6 kHz against a 180 kHz drive).

**Two-fluid and multi-fluid precedent.** Belova et al. (HYM, hybrid): "lower plasma density and
larger RMF amplitudes result in faster RMF field penetration, in agreement with previous
two-fluid studies." Sousa (APS DPP 2016, "Rotating Magnetic Field FRC Formation Studies using
the Multi-Fluid Plasma Model"): "aspects of the FRC formation physics using a rotating
magnetic field (RMF) at low power are simulated using a multi-fluid plasma model, with results
compared with experimental observations with emphasis on the development of instabilities and
robustness of the field reversal" — the direct precedent for this example.

**Thrusters — the xenon deck's context.** The Michigan RMF thruster (Woods, Gill, Sercel, Jorns;
2021–2024): radius R = 10 cm with a measured classical skin depth of 1 cm, T_e ≈ 9 eV, xenon,
~500 kHz, induced currents to 2500 A "sufficient to form a FRC plasmoid"; the Lorentz force
"contributes ∼25% of measured thrust"; and "the RMF may not penetrate the plasma as expected due
to screening caused by a combination of collisionality and classical skin depth." Thrust is an
axial effect: an r-θ cross-section can model penetration and current drive but not the
acceleration.

---

## 2. What Apollo implements

Described in full in the conversation that produced this document; in brief: 18 unknowns per
node (electron and ion Euler fluids, perfectly-hyperbolic Maxwell), coupled by Lorentz-force,
current and charge sources plus `constantResistivity` friction; speed of light reduced 100×;
the mesh is the plasma disc r ≤ a = 3 cm; the RMF enters through `twoFluidSimplifiedRMFBC`
at r = a, which imposes the vacuum rotating field, an axial E, a conducting-wall in-plane E, and
B_z from an axial flux conserver at b = 3.5 cm evaluated from ∫B_z dA over the plasma.

---

## 3. Assessment

### 3.1 Model class — appropriate ✓

A two-fluid model with full Maxwell contains the Hall physics that the literature identifies as
zeroth-order for RMF drive, without an Ohm's-law closure; electron inertia and displacement
current are retained. This is the model class Sousa 2016 used for the same problem and the one
Belova cites two-fluid results against. Nothing to change.

### 3.2 Regime — right window ✓

From the hydrogen deck (n = 10²⁰ m⁻³, T_e = 30 eV, B_bias = 60 G, B_ω = 50 G, f = 795 kHz,
a = 3 cm, η = 5 × 10⁻⁶ Ω m):

| quantity | value | reading |
| --- | --- | --- |
| ω_ci, ω, ω_ce | 5.7 × 10⁵, 5.0 × 10⁶, 1.1 × 10⁹ rad/s | ω_ci ≪ ω ≪ ω_ce holds (×8.7, ×211) |
| γ = ω_ce/ν_ei | 75 | ν_ei from η: 1.4 × 10⁷ s⁻¹ |
| δ, λ = a/δ | 1.26 mm, 24 | γ/λ ≈ 3: penetration expected |
| η | 5 × 10⁻⁶ Ω m | Spitzer at 30 eV, lnΛ = 10: ~3 × 10⁻⁶; plausible |
| μ₀ n e ω a²/2 | 450 G | synchronous reversal would be 7.5× the bias: strong drive |
| B_ω/B_bias | 0.83 | TCS-like ratio |
| β (T_e + T_i vs bias) | 50 | the bias is weak; the column is essentially β ≫ 1 |
| ion gyroperiod | 11 μs | longer than the 2 μs run: ions barely respond |
| run length | 1.6 RMF periods | rise time 30 ns ≪ period: an abrupt switch-on |

Two remarks. β ≈ 50 against the bias field means the initial column is far from any
equilibrium with that field; that is normal for the formation problem, but the early
transient is violent and the run is short enough that it is mostly transient. And 1.6 periods
is too short to see the torque–friction balance the literature describes; this deck exercises
penetration and initial reversal, not sustainment.

### 3.3 The formation scenario — matches the literature ✓

Uniform column, uniform bias, RMF applied, flux conserver: Milroy 2000's formation problem and
Guo 2002's description, with the flux conserver providing the compression feedback.

### 3.4 The RMF boundary condition — not physically consistent ✗

This is the central finding. `twoFluidSimplifiedRMFBC` sets, at the plasma edge r = a, the
ghost state (B_x, B_y) = −B_t(t)(sin ωt, cos(ωt+φ)) — the **vacuum** rotating field.

1. **It prescribes the field the plasma is supposed to be screening.** The TCS observation
   above is that the field at the plasma edge is "mostly azimuthal … with a very small radial
   component" because the plasma's own induced currents cancel much of the applied field. A
   boundary condition that pins B⊥ at the edge to the applied value forbids that response
   from reaching the edge: it is equivalent to an antenna of infinite stiffness located *at* the
   plasma surface, one that cannot be loaded by the plasma. The consequences are that
   penetration is over-driven relative to a real antenna at 3.6 cm, the penetration threshold
   in B_ω is not the literature's threshold, and the antenna power (torque × ω) is not
   computable from the simulation because the antenna current is never represented. Every
   published model that aims to be quantitative puts the drive at the coils and lets the field
   at the edge be an outcome (Milroy 2000; NIMROD's antenna and discrete-coil boundary
   conditions).

2. **The axial field E_z has no azimuthal structure.** For a spatially uniform rotating
   transverse field, Faraday's law with E = E_z ẑ gives ∂E_z/∂x = ∂B_y/∂t and
   ∂E_z/∂y = −∂B_x/∂t, so E_z = x Ḃ_y − y Ḃ_x (+ const): a field that rotates with B⊥. The
   code imposes E_z = r·(−Ḃ_t cos ωt − B_t ω sin(ωt+φ)), which is that expression's magnitude
   with the angular dependence dropped (the unused `theta = atan(y/x)` in the same function is
   the residue). The boundary E_z drives the oscillating axial currents that — per Hoffman —
   "create the azimuthal torque", so their θ-structure is the torque's structure. This one is
   cheap to make exact.

3. **The hyperbolic system is over-specified at the boundary.** At a face of the Maxwell
   system, only the incoming characteristics may be imposed; the outgoing ones must come from
   the interior. The ghost state sets E_z, B_x, B_y, B_z all to prescribed values, and the
   Lax–Friedrichs flux then averages ghost and interior. The net effect is a partially
   reflecting, non-characteristic injection whose effective imposed field is neither the
   prescribed one nor the interior one. This is a numerics issue rather than a physics one, but
   it means the "prescribed" B_ω is not what the plasma sees even by the code's own logic.

4. **The geometry is inconsistent unless an assumption is stated.** The coil sits at 3.6 cm,
   *outside* the flux conserver at 3.5 cm. A conducting shell at 795 kHz screens a transverse
   field completely (copper skin depth ~70 μm), so a coil outside it could drive nothing
   inside. The configuration is only sensible if the flux conserver is slotted — as TCS's
   segmented conservers are — so that it passes B⊥ while conserving axial flux. The code
   silently assumes that; it should say so.

5. **The sibling boundary conditions do not offer a corrected alternative.** `twoFluidRMFBC`
   uses 0.5·B_ω where the simplified one uses B_ω (an undocumented factor of two between
   two BCs that claim to impose the same drive); `twoFluidRMFHarmonicsBC` adds a harmonic
   series but keeps the same prescribed-field philosophy; `twoFluidOMFBC` is the oscillating
   variant with the flux term commented out; and `twoFluidRMFAntennaBC` — the one whose name
   suggests the right approach — imposes a purely azimuthal oscillating B and, in `setup()`,
   never reads `_baxial`, so its flux-conserver formula uses an uninitialised member.

### 3.5 The flux conserver — correct ✓ (with the assumption above)

B_z(a) = [π b² B_bias − Φ] / [π (b² − a²)] is exactly conservation of the total axial flux
inside an ideal shell at r = b, with the plasma's flux Φ measured each step and the remainder
assigned to the vacuum annulus. It reproduces the uniform field at t = 0 and the one-step lag is
harmless. It is what makes this an FRC problem rather than a driven column, and it is the part
of the boundary treatment that is right. It presumes the annulus a < r < b carries a uniform
B_z, which is exact only if the annulus is field-free of plasma currents — true here because
the annulus is not in the domain.

### 3.6 Friction energetics — energy is not conserved ✗ (small on this run)

All four friction sources (`constantResistivity`, `anisotropicResistivity`, `braginskiiFriction`,
`MomentumXfer2D`) have the same energy terms. With R = −η n² e² (u_e − u_i) the friction on the
electrons:

```
electron energy source   s_e = -(R . u_i) - Q_delta
ion energy source        s_i =              + Q_delta
```

For total-energy variables, the correct sources are s_e = R·u_e + Q_e and s_i = −R·u_i + Q_i,
with Q_e + Q_i = −R·(u_e − u_i) = ηJ² ≥ 0 the frictional heat, which Braginskii assigns to the
electrons (Q_e = ηJ² − Q_Δ, Q_i = Q_Δ). Substituting, s_e = R·u_i − Q_Δ and s_i = −R·u_i + Q_Δ,
which sum to zero.

So the code has the right structure and the wrong sign on the electron work term, and omits
the ion work term; the total energy has a spurious source −R·u_i. Two things follow. First, an
earlier note in `docs/known-issues.md` said Ohmic heating was absent altogether; that was
wrong for a total-energy formulation — with ions at rest the code and the correct form agree
(both −Q_Δ), because the friction thermalises the electron drift *within* the electron fluid
and the ηJ² heating is exactly the kinetic energy the momentum source removes. Second, the
error is proportional to u_i, so it is negligible on this 2 μs run (ion gyroperiod 11 μs) and
grows as the ions spin up — which the literature says they do, quickly. Any run long enough to
reach the torque–friction balance carries it.

Separately, `braginskiiFriction` and `MomentumXfer2D` computed Q_Δ with 3·(m_i/m_e) where the
Braginskii coefficient is 3·(m_e/m_i) — an overstatement by (m_i/m_e)², about 3.4e6 for hydrogen.
`constantResistivity`, the source this deck actually uses, always had the right mass scaling.
Both were corrected in Phase 0 and are covered by `test/cxx/test_friction_sources.cc`.

### 3.7 Resistivity — constant η, no feedback (limitation)

η is a deck constant. Hoffman 2006 notes that "higher temperatures have been noted to reduce the
effective plasma resistivity", and Milroy 1999 that anisotropic resistivity changes penetration
qualitatively. With a fixed η, γ and λ are fixed for the whole run, so the simulation cannot
show the self-consistent penetration the experiments show, and the anomalous absorption ∝ B_ω²
that Guo 2007 finds dominant at low current is not represented at all. Reasonable for a first
study; a limitation for comparison with data.

### 3.8 Ions and neutrals (limitation)

No neutral species, so no ion-neutral drag. In the literature that drag (with end-shorting) is
what limits ion spin-up and preserves the current drive; without it a two-fluid model will spin
the ions up on the electron-ion collision timescale and lose the drive. Irrelevant at 2 μs;
decisive for sustainment runs. The multi-fluid machinery for a neutral fluid exists in the code.

### 3.9 Dimensionality (limitation)

r-θ is the right plane for penetration and drive — the classical theory is an infinite cylinder
— but cannot represent finite antenna length, end-shorting, axial flux compression, or thrust.
For the xenon deck that means the r-θ model addresses whether the RMF penetrates and drives
current, not whether the thruster produces thrust.

### 3.10 Reduced speed of light — acceptable ✓

c/100 with ε₀ rescaled keeps c² = 1/(μ₀ε₀). Light transit across the column is 10 ns against a
1258 ns RMF period, so the displacement current remains negligible for the RMF; ω_pe drops to
5.6 × 10⁹ s⁻¹, still 5× ω_ce and resolved (ω_pe Δt = 0.15); the Debye length grows to 0.4 mm,
comparable to the finest cells. Reduced c is standard in five-moment work; a two-fluid study in
the search results ran c from 3 to 12 v_Te and found "only a minor effect on the modeled
instabilities." The price is 73,000 steps for 2 μs.

### 3.11 Smaller items

- The initial condition computes a Gaussian density profile (`gaus`, `mval`, `value`, `rds`,
  `beta`) and never uses it; every component is uniform.
- `coilRadius` and `numberOfHarmonics` are read only by the harmonics BC; in this deck they are
  ignored.
- The heavy-ion deck sets `MI = 131*MP` — xenon — under a comment saying "hydrogen".
- No limiter is active (`#Limiter = [twoFluidLimiter]`); the run relies on Lax–Friedrichs
  dissipation, which is the most diffusive flux available.
- Both decks stop at ~1.6 periods; the literature's diagnostics (ζ, torque balance, density
  compression) need tens of periods.

---

## 4. Plan

Ordered so that each step has a verification and the early ones are cheap.

### Phase 0 — correctness fixes with tests — **DONE**

Implemented; see `test/cxx/` for the tests and the git history for the changes.
The items below are kept for the record, each annotated with what was done.

#### The items

1. **Friction energy terms** — done, all four sources. `s_e = +R·u_i − Q_Δ`,
   `s_i = −R·u_i + Q_Δ`. `constantResistivity` additionally had a three-dimensional friction
   with a two-dimensional work term, which broke conservation on its own; its work term is now
   3-D too. `test/cxx/test_friction_sources.cc` configures each real class from a real deck
   fragment and checks: the ion momentum source is exactly −R; the two energy sources sum to
   zero; each equals the closed form; the frictional heat −R·w is non-negative (which holds
   even with the anisotropic cross term, since that term does no work); and, for the isotropic
   law where it can be computed independently, that −R·w equals ηJ². 31 checks.
2. **Braginskii Q_Δ mass ratio** — done. Both now use `3*_me/_mi`. The cross-check
   confirmed `constantResistivity`'s `(3/m_i) n² η e²` reduces exactly to `3 (m_e/m_i) n ν_ei`
   under η = m_e ν_ei/(n e²), so that form was always right; and that the `nue` these two
   classes compute is an electron–ion collision frequency (within 1.4× of the NRL value at the
   deck's conditions), so the substitution is like-for-like. Tested by mass-ratio scaling
   (quadrupling m_i must quarter Q_Δ) and by a physical bound that an inverted ratio overshoots
   by ~10⁶.
3. **E_z at the boundary** — done. Now `E_z = x Ḃ_y − y Ḃ_x`, computed from the same
   B_x(t), B_y(t) the boundary imposes, so it stays consistent with whatever field convention
   is used. `test/cxx/test_rmf_boundary.cc` checks both Faraday identities by central
   differences against the shipped class (agreement to 1e-10), and that E_z now varies around
   the boundary circle where the old form was axisymmetric.

   **Also found and fixed while here:** `phase` was applied to only one of the two transverse
   components, so it changed the *polarisation* rather than the phase. At `PHASE = PI/2` — the
   value in the heavyIons deck — the applied field was linearly polarised along a fixed axis
   with a magnitude swinging between 0 and √2·B_ω: an oscillating field, not a rotating one,
   and a different experiment (Apollo has `twoFluidOMFBC` for that). The phase is now applied
   to both components. `PHASE = 0`, the hydrogen deck, is unaffected. The test checks the
   applied field has constant magnitude over a period and sweeps exactly one full turn, at
   both phases.
4. **`twoFluidRMFAntennaBC`** — done; `B_axial` is read in `setup()`, so its flux-conserver
   expression no longer runs on an uninitialised member.
5. **Documentation** — done, in the header of `aptwofluidsimplifiedrmfbc.h`: what the BC
   imposes, the three assumptions (edge-prescribed field, slotted flux conserver,
   over-specified hyperbolic boundary), and what each sibling BC actually applies. The
   "0.5·B_ω discrepancy" turned out not to be one: `twoFluidRMFBC` imposes an *azimuthal
   oscillating* field, not a uniform transverse rotating one, so the two are not alternative
   spellings of the same drive. The simplified BC is the only one of the group that imposes
   the classical RMF.
6. **Decks** — done. The heavy-ion deck says xenon. The dead Gaussian profile (`value` was
   computed and never referenced by any `exprList` entry) is removed from both decks, with a
   comment recording that the initial column is deliberately uniform — the literature's
   formation problem — and how to switch a profile back on.

### Phase 1 — put the drive at the antenna, not the plasma edge — **DONE, with a finding**

Implemented, and it changed the plan. The drive is now a current; the vacuum annulus turned
out not to be representable in this model, and the reason is quantitative and worth keeping.

**What was built.** `examples/unstructuredDG/multifluid/rmf_frc/antenna/`:

- **The antenna** — `maxwellRMFAntenna` (`src/hyperapps/maxwell/aprmfantennasrc.h`), a rotating
  m = 1 axial current sheet, which is what the saddle coils of an RMF machine are. It sources
  `dE_z/dt` with `-J_z/eps0`, the same normalisation the plasma's own current uses
  (`WxCurrentSrc`), so the field it produces is a solution of Maxwell's equations rather than a
  prescription. The deck states the field it wants and the class solves for the current,
  reporting what it arrived at: for the shipped deck, 1190 A peak per coil.
- **The flux conserver as an actual boundary** — `twoFluidConductingWallBC`, which already
  existed. It zeroes *both* tangential components of E, so by Faraday's law the axial flux it
  encloses is conserved exactly. That replaces the area-integral feedback of §3.5 with a local
  condition that needs no assumption about how the expelled flux distributes itself.
- **A mesh generator** — `scripts/mkdiscmesh.py`. The `rmf_frc` meshes ship with no `.geo`
  alongside them (unlike, say, `euler/backwardFacingStep`), so their domain could not be changed
  at all. It writes graded polar-triangulated discs
  in the gmsh 2.2 ASCII subset the reader accepts, with a node ring placed exactly on a chosen
  radius. Extending the domain from r = 0.03 to r = 0.05 costs 1.54× the cells (11 989 triangles
  against `optimizedCircle2.msh`'s 7 792) and *nothing* in time step: dt is set by the smallest
  cell, which is still in the plasma, and is in fact larger on the new mesh — 3.57e-11 s against
  the shipped mesh's 2.74e-11 s, both measured.
- **Verification with an exact answer** — `antenna/vacuum.pin` and
  `test/test_rmf_antenna_field.py`. With the plasma removed, a winding inside a perfectly
  conducting cylinder produces a uniform transverse field `B_rmf (1 - exp(-t/rise))` rotating at
  ω, in closed form. The test checks all three of magnitude, uniformity and rotation rate
  against it, at 5%, 6% and 10%; those tolerances are set by the O((ωb/c)²) ≈ 1% quasi-static
  correction and by the coarse mesh, not chosen to fit. Measured on the shipped deck, the error
  in |B| against the closed form falls from 9.8% a quarter of a rise time in to 0.1–1.9% over
  the seven frames from 1.5 rise times on, the field's variation across the interior from 3.8%
  to 0.6–3.1%, and the rotation rate over that window matches ω to 0.03%. The test looks only at
  t ≥ 1.5·RISE; the transient before that is the conducting cylinder's own undamped modes, which
  nothing absorbs when there is no plasma. This is the Phase 2 test
  applied to Phase 1's geometry, and it is what makes the antenna trustworthy in the deck where
  nothing is exact.
- **The antenna power** — `scripts/rmf_antenna_power.py`, which evaluates
  P = −∫E·J over the winding from a run's output and the deck's own antenna parameters. This
  is the second of the two numbers the literature is built on, and it exists only because the
  drive is a current: a prescribed edge field has no power. It is post-processing rather than an
  in-solver diagnostic because Apollo's area integrals are not observable — `WxNodalDG2dMethod`
  reduces them into `_AgregateAreaIntegral` and passes that to boundary conditions only, never
  logging them and never writing them out. The script cross-checks itself against the solver:
  it re-derives the antenna current amplitude and gets 16534.8 A/m, the figure the solver
  prints at setup. It reports a cycle average only when the frames span a whole RMF period, and
  says so otherwise — most of the instantaneous power is reactive.
- **Unit tests** — `test/cxx/test_rmf_antenna.cc`. The strongest of them inverts the class:
  sample `J_z` from `src()`, extract its m = 1 Fourier component shell by shell, superpose the
  shell fields analytically, and check the result equals the requested `B_rmf`. It agrees to
  2.6e-11, along a path that shares no arithmetic with `setup()` or `src()`.

**The finding: the plasma–vacuum interface cannot be carried, so there is no annulus.** The plan
this section replaces proposed the physically right structure — plasma disc r < a, vacuum annulus
a < r < b, coils, conserver at b — with the annulus as a very-low-density two-fluid region, on
the grounds that Apollo's density and pressure floors already support one. The floors are not
the problem. The *edge* is.

Apollo carries explicit charge separation, so the electron dynamics have two scales that must be
resolved. At the column's density, and with the deck's numbers:

| | value | against the mesh / step |
|---|---|---|
| Debye length λ_D = √(ε₀ k T_e / n e²) | 4.07e-4 m | 0.40–0.45 of a cell |
| plasma period 2π/ω_pe, ω_pe = 5.64e9 rad/s | 1.11e-9 s | 31 time steps |

Both are marginal, and the shipped deck gets away with them only because its state is uniform:
nothing excites either mode. A density ramp at the column edge excites exactly them, and the
electron fluid there runs away — the runs diverge at the ramp, at r ≈ 0.035, with electron
speeds reaching 7e7 m/s and the pressure going negative, with or without the antenna running.
A hard step fails on the first residual evaluation; a 1.5 mm tanh ramp lasts about 30 steps; a
3 mm ramp about 150; a 5 mm ramp survived at least 244 before it was stopped. Widening the ramp
buys steps, not stability, and it cannot be widened far — past a few millimetres the "column
edge" is the whole gap to the winding. There is no slope limiter to fall back on, because the
only one in the tree is broken (below).

Two things are worth stating precisely, because the obvious intuitions about them are backwards:

- **The reduced speed of light helps here, it does not hurt.** Rescaling ε₀ up by 10⁴ raises λ_D
  by 100: at the true speed of light λ_D would be 4.07e-6 m, about 1/200 of a cell, and the mesh
  would have to shrink by that factor to resolve it. The plasma *period* is a different matter
  and does not improve — dt is CFL-limited by c, so raising c by 100 shortens the step by 100
  too and the period stays at about 31 steps either way. It is the Debye length, whose yardstick
  is the fixed mesh, that the reduced c rescues.
- **A tenuous region is better resolved, not worse.** λ_D ∝ 1/√n, so the annulus would be the
  best-resolved part of the domain: at 1e-4 of the column density λ_D is 41 mm against 1.6 mm
  cells. The barrier is not the annulus; it is the gradient between it and the column, which
  necessarily passes through the column's own under-resolved density.

For completeness, the transparency requirement the annulus would have had to meet: the
collisionless skin depth c/ω_pe = 5.32e-4/√f must exceed the 3.5 mm from the column edge to the
inner edge of the winding, which needs f < 2.3e-2 — or f < 2.3e-4 for a comfortable tenfold
margin. That is easily satisfiable and was never the binding constraint. (c/ω_pe is independent
of the c rescaling: c²/ω_pe² = m_e/(μ₀ n e²), in which neither c nor ε₀ appears.)

**What was built instead.** The plasma fills the domain out to the conducting wall and the
antenna is embedded in it. That keeps everything Phase 1 was for — the drive is a current, the
field everywhere is an outcome, the plasma screens it, the antenna is loaded, and the absorbed
power is the ∫E·J over the winding — and gives up only the vacuum gap. Whether the field reaches
r < 0.03 is now a result of the run rather than a boundary condition, which is the point. It is
also not an unreasonable machine: RMF thrusters run with plasma filling the tube.

**What would restore the vacuum gap**, in increasing order of work: a working slope limiter; a
finer mesh in the column, which is what actually buys Debye resolution (λ_D there is fixed by the
plasma, so the cell has to come down to meet it — roughly a factor of three, at a factor of nine
in cells and three in time step); or a quasi-neutral or Hall formulation that does not carry
charge separation at all. The third is the real answer and is what the NIMROD work implies when
it says "the Hall term is a zeroth order effect".

**What a short run shows.** The transverse field peaks at the winding (2.9e-4 T at t = 3e-8 s)
and falls about two and a half decades over the centimetre on either side of it — the
non-penetrated skin state, which is what to expect at t ≪ one RMF period. The peak field at the winding is about 11× below the free-space drive the
same current would produce, which is the plasma loading the antenna: precisely the effect the
edge-driven boundary condition forbids by construction. Whether it penetrates over tens of
periods is Phase 3.

Two resolution caveats before anyone quotes a number from it. The collisionless skin depth,
c/ω_pe = 0.53 mm, is *below* the 0.9–1.5 mm cells, so the screening layer is only marginally
represented; and the resistive skin depth δ = 1.26 mm — the length RMF penetration theory is
written in terms of, through λ = r_s/δ — is resolved at about 1.4 cells. Both want refining
before a threshold is quoted.

**Bugs found and fixed on the way.** None of these are Phase 1 features; all were blocking it.

1. **Every RMF boundary condition dereferenced a null pointer.** `WxTuAliabadiLimiter::applyBc`
   passes no area integrals, and all five of `twoFluidSimplifiedRMFBC`, `twoFluidRMFBC`,
   `twoFluidRMFHarmonicsBC`, `twoFluidRMFAntennaBC` and `twoFluidOMFBC` read `AreaInts[15]`
   unguarded. Enabling the limiter on any RMF deck segfaulted on the first step.
2. **A heap buffer overflow in the DG geometry, present in every run.**
   `wxNodalDGgeometry2D` allocated `_ETETF` with `_Ktotal` rows while indexing it by local cell
   id up to `_kLocalInt - 1`. `_Ktotal` counts only the cells with three vertices *and* is an
   `MPI_Allreduce` **sum across ranks**, so on more than one rank it is the global element count
   rather than a local size at all. On the shipped `rmf_frc` mesh the stratum holds 7984 cells
   and `_Ktotal` is 7792: the last 192 rows were past the end of the pointer array.
3. **Ghost cells in the height-0 stratum.** `ApSolver::createMesh` calls
   `DMPlexConstructGhostCells`, which appends one ghost cell per boundary facet — the whole of
   the 7984-against-7792 discrepancy, and matching on every mesh checked (3929 triangles + 157
   boundary lines = 4086 for `vacuumDisc.msh`, 11 989 + 211 = 12 200 for `disc.msh`). They have
   no geometry, so the DG scheme tolerates them but the limiter asked them for their normals and
   stopped with "Edge length of 0 found in element 7792".
   `wxNodalDGgeometry2D::isRealCell()` now reports which cells are real.
4. **`maxwellRMFSrc` would leak stale source values across the whole plasma interior.**
   `WxHyperbolicSrc::compSource` keeps `_outValues` as a member, never clears it, and does
   `sfull[idx] += _outValues[i]`. `ApRMFSrc::src` wrote its three outputs only inside
   `if (r > _r0)`, so every point inside r₀ would receive the previous quadrature point's values
   in E_z, B_x and B_y. **Latent, not active:** the `heavyIons` deck is the only one carrying an
   `RMFsrc` block and it is commented out of that deck's `Sources` list, and
   `WxHyperbolicSrcSet::setup` instantiates only what `Sources` names — so the class is
   constructed by no shipped deck and nothing was ever corrupted. Fixed anyway, since Phase 1
   builds on this file; the class is also now documented as *not* an antenna (it has no current
   and cannot be loaded) and superseded by `maxwellRMFAntenna`.

**Still broken: the two-fluid slope limiter.** With those three fixes `tuAliabadiLimiter` runs
instead of crashing, but still produces NaN energy after a few tens of steps — on the shipped
`rmf_frc` deck as much as on the new one. The remaining fault is in its own numerics and was not
chased down. It is the only limiter that handles the eighteen-component two-fluid state
(`eulerLimiterHW` limits the five-component Euler state and is enabled in several shipped Euler
decks), so with it broken no multifluid case can be limited at all.

### Phase 2 — characteristic-consistent boundary injection (with Phase 1)

For the Maxwell subsystem, set only the incoming characteristic variables at a boundary face
from the prescribed field and take the outgoing ones from the interior state (the standard
Riemann-invariant / ghost-state construction: for a prescribed field F, use ghost = 2F − F_int
for the prescribed components so the face average is F, or inject through the eigenvectors of
the flux Jacobian, which `eigenSystem()` in `wxphmaxwelleqn.cc` already provides). Verification
that costs nothing: **a plasma-free run must reproduce the applied field exactly.** With n → 0
(or the plasma sources switched off) in the Phase 1 geometry, the interior field must equal the
uniform rotating field of the coils at every point and time — an exact solution, so the test
tolerance is discretisation error. Add it to `test/`.

### Phase 3 — validate against the literature (weeks, mostly compute)

0. **Which way the drive turns.** Confirm, over tens of RMF periods, that the driven B_z on
   axis moves *against* the bias field, for each of the two decks. The antenna deck rotates
   counter-clockwise and the edge-driven deck clockwise (see above); with the same +60 G bias
   they cannot both be forming an FRC. Everything below assumes this is settled.
1. **Penetration threshold.** At fixed γ and λ, scan B_ω and locate the transition between
   skin-depth-limited (interior B⊥ ≈ 0, current in a layer of thickness δ) and penetrated
   (near-synchronous rotation) states. Compare with the Hugrass–Grimm 1981 threshold and
   Milroy 1999's empirical expression. This is the test the whole subject rests on.
2. **Penetrated limit.** Check the driven field reversal against μ₀ n e ω a²/2 · ζ with ζ
   measured from the simulation's electron rotation.
3. **Skin-depth limit.** At low B_ω, check the current layer thickness against δ.
4. **Formation dynamics.** Reproduce Guo 2002's sequence — reversal, radial expansion, bias-flux
   compression, density rise, torque–friction balance — which needs tens of periods and the
   Phase 0 energy fix.
5. **Sousa 2016.** Re-run the cases behind the APS 2016 abstract with the corrected boundary and
   compare instability onset and reversal robustness.
6. **Thruster regime.** For the xenon deck, compare penetration against the Michigan
   measurements (R = 10 cm, δ = 1 cm, T_e ≈ 9 eV, currents to 2500 A), restricted to what r-θ
   can say.

### Phase 4 — physics the literature says matters next

- **η(T_e)** (Spitzer, with lnΛ) so that heating feeds back on penetration; the Phase 0 energy
  fix is a prerequisite or the feedback is fed wrong numbers.
- **Anomalous absorption ∝ B_ω²** as an optional resistivity enhancement, to compare with the
  Guo 2007 power balance.
- **Neutral fluid with ion-neutral drag**, using the existing multi-fluid machinery, for any run
  intended to reach sustainment.
- **r-z or 3-D** for finite antenna length, end-shorting and thrust; out of scope for this
  example but the reason an r-θ result on the xenon deck is not a thrust prediction.

### What can be compared today

With the antenna deck, the field the plasma sees is an outcome and the antenna current is an
input, so the two numbers the RMF literature is built on — the penetration threshold in B_ω and
the absorbed antenna power — are now *measurable* from a run rather than pinned by a boundary
condition. That is the door Phase 1 opens; walking through it is Phase 3.

Two caveats on any comparison made through that door. There is no vacuum gap between plasma and
winding (see Phase 1), so the antenna is more heavily loaded than a real one and the threshold
will not be the experiment's number. And the geometry is r-θ, so nothing about finite antenna
length, end-shorting or thrust follows.

The original edge-driven deck, `rmf_frc/frc2d.pin`, is unchanged and remains what it was: a
demonstration of what a two-fluid model with these parameters does when a uniform rotating field
is imposed at its edge. Both decks are worth keeping.

They differ in **two** things, not one, and the second was not intended. The antenna's current
goes as cos(θ − ωt − φ), whose field turns counter-clockwise; `twoFluidSimplifiedRMFBC` writes
B = (−B_t sin(ωt+φ), −B_t cos(ωt+φ)), which turns clockwise. In the synchronous limit the
electrons are dragged in the sense of the rotation, so J_θ = −e n u_θ takes the opposite sign to
it: a counter-clockwise RMF drives J_θ < 0, whose axial field opposes a +ẑ bias — field reversal
— while a clockwise one reinforces it. The antenna deck turns counter-clockwise for that reason.
If the argument holds, the shipped deck's pairing of a clockwise drive with a +60 G bias is the
wrong way round for formation. It is an analytic expectation only: over the fraction of a period
a short run covers, the driven ΔB_z is about 2e-6 of the bias and oscillates in sign, so it
settles nothing. **Confirming it is the first thing Phase 3 should do**, before any threshold
scan, because the sign decides whether the scan is measuring formation or its opposite.


---

## References (as located by search; publisher pages were not fetchable from this environment)

- Jones, I. R. & Hugrass, W. N., "Steady-state solutions for the penetration of a rotating magnetic field into a plasma column", *J. Plasma Phys.* 26, 441 (1981). https://www.cambridge.org/core/journals/journal-of-plasma-physics/article/abs/steadystate-solutions-for-the-penetration-of-a-rotating-magnetic-field-into-a-plasma-column/27A1A38DCB4E3819FF1DD153C6FBD641
- Hugrass, W. N. & Grimm, R. C., "A numerical study of the generation of an azimuthal current in a plasma cylinder using a transverse rotating magnetic field", *J. Plasma Phys.* 26, 455 (1981).
- Milroy, R. D., "A numerical study of rotating magnetic fields as a current drive for field reversed configurations", *Phys. Plasmas* 6, 2771 (1999). https://pubs.aip.org/aip/pop/article-abstract/6/7/2771/464976
- Milroy, R. D., "A magnetohydrodynamic model of rotating magnetic field current drive in a field-reversed configuration", *Phys. Plasmas* 7, 4135 (2000). https://pubs.aip.org/aip/pop/article-abstract/7/10/4135/264692
- Guo, H. Y., Hoffman, A. L., Brooks, R. D., et al., "Formation and steady-state maintenance of field reversed configuration using rotating magnetic field current drive", *Phys. Plasmas* 9, 185 (2002). https://pubs.aip.org/aip/pop/article/9/1/185/264956
- Hoffman, A. L., Guo, H. Y., Miller, K. E. & Milroy, R. D., "Principal physics of rotating magnetic-field current drive of field reversed configurations", *Phys. Plasmas* 13, 012507 (2006). https://pubs.aip.org/aip/pop/article-abstract/13/1/012507/896633
- Guo, H. Y., Hoffman, A. L. & Milroy, R. D., "Rotating magnetic field current drive of high-temperature field reversed configurations with high ζ scaling", *Phys. Plasmas* 14, 112502 (2007). https://pubs.aip.org/aip/pop/article-abstract/14/11/112502/936710
- Milroy, R. D., Kim, C. C. & Sovinec, C. R., NIMROD simulations of RMF-driven FRCs, *J. Comput. Phys.* 195, 355 (2004); IAEA FEC 2008 IC/P4-3. https://www-pub.iaea.org/mtcd/meetings/fec2008/ic_p4-3.pdf
- Hoffman, A. L. et al., "The TCS rotating magnetic field FRC current-drive experiment"; "The TCS upgrade: Design, construction, conditioning, and enhanced RMF FRC performance". https://www.researchgate.net/publication/224002692_The_TCS_rotating_magnetic_field_FRC_current-drive_experiment
- Belova, E. V. et al., "Hybrid magneto-hydrodynamic simulation of a driven FRC", *Phys. Plasmas* 21, 032507 (2014). https://pubs.aip.org/aip/pop/article-abstract/21/3/032507/1032494
- Sousa, E. M. et al., "Rotating Magnetic Field FRC Formation Studies using the Multi-Fluid Plasma Model", APS DPP 2016, J10.160. https://ui.adsabs.harvard.edu/abs/2016APS..DPPJ10160S/abstract
- Woods, J. M., Jorns, B. A. et al., "Equivalent Circuit Model for a Rotating Magnetic Field Thruster", AIAA 2021. https://pepl.engin.umich.edu/pdf/2021_AIAA_PE_Woods.pdf
- Gill, T. M., Woods, J. M., Sercel, C. L., Jorns, B. A., "Experimental investigation into efficiency loss in rotating magnetic field thrusters", *Plasma Sources Sci. Technol.* (2024). https://iopscience.iop.org/article/10.1088/1361-6595/ad107a
- Sercel, C. L. et al., "Inductive probe measurements in a rotating magnetic field thruster", *Plasma Sources Sci. Technol.* (2023). https://iopscience.iop.org/article/10.1088/1361-6595/acfd5a
- Hakim, A. & Shumlak, U., "Two-fluid physics and field-reversed configurations", *Phys. Plasmas* 14, 055911 (2007). https://www.aa.washington.edu/sites/aa/files/research/cpdlab/docs/Hakim_PoP2007.pdf
