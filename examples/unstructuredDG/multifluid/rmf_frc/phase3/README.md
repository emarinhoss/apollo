# Phase 3 runs: what to run, in what order, and what it costs

These folders hold the runs that Phase 3 of
[`docs/rmf-frc-model-assessment.md`](../../../../../docs/rmf-frc-model-assessment.md)
asks for. Each is self-contained — decks and mesh, nothing to resolve at run
time — so a folder can be copied to a cluster or a Colab session and started.

Read [Before you start](#before-you-start) first. One of the four runs is a gate
that decides whether the others are worth their compute, and it takes a minute.

---

## The runs

| folder | what it answers | runs | wall clock, 1 rank |
| --- | --- | --- | --- |
| [`00-divergence-gate`](00-divergence-gate) | Is the answer a property of the plasma or of the divergence-cleaning scheme? | 7 | 8 min, then **5 h** |
| [`01-c-sensitivity`](01-c-sensitivity) | Does the reduced speed of light being too low actually change anything? | 2 | 70 h + 160 h |
| [`02-threshold-scan`](02-threshold-scan) | Where is the RMF penetration threshold? (items 1 and 3) | 9 | 35 h each, 13 days total |
| [`03-formation`](03-formation) | Does the FRC formation sequence reproduce? (item 4) | 1 | 140 h |

`02` is nine independent runs and `01` is two, so on any machine with spare
cores the campaign is far shorter in wall clock than in core-hours — see the
note about MPI below, which is the reason to spend those cores on separate runs
rather than on one.

Those are one-rank figures on the machine the cost model was measured on
(0.464 s/step for 7792 triangles). Get your own with

```bash
python3 scripts/rmf_scan.py <deck.pin> --periods <N> --dry-run
```

and pass `--sec-per-step` once you have measured your node. **Treat the step
counts as exact and the hours as an order of magnitude**: the timestep model
reproduces the solver's `dt` to six significant figures, but the same binary on
the same mesh measured 0.464, 0.569 and 0.720 s/step under different load in one
afternoon.

---

## Before you start

### 1. Run the gate, in three stages. Only the third one answers the question.

`twoFluidSimplifiedRMFBC` writes the applied field into the ghost state at every
boundary face, and the resulting jump in the normal component of **B** feeds the
divergence-cleaning potential ψ faster than it can be carried away. On the
shipped edge-driven deck ψ passes 12× the electric field within 210 steps, and
two runs differing in nothing but `DIVB_SPEED` end up 92% apart in B_x. See
[`docs/known-issues.md`](../../../../../docs/known-issues.md) §13. That is why
**every long run here uses the antenna deck**.

**Stage 1, four runs, about eight minutes.** `edge-cleaning-{on,off}` and
`antenna-cleaning-{on,off}`. Cheap, and it separates a pathological deck from a
plausible one: the edge deck fails at 92%, the antenna deck passes at 2.4%.

**Stage 2, two runs, about 1.7 h each.** `antenna-long-cleaning-{on,off}`, 0.25
periods — 10.5 × RISE, so the drive actually reaches full amplitude, which the
stage-1 window never does (it spans 0.17 × RISE, entirely inside the ramp).
Measured result: **11.1%**, still growing at the end of the window, with ψ/E
peaking at 0.122 and then *falling* to 0.071. That last combination is worth
keeping: ψ bounded and improving while the sensitivity worsens, which is the
proof that the ψ check is necessary and not sufficient.

**But stages 1 and 2 both compare against `DIVB_SPEED = DIVE_SPEED = 0`, and
that is not a slower cleaning scheme — it is no cleaning at all.** γ and χ are
bare multiplicative factors on every term coupling φ and ψ to **E** and **B**
(`wxphmaxwelleqn.cc:490-497`), and χ also scales the charge source that
generates φ (`wxchargesrc.h:30`). At zero, all of them vanish and the potentials
are inert. So those two stages measure **how much divergence error the deck
carries** — real and worth knowing — but they cannot tell you whether the answer
depends on the cleaning *speed*, which is the question that decides the campaign.

**Stage 3, one run, about 1.7 h.** `antenna-long-cleaning-half` at
γ = χ = 0.5, compared against the stage-2 `antenna-long-cleaning-on` you already
have. Both speeds are non-zero, and because the Maxwell wave speed is
`dmax(χc₀, γc₀, c₀)` both give exactly c₀ — so the timestep is **identical** and
nothing is confounded. (2.0 would double the wave speed, halve dt, and confound
the comparison with a resolution change.) **This is the run that licenses, or
refuses, `01`, `02` and `03`.**

```bash
./run_one.sh 00-divergence-gate/antenna-long-cleaning-half.pin
python3 ../../../../../scripts/rmf_cleaning_check.py \
    results/00-divergence-gate/antenna-long-cleaning-on \
    results/00-divergence-gate/antenna-long-cleaning-half
```

The script now reads both decks' cleaning speeds and says so, and refuses to
call a comparison against zero a sensitivity test.

### 1a. The gate has been run. Here is what it decided.

Stage 3 (γ=χ=1.0 against 0.5, 0.25 periods, identical timestep) reports
**PASS at 1.1%**. Three things follow, and the third is the one that governs
the campaign.

**The cleaning scheme is doing its job.** The late-window growth rate is
0.016 of rms B_x per RMF period against 0.189 for the on/off comparison — the
cleaning *speed* matters about **12× less** than whether you clean at all. That
is the reassurance stage 3 existed to provide, and stages 1 and 2 could not.

**But it does not saturate.** From frame 9 onward the difference grows linearly
at a steady 0.0002 per frame. Extrapolating that late trend:

| | B_x (drives `penetration`) | B_z (drives field reversal) |
| --- | --- | --- |
| 0.25 periods (measured) | **1.1%** | 0.025% of bias |
| 5 periods — `02-threshold-scan` | ~8.7% | ~0.45% |
| 10 periods — `01-c-sensitivity` | ~17% | ~0.92% |
| 20 periods — `03-formation` | ~33% | ~1.8% |

B_x crosses the 5% tolerance at about **2.7 periods**. Every long folder here is
longer than that.

**The two observables are not equally exposed, and that is the useful part.**
B_z stays under 2% out to twenty periods, so the items that read the axial field
— field reversal, the penetrated limit, formation — are within tolerance for the
runs as specified. The penetration ratio is built from B⊥, so `02` carries a
systematic of order 8.7%. That is above the tolerance but it is **not** fatal to
what `02` is for: the scan brackets a transition across which `penetration`
moves by two orders of magnitude (0.45 to 0.003 within this very run), so an 8.7%
systematic on B_x does not blur the bracket. Quote it as a systematic rather than
treating the scan as blocked.

**Treat the table as a projection, not a measurement.** It extrapolates a
0.25-period window out to 20 — a factor of 80 — and an earlier extrapolation of
exactly this kind, from the 0.004-period stage-1 window, got the shape wrong even
though it got the order of magnitude right. If cluster time allows, a 1-period
pair (about 7 h per run) tests the linear trend at 4× the window and costs a day
against the thirteen `02` asks for.

### 2. Do not use MPI for a run whose output you intend to analyse.

The conclusion is unchanged but the reason has changed, and the new one is worse.

This section used to say a frame held only about 1/N of the cells on N ranks.
That was a misdiagnosis: PETSc writes one `<Piece>` per rank and every piece is
in the file, but `test/vtu.py` flattened them so the last rank's copy of each
array won. The counts that looked like a truncated writer — 7792, 3896, 1961 on
1, 2 and 4 ranks — were the size of the last piece. The reader is fixed and a
two-rank frame now reads back 4304 cells against the one-rank file's 4304.

Fixing it made the runs comparable for the first time, and **they do not agree**.
On this gate deck the penetration ratio is 0.0000 at every frame on one rank and
climbs to 0.71 on two, over the same 210 steps; B_z on axis moves in the fifth
digit. The same divergence appears on the Maxwell circular pulse — a different
module and a different boundary condition — from an identical frame 0. See
[`known-issues.md`](../../../../../docs/known-issues.md) §15.

> **Run one rank per run, and get throughput from running many runs at once.**
> That still suits this campaign: `02-threshold-scan` is nine independent points,
> which is a SLURM job array, not an MPI job.

The solver does scale — 2.06× on 2 ranks, 3.88× on 4 — so there is real time to
be had here once §15 is understood. It is not available yet.

`scripts/rmf_diagnostics.py` prints a banner when a frame does not cover the
domain. With the reader fixed it should not fire on a multi-rank run; if it does,
that is a genuine gap in the output.

### 3. There is no checkpoint or restart.

A run that is interrupted is a run that is lost
([`known-issues.md`](../../../../../docs/known-issues.md) §6). Size your job's
wall-clock limit against the estimate **plus a wide margin**, and prefer several
shorter runs to one long one where the physics allows it.

---

## Running a folder

Every folder works the same way. From inside it:

```bash
# 1. Expand the .pin macro block into the .inp the solver reads.
PYTHONPATH=$APOLLO/scripts python3 $APOLLO/scripts/wxinpparse.py -i <name>.pin

# 2. Run. One rank; see "Do not use MPI" above.
$APOLLO/src/build-opt/apollo -i <name>.inp

# 3. Analyse. The deck is read for its own parameters, so this cannot drift out
#    of step with the run it is describing.
python3 $APOLLO/scripts/rmf_diagnostics.py <name>.pin <name>_*.vtu
```

with `APOLLO` set to your checkout. Each run writes its frames into the working
directory, so **give each run its own directory** or they will overwrite each
other — the frame names depend only on the deck's `Simulation` name.

`run_one.sh` in this folder does the three steps with that isolation, and
`slurm_array.sh` is a job-array template for the scan.

### Output volume

One frame is 3.5 MB on the edge deck's mesh and 5.4 MB on the antenna deck's.
The folders ask for 12 frames per RMF period, which is what
`rmf_diagnostics.cycle_average` needs to call an average a cycle average:

| folder | frames | disk |
| --- | --- | --- |
| `00-divergence-gate` | 21 | 90 MB |
| `01-c-sensitivity` | 242 | 1.3 GB |
| `02-threshold-scan` | 549 | 3.0 GB |
| `03-formation` | 241 | 1.3 GB |

Reduce `OUT` in a deck if that is a problem, but not below 8 frames per period
or the cycle averages stop being cycle averages and the diagnostics will say so.

---

## What to do with the results

The diagnostics print a per-frame table and a cycle average. The numbers to take
to the literature are the **cycle averages**, and only from frames spanning at
least one whole period — everything before that is the switch-on transient, and
the tool says so in a banner.

| quantity | what it is | compare against |
| --- | --- | --- |
| `B_z axis` | driven axial field | the bias: opposing it is field reversal |
| `zeta` | u_eθ/(ωr), the electron rotation parameter | 1 is synchronous |
| `penetration` | interior \|B⊥\| over edge \|B⊥\| | 1 penetrated, ~0 screened |
| `layer` | e-folding length of J_θ inward from the edge | δ = 1.26 mm |

Two things the assessment asks you to keep in view while reading them:

- **The penetration threshold's numeric coefficient could not be sourced.** The
  scan is designed to *bracket* the transition rather than confirm a predicted
  one, which measures the coefficient instead of assuming it. γ/λ = 2.63 for
  this deck, so it sits near the threshold, not far above it.
- **The penetrated state is the one the model represents least well.** In it the
  electron fluid's fastest characteristic exceeds the model's reduced speed of
  light by 3.9%. That is what `01-c-sensitivity` exists to quantify, and it is
  the reason to run that one first among the long ones.

---

## Regenerating these folders

The decks are generated, not hand-edited:

```bash
python3 scripts/rmf_make_phase3_runs.py           # rewrite them
python3 scripts/rmf_make_phase3_runs.py --check   # verify they are current
```

Generate rather than edit by hand, and if you do edit one by hand, know the trap:
Apollo's parser reads `3000000` as an integer and `get<REAL>` on it aborts the
run at setup with `std::bad_cast` naming no key — and the mirror is just as
fatal, because `Output_files` is read with `get<int>` and `OUT = 6.0` dies the
same way. `known-issues.md` §12 has both. The generator's `deck_number()` and
`deck_integer()` exist for exactly this.
