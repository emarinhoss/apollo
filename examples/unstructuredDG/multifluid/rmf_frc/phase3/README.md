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
| [`00-divergence-gate`](00-divergence-gate) | Is the answer a property of the plasma or of the divergence-cleaning scheme? | 3 | **6 min total** |
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

### 1. Run the gate. It takes a minute and it can save you a week.

`twoFluidSimplifiedRMFBC` writes the applied field into the ghost state at every
boundary face, and the resulting jump in the normal component of **B** feeds the
divergence-cleaning potential ψ faster than it can be carried away. On the
shipped edge-driven deck ψ passes 12× the electric field within 210 steps, and
two runs differing in nothing but `DIVB_SPEED` end up 92% apart in B_x — over
0.005 of one RMF period. See [`docs/known-issues.md`](../../../../../docs/known-issues.md)
§13.

That is why **every long run here uses the antenna deck**, whose ψ/E stays near
0.05. `00-divergence-gate` demonstrates both, so you can confirm it on your build
rather than taking it from this file.

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
