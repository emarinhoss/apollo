# 02 — the RMF penetration threshold (Phase 3 items 1 and 3)

## The question

The test the whole subject rests on: at fixed γ and λ, where does the plasma
stop screening the rotating field and start rotating with it? Nine runs sweep
B_ω from 5 to 150 G on the antenna deck, five RMF periods each.

The low end doubles as item 3 (the skin-depth limit): where the field is
screened, the driven current sits in a layer whose thickness should be the
resistive skin depth δ = 1.26 mm.

## Why the range is so wide

The numeric coefficient in the Hugrass–Grimm and Milroy threshold conditions
**could not be sourced** from the environment this was assessed in, and the
assessment deliberately quotes none. So the scan is designed to *bracket* the
transition rather than to confirm a predicted one — which measures the
coefficient instead of assuming it, and is the better experiment.

γ/λ = 2.63 for this deck, which is near the threshold rather than far above it,
so the transition should fall inside 5–150 G. If it does not, the scan has told
you something too.

## Running it

Nine independent runs, 35 h each on one core. This is a job array, not an MPI
job — see the note about MPI in the parent README.

```bash
sbatch --array=0-8 ../slurm_array.sh 02-threshold-scan     # cluster
./../run_one.sh Bomega_50G.pin                             # one point by hand
```

## Reading the result

Plot the cycle-averaged `penetration` against B_ω. The transition is where it
moves from near 0 to near 1. Then:

- **Below it**, `layer` should be measurable and close to δ = 1.26 mm. Beware:
  the mesh is essentially uniform at 1.035 mm per cell, so δ is 1.2 cells thick
  and the layer is barely resolved. `current_layer_thickness` refuses profiles
  it cannot honestly fit, so expect "no layer" rather than a wrong number, and
  read that as a statement about the mesh.
- **Above it**, ζ should approach 1 and B_z on axis should move against the
  bias. Compare its magnitude with μ₀ n e ω a²/2 · ζ = 452 G × ζ.
