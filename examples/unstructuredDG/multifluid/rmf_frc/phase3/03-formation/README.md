# 03 — formation dynamics (Phase 3 item 4)

## The question

Guo 2002's sequence: field reversal, radial expansion, bias-flux compression,
density rise, and a torque–friction balance that settles. Twenty RMF periods on
the antenna deck, which is 140 h on one core and the longest single run in the
campaign.

## Run it last

Not because it is least interesting — it is the one that would show an FRC
actually forming — but because it is the one most exposed to everything else:

- It is the only run here that depends on **β**, so the `--speedup` scaling that
  makes the others cheaper cannot be used on it.
- It spends most of its length in the **penetrated** state, which is the state
  §3.10 says the model represents least well. Run `01-c-sensitivity` first; if
  those two disagree, this run inherits the disagreement over twenty periods
  rather than ten.
- At 140 h it is the longest run in the campaign, so size the job's wall-clock
  limit generously. It now checkpoints every frame, so an interruption costs one
  output interval (~7 minutes) rather than everything — resume with
  `APOLLO_RESUME=1 ./run_one.sh 03-formation/formation.pin`, or with
  `apollo -i formation.inp -r formation.checkpoint` from inside the results
  directory.

```bash
./../run_one.sh formation.pin
```

## Reading the result

The sequence is a story over time rather than a single number, so read the
per-frame table, not just the cycle average at the end:

| watch | for |
| --- | --- |
| `B_z axis` | crossing zero and going negative: that is the reversal |
| `zeta` | rising toward 1 and then holding: the torque–friction balance |
| `penetration` | staying near 1 once the field is in |
| `layer` | "no layer" once penetrated, which is the correct answer there |

241 frames at 12 per period, which is what makes a cycle average meaningful and
also 1.3 GB of output.
