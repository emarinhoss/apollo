# 01 — does the reduced speed of light matter?

**The most informative pair of runs in Phase 3. Buy this one first.**

## The question

The decks run with a reduced speed of light, c = 3.0e6 m/s. §3.10 of the
assessment establishes that this is not comfortably above the electron speeds
but level with them: the electron sound speed is 2.966e6, and in a *penetrated*
state — the one items 1, 2 and 4 exist to measure — the fastest electron
characteristic |u| + c_se reaches 3.116e6, which is 3.9% **above** the model's
own speed of light.

Whether that damages the answer or merely violates a principle is not something
argument settles. These two runs settle it.

| deck | c | c/v_Te |
| --- | --- | --- |
| `c-as-shipped.pin` | 3.0e6 m/s | 1.31 |
| `c-3vTe.pin` | 6.892e6 m/s | 3.00 — the bottom of the range §3.10's citation covers |

Ten RMF periods each, on the antenna deck, identical in everything else.

## Cost

70 h and 160 h on one core. The second is longer because dt scales as 1/c: a
defensible speed of light costs 2.3× the compute, which is itself part of the
finding.

```bash
./../run_one.sh c-as-shipped.pin
./../run_one.sh c-3vTe.pin
```

## Reading the result

Compare the **cycle averages** of ζ, penetration and B_z on axis between the two.

- **They agree** — the shipped setting is defensible after all, `02` and `03`
  can be run on it, and the rest of §3.10 is a caveat rather than a blocker.
- **They disagree** — the size of the disagreement is the measurement of how
  much this matters, and every other Phase 3 number carries it as a systematic.

Either way the answer is worth more than the runs it took, because it applies to
every subsequent comparison. Note what it does *not* settle: raising c to
3 v_Te leaves the Debye length at 0.17 of a typical cell, so a disagreement
could be either effect, and separating them needs a finer mesh.
