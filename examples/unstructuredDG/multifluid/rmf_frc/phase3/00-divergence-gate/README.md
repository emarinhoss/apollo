# 00 — the divergence-cleaning gate

**Run this first. It takes about six minutes and it decides whether the other
three folders are worth their compute.**

## The question

Apollo's perfectly-hyperbolic Maxwell system carries a potential ψ whose job is
to transport away any divergence in **B**. `twoFluidSimplifiedRMFBC` writes the
applied transverse field into the ghost state at every boundary face, against
whatever interior field the run has produced, and the jump that creates feeds ψ
faster than it can be carried. ψ then enters the induction equation in the slot
**E** occupies, so it is not a bookkeeping quantity sitting off to one side.

The test is not "is ψ large" but "does the answer depend on the cleaning speed".
Two runs identical except for `DIVB_SPEED` must agree. If they do not, whatever
they produce is partly a property of the scheme.

## The runs

| deck | what it is |
| --- | --- |
| `edge-cleaning-on.pin` | the shipped edge-driven deck, `DIVB_SPEED = DIVE_SPEED = 1.0` |
| `edge-cleaning-off.pin` | the same deck with cleaning off |
| `antenna-cleaning-on.pin` | the antenna deck, for contrast |

All three are about 200 steps — 0.004 of an RMF period.

```bash
./../run_one.sh edge-cleaning-on.pin
./../run_one.sh edge-cleaning-off.pin
./../run_one.sh antenna-cleaning-on.pin

python3 $APOLLO/scripts/rmf_cleaning_check.py \
    ../results/00-divergence-gate/edge-cleaning-on \
    ../results/00-divergence-gate/edge-cleaning-off
```

## What to expect

On the machine this was written on, the edge-driven pair gives ψ/E climbing past
10 within 210 steps and B_x differing by over 90% of its own rms between the two
cleaning speeds — the check exits non-zero and says so. The antenna deck's ψ/E
stays near 0.05.

If you see that, it is reproducing, and the reason the long runs in `01`, `02`
and `03` all use the antenna deck. If you see something else, that is worth
knowing before spending a fortnight of compute, and it means this file is wrong
about your build.

`docs/known-issues.md` §13 has the measured numbers and the mechanism.
