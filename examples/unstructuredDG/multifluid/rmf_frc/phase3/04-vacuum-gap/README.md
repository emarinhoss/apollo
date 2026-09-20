# 04 — the vacuum-gap arm (known-issues §23)

## The question

Can the antenna sit in a tenuous annulus around the plasma column, the way
Phase 1 of the assessment first proposed, instead of being embedded in the
plasma as `03-formation` has it? Known-issues §22 removed the reason the
shipped deck gives for saying no (the Debye length), and measured that a
background of 1e-3 of the column density survives the column-edge expansion
that kills 1e-4. This deck is that annulus with the antenna switched on, asked
to run two RMF periods.

## The answer, in 18 minutes on one core

It dies at t = 3.53e-8 s, 993 steps, five nanoseconds after the antenna's 30 ns
rise ends — and not for §22's reason or §19's. An electron in a region too thin
to screen the winding picks up e·ΔA_z/m_e when the winding's vector potential
rises, 1.5e7 m/s here, which is five times the reduced speed of light. The
deck header and §23 have the numbers; the short form is that the drift a
background of density n needs in order to screen the antenna goes as n^-1/2,
and crosses c₀ between a tenth and a hundredth of the column density.

```bash
./../run_one.sh frc2d.pin          # dies at step 993; four frames
python3 ../../../../../../scripts/winding_anatomy.py frc2d.pin frc2d_*.vtu
```

The row to read in `winding_anatomy.py`'s table is `max(lam_e)/c0`: it is 2.6
by the first frame (10 ns), while the Gauss residual and c₀|φ|/|E⊥| — the two
numbers that grow in the embedded run — stay at 0.05 and 0.02.

## What it is for now

The negative result closes a design option: at LIGHT = 3e6 and B_ω = 50 G the
winding has to sit in plasma of at least about a tenth of the column density,
so "move the antenna out of the plasma" is not available as a cure for §19
without raising `LIGHT` (which costs ∝ c₀ in steps and shrinks λ_D ∝ 1/c₀)
or lowering `Bomega` (which is no longer the formation drive). The two
companions, `VAC_FRAC = 1.e-2` and `1.e-1`, are the same deck with that one
line changed; §23 records how far each got. The 1e-1 arm survives the rise and
then shows the other reason a gap deck cannot carry the formation question:
at β = 50 nothing confines the column, and it expands into the gap at the ion
sound speed, filling it in a third of a period.
