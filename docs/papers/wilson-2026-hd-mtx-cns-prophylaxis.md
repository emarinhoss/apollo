# High-dose methotrexate as CNS prophylaxis in ultra high-risk LBCL: an appraisal

> **Scope.** This appraises a clinical oncology paper. It has nothing to do with
> Apollo and sits in this repository at the owner's request. It is a reading of
> the evidence, not advice about any patient.

**Paper.** Wilson MR, Lewis KL, Kirkwood AA, et al. *High-Dose Methotrexate as
CNS Prophylaxis in Ultra High-Risk Large B-Cell Lymphoma: An International
Multicenter Analysis.* J Clin Oncol 44:2290–2302, 2026.
[doi:10.1200/JCO-26-00947](https://doi.org/10.1200/JCO-26-00947). Accepted
7 May 2026, published 11 June 2026, CC BY-NC-ND 4.0.

**What was read.** The 16-page article: text, Tables 1–3, Figs 1–3, references
and disclosures. **Not read:** the Data Supplement (Tables S1–S7, Fig S1) and
the accompanying editorial (p 2263), which were not in the copy provided.
Anything that rests only on them is marked as unverified.

**Where numbers come from.** A number cited to a page, table or figure is
transcribed from the paper. A number marked *(check X)* is derived here and is
printed by section X of `python3 docs/papers/wilson2026_checks.py`, which uses
the standard library only and names the source of every input. Nothing is
quoted from memory.

**Abbreviations.** LBCL, large B-cell lymphoma. HD-MTX, high-dose
methotrexate. UHR, ultra high-risk. CNS-IPI, CNS international prognostic
index. IT, intrathecal. EN, extranodal. HR, hazard ratio. NNT, number needed
to treat. Throughout, **HR is the paper's: the hazard without HD-MTX relative
to the hazard with it, so HR > 1 favours HD-MTX.**

## Bottom line

- **The null finding holds up.** Across 1,923 UHR patients, HD-MTX shows no
  reduction in CNS relapse. That holds from diagnosis (3-year rate 9.3% v 8.1%,
  adjusted HR 1.13, 0.82–1.57) and from a 6-month landmark (6.7% v 6.6%, HR
  0.95, 0.62–1.44), for any and for isolated relapse, with no dose gradient.
  The measured imbalances left unadjusted mostly favour HD-MTX, so correcting
  them would not create a benefit.
- **It rules out a large benefit, not a small one.** The pooled analyses had
  80% power only for relative reductions of 34–46% *(check C)*. In absolute
  terms the adjusted intervals still allow up to 2.0–3.3 fewer CNS relapses
  per 100 patients at 3 years, an NNT of 31–50 *(check D)*. Whether that is
  "meaningful" depends on HD-MTX's harms, which the paper does not report for
  this cohort.
- **The headline comparison is the weaker one.** The abstract leads with the
  from-diagnosis HR. Survival shows that comparison is still confounded after
  adjustment and after matching: a 13–15 point overall-survival gap, of which
  CNS prophylaxis could explain less than 3 points *(check I)*. Only the matched
  landmark cohort balances survival (82.7% v 82.6%). Its CNS-relapse curves
  are also null, so the conclusion stands, but its hazard ratio appears only
  in the supplement.
- **It is mostly a between-study comparison, and a re-analysis.** 85% of the
  treated arm comes from one study and all of the untreated arm from the
  other *(check F)*. No adjustment for study, site or year of diagnosis is
  reported, though median diagnosis dates differ by four years. These are also
  the same patients as the "previous studies" the paper says it supports.
- **The subgroup results cannot show "no benefit".** Each subgroup has 13–52
  landmark events, enough to detect only 54–79% reductions *(check C)*, and the
  main text reports no interaction tests.
- **The published version contains ten internal inconsistencies.** The two
  that matter most are quantitative. Fig 2's absolute differences are
  inflated by a factor of the hazard ratio, which overstates every benefit
  bound. Two rows of Table 3 appear transposed. None of them reverses the
  conclusion, and correcting Fig 2 strengthens it.

## The question

CNS relapse after first-line treatment of LBCL is uncommon, about 2–6%, and
usually fatal, with median survival under 6 months (pp 2290–2291).
Intravenous HD-MTX has been used for about 20 years to prevent it, on the
strength of retrospective series only. Two large analyses had already
questioned the benefit (refs 19 and 20, summarised below). That left the
highest-risk subgroups, for which guidelines still recommend considering
HD-MTX, as the open question. This paper pools the two datasets to address
that group.

UHR means at least one of: CNS-IPI 5–6; testicular, renal/adrenal or breast
involvement; three or more extranodal sites (abstract).

## What was done

| | |
| --- | --- |
| Sources | **Study 1**, Wilson et al. 2022 (ref 19): 1,384 patients from 47 sites who all received HD-MTX with R-CHOP-like therapy. It had no comparator arm. **Study 2**, Lewis et al. 2023 (ref 20): 2,418 patients from 23 sites with prespecified high-risk features, 425 of them given HD-MTX. Its features included double/triple-hit lymphoma and stage I–II primary breast or testicular LBCL. |
| Amalgamation (p 2292) | All of study 1, plus study 2 where data-transfer agreements allowed. Study 2 contributed both arms from 9 sites, and untreated patients only from 6 further sites that had also contributed to study 1. The result was 3,491 patients: 1,866 without HD-MTX, 1,625 with. |
| Analysed | 1,923 UHR patients, 1,051 without HD-MTX and 872 with (Table 1). The landmark cohort had 1,555 (782 v 773). |
| HD-MTX | Study 1: any IV MTX dose intended to cross the blood–brain barrier, at least one cycle. Study 2: at least one cycle of IV MTX at any dose. Median 2 cycles, median cumulative dose 6 g/m² (p 2297). |
| Outcome | First CNS relapse, either isolated or synchronous (CNS and systemic relapse within 30 days). Death, systemic relapse and stable or progressive disease at the end of chemotherapy are competing events. Measured from diagnosis, and from a 6-month landmark limited to patients alive, progression-free and in complete or partial response. |
| Analysis | Fine–Gray competing-risks regression is primary, with cause-specific Cox as a sensitivity analysis. Adjustment covers age, sex, the CNS-IPI factors, and breast and testicular involvement, plus IT prophylaxis in the landmark cohort. Propensity-score matching used a caliper of 0.25 SD, giving 618 pairs from diagnosis and 501 in the landmark cohort. |

## Results

| 3-year CNS relapse, % (95% CI) | without HD-MTX | with HD-MTX | adjusted HR (95% CI) |
| --- | --- | --- | --- |
| from diagnosis, any | 9.3 (7.7–11.3), 111 events | 8.1 (6.4–10.3), 71 events | 1.13 (0.82–1.57) |
| from diagnosis, isolated | 5.9 (4.6–7.6) | 5.7 (4.3–7.6) | 1.03 (0.69–1.53) |
| landmark, any | 6.7 (5.1–8.7), 65 events | 6.6 (4.9–8.9), 50 events | 0.95 (0.62–1.44) |
| landmark, isolated | 4.6 (3.3–6.4) | 4.8 (3.4–6.8) | supplement only |

Sources: abstract, Table 2, p 2297.

- **Dose.** No gradient by cumulative dose (under 6 v at least 6 g/m²) or by
  number of doses of at least 3 g/m². Trend P = .86 and .63 for any relapse,
  .95 and .51 for isolated relapse (Table 3).
- **Subgroups.** No significant difference in any of CNS-IPI 5–6, CNS-IPI 5,
  CNS-IPI 6, three or more extranodal sites, testicular, breast, or
  renal/adrenal (Fig 2). The largest numerical gap is isolated relapse in
  CNS-IPI 5: 6.2% v 2.5%, or 16/262 v 4/189. The authors judge it probably
  spurious because CNS-IPI 6 shows nothing similar.
- **Survival.** 3-year overall survival without and with HD-MTX was 64.6% v
  79.2% from diagnosis (adjusted HR 1.65), still 65.7% v 79.1% after matching
  (HR 1.69), and 82.7% v 82.6% in the matched landmark cohort (HR 0.99)
  (Table 2).
- **After a CNS relapse.** The whole amalgamated cohort had 264 CNS relapses
  and 221 deaths, with median survival 3.8 months. Survival did not differ by
  prior HD-MTX (P = .95), by CNS compartment, or between isolated and
  synchronous relapse (Fig 3).

## Appraisal

### 1. The matched landmark comparison is the fair test, and it is null

HD-MTX can change survival only through the CNS relapses it prevents. From
diagnosis, the largest reduction in CNS relapse the adjusted interval allows
is 3.3 points. With 84% of CNS relapses fatal, that could explain about
2.7 points of the 3-year survival gap. The observed gap is 14.6 points, and
13.4 after matching *(check I)*. So the from-diagnosis comparison is still
dominated by selection after both adjustment and matching.

The source of that selection is visible. Patients who died, progressed or
failed to respond early could not go on to receive end-of-treatment HD-MTX.
A quarter of the untreated arm (25.6%) leaves before the landmark, against
11.4% of the treated arm. Those patients take 41% and 30% of each arm's CNS
events with them *(check F)*. The landmark removes this, and matching then
balances survival almost exactly.

That makes the matched landmark cohort the fair test. Its CNS-relapse curves
(Fig 1C–D) show no benefit, with the HD-MTX curve running slightly above after
about two years. Its hazard ratio, however, is only in the supplement. The
abstract instead leads with the from-diagnosis HR of 1.13, the estimate the
survival data show to be least trustworthy. The conclusion survives because
both estimates are null, but the emphasis is the wrong way round.

### 2. The imbalances left over mostly favour HD-MTX

| imbalance | evidence in the paper | pushes the HD-MTX estimate toward | handled by |
| --- | --- | --- | --- |
| immortal time and early drop-out | 25.6% v 11.4% leave before 6 months *(check F)* | benefit | the landmark |
| fitter patients, more complete chemotherapy | ECOG ≥ 2 in 51.0% v 30.1%; six cycles in 75.5% v 91.8% (Tables 1–2) | benefit, through fewer synchronous relapses | age and ECOG adjusted; cycles not |
| more IT prophylaxis | 32.4% v 50.3% (Table 2) | benefit, if IT works | adjustment in the landmark only |
| more baseline CNS staging | CSF analysis in 24.7% v 42.1% of each arm, imaging in 6.2% v 14.0% *(check H)* | benefit, since occult CNS disease is excluded more often | not adjusted |
| later era | median diagnosis Nov 2012 v Dec 2016 (p 2292) | uncertain; PET upstaging, which the authors invoke, would favour the later treated arm | not adjusted |
| confounding by indication | treated arm has more ≥ 3 extranodal sites, testicular, renal/adrenal and stage III–IV disease, but less CNS-IPI 5–6 (31.2% v 47.5% of those scored; Table 1, *check H*) | mixed; selection on unrecorded CNS risk would hide a benefit | measured factors adjusted |
| competing risks | treated patients live longer, so they spend longer at risk | harm, in the Fine–Gray model | cause-specific Cox, reported unchanged |
| diluted treatment | 13.6% had one cycle, 29.8% under 6 g/m² in total; infusion times unknown (Table 2, p 2299) | the null | dose analysis, no gradient |
| cross-study structure | see section 3 | unknown | not reported |

Most of what was left unadjusted favours HD-MTX, and the estimate is still
null. That strengthens the null result. What could hide a real benefit is
unrecorded confounding by indication: clinicians reserving HD-MTX for patients
they judged at highest risk on grounds not in the data. To hide a true 30%
reduction, such a confounder would need a risk ratio of about 2.4 with both
receipt of HD-MTX and CNS relapse in the landmark analysis, or 1.8 from
diagnosis. To hide a 20% reduction it would need 2.0 and 1.5 *(check E,
E-values)*. That is large but not implausible for a treatment given on
clinical judgement.

### 3. It compares two studies more than two treatments

85.2% of the treated patients in the amalgamated cohort come from study 1,
and every untreated patient comes from study 2 *(check F)*. The studies
differ in several ways:

- **Eligibility.** Study 2 selected on high-risk features. Among patients
  tested in the full cohort, double/triple-hit lymphoma is 21.0% of the
  untreated and 6.7% of the treated (p 2297).
- **Organ criteria.** Study 2's criteria were stage I–II *primary* breast or
  testicular LBCL, while study 1 took any HD-MTX recipient (p 2292). So
  "testicular" and "breast" do not select the same patients in the two arms.
- **Era and follow-up.** Median diagnosis dates are four years apart, and
  median follow-up is 76.0 v 40.6 months (p 2292).
- **Data completeness.** Baseline CNS assessment is unknown for 25.6% v 3.1%,
  and double/triple-hit status for 77.7% v 30.7% *(check H)*.

Any difference in how the two studies ascertained CNS relapse is therefore
inseparable from treatment. The main text reports no adjustment or
stratification by study, site or year, and does not say how missing
covariates were handled. It also does not report the within-study comparison
from the nine study-2 sites that contributed both arms.

### 4. What the intervals exclude

| analysis | events | smallest relative reduction detectable with 80% power |
| --- | --- | --- |
| from diagnosis, any | 182 | 34% |
| from diagnosis, isolated | 125 | 40% |
| landmark, any | 115 | 41% |
| landmark, isolated | 84 | 46% |
| landmark subgroups, any | 13–52 | 54–79% |

*(check C, Schoenfeld's formula at two-sided 0.05)*

| adjusted analysis | HR (95% CI) | 3-year difference, points (95% CI) | NNT at the most favourable bound |
| --- | --- | --- | --- |
| from diagnosis, any | 1.13 (0.82–1.57) | +1.0 (−1.9 to +3.3) | 31 |
| from diagnosis, isolated | 1.03 (0.69–1.53) | +0.2 (−2.5 to +2.0) | 50 |
| landmark, any | 0.95 (0.62–1.44) | −0.3 (−3.9 to +2.0) | 50 |

*(check D; positive means fewer relapses with HD-MTX; the difference is
anchored at the untreated rate using the subdistribution relation
1 − r₁ = (1 − r₀)^(1/HR))*

The data rule out the large effects historically hoped for, reductions of
roughly 35–45% or more. They remain compatible with preventing up to 2–3 CNS
relapses per 100 UHR patients treated. Whether that is "meaningful" depends on
HD-MTX's harms: treatment-related deaths, renal toxicity, and delays to
R-CHOP. The paper invokes those harms (p 2298, refs 19, 21, 30) but does not
quantify them in this cohort. Its conclusion that HD-MTX has "no meaningful
benefit for most patients" is a benefit–harm judgement, and the paper
supplies only the benefit side.

### 5. Subgroups

Each subgroup has 13–52 landmark events, so a finding of no significant
reduction in any individual subgroup was expected even under a moderate true
benefit. The main text reports a test within each subgroup but no
subgroup-by-treatment interaction tests. Fig 2's intervals are also distorted
(Inconsistencies, item 6).

Two subgroups deserve a note:

- **CNS-IPI 5, isolated relapse.** This signal (6.2% v 2.5%) is best read the
  way the authors read it. Corrected, the difference is +3.8 points (−0.4 to
  +5.2), not the +10.1 (−0.4 to 35) printed.
- **Testicular.** The authors report a non-significant trend in favour of
  HD-MTX. Fig 2 prints +4.7 points for any relapse, which corrects to +3.0
  (−2.4 to +5.9). Isolated relapse shows +0.5 points, which the correction
  leaves at +0.4, and the trend disappears on matching (p 2299). The
  IELSG30 regimen (ref 17) remains the only prospective data for this
  subgroup.

### 6. A re-analysis, not a replication

The abstract concludes that the data "strongly support previous studies".
Those previous studies are refs 19 and 20: the same two datasets and largely
the same patients. What is new is pairing study 1's treated patients with
study 2's untreated ones. That adds numbers but is not an independent test.
Lewis 2023 compared the two arms within one study's data collection; this
design mostly replaces that with a comparison between studies.

### 7. Smaller points

- **Dose analysis.** Its contrasts rest on few events. The reference group
  of three or more doses has 10 events, so "no dose effect" is weak evidence
  on its own.
- **Early relapse.** The landmark cannot test whether prophylaxis given
  during R-CHOP prevents early relapse. Study 1 addressed timing and found no
  difference.
- **Exploratory: relapse site.** Among landmark relapses with a known site,
  leptomeningeal-only relapse was 4/47 without HD-MTX and 12/46 with it
  (Fisher exact P = 0.030, uncorrected; *check G*). That fits intravenous MTX
  protecting the parenchyma better than the meninges. But the site is unknown
  for 28% v 8% of events, and this is one of many possible comparisons, so it
  is a hypothesis, not a finding.

## Inconsistencies in the published version

| # | where | printed | consistent reading | evidence |
| --- | --- | --- | --- | --- |
| 1 | Results, p 2292 | "no HD-MTX n = 872, HD-MTX n = 1,051" | 1,051 without, 872 with | abstract, Tables 1–2; Table 2's HD-MTX cycle rows sum to 872 |
| 2 | Survival Outcomes, p 2297 | landmark survival differences "remained significant" | event-free survival did (HR 1.22, P = .043); overall survival did not (HR 1.18, 0.94–1.49, P = .16) | Table 2 |
| 3 | Table 2, matched landmark EFS | HR 1.05 (1.84–1.33) | lower bound about 0.83; "1.84" is most likely 0.84 | the interval excludes its own estimate; *check J* |
| 4 | Table 2, per-cycle dose rows | "≥3 mg/m²" | ≥ 3 g/m² | the row heading says g/m², and 3 mg/m² is not a prophylactic dose |
| 5 | Table 3, cumulative-dose blocks | no HD-MTX 1.02 (0.57–1.82) and < 6 g/m² 0.96 (0.63–1.47); isolated 1.36 (0.73–2.53) and 0.96 (0.58–1.58) | the two rows' results appear transposed | each interval's width fits the *other* row's events, while the dose-count blocks fit to two decimals *(check B)* |
| 6 | Fig 2, "3-year difference" | r₀ × (HR − 1) | r₀ − r₁, about r₀ × (1 − 1/HR) | see below *(check A)* |
| 7 | p 2298 | UHR v non-UHR survival after relapse cited to "Fig 1D" | Fig 3D | Fig 1D is matched isolated CNS relapse |
| 8 | p 2297 | from-diagnosis rates 9.3% v 8.1% cited to Fig 1 | not shown in Fig 1 | Fig 1 covers only the landmark cohort (caption, axes) |
| 9 | Fig 1 caption | matching list omits breast involvement | Methods say matching used every adjusted factor, which includes breast | p 2292; Fig 2 footnote |
| 10 | Fig 2 against the text | adjusted landmark row implies HR 0.69–1.36 | text gives 0.95 (0.62–1.44) | *check A*; cannot be resolved without the supplement |

### Item 6 in detail

The Methods say the differences "were calculated by applying the hazard ratio
(HR) to the 3-year rate in the no HD-MTX group and taking the difference"
(p 2292). With the paper's HR, the hazard without HD-MTX over the hazard with
it, multiplying the untreated rate by HR does not give the treated rate;
dividing does. Read literally, the sentence computes the untreated rate times
(HR − 1), and the printed values fit that formula. Solving each row for its
implied HR interval gives an interval symmetric on the log scale, as a Wald
interval must be. The ratio of its two log-distances is 0.90–1.07 across all
18 rows, against 0.61–1.57 if the same differences are anchored at the
treated rate instead *(check A)*. Taking HR the other way round is not
possible at all: it would need a negative hazard ratio for the CNS-IPI 5 row.

The model-consistent difference is r₀ − r₁ with 1 − r₁ = (1 − r₀)^(1/HR),
close to r₀(1 − 1/HR). Each printed difference is therefore the correct one
multiplied by HR, exactly so for small rates. The error is always in the
same direction: benefit bounds are inflated and harm bounds shrunk. It is
small in the pooled point estimates but not in their bounds:

| landmark row | printed | corrected | printed rates, r₀ − r₁ |
| --- | --- | --- | --- |
| All UHR, any | +0.3 (−1.8 to +3.4) | +0.3 (−2.3 to +2.2) | +0.1 |
| All UHR adjusted, any | −0.3 (−2.1 to +2.4) | −0.3 (−2.9 to +1.7) | – |
| CNS-IPI 5–6, isolated | +4.7 (−1.4 to +17) | +2.7 (−1.7 to +4.7) | +2.5 |
| CNS-IPI 5, isolated | +10.1 (−0.4 to +35) | +3.8 (−0.4 to +5.2) | +3.7 |
| CNS-IPI 6, any | +0.4 (−7.5 to +20.7) | +0.4 (−17.8 to +7.3) | −2.8 |
| Testicular, any | +4.7 (−2.0 to +16.5) | +3.0 (−2.4 to +5.9) | +2.6 |

*(all 18 rows: check A)*

The distortion is largest where the hazard ratio is far from 1. In the three
rows with HR of 1.5 or more, the corrected differences (+2.7, +3.8, +3.0) sit
close to the rate differences printed in the same rows (+2.5, +3.7, +2.6),
and the published ones (+4.7, +10.1, +4.7) do not. In rows with HR near 1
the two versions nearly coincide, and their gaps from the printed rates are
within the noise of those rows' events, so they settle nothing. The CNS-IPI 5
row prints a 10-point difference against an untreated rate of 6.2%, which no
absolute difference can exceed. Its implied HR of 2.63 reproduces the printed
treated rate (2.40% against 2.5%). Because the error inflates apparent
benefit, correcting it makes the paper's own conclusion stronger.

## Where it sits in the evidence the paper cites

- **Wilson 2022 (ref 19).** 1,384 HD-MTX recipients. Giving HD-MTX between
  R-CHOP cycles was no better than giving it after R-CHOP. The crude CNS
  relapse rate at high CNS-IPI was 9.1% despite HD-MTX in every patient.
- **Lewis 2023 (ref 20).** 2,418 patients. It found a small significant
  reduction overall that did not survive a landmark restricted to complete
  responders. The 5-year adjusted risk difference was 1.4% (−1.5 to 4.1).
- **Two randomised trials of IT v HD-MTX (refs 31 and 32).** 142 and 100
  patients, reported as conference abstracts. Both stopped early and neither
  showed a difference.
- **IELSG30 (ref 17).** Single-arm, testicular LBCL: R-CHOP with IT therapy,
  contralateral testicular radiotherapy and two cycles of 1.5 g/m² IV MTX. No
  CNS relapse in 54 patients at 5 years, but no way to attribute that to any
  one component.
- **Ferreri 2026 (ref 34).** Single centre: 69 HD-MTX patients against a
  historical cohort, 0% v 14% CNS relapse. The historical control limits what
  it can show.

No randomised evidence shows a benefit. This paper narrows how large a benefit
could be in the highest-risk patients. The editor-in-chief's Relevance note
says these data "should definitively end the practice of CNS prophylaxis with
methotrexate in LBCL, even in patients with highest risk". The direction is
well supported. "Definitively" asks more than an observational, largely
between-study comparison with 115 landmark events can give: it rules out
large benefits but not a 2-point one, and the harm side of the comparison is
not in the paper.

## What would close the gap

1. Report the matched landmark CNS-relapse HR and absolute difference in the
   main text, and lead with them.
2. Stratify by study and site, or restrict to the nine sites that contributed
   both arms, and adjust for year of diagnosis.
3. Test subgroup-by-treatment interactions instead of significance within
   each subgroup.
4. Report harms in this cohort (treatment-related deaths, renal toxicity,
   R-CHOP delays), so that "no meaningful benefit" is a benefit–harm
   statement.
5. State how missing covariates were handled. Baseline CNS assessment is
   missing for 25.6% v 3.1%.
6. Correct Fig 2 and Table 3 (items 5 and 6).

## Reproducing the derived numbers

```bash
python3 docs/papers/wilson2026_checks.py
```

| section | what it computes |
| --- | --- |
| A | the Fig 2 formula and the corrected differences |
| B | Table 3 interval widths against event counts |
| C | power |
| D | absolute bounds and NNT |
| E | E-values |
| F | selection before the landmark and matching coverage |
| G | the exploratory relapse-site test |
| H | CNS staging and missingness |
| I | survival as a negative control |
| J | the misprinted EFS interval |

Inputs are transcribed in the script with page, table and figure references.
The PDF itself is not in the repository.
