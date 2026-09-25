#!/usr/bin/env python3
"""Recompute every derived number in wilson-2026-hd-mtx-cns-prophylaxis.md.

Wilson MR, Lewis KL, Kirkwood AA, et al. High-dose methotrexate as CNS
prophylaxis in ultra high-risk large B-cell lymphoma: an international
multicenter analysis. J Clin Oncol 44:2290-2302, 2026. doi:10.1200/JCO-26-00947

Inputs are transcribed from the paper's main text, Tables 1-3 and Fig 2; each
block says where. Nothing here comes from the Data Supplement, which was not
available. Standard library only:

    python3 docs/papers/wilson2026_checks.py

Conventions. HR is the paper's: hazard WITHOUT HD-MTX relative to WITH HD-MTX
(the abstract's "9.3% v 8.1%; adjusted HR 1.13"), so HR > 1 favours HD-MTX.
r0 and r1 are 3-year cumulative incidences in % without / with HD-MTX. An
absolute difference is r0 - r1, positive favouring HD-MTX, as on Fig 2's axis.
"""

import math

Z975 = 1.959964   # two-sided alpha = 0.05
Z80 = 0.841621    # power = 0.80


def r1_from_hr(r0, hr):
    """3-year incidence WITH HD-MTX implied by r0 and a subdistribution HR.

    Under proportional subdistribution hazards 1 - F1 = (1 - F0)^(1/HR) when
    HR = h0/h1. For incidences of a few percent this is r0/HR to within a
    few hundredths of a point.
    """
    return 100.0 * (1.0 - (1.0 - r0 / 100.0) ** (1.0 / hr))


def section(title):
    print('\n' + title)
    print('-' * len(title))


# ---------------------------------------------------------------------------
# A. Fig 2 (p 2297): the "3-year difference" column.
#    (row, outcome, reported diff, lo, hi, r0, r1, events0/n0, events1/n1)
# ---------------------------------------------------------------------------
FIG2 = [
    ('All UHR', 'any', +0.3, -1.8, 3.4, 6.7, 6.6, '65/782', '50/773'),
    ('All UHR', 'isolated', 0.0, -1.6, 2.4, 4.6, 4.8, '47/782', '37/773'),
    ('All UHR adj.', 'any', -0.3, -2.1, 2.4, 6.7, None, '', ''),
    ('All UHR adj.', 'isolated', -0.8, -2.1, 1.3, 4.6, None, '', ''),
    ('CNS-IPI 5-6', 'any', -0.1, -3.8, 6.2, 8.8, 9.1, '29/315', '19/226'),
    ('CNS-IPI 5-6', 'isolated', +4.7, -1.4, 17.0, 6.5, 4.0, '21/315', '8/226'),
    ('CNS-IPI 5', 'any', +0.1, -3.9, 7.4, 8.2, 8.0, '21/262', '14/189'),
    ('CNS-IPI 5', 'isolated', +10.1, -0.4, 35.0, 6.2, 2.5, '16/262', '4/189'),
    ('CNS-IPI 6', 'any', +0.4, -7.5, 20.7, 11.6, 14.4, '8/53', '5/37'),
    ('CNS-IPI 6', 'isolated', -1.4, -6.0, 13.9, 7.8, 11.2, '5/53', '4/37'),
    ('>=3 EN sites', 'any', +2.0, -2.2, 8.9, 8.3, 6.2, '27/317', '25/398'),
    ('>=3 EN sites', 'isolated', +1.5, -1.9, 7.8, 5.5, 3.9, '18/317', '16/398'),
    ('Testicular', 'any', +4.7, -2.0, 16.5, 9.4, 6.8, '19/153', '15/221'),
    ('Testicular', 'isolated', +0.5, -3.0, 7.5, 6.1, 6.3, '13/153', '14/221'),
    ('Breast', 'any', -0.9, -3.3, 6.5, 4.3, 6.9, '11/117', '4/61'),
    ('Breast', 'isolated', -0.4, -2.7, 7.7, 3.5, 5.3, '10/117', '3/61'),
    ('Renal/adrenal', 'any', +0.9, -3.6, 8.5, 10.1, 9.8, '25/221', '24/265'),
    ('Renal/adrenal', 'isolated', +1.4, -2.7, 9.1, 6.9, 6.0, '17/221', '15/265'),
]


def fig2():
    section('A. Fig 2: the reported differences are r0*(HR - 1)')
    print('If each reported difference is d = r0*(HR - 1), then HR = 1 + d/r0, and the\n'
          'implied HR interval must be symmetric on the log scale (a Wald interval).\n'
          '"sym" = log-distance above / log-distance below the point; 1.00 is exact.\n'
          'The subdistribution-consistent difference is r0 - r1_from_hr(r0, HR).\n')
    hdr = ('%-14s %-8s %-22s %-19s %5s %-24s %s'
           % ('row', 'outcome', 'reported d (95% CI)', 'implied HR (CI)', 'sym',
              'model-consistent d (CI)', 'displayed r0-r1'))
    print(hdr)
    print('-' * len(hdr))
    for name, outcome, d, lo, hi, r0, r1, _, _ in FIG2:
        hr, hr_lo, hr_hi = (1 + x / r0 for x in (d, lo, hi))
        up = math.log(hr_hi) - math.log(hr)
        dn = math.log(hr) - math.log(hr_lo)
        sym = up / dn if abs(dn) > 1e-12 else float('nan')
        c, c_lo, c_hi = (r0 - r1_from_hr(r0, h) for h in (hr, hr_lo, hr_hi))
        shown = '%+.1f' % (r0 - r1) if r1 is not None else '-'
        print('%-14s %-8s %+5.1f (%+5.1f, %+5.1f)   %4.2f (%4.2f, %4.2f)  %5.2f  '
              '%+5.1f (%+6.1f, %+5.1f)    %s'
              % (name, outcome, d, lo, hi, hr, hr_lo, hr_hi, sym, c, c_lo, c_hi, shown))
    syms0, syms1 = [], []
    for name, outcome, d, lo, hi, r0, r1, _, _ in FIG2:
        for anchor, out in ((r0, syms0), (r1, syms1)):
            if anchor is None:
                continue
            h, h_lo, h_hi = (1 + x / anchor for x in (d, lo, hi))
            if min(h, h_lo, h_hi) <= 0 or abs(math.log(h) - math.log(h_lo)) < 1e-12:
                out.append(float('nan'))
                continue
            out.append((math.log(h_hi) - math.log(h)) / (math.log(h) - math.log(h_lo)))
    ok0 = [x for x in syms0 if x == x]
    ok1 = [x for x in syms1 if x == x]
    print('\nsym anchored at the untreated rate r0: %.2f-%.2f over %d rows; anchored at the\n'
          'treated rate r1 instead: %.2f-%.2f over %d rows. Anchoring at r0 with HR = h0/h1\n'
          'is the only reading that gives Wald-shaped intervals throughout.'
          % (min(ok0), max(ok0), len(ok0), min(ok1), max(ok1), len(ok1)))
    print('\nReported = model-consistent x HR, exactly in the linear limit: r0*(HR-1)\n'
          '= HR * r0*(1 - 1/HR). The CNS-IPI 5 isolated row reports +10.1 (-0.4 to 35)\n'
          'percentage points where the untreated rate is 6.2%%, so no absolute\n'
          'difference in that row can exceed 6.2. Its implied HR of 2.63 reproduces\n'
          'the displayed treated rate: r1_from_hr(6.2, 2.63) = %.2f%% against 2.5%%.'
          % r1_from_hr(6.2, 2.63))
    # The adjusted row against the text's adjusted landmark HR, 0.95 (0.62-1.44).
    lo_t, hi_t = 6.7 * (0.62 - 1), 6.7 * (1.44 - 1)
    print('\nText (p 2297): landmark adjusted HR 0.95 (0.62-1.44). Under Fig 2\'s own\n'
          'formula that interval would print as (%+.1f to %+.1f); Fig 2 prints\n'
          '(-2.1 to +2.4), i.e. HR (0.69-1.36). The two adjusted intervals differ.'
          % (lo_t, hi_t))


# ---------------------------------------------------------------------------
# B. Table 3 (p 2298): interval widths against event counts.
#    (block, row, events, HR, lo, hi, reference events)
# ---------------------------------------------------------------------------
TABLE3 = [
    ('any, cumulative dose', 'No HD-MTX', 65, 1.02, 0.57, 1.82, 33),
    ('any, cumulative dose', '<6 g/m2', 14, 0.96, 0.63, 1.47, 33),
    ('any, doses >=3 g/m2', 'No HD-MTX', 65, 1.06, 0.54, 2.06, 10),
    ('any, doses >=3 g/m2', '1 dose', 5, 1.23, 0.42, 3.58, 10),
    ('any, doses >=3 g/m2', '2 doses', 19, 0.89, 0.42, 1.92, 10),
    ('isolated, cumulative', 'No HD-MTX', 47, 1.36, 0.73, 2.53, 23),
    ('isolated, cumulative', '<6 g/m2', 13, 0.96, 0.58, 1.58, 23),
    ('isolated, doses', 'No HD-MTX', 47, 1.50, 0.59, 3.77, 5),
    ('isolated, doses', '1 dose', 4, 2.00, 0.54, 7.43, 5),
    ('isolated, doses', '2 doses', 15, 1.43, 0.52, 3.95, 5),
]


def table3():
    section('B. Table 3: are the interval widths those the event counts allow?')
    print('For a two-group hazard ratio, SE(log HR) is close to sqrt(1/d1 + 1/d2).\n'
          'Reported SE = (ln hi - ln lo)/(2*1.96). A ratio near 1 means the interval\n'
          'fits the row\'s events; the reference is >=6 g/m2 or >=3 doses >=3 g/m2.\n')
    hdr = '%-22s %-10s %6s %5s %-13s %8s %8s %6s' % (
        'block', 'row', 'events', 'HR', '95% CI', 'SE rep', 'SE exp', 'ratio')
    print(hdr)
    print('-' * len(hdr))
    for block, row, d, hr, lo, hi, dref in TABLE3:
        se_rep = (math.log(hi) - math.log(lo)) / (2 * Z975)
        se_exp = math.sqrt(1.0 / d + 1.0 / dref)
        print('%-22s %-10s %3d/%-2d %5.2f %4.2f-%-8.2f %8.3f %8.3f %6.2f'
              % (block, row, d, dref, hr, lo, hi, se_rep, se_exp, se_rep / se_exp))
    a = math.sqrt(1 / 65 + 1 / 33)
    b = math.sqrt(1 / 14 + 1 / 33)
    print('\nIn both cumulative-dose blocks the "No HD-MTX" interval has the width of a\n'
          '14-event comparison and the "<6" interval that of a 65-event one (any:\n'
          'expected SE %.3f and %.3f; reported %.3f and %.3f). The dose-count blocks\n'
          'match to the second decimal. The two rows\' results look transposed.'
          % (a, b, (math.log(1.82) - math.log(0.57)) / (2 * Z975),
             (math.log(1.47) - math.log(0.63)) / (2 * Z975)))


# ---------------------------------------------------------------------------
# C. Power (Schoenfeld): smallest hazard ratio detectable with 80% power.
#    (analysis, events without, events with, n without, n with) from Table 2
#    (p 2295) and Fig 2 (p 2297).
# ---------------------------------------------------------------------------
POWER = [
    ('all UHR from diagnosis, any', 111, 71, 1051, 872),
    ('all UHR from diagnosis, isolated', 74, 51, 1051, 872),
    ('landmark, any', 65, 50, 782, 773),
    ('landmark, isolated', 47, 37, 782, 773),
    ('landmark CNS-IPI 5-6, any', 29, 19, 315, 226),
    ('landmark >=3 EN sites, any', 27, 25, 317, 398),
    ('landmark renal/adrenal, any', 25, 24, 221, 265),
    ('landmark testicular, any', 19, 15, 153, 221),
    ('landmark breast, any', 11, 4, 117, 61),
    ('landmark CNS-IPI 6, any', 8, 5, 53, 37),
]


def power():
    section('C. Power: the smallest effect each analysis could reliably detect')
    print('Schoenfeld: |ln HR| = (z_0.975 + z_0.80) / sqrt(D p (1 - p)), D = events,\n'
          'p = share of patients given HD-MTX. Reported as the relative reduction in\n'
          'hazard that HD-MTX would need to have for 80% power at two-sided 0.05.\n')
    hdr = '%-36s %7s %6s %11s' % ('analysis', 'events', 'p', 'reduction')
    print(hdr)
    print('-' * len(hdr))
    for name, d0, d1, n0, n1 in POWER:
        dd, p = d0 + d1, n1 / (n0 + n1)
        hr = math.exp(-(Z975 + Z80) / math.sqrt(dd * p * (1 - p)))
        print('%-36s %7d %6.3f %10.0f%%' % (name, dd, p, 100 * (1 - hr)))


# ---------------------------------------------------------------------------
# D. What the adjusted hazard ratios allow in absolute terms (text, p 2297).
# ---------------------------------------------------------------------------
BOUNDS = [
    ('from diagnosis, any (adjusted)', 1.13, 0.82, 1.57, 9.3),
    ('from diagnosis, isolated (adjusted)', 1.03, 0.69, 1.53, 5.9),
    ('landmark, any (adjusted)', 0.95, 0.62, 1.44, 6.7),
]


def bounds():
    section('D. Absolute effect the adjusted hazard ratios are compatible with')
    print('Anchored at the untreated 3-year rate r0: difference = r0 - r1_from_hr.\n'
          'The upper HR bound gives the largest benefit the interval allows.\n')
    hdr = '%-38s %-17s %5s %-22s %s' % (
        'analysis', 'HR (95% CI)', 'r0', 'difference, pp (CI)', 'NNT at max benefit')
    print(hdr)
    print('-' * len(hdr))
    for name, hr, lo, hi, r0 in BOUNDS:
        d, d_lo, d_hi = (r0 - r1_from_hr(r0, h) for h in (hr, lo, hi))
        print('%-38s %4.2f (%4.2f-%4.2f)  %4.1f  %+4.1f (%+4.1f, %+4.1f)      %.0f'
              % (name, hr, lo, hi, r0, d, d_lo, d_hi, 100 / d_hi))


# ---------------------------------------------------------------------------
# E. E-values: how strong unmeasured confounding would have to be to hide a
#    true benefit (VanderWeele & Ding 2017; HR treated as a risk ratio, which
#    is reasonable for outcomes under ~15%).
# ---------------------------------------------------------------------------
def evalues():
    section('E. Confounding strength that could hide a true benefit')
    print('Observed hazard ratio of HD-MTX against none is 1/HR. To be hiding a true\n'
          'ratio T < 1, confounding must bias it by B = (1/HR)/T; the E-value\n'
          'B + sqrt(B (B - 1)) is the risk ratio an unmeasured confounder would need\n'
          'with BOTH receipt of HD-MTX and CNS relapse.\n')
    hdr = '%-30s %9s %-24s %-24s' % ('analysis', 'observed', 'true 20% reduction', 'true 30% reduction')
    print(hdr)
    print('-' * len(hdr))
    for name, hr in (('from diagnosis, any (1.13)', 1.13), ('landmark, any (0.95)', 0.95)):
        obs = 1 / hr
        cells = []
        for true in (0.8, 0.7):
            b = obs / true
            cells.append('B %.2f, E-value %.2f' % (b, b + math.sqrt(b * (b - 1))))
        print('%-30s %9.3f %-24s %-24s' % (name, obs, cells[0], cells[1]))


# ---------------------------------------------------------------------------
# F. Selection before the landmark, matching coverage, source of each arm.
# ---------------------------------------------------------------------------
def selection():
    section('F. Who leaves before the landmark, and who gets matched')
    for name, n, nl, ev, evl in (('no HD-MTX', 1051, 782, 111, 65),
                                 ('HD-MTX', 872, 773, 71, 50)):
        print('%-10s patients %4d -> %3d at 6 months (%4.1f%% gone); CNS events %3d -> %2d '
              '(%2d, %4.1f%%, before the landmark or in excluded patients)'
              % (name, n, nl, 100 * (n - nl) / n, ev, evl, ev - evl, 100 * (ev - evl) / ev))
    print('propensity-matched pairs: %d of %d HD-MTX patients (%.1f%%) from diagnosis, '
          '%d of %d (%.1f%%) in the landmark cohort'
          % (618, 872, 100 * 618 / 872, 501, 773, 100 * 501 / 773))
    print('HD-MTX patients in the amalgamated cohort from study 1: 1,384 of 1,625 (%.1f%%);'
          ' no-HD-MTX patients from study 2: all 1,866' % (100 * 1384 / 1625))


# ---------------------------------------------------------------------------
# G. Exploratory: site of CNS relapse in the landmark cohort (Table 2, p 2295).
# ---------------------------------------------------------------------------
def fisher_two_sided(a, b, c, d):
    """Exact two-sided P for the 2x2 table [[a, b], [c, d]]."""
    r1, r2, c1, n = a + b, c + d, a + c, a + b + c + d

    def prob(x):
        return math.comb(r1, x) * math.comb(r2, c1 - x) / math.comb(n, c1)
    p_obs = prob(a)
    lo, hi = max(0, c1 - r2), min(r1, c1)
    return sum(prob(x) for x in range(lo, hi + 1) if prob(x) <= p_obs * (1 + 1e-9))


def relapse_site():
    section('G. Exploratory: leptomeningeal-only relapse in the landmark cohort')
    # known sites: parenchymal only, leptomeningeal only, both; unknown separately
    no = dict(par=34, lep=4, both=9, unk=18)
    yes = dict(par=26, lep=12, both=8, unk=4)
    kn0 = no['par'] + no['lep'] + no['both']
    kn1 = yes['par'] + yes['lep'] + yes['both']
    p = fisher_two_sided(no['lep'], kn0 - no['lep'], yes['lep'], kn1 - yes['lep'])
    print('leptomeningeal only: %d/%d known sites without HD-MTX, %d/%d with; Fisher exact P = %.3f'
          % (no['lep'], kn0, yes['lep'], kn1, p))
    print('site unknown: %d of %d events without HD-MTX (%.0f%%), %d of %d with (%.0f%%)'
          % (no['unk'], kn0 + no['unk'], 100 * no['unk'] / (kn0 + no['unk']),
             yes['unk'], kn1 + yes['unk'], 100 * yes['unk'] / (kn1 + yes['unk'])))



# ---------------------------------------------------------------------------
# H. Baseline CNS staging and missing data (Table 1, p 2293), as shares of
#    ALL patients in each arm rather than of those with a known answer.
# ---------------------------------------------------------------------------
def staging():
    section('H. Baseline CNS staging and missingness, as a share of each whole arm')
    rows = (
        ('CNS imaging (with or without CSF)', 37 + 28, 75 + 47),
        ('CSF analysis (with or without imaging)', 223 + 37, 292 + 75),
        ('baseline CNS assessment: unknown', 269, 27),
        ('double/triple-hit status: unknown', 817, 268),
    )
    for name, a, b in rows:
        print('%-40s no HD-MTX %4d/1051 (%4.1f%%)   HD-MTX %4d/872 (%4.1f%%)'
              % (name, a, 100 * a / 1051, b, 100 * b / 872))
    known0, known1 = 1051 - 29, 872 - 13
    print('%-40s no HD-MTX %4d/%d (%4.1f%%)   HD-MTX %4d/%d (%4.1f%%)'
          % ('CNS-IPI 5-6, of those with a score', 399 + 86, known0,
             100 * (399 + 86) / known0, 218 + 50, known1, 100 * (218 + 50) / known1))


# ---------------------------------------------------------------------------
# I. Survival as a negative control (Table 2, p 2295). HD-MTX can change
#    survival only through the CNS relapses it prevents.
# ---------------------------------------------------------------------------
def survival_control():
    section('I. Survival as a negative control')
    for name, os0, os1, hr in (('from diagnosis, unmatched', 64.6, 79.2, '1.65 (1.38-1.98) adj.'),
                               ('from diagnosis, matched', 65.7, 79.1, '1.69 (1.37-2.08)'),
                               ('landmark, unmatched', 79.7, 83.8, '1.18 (0.94-1.49) adj.'),
                               ('landmark, matched', 82.7, 82.6, '0.99 (0.75-1.30)')):
        print('%-26s 3-year OS %4.1f%% v %4.1f%%  gap %+5.1f points  OS HR %s'
              % (name, os0, os1, os1 - os0, hr))
    fatal = 221 / 264
    best = 9.3 - r1_from_hr(9.3, 1.57)
    print('\nLargest CNS-relapse reduction the adjusted interval allows from diagnosis:\n'
          '%.1f points (section D). With %d of %d CNS relapses fatal (%.0f%%, median OS\n'
          '3.8 months), that could account for about %.1f points of the 3-year OS gap,\n'
          'against the 13.4-14.6 points observed from diagnosis.'
          % (best, 221, 264, 100 * fatal, best * fatal))



# ---------------------------------------------------------------------------
# J. Table 2 (p 2295): matched landmark EFS is printed HR 1.05 (1.84-1.33).
# ---------------------------------------------------------------------------
def efs_interval():
    section('J. The matched landmark EFS interval printed as 1.05 (1.84-1.33)')
    lo = 1.05 ** 2 / 1.33
    lo_min, lo_max = 1.045 ** 2 / 1.335, 1.055 ** 2 / 1.325
    print('A Wald interval is symmetric on the log scale, so the lower bound that goes\n'
          'with 1.05 and 1.33 is 1.05^2/1.33 = %.3f (%.3f-%.3f allowing for rounding):\n'
          'the printed "1.84" is most likely 0.84.' % (lo, lo_min, lo_max))


if __name__ == '__main__':
    fig2()
    table3()
    power()
    bounds()
    evalues()
    selection()
    relapse_site()
    staging()
    survival_control()
    efs_interval()
