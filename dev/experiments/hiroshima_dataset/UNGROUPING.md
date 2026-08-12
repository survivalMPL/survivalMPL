# Ungrouping LSS14 into individual records

How `data(hiroshima)` is built from `lss14/lss14.csv`, and what it is and isn't
good for. Canonical code: `dev/data-raw/hiroshima.R`. Verification:
`dev/experiments/hiroshima_reconstruction_checks.R`.

## The problem

`lss14.csv` is not individual records. It is a cross-tabulation: 53,782 cells
over city x sex x ground distance x AHS participation x age-at-exposure band x
attained-age band x calendar-time band x dose band, giving counts and
person-years per cell. Nothing in it says "here is a person".

## What makes the cohort recoverable

`subjects` is an **entry census**, not a risk-set count — each subject is
counted exactly once, in their first calendar-time cell:

```
tapply(raw$subjects, raw$ctime, sum)
    1     2     3  ...  11
86611     0     0  ...   0        sum = 86,611 = the cohort size
```

Had it been a risk-set count, one person would appear in ~11 calendar cells and
the cohort size would be unrecoverable.

The attained-age bands are given only as an integer code `agecat`. Recovered
empirically from the min/max of each code's mean age: `k -> [5(k-1), 5k)` for
k = 2..21, with band 21 open-ended (observed max 113.3).

## The reconstruction

The part that matters is `entry`. LSS follow-up opens 1950-10-01, 5.15 years
after the bombings, so a subject enters the risk set at the age they had
reached then — anyone who died in those five years is absent from the cohort
entirely. Hence `entry = agex + 5.15`, and then:

- **deaths** — one row per death in each cell (`rep.int` on the cell index,
  expanded by the `death` count), age at death drawn uniformly within the
  cell's five-year band, floored at `entry`. Covariates and dose from the cell.
  Within a cell the first `colon` rows are flagged as colon-cancer deaths.
- **survivors** — per subject stratum (the non-time-varying part of the cell
  key), `subjects - deaths` rows, censored at `agex + 58.4` (age on
  2003-12-31). Covariates and mean dose from the `ctime == 1` rows, where the
  subject counts live, subject-count weighted.

## Validation

Counts are easy to hit and prove little. The check that tests the reconstructed
*ages* rather than the row count is total person-years, since it equals
`sum(time - entry)` and the table reports it independently. All four are
`stopifnot()` assertions in the builder:

| quantity | source | rebuilt |
|---|---|---|
| subjects | 86,611 | 86,611 (exact) |
| deaths | 50,620 | 50,620 (exact) |
| colon-cancer deaths | 621 | 621 (exact) |
| person-years | 3,294,282 | 3,305,240 (+0.33%) |

## Known approximations

All three leave the risk sets *by attained age* intact, which is what a
left-truncation analysis depends on:

1. within-cell variation in `agex` and in dose collapses to the cell mean;
2. which subject in a stratum gets which death age is arbitrary (they are
   exchangeable within stratum given the covariates);
3. other-cause deaths are censoring for `status_colon` (cause-specific hazard).

---

# Open questions — resolved

## 1. How did Ozasa et al. do it? Is there a better approach?

**They did not ungroup at all.** Their own analysis script is in this directory
(`lss14/lss14.scr`, with output in `lss14.log`) — it is Epicure/AMFIT, fitting
Poisson regression directly to the grouped table:

```
strata city sex agexcat agecat@     ! baseline stratified, 1200 strata
pyr py1@  cases death@              ! person-years as offset
null@line 1 colon10@fit@            ! ERR linear in dose
logl 1 ew30 lage70@                 ! ERR modified by age at exposure
pline 1 %con = 1 msex@fit@          !   attained age, and sex
```

So their model is `lambda = lambda0(strata) * (1 + beta*d * modifiers)` — a
**product-additive excess relative risk** model, not proportional hazards, with
`beta` an ERR per Gy rather than a log hazard ratio.

For analysing *this table* their approach is strictly better: no
reconstruction, no arbitrary choices, full person-year weighting, and a much
richer baseline (1200 strata vs our single attained-age baseline). The
reconstruction exists only because `coxph_mpl()` needs one row per subject; it
buys individual-level records at the cost of approximation. It is not an
improvement on Ozasa et al. and should not be presented as one.

Worth doing: check whether the published Report 14 paper documents any
individual-level reconstruction, and whether RERF distributes a subject-level
file under any data agreement (which would retire this reconstruction).

## 2. "Ungrouping is basically pyr/subjects" — does that match what we did?

Partly — the ratio is the right *check*, but it is not a per-cell recipe, and
there are two distinct ideas hiding under it.

- **As a per-cell recipe it does not work.** `subjects` is nonzero in only
  3,420 of 53,782 cells (all at `ctime == 1`), so the ratio is undefined for
  47,378 cells. Worse, where it is defined it is not interpretable: it ranges
  0.01 to **536.9**, and 537 years cannot be a duration inside a 5-year band.
  The reason is that a subject is counted in the cell they *enter* but ages out
  of that attained-age band during the calendar band, contributing person-years
  to cells where they were never counted — numerator and denominator refer to
  different sets of people.

- **As an aggregate it is exactly right, and we already match it.** Over the
  whole table `sum(pyr)/sum(subjects) = 38.04` years is the mean follow-up per
  subject, and the reconstruction gives `mean(time - entry) = 38.16` (+0.33%) —
  the same 0.33% that appears in the validation table above. Restricted to
  `ctime == 1` the ratio is 5.12 years, i.e. one calendar band's width, which
  is a useful sanity check on the band structure.

- **The other reading is a different object.** Taking "at risk in a cell =
  pyr / band width" gives the piecewise-exponential **episode expansion**:
  658,807 person-interval rows rather than 86,611 subjects. That expansion is
  *lossless* for the table — it reproduces the grouped Poisson likelihood
  exactly — so it is the better choice if the goal is to match Ozasa et al.
  But every episode's entry time is a band boundary, so individual-level left
  truncation is not represented, which defeats the purpose here.

Summary: same underlying identity, different granularity. We reconstruct
subjects (86,611 rows, real entry ages); `pyr`-based expansion reconstructs
episodes (658,807 rows, boundary entry times).

## 3. "Illustrative only" — prove it

The earlier "treat the coefficients as illustrative" claim was **too
dismissive, and the scepticism was justified.** Two separate things need
separating, and they behave very differently.

**The truncation effect is real, not an artefact.** Benchmarked against the
grouped-data Poisson fit — `glm(death ~ dose + sex + city + factor(agecat),
offset = log(pyr), family = poisson)`, which *is* the piecewise-exponential
proportional-hazards likelihood for this table, uses the full person-year
weighting, and needs no reconstruction:

| | dose | sexFemale | cityNagasaki |
|---|---|---|---|
| (a) Poisson benchmark | 0.1618 | −0.5687 | 0.0565 |
| (b) reconstruction, LTRC Cox | 0.1674 | −0.5623 | 0.0580 |
| (c) reconstruction, naive Cox | 0.1862 | −0.5520 | 0.0757 |
| Poisson SE, `s` | 0.0139 | 0.0091 | 0.0098 |

The comparison is just the two differences from row (a), divided by `s`:

| | dose | sexFemale | cityNagasaki |
|---|---|---|---|
| (b) − (a) | +0.0056 | +0.0064 | +0.0015 |
| (c) − (a) | +0.0244 | +0.0167 | +0.0192 |
| ((b)−(a)) / `s` | **0.40** | **0.70** | **0.15** |
| ((c)−(a)) / `s` | **1.76** | **1.84** | **1.96** |
| how much further off (c) is | 4.4x | 2.6x | 12.8x |

The LTRC fit is nearer the benchmark on all three coefficients, by factors of
2.6 to 12.8.

Two things this does **not** say. `s` is the precision of the benchmark, used
as a ruler; it is *not* the standard error of the difference, because both fits
use the same underlying data and are strongly correlated. So "1.8 SE" means "the
naive estimate is off by about twice the noise level of the correct analysis",
not "significant at 5%". And the benchmark shares the cell-mean-dose collapse
with the reconstruction, so it cannot validate that particular approximation —
what it validates is the entry-age, death-placement and censoring assumptions.

### Why the Poisson fit is the right benchmark for a truncation question

Because the grouped person-time table is left-truncation-correct **by
construction**: `pyr` accumulates only time actually under observation, and
follow-up starts in 1950, so a subject who was 30 in 1950 contributes zero
person-years to every attained-age band below 30. The table cannot represent
pre-1950 risk time at all. A naive Cox fit invents risk time from age 0 to
entry for every subject.

That mechanism is testable: rebuild the person-time table *with* the fictitious
pre-entry time and the Poisson fit should land on the naive Cox estimate.
Splitting the same 20,000 subjects at 5-year attained-age boundaries and fitting
both ways (`Q3a-bis` in the checks script):

| | dose | sexFemale | cityNagasaki |
|---|---|---|---|
| Poisson, person-time from entry | 0.1813 | −0.5511 | 0.0342 |
| Cox LTRC, same subjects | 0.1861 | −0.5569 | 0.0346 |
| Poisson, person-time from age 0 | 0.1935 | −0.5404 | 0.0528 |
| Cox naive, same subjects | 0.1983 | −0.5465 | 0.0531 |

Within each pair the two estimators agree to 0.0058 and 0.0061, while the pairs
differ from each other by up to 0.019. Starting the clock at age 0 adds 675,000
fictitious person-years to 766,000 real ones — **+88%**. So the gap tracks the
person-time bookkeeping, not the estimator and not the reconstruction: swapping
Cox for Poisson changes nothing, swapping the risk-set definition changes
everything.

It is also **insensitive to every arbitrary choice** in the reconstruction. Over
12 variants (5 RNG seeds, 3 within-band placement rules, cell vs stratum dose,
entry offset 4.5/5.8y, follow-up end 56/60y) the naive−LTRC difference stays in
[0.0165, 0.0200] for dose and is **positive in all 12 for all three
coefficients**. Across seeds it is identical to four decimal places, because it
depends only on the entry-age distribution, which is fixed by `agex` — a column
the table reports directly and the reconstruction does not invent.

**The absolute coefficients are the part to hedge on.** The same 12 variants
move them by 1.8 SE (dose), 2.9 SE (city) and 4.1 SE (sex) — mostly driven by
`placement = cellmean` and by the assumed end of follow-up. So the absolute
values carry real reconstruction uncertainty beyond their nominal SEs.

Even so they are better than "illustrative". Matching Ozasa et al.'s
stratification as closely as the reconstruction allows (attained age as time
scale, stratified on city and 5-year `agex` bands) and converting a log hazard
ratio to their ERR scale via `exp(b) − 1`:

| | reconstruction ERR/Gy | Ozasa et al. (`lss14.log`) |
|---|---|---|
| males | 0.1450 [0.1015, 0.1902] | 0.1500 [0.1035, 0.1986] |
| females | 0.2764 [0.2281, 0.3265] | 0.2978 [0.2438, 0.3539] |

Both point estimates fall inside the published profile-likelihood intervals,
despite the differing model form (linear-in-dose excess vs log-linear) and much
coarser baseline. Both are ~3–7% low, consistent with `exp(b) − 1` understating
a linear ERR plus the coarser stratification.

**Revised wording:** the dataset supports conclusions about *left truncation*
robustly, and recovers published dose-response estimates to within their
confidence intervals. What it does not support is quoting a coefficient to
three decimals as an LSS14 estimate — the reconstruction adds 2–4 SE of
uncertainty that no standard error reports.
