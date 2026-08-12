# Constructs the `hiroshima` package dataset from RERF's Life Span Study
# Report 14 grouped person-time table (dev/experiments/hiroshima_dataset/), by
# reconstructing approximate INDIVIDUAL-level left-truncated / right-censored
# records - one row per subject, 86,611 subjects.
#
# HISTORY: an earlier version of this dataset treated each of the 53,782
# grouped-table cells as a single pseudo-subject (entry = cell mean age at
# exposure, time = cell mean attained age). That is NOT usable as an analysis
# dataset: it inverts the dose-response (it gives a *negative* colon-dose
# coefficient for all-cause mortality), because a cell with high dose is a
# cell with many person-years and therefore a long mean attained age, and that
# artefact swamps the real signal. It has been replaced by the reconstruction
# below. See dev/experiments/demo_hiroshima_ungrouped_ltrc.R, which refits the
# old pseudo-subject version alongside this one to show the difference.
#
# ---------------------------------------------------------------------------
# STRUCTURE OF THE SOURCE TABLE (dev/experiments/hiroshima_dataset/lss14/)
# ---------------------------------------------------------------------------
#   cell     : city x sex x gd3 x ahs x agexcat x agecat x ctime x dosecat
#   subjects : counted ONCE per subject, in that subject's ctime == 1 cell.
#              Verified: sum(subjects) over the whole table = 86,611 = the
#              cohort size, and tapply(subjects, ctime, sum) is 86,611 for
#              ctime == 1 and 0 for ctime 2..11. This is what makes the cohort
#              recoverable from the grouped table at all.
#   pyr      : person-years contributed to the cell (total 3,294,282)
#   death    : deaths in the cell (total 50,620); colon = colon-cancer deaths
#   agex     : person-year weighted mean age at exposure in the cell
#   agecat   : attained-age band, k -> [5*(k-1), 5*k) years, for k = 2..21
#              (verified against the observed min/max of `age` per band; band
#              21 is open-ended, observed max age 113.3)
#
# ---------------------------------------------------------------------------
# RECONSTRUCTION
# ---------------------------------------------------------------------------
# The scientifically important part is `entry`. LSS follow-up begins
# 1950-10-01, five years after the 1945 bombings, so a subject enters the risk
# set only at the attained age they had reached on that date - anyone exposed
# in 1945 who died before follow-up began is not in the cohort at all. Hence:
#
#   entry  = agex + YRS_TO_STUDY_START        (attained age on 1950-10-01)
#   deaths : one record per death in a cell; the attained age at death is drawn
#            uniformly inside that cell's agecat band, truncated below at
#            entry. Covariates and dose come from the cell.
#   others : subjects - deaths per subject stratum, censored administratively
#            at age agex + YRS_TO_STUDY_END   (attained age on 2003-12-31)
#
# APPROXIMATIONS (all of them leave the risk sets by attained age intact,
# which is what a left-truncation analysis turns on):
#   - within-cell variation in agex and in dose is collapsed to the cell mean;
#   - which subject in a stratum gets which death age is arbitrary;
#   - deaths from other causes are treated as censoring for `status_colon`
#     (cause-specific hazard, the usual convention).
#
# VALIDATION (asserted below): subjects, deaths and colon deaths reproduce the
# source table exactly, and the total person-years - which is the check that
# actually tests the reconstructed entry/exit ages rather than just the counts
# - comes within 0.5%.
#
# This remains a demonstration dataset for left truncation, not a reanalysis
# of LSS Report 14: Ozasa et al. fit grouped-data Poisson excess-relative-risk
# models to the person-years, which is a different estimand.

# ---------------------------------------------------------------------------
# OPTIONS. Defaults produce the packaged dataset. Each can be overridden by
# defining it before source()-ing this file, which is how
# dev/experiments/hiroshima_reconstruction_checks.R probes the sensitivity of
# the results to every arbitrary choice made here.
# ---------------------------------------------------------------------------
opt <- function(name, default) if (exists(name)) get(name) else default

YRS_TO_STUDY_START <- opt("HIROSHIMA_YRS_START", 5.15) # 1945-08 -> 1950-10-01
YRS_TO_STUDY_END <- opt("HIROSHIMA_YRS_END", 58.40) # 1945-08 -> 2003-12-31
SEED <- opt("HIROSHIMA_SEED", 20260811)
# where in its attained-age band a death is placed: "uniform" draw (default),
# band "midpoint", or the cell's own reported mean attained age ("cellmean")
PLACEMENT <- opt("HIROSHIMA_PLACEMENT", "uniform")
# dose assigned to a death: its own "cell" mean (default) or its subject
# "stratum" mean, i.e. the same value the stratum's survivors receive
DOSE_LEVEL <- opt("HIROSHIMA_DOSE", "cell")
SAVE <- opt("HIROSHIMA_SAVE", TRUE)

set.seed(SEED)
raw <- read.csv("dev/experiments/hiroshima_dataset/lss14/lss14.csv")

# subject strata = the fixed (non-time-varying) part of the cell definition
strat_all <- with(raw, paste(city, sex, gd3, ahs, agexcat, dosecat, sep = "|"))

## ---- 1. one record per death ---------------------------------------------
# `idx` maps each expanded row back to its cell
d <- raw[raw$death > 0, , drop = FALSE]
idx <- rep.int(seq_len(nrow(d)), d$death)

lo <- 5 * (d$agecat[idx] - 1)
hi <- ifelse(d$agecat[idx] == 21, 115, 5 * d$agecat[idx])
ent <- d$agex[idx] + YRS_TO_STUDY_START

# place the age at death in the band, but never at or below entry
a <- pmax(lo, ent + 1e-3)
b <- pmax(hi, a + 1e-3)

# inside each cell, the first `colon` of the `death` records are colon deaths
colon_flag <- as.integer(sequence(d$death) <= d$colon[idx])

t_death <- switch(PLACEMENT,
  uniform  = a + runif(length(a)) * (b - a),
  midpoint = pmin(pmax((lo + hi) / 2, a), b),
  cellmean = pmin(pmax(d$age[idx], a), b),
  stop("unknown HIROSHIMA_PLACEMENT: ", PLACEMENT)
)

deaths <- data.frame(
  city = d$city[idx],
  sex = d$sex[idx],
  gd3 = d$gd3[idx],
  agex = d$agex[idx],
  entry = ent,
  time = t_death,
  status = 1L,
  status_colon = colon_flag,
  dose = d$colon10[idx],
  strat = strat_all[raw$death > 0][idx]
)

## ---- 2. one record per surviving subject ---------------------------------
strat <- strat_all
n_sub <- tapply(raw$subjects, strat, sum)
n_die <- tapply(raw$death, strat, sum)[names(n_sub)]

# the ctime == 1 rows carry the subject counts, so each stratum's covariates
# and its mean agex / mean dose are taken from there, subject-count weighted
keep <- raw$ctime == 1 & raw$subjects > 0
e1 <- raw[keep, , drop = FALSE]
s1 <- strat[keep]
wmean <- function(x) tapply(x * e1$subjects, s1, sum) / tapply(e1$subjects, s1, sum)
first <- function(x) tapply(x, s1, function(z) z[1])

ref <- data.frame(
  strat = names(tapply(e1$subjects, s1, sum)),
  city = first(e1$city), sex = first(e1$sex), gd3 = first(e1$gd3),
  agex = wmean(e1$agex), dose = wmean(e1$colon10),
  stringsAsFactors = FALSE
)
ref$n_cens <- pmax(as.integer(n_sub[ref$strat]) - as.integer(n_die[ref$strat]), 0L)

j <- rep.int(seq_len(nrow(ref)), ref$n_cens)
cens <- data.frame(
  city = ref$city[j],
  sex = ref$sex[j],
  gd3 = ref$gd3[j],
  agex = ref$agex[j],
  entry = ref$agex[j] + YRS_TO_STUDY_START,
  time = ref$agex[j] + YRS_TO_STUDY_END,
  status = 0L,
  status_colon = 0L,
  dose = ref$dose[j],
  strat = ref$strat[j]
)

## ---- 3. combine, label, shuffle ------------------------------------------
hiroshima <- rbind(deaths, cens)

# Optional diagnostic variant: give each death its subject stratum's mean dose
# instead of its own cell's mean dose. Within-dosecat variation in the cell
# mean dose is what wrecked the old pseudo-subject dataset (high-dose cells are
# cells with many person-years, hence long mean attained age), so if the dose
# coefficient is sensitive to this switch, the reconstruction is still leaking
# that artefact. See hiroshima_reconstruction_checks.R.
if (identical(DOSE_LEVEL, "stratum")) {
  m <- match(hiroshima$strat, ref$strat)
  hiroshima$dose <- ifelse(is.na(m), hiroshima$dose, ref$dose[m])
}
hiroshima$strat <- NULL
hiroshima$city <- factor(hiroshima$city, levels = c(1, 2), labels = c("Hiroshima", "Nagasaki"))
hiroshima$sex <- factor(hiroshima$sex, levels = c(1, 2), labels = c("Male", "Female"))
hiroshima$gd3 <- factor(hiroshima$gd3, levels = c(1, 2), labels = c("<3km", "3-10km"))
# NOTE: dose is stored in Gy, not mGy as in the source column `colon10`. Gy is
# the LSS14 reporting scale and makes a coefficient read as "per Gy"; per-mGy
# coefficients are ~1e-4 and awkward to report.
hiroshima$dose <- hiroshima$dose / 1000
hiroshima <- hiroshima[sample.int(nrow(hiroshima)), , drop = FALSE]
rownames(hiroshima) <- NULL

## ---- 4. validate against the source table --------------------------------
pyr_reconstructed <- sum(hiroshima$time - hiroshima$entry)
pyr_source <- sum(raw$pyr)

stopifnot(
  all(hiroshima$time > hiroshima$entry),
  nrow(hiroshima) == sum(raw$subjects), # 86,611
  sum(hiroshima$status) == sum(raw$death), # 50,620
  sum(hiroshima$status_colon) == sum(raw$colon) # 621
)
# The person-year total is the check that actually tests the reconstructed
# entry/exit AGES rather than just the row counts, so it is an assertion for
# the packaged configuration. The diagnostic placement rules deliberately move
# the death ages, so for those it is only reported.
if (identical(PLACEMENT, "uniform") &&
  isTRUE(all.equal(YRS_TO_STUDY_START, 5.15)) &&
  isTRUE(all.equal(YRS_TO_STUDY_END, 58.40))) {
  stopifnot(abs(pyr_reconstructed / pyr_source - 1) < 0.005)
}

if (!SAVE) {
  message(sprintf(
    "Built hiroshima in memory (seed %s, placement %s, dose %s): %+.2f%% person-years",
    SEED, PLACEMENT, DOSE_LEVEL, 100 * (pyr_reconstructed / pyr_source - 1)
  ))
} else {
save(hiroshima, file = "data/hiroshima.rda", compress = "xz")
message(sprintf(
  paste0(
    "Saved data/hiroshima.rda\n",
    "  subjects     : %d (source %d)\n",
    "  deaths       : %d (source %d)\n",
    "  colon deaths : %d (source %d)\n",
    "  person-years : %.0f (source %.0f, %+.2f%%)\n",
    "  entry age    : %.1f - %.1f\n",
    "  exit age     : %.1f - %.1f\n",
    "  file size    : %.2f MB"
  ),
  nrow(hiroshima), sum(raw$subjects),
  sum(hiroshima$status), sum(raw$death),
  sum(hiroshima$status_colon), sum(raw$colon),
  pyr_reconstructed, pyr_source, 100 * (pyr_reconstructed / pyr_source - 1),
  min(hiroshima$entry), max(hiroshima$entry),
  min(hiroshima$time), max(hiroshima$time),
  file.size("data/hiroshima.rda") / 1024^2
))
}
