# Does the `hiroshima` reconstruction support the conclusions drawn from it?
#
# Answers the three open questions in
# dev/experiments/hiroshima_dataset/UNGROUPING.md:
#
#   Q2. Is "ungrouping = pyr/subjects" the same idea as what we implemented?
#   Q3. Are the coefficients merely illustrative, or do they estimate
#       something real? -> benchmark against the grouped-data Poisson fit,
#       which is the CORRECT analysis of this table, plus a sensitivity sweep
#       over every arbitrary choice the reconstruction makes.
#
# (Q1, what Ozasa et al. actually did, is answered from their own analysis
# script, dev/experiments/hiroshima_dataset/lss14/lss14.scr - no code needed.)

devtools::load_all()
library(survival)

RAW <- "dev/experiments/hiroshima_dataset/lss14/lss14.csv"
raw <- read.csv(RAW)
COVS <- c("dose", "sex", "city")

# Build a reconstruction in memory under a given configuration, by sourcing the
# canonical builder rather than reimplementing it.
build <- function(seed = 20260811, placement = "uniform", dose = "cell",
                  yrs_start = 5.15, yrs_end = 58.40) {
  e <- new.env()
  assign("HIROSHIMA_SEED", seed, envir = e)
  assign("HIROSHIMA_PLACEMENT", placement, envir = e)
  assign("HIROSHIMA_DOSE", dose, envir = e)
  assign("HIROSHIMA_YRS_START", yrs_start, envir = e)
  assign("HIROSHIMA_YRS_END", yrs_end, envir = e)
  assign("HIROSHIMA_SAVE", FALSE, envir = e)
  suppressMessages(sys.source("dev/data-raw/hiroshima.R", envir = e))
  get("hiroshima", envir = e)
}

cox_ltrc <- function(data) {
  coef(coxph(Surv(entry, time, status) ~ dose + sex + city, data = data))
}
cox_naive <- function(data) {
  coef(coxph(Surv(time, status) ~ dose + sex + city, data = data))
}

## =========================================================================
## Q2. "ungrouping is basically pyr/subjects"
## =========================================================================
cat("\n########## Q2: what is pyr/subjects? ##########\n\n")

cat("(a) Per cell, the ratio is only defined where subjects > 0, and subjects\n")
cat("    is nonzero ONLY in ctime == 1 (the entry census):\n")
nz <- tapply(raw$subjects > 0, raw$ctime, sum)
cat("    cells with subjects > 0, by ctime:", paste(nz, collapse = " "), "\n")
cat("    -> the ratio is unavailable for", sum(raw$ctime != 1), "of", nrow(raw), "cells.\n\n")

e1 <- raw[raw$ctime == 1 & raw$subjects > 0, ]
r1 <- e1$pyr / e1$subjects
cat("(b) Even where it IS defined, the per-cell ratio is not interpretable:\n")
cat(sprintf(
  "    pyr/subjects over ctime==1 cells: median %.2f, mean %.2f, range %.2f-%.2f\n",
  median(r1), mean(r1), min(r1), max(r1)
))
cat("    A ratio of 537 is impossible as a duration inside a 5-year band. The\n")
cat("    reason: a subject is counted in the cell they ENTER, but ages out of\n")
cat("    that attained-age band during the calendar band, contributing person-\n")
cat("    years to cells where they were never counted as a subject. So the\n")
cat("    numerator and denominator of the per-cell ratio refer to different\n")
cat("    sets of people.\n")
cat(sprintf(
  "    Aggregated over ctime == 1 it does behave: %.0f pyr / %.0f subjects\n",
  sum(raw$pyr[raw$ctime == 1]), sum(raw$subjects)
))
cat(sprintf(
  "    = %.2f years, i.e. the ~5-year width of one calendar band.\n\n",
  sum(raw$pyr[raw$ctime == 1]) / sum(raw$subjects)
))

cat("(c) Over the WHOLE table the ratio is the mean follow-up per subject,\n")
cat("    and that is exactly the quantity the reconstruction must reproduce:\n")
hiro <- build()
cat(sprintf(
  "    source      sum(pyr)/sum(subjects) = %.4f years\n", sum(raw$pyr) / sum(raw$subjects)
))
cat(sprintf(
  "    rebuilt     mean(time - entry)     = %.4f years   (%+.2f%%)\n",
  mean(hiro$time - hiro$entry),
  100 * (mean(hiro$time - hiro$entry) / (sum(raw$pyr) / sum(raw$subjects)) - 1)
))

cat("\n(d) The OTHER reading - n at risk in a cell = pyr / band width - is the\n")
cat("    piecewise-exponential episode expansion. It is a different object:\n")
w <- ifelse(raw$agecat == 21, 15, 5)
cat(sprintf(
  "    episodes implied by sum(pyr/width) = %s  vs  %s subjects\n",
  format(round(sum(raw$pyr / w)), big.mark = ","), format(sum(raw$subjects), big.mark = ",")
))
cat("    That expansion is lossless for the table but gives one row per\n")
cat("    person-INTERVAL, with entry at a band boundary - so individual-level\n")
cat("    left truncation is not represented. See UNGROUPING.md.\n")

## =========================================================================
## Q3a. Benchmark against the correct analysis of this table
## =========================================================================
cat("\n\n########## Q3a: vs grouped-data Poisson (the correct analysis) ##########\n\n")

# A Poisson/log-linear model for the cell death counts with log person-years as
# offset and a free effect per attained-age band IS the piecewise-exponential
# proportional-hazards likelihood for this table, with attained age as the time
# scale. Its dose/sex/city coefficients are log hazard ratios, directly
# comparable to a Cox coefficient. It uses the full person-year weighting and
# needs no reconstruction, so it is the benchmark.
raw$dose <- raw$colon10 / 1000
raw$sex <- factor(raw$sex, levels = c(1, 2), labels = c("Male", "Female"))
raw$city <- factor(raw$city, levels = c(1, 2), labels = c("Hiroshima", "Nagasaki"))
pos <- raw[raw$pyr > 0, ]

p1 <- glm(death ~ dose + sex + city + factor(agecat),
  offset = log(pyr), family = poisson, data = pos
)
p2 <- glm(death ~ dose + sex + city + factor(agecat) + factor(agexcat),
  offset = log(pyr), family = poisson, data = pos
)

tab <- rbind(
  `Poisson, age bands` = coef(p1)[c("dose", "sexFemale", "cityNagasaki")],
  `Poisson, + agexcat` = coef(p2)[c("dose", "sexFemale", "cityNagasaki")],
  `reconstruction, LTRC Cox` = cox_ltrc(hiro),
  `reconstruction, naive Cox` = cox_naive(hiro)
)
print(round(tab, 4))
se_p1 <- summary(p1)$coefficients[c("dose", "sexFemale", "cityNagasaki"), 2]
cat("\nPoisson SEs:", paste(sprintf("%.4f", se_p1), collapse = "  "), "\n")
cat("\nLTRC Cox - Poisson, in Poisson SE units:\n")
print(round((tab[3, ] - tab[1, ]) / se_p1, 2))
cat("naive Cox - Poisson, in Poisson SE units:\n")
print(round((tab[4, ] - tab[1, ]) / se_p1, 2))

## =========================================================================
## Q3a-bis. WHY the Poisson fit is the right benchmark for a truncation
##          question, demonstrated rather than asserted.
## =========================================================================
# The claim: the grouped person-time table is left-truncation-correct BY
# CONSTRUCTION. `pyr` accumulates only time actually under observation, and
# follow-up starts in 1950, so a subject who was 30 in 1950 contributes zero
# person-years to every attained-age band below 30. The table simply cannot
# represent pre-1950 risk time. A naive Cox fit, by contrast, invents risk time
# from age 0 up to entry for every subject.
#
# If that is really the mechanism, then rebuilding the person-time table WITH
# the fictitious pre-entry time should move the Poisson fit onto the naive Cox
# estimate. Both Poisson fits and both Cox fits below use exactly the same
# subjects, so nothing but the person-time bookkeeping differs.
cat("\n\n########## Q3a-bis: is the Poisson benchmark truncation-correct? ##########\n\n")

set.seed(7)
sm <- hiro[sort(sample.int(nrow(hiro), 20000)), ]
cuts <- seq(0, 115, by = 5)

pois_from_episodes <- function(data, truncated) {
  # split each subject's follow-up at 5-year attained-age boundaries; entry is
  # either the real entry age (truncated) or 0 (fictitious pre-entry risk time)
  d <- data
  d$start <- if (truncated) d$entry else 0
  # name the output columns explicitly: survSplit otherwise inherits the names
  # used inside Surv(), which differ between the two calls
  sp <- survSplit(Surv(start, time, status) ~ .,
    data = d, cut = cuts, episode = "band",
    start = "tstart", end = "tstop", event = "ev"
  )
  sp$pyr <- sp$tstop - sp$tstart
  sp <- sp[sp$pyr > 0, ]
  fit <- glm(ev ~ dose + sex + city + factor(band),
    offset = log(pyr), family = poisson, data = sp
  )
  c(coef(fit)[c("dose", "sexFemale", "cityNagasaki")],
    episodes = nrow(sp), pyr = sum(sp$pyr))
}

pe_t <- pois_from_episodes(sm, TRUE)
pe_n <- pois_from_episodes(sm, FALSE)
cmp <- rbind(
  `Poisson, person-time from entry` = pe_t[1:3],
  `Cox LTRC (same subjects)` = cox_ltrc(sm),
  `Poisson, person-time from age 0` = pe_n[1:3],
  `Cox naive (same subjects)` = cox_naive(sm)
)
print(round(cmp, 4))
cat(sprintf(
  "\nperson-years: from entry %s over %s episodes\n              from age 0 %s over %s episodes (+%.0f%% fictitious)\n",
  format(round(pe_t["pyr"]), big.mark = ","), format(pe_t["episodes"], big.mark = ","),
  format(round(pe_n["pyr"]), big.mark = ","), format(pe_n["episodes"], big.mark = ","),
  100 * (pe_n["pyr"] / pe_t["pyr"] - 1)
))
cat("\nmax |Poisson - Cox| within each pair:\n")
cat(sprintf(
  "  truncated person-time %.4f\n  age-0 person-time     %.4f\n",
  max(abs(cmp[1, ] - cmp[2, ])), max(abs(cmp[3, ] - cmp[4, ]))
))
cat("\nIf the two pairs agree internally while differing from each other, the\n")
cat("naive-vs-LTRC gap is a property of the person-time bookkeeping - i.e. of\n")
cat("the truncation - and not of the reconstruction or of either estimator.\n")

## =========================================================================
## Q3c. Against Ozasa et al.'s OWN published estimates
## =========================================================================
# lss14.log records the fitted all-cause ERR per Gy from Ozasa et al.'s AMFIT
# script, with baseline stratified on city x sex x agexcat x agecat:
#   males   0.1500  (97.5% profile bounds 0.1035 - 0.1986)
#   females 0.2978  (97.5% profile bounds 0.2438 - 0.3539)
# Their model is  lambda = lambda0(strata) * (1 + beta*d)  so beta is an ERR
# per Gy, whereas a Cox coefficient b is a log hazard ratio; the comparable
# quantity is exp(b) - 1. Matching their stratification as closely as the
# reconstruction allows: attained age is the Cox time scale (= their agecat),
# and we stratify on city and on 5-year bands of age at exposure (= agexcat).
cat("\n\n########## Q3c: vs Ozasa et al.'s published all-cause ERR/Gy ##########\n\n")

published <- list(Male = c(0.1500, 0.1035, 0.1986), Female = c(0.2978, 0.2438, 0.3539))
hiro$agexband <- cut(hiro$agex, breaks = seq(0, 100, 5))

for (s in c("Male", "Female")) {
  sub <- hiro[hiro$sex == s, ]
  f <- coxph(Surv(entry, time, status) ~ dose + strata(city, agexband), data = sub)
  b <- coef(f)[["dose"]]
  ci <- confint(f)["dose", ]
  cat(sprintf(
    "%-7s  reconstruction ERR/Gy = exp(b)-1 = %.4f [%.4f, %.4f]   (b = %.4f)\n",
    s, exp(b) - 1, exp(ci[1]) - 1, exp(ci[2]) - 1, b
  ))
  p <- published[[s]]
  cat(sprintf(
    "%-7s  Ozasa et al. AMFIT   ERR/Gy = %.4f [%.4f, %.4f]   -> %s\n\n",
    "", p[1], p[2], p[3],
    if (exp(b) - 1 >= p[2] && exp(b) - 1 <= p[3]) {
      "inside their interval"
    } else {
      "OUTSIDE their interval"
    }
  ))
}
cat("Caveats on this comparison: their dose term is linear inside (1 + beta*d)\n")
cat("while ours is log-linear (the two agree only to first order in beta*d);\n")
cat("their baseline is stratified on 1200 cells vs our attained-age time scale\n")
cat("plus city and 5-year agex bands; and they use the person-year weighting\n")
cat("directly rather than a reconstruction.\n")

## =========================================================================
## Q3b. Sensitivity to every arbitrary choice in the reconstruction
## =========================================================================
cat("\n\n########## Q3b: sensitivity of the coefficients ##########\n\n")

variants <- list(
  `default (seed 20260811)` = list(),
  `seed 2` = list(seed = 2),
  `seed 3` = list(seed = 3),
  `seed 4` = list(seed = 4),
  `seed 5` = list(seed = 5),
  `placement = midpoint` = list(placement = "midpoint"),
  `placement = cellmean` = list(placement = "cellmean"),
  `dose = stratum mean` = list(dose = "stratum"),
  `entry offset 4.5y` = list(yrs_start = 4.5),
  `entry offset 5.8y` = list(yrs_start = 5.8),
  `end of f-up 56y` = list(yrs_end = 56.0),
  `end of f-up 60y` = list(yrs_end = 60.0)
)

rows <- lapply(names(variants), function(nm) {
  x <- do.call(build, variants[[nm]])
  c(ltrc = cox_ltrc(x), naive = cox_naive(x))
})
S <- do.call(rbind, rows)
rownames(S) <- names(variants)

pars <- c("dose", "sexFemale", "cityNagasaki")
cat("LTRC Cox coefficients:\n")
print(round(S[, paste0("ltrc.", pars)], 4))
cat("\nnaive - LTRC (the truncation effect) under each variant:\n")
print(round(S[, paste0("naive.", pars)] - S[, paste0("ltrc.", pars)], 4))

cat("\nspread across variants, vs the Poisson SE:\n")
for (p in pars) {
  v <- S[, paste0("ltrc.", p)]
  dd <- S[, paste0("naive.", p)] - S[, paste0("ltrc.", p)]
  cat(sprintf(
    "  %-13s coef range %.4f (%.1f SE)   |   truncation effect range %.4f, sign %s\n",
    p, diff(range(v)), diff(range(v)) / se_p1[p], diff(range(dd)),
    if (all(dd > 0)) "always +" else if (all(dd < 0)) "always -" else "MIXED"
  ))
}
