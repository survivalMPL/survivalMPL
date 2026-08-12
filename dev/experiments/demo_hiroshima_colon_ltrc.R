# Demonstration: colon-cancer-specific left-truncated Cox regression on the
# hiroshima dataset, based on a colleague's MPLE-LTRC analysis of LSS14.
#
# The colleague's script fits a custom mpleltrc() function (source not
# available here - likely a sibling of dev/experiments/old_sources/
# LTRC_codes/Cox_LT_RC) and cross-checks against
# survival::coxph(Surv(agex, age, colon) ~ city + gd3 + colon10).
# Same cross-check strategy here, using coxph_mpl(entry = ...) in place of
# mpleltrc(), plus a naive (no entry=) fit to show the effect of ignoring
# truncation on this real, colon-cancer-specific outcome.
#
# Column mapping vs the colleague's script: agex -> entry (note: `entry` is the
# attained age at the 1950-10-01 start of follow-up, which is the correct
# risk-set entry age; the dataset's `agex` column is the age at exposure in
# 1945, i.e. entry - 5.15), age -> time, colon -> status_colon, colon10 ->
# dose (stored in Gy here, not mGy; a standardized version is also shown).
# city/gd3 retained as covariates, sex dropped as the colleague found it
# non-significant and removed it from their final model.
#
# entry= is only correct for basis = "uniform", so that is used throughout.

devtools::load_all()
library(survival)

data(hiroshima)
hiroshima$dose_stan <- as.numeric(scale(hiroshima$dose))

# WITH entry= (correct handling of left truncation)
t0 <- proc.time()[3]
fit_trunc <- coxph_mpl(Surv(time, status_colon) ~ city + gd3 + dose_stan,
  data = hiroshima, entry = entry,
  basis = "uniform", smooth = 0
)
cat(sprintf("entry= : %.1fs\n", proc.time()[3] - t0))

# WITHOUT entry= (naive, ignores truncation)
fit_naive <- coxph_mpl(Surv(time, status_colon) ~ city + gd3 + dose_stan,
  data = hiroshima, basis = "uniform", smooth = 0
)

print(round(rbind(
  `MPL entry=` = coef(fit_trunc, "Beta"),
  `coxph LTRC` = coef(coxph(
    Surv(entry, time, status_colon) ~ city + gd3 + dose_stan,
    data = hiroshima
  )),
  `MPL naive` = coef(fit_naive, "Beta"),
  `coxph naive` = coef(coxph(
    Surv(time, status_colon) ~ city + gd3 + dose_stan,
    data = hiroshima
  ))
), 4))

summary(fit_trunc)
