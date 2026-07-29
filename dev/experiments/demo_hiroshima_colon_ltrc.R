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
# Column mapping vs the colleague's script: agex -> entry, age -> time,
# colon -> status_colon, colon10 -> dose (both raw and standardized versions
# shown), city/gd3 retained as covariates (sex dropped, as the colleague
# found it non-significant and removed it from their final model).

devtools::load_all()
library(survival)

data(hiroshima)
hiroshima$dose_stan <- as.numeric(scale(hiroshima$dose))


# coxph_mpl WITH entry= (correct handling of left truncation).

tic()
fit_trunc <- coxph_mpl(
  Surv(time, status_colon) ~ city + gd3 + dose_stan,
  data   = hiroshima,
  entry  = entry,
  # basis  = "msplines",
  # smooth = 0,
)
toc()

tic()
fit_trunc_0 <- coxph_mpl(
  Surv(time, status_colon) ~ city + gd3 + dose_stan,
  data   = hiroshima,
  entry  = entry,
  basis  = "msplines",
  smooth = 5,
)
toc()

# coxph_mpl WITHOUT entry= (naive, ignores truncation).
fit_naive <- coxph_mpl(
  Surv(time, status_colon) ~ city + gd3 + dose_stan,
  data   = hiroshima,
  # basis  = "msplines",
  # smooth = 0,
)

coef(fit_trunc_0)
coef(fit_trunc)
coef(fit_naive)
