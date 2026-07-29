# Demonstration: left truncation (entry=) on the real hiroshima dataset
# (RERF Life Span Study Report 14; see R/data-hiroshima.R for provenance and
# the required RERF/MHLW/DOE acknowledgment).
#
# Each row is a pseudo-subject: entry = age at exposure, time = attained age,
# status = any death in that (grouped-table) stratum, dose = colon dose.

devtools::load_all()
data(hiroshima)

fit_entry <- coxph_mpl(
  Surv(time, status) ~ dose + sex + city,
  data   = hiroshima,
  entry  = entry,
  basis  = "msplines",
  # smooth = 0
)

fit_naive <- coxph_mpl(
  Surv(time, status) ~ dose + sex + city,
  data   = hiroshima,
  # basis  = "msplines",
  # smooth = 0
)

args(coxph_mpl)
print(coef(fit_entry, "Beta"))

print(coef(fit_naive, "Beta"))
 
# NOTE: unlike the synthetic demo (demo_left_truncation.R), the two fits are
# fairly close here because age-at-exposure isn't strongly correlated with
# dose/sex/city in this data - the point of this script is to show entry=
# running end-to-end on a real dataset, not to reproduce the informative-
# truncation bias the synthetic demo illustrates. Also note this is NOT a
# reanalysis of LSS Report 14's actual excess-relative-risk model (which uses
# grouped-data Poisson regression on person-years, not per-stratum pseudo-
# subjects) - treat the coefficients here as illustrative only.
