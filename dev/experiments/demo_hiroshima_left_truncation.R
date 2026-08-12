# Demonstration: left truncation (entry=) on the real hiroshima dataset
# (RERF Life Span Study Report 14; see R/data-hiroshima.R for provenance and
# the required RERF/MHLW/DOE acknowledgment).
#
# One row per subject. entry = attained age on 1950-10-01, when LSS follow-up
# starts, so entry ages run from 5 to 94 years.
#
# entry= is only correct for basis = "uniform" (piecewise-constant baseline).

devtools::load_all()
data(hiroshima)

fit_entry <- coxph_mpl(
  Surv(time, status) ~ dose + sex + city,
  data   = hiroshima,
  entry  = entry,
  basis  = "mspline"
)

fit_naive <- coxph_mpl(
  Surv(time, status) ~ dose + sex + city,
  data   = hiroshima,
  basis  = "uniform"
)

print(coef(fit_entry, "Beta"))
print(coef(fit_naive, "Beta"))

# NOTE: ignoring truncation inflates dose by ~1.3 SE and city by ~2.1 SE, both
# upwards; sex is unaffected, entry age being unrelated to sex. That the shift
# is systematic rather than noise, and that entry= moves towards the correct
# grouped-data analysis rather than away from it, is established in
# dev/experiments/hiroshima_reconstruction_checks.R.
