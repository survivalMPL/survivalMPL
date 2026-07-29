# Demonstration: run the LTRC_codes example data generator through the
# updated coxph_mpl(entry = ...) instead of the old piecewise Cox_LT_RC().
#
# See dev/experiments/old_sources/LTRC_codes/runPiecewiseSimul.R for the
# original example this is based on (same n, beta, dist).
#
# NOTE: the original example's c_max = 0.95 keeps every observed time under
# 1, which trips an unrelated, pre-existing quirk: survival::Surv()'s
# right -> interval2 conversion sets a placeholder time2 = 1 for every
# right-censored/event row (R/coxph.r), and coxph_mpl() evaluates its basis
# on that placeholder too, before subsetting it away. If the knot range
# doesn't reach 1, that raises an error - unrelated to entry=. c_max is
# bumped here so real times span past 1 and the placeholder falls safely
# inside the basis support; nothing else about the original example changed.

devtools::load_all()

source("dev/experiments/old_sources/LTRC_codes/dataGenLeftTruncation.R")

set.seed(1)
sim <- datagen_LT_RC(
  n = 1000, beta = c(1, -1), dist = "weibull",
  c_max = 1.6, alpha = 3, psi = 1
)


df <- data.frame(
  y     = sim$Y,
  x1    = sim$X[, 1],
  x2    = sim$X[, 2],
  entry = sim$L,
  event = sim$del
)

fit <- coxph_mpl(
  Surv(y, event) ~ x1 + x2,
  data   = df,
  entry  = entry,
  basis  = "msplines",
  smooth = 0
)

cat("True beta:            ", c(1, -1), "\n")
cat("coxph_mpl(entry = L): ", coef(fit, "Beta"), "\n")
