library(survivalMPL)

# Verification: coxph_mpl(entry=...) against the reference LTRC implementation
# in dev/experiments/old_sources/LTRC_codes/. Confirms the which=2 basis-matrix
# differencing reproduces the old piecewise-hazard profile-likelihood estimator
# under a matched configuration (basis="uniform", smooth=0, matched knots).

source(testthat::test_path("../../dev/experiments/old_sources/LTRC_codes/dataGenLeftTruncation.R"))
source(testthat::test_path("../../dev/experiments/old_sources/LTRC_codes/LT_RC_optimization_updated.R"))

test_that("entry= reproduces the reference LTRC piecewise estimator", {
  set.seed(1)
  # NOTE on parameter choice: datagen_LT_RC()'s Surv(right-type) -> interval2
  # conversion (R/coxph.r) sets a placeholder time2 = 1 for every non-interval
  # censored row (a quirk of survival::Surv(), not of this package). The
  # simulated event/censoring times must span past 1 or basis = "uniform"
  # errors with "subscript out of bounds" when it tries to place that
  # placeholder outside the knot range. alpha = 1, c_max = 3 gives an event
  # time range of roughly [0, 3], safely covering the placeholder.
  sim <- datagen_LT_RC(
    n = 500, beta = c(1, -1), dist = "weibull",
    c_max = 3, alpha = 1, psi = 0.7
  )

  df <- data.frame(
    y = sim$Y, x1 = sim$X[, 1], x2 = sim$X[, 2],
    l = sim$L, event = sim$del
  )

  ref <- Cox_LT_RC(
    Y = sim$Y, X = sim$X, L = sim$L, delVec = sim$del,
    maxIter = 1000, tol = 1e-6, bin_option = "quantile"
  )

  fit <- coxph_mpl(
    Surv(y, event) ~ x1 + x2,
    data    = df,
    entry   = l,
    basis   = "uniform",
    smooth  = 0,
    control = coxph_mpl.control(
      n.obs = sum(sim$del),
      n.knots = c(length(ref$bin_edges) - 2, 0)
    )
  )

  # Tolerance is intentionally loose (not exact-equality): coxph_mpl()'s
  # compute_knots() places quantile knots from a different event sample
  # (event/censoring/entry times combined) than the reference script's
  # bin_edges (quantiles of observed event times only, floored at min(L)), so
  # the two estimators are not expected to match exactly even though both
  # target the same likelihood.
  expect_equal(unname(coef(fit, "Beta")), as.vector(ref$bVal), tolerance = 0.1)

  # The entry= correction should move the estimate closer to the reference
  # than ignoring truncation entirely (fitting the same model without entry=).
  fit_ignoring_entry <- coxph_mpl(
    Surv(y, event) ~ x1 + x2,
    data    = df,
    basis   = "uniform",
    smooth  = 0,
    control = coxph_mpl.control(
      n.obs = sum(sim$del),
      n.knots = c(length(ref$bin_edges) - 2, 0)
    )
  )
  dist_with_entry <- sum((unname(coef(fit, "Beta")) - as.vector(ref$bVal))^2)
  dist_ignoring_entry <- sum((unname(coef(fit_ignoring_entry, "Beta")) - as.vector(ref$bVal))^2)
  expect_lt(dist_with_entry, dist_ignoring_entry)
})
