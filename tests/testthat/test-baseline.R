library(survivalMPL)

# Golden-value regression tests for the four basis types.
#
# Coefficients are checked against reference values with a relative tolerance
# rather than snapshotted. The gaussian coefficient is 0.9494945086, which sits
# 8.6e-9 above the rounding boundary of round(x, 6), so an exact snapshot of the
# rounded value flips between 0.949494 and 0.949495 across R versions and BLAS
# implementations while the fit itself is unchanged. A 1e-6 tolerance still
# catches any real regression.
#
# The penalised log-likelihood and dimension snapshots are kept as snapshots:
# run devtools::test() once on a clean state to create them, and thereafter any
# difference signals a regression.

test_that("bcos2 msplines: Beta coef stable", {
  data(bcos2)
  fit <- coxph_mpl(
    Surv(left, right, type = "interval2") ~ treatment,
    data   = bcos2,
    basis  = "msplines",
    smooth = 0
  )
  expect_equal(unname(coef(fit, "Beta")), 0.8731076, tolerance = 1e-6)
  expect_snapshot(round(fit$ploglik, 4))
  expect_snapshot(fit$dim[c("n", "n.obs", "p", "m")])
})

test_that("bcos2 uniform: Beta coef stable", {
  data(bcos2)
  fit <- coxph_mpl(
    Surv(left, right, type = "interval2") ~ treatment,
    data   = bcos2,
    basis  = "uniform",
    smooth = 0
  )
  expect_equal(unname(coef(fit, "Beta")), 0.9073718, tolerance = 1e-6)
  expect_snapshot(round(fit$ploglik, 4))
})

test_that("bcos2 gaussian: Beta coef stable", {
  data(bcos2)
  fit <- coxph_mpl(
    Surv(left, right, type = "interval2") ~ treatment,
    data   = bcos2,
    basis  = "gaussian",
    smooth = 0
  )
  expect_equal(unname(coef(fit, "Beta")), 0.9494945, tolerance = 1e-6)
  expect_snapshot(round(fit$ploglik, 4))
})

test_that("bcos2 epanechikov: Beta coef stable", {
  data(bcos2)
  fit <- coxph_mpl(
    Surv(left, right, type = "interval2") ~ treatment,
    data   = bcos2,
    basis  = "epanechikov",
    smooth = 0
  )
  expect_equal(unname(coef(fit, "Beta")), 0.7462842, tolerance = 1e-6)
  expect_snapshot(round(fit$ploglik, 4))
})
