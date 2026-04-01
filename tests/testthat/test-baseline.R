library(survivalMPL)

# Golden-value snapshot tests.
# Run devtools::test() ONCE on a clean state to CREATE the snapshots.
# Thereafter, any difference signals a regression.

test_that("bcos2 msplines: Beta coef stable", {
  data(bcos2)
  fit <- coxph_mpl(
    Surv(left, right, type = "interval2") ~ treatment,
    data   = bcos2,
    basis  = "msplines",
    smooth = 0
  )
  expect_snapshot(round(coef(fit, "Beta"), 6))
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
  expect_snapshot(round(coef(fit, "Beta"), 6))
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
  expect_snapshot(round(coef(fit, "Beta"), 6))
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
  expect_snapshot(round(coef(fit, "Beta"), 6))
  expect_snapshot(round(fit$ploglik, 4))
})
