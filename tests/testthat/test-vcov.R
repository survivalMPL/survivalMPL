library(survivalMPL)

# vcov() / confint() for coxph_mpl objects.
#
# The fit already stores standard errors in fit$se, so these tests mostly
# check that vcov() exposes the same numbers from the covariance matrix the
# errors were derived from - if the two ever disagree, one of them is wrong.

data(bcos2)

fit_bc <- coxph_mpl(
  Surv(left, right, type = "interval2") ~ treatment,
  data   = bcos2,
  basis  = "uniform",
  smooth = 0
)

# ---------------------------------------------------------------------------
# vcov()
# ---------------------------------------------------------------------------

test_that("vcov returns the Beta block, named and symmetric", {
  V <- vcov(fit_bc)
  expect_true(is.matrix(V))
  expect_identical(dim(V), c(fit_bc$dim$p, fit_bc$dim$p))
  expect_identical(rownames(V), colnames(fit_bc$data$X))
  expect_identical(rownames(V), colnames(V))
  expect_equal(V, t(V))
})

test_that("vcov standard errors match the ones summary() reports", {
  expect_equal(
    unname(sqrt(diag(vcov(fit_bc)))),
    unname(coef(summary(fit_bc))[, "Std. Error"])
  )
})

test_that("the Theta block is m x m and matches fit$se", {
  V <- vcov(fit_bc, parameters = "Theta")
  expect_identical(dim(V), c(fit_bc$dim$m, fit_bc$dim$m))
  expect_equal(
    unname(sqrt(diag(V))),
    unname(fit_bc$se$Theta$M2QM2)
  )
})

test_that("se= selects among the stored covariance matrices", {
  expect_false(identical(vcov(fit_bc, se = "M2QM2"), vcov(fit_bc, se = "H")))
  for (s in c("M2QM2", "M2HM2", "H")) {
    expect_equal(
      unname(sqrt(diag(vcov(fit_bc, se = s)))),
      unname(fit_bc$se$Beta[[s]])
    )
  }
  expect_error(vcov(fit_bc, se = "nonsense"))
})

test_that("parameters= is matched, not taken on trust", {
  expect_error(vcov(fit_bc, parameters = "Gamma"))
})

# ---------------------------------------------------------------------------
# confint()
# ---------------------------------------------------------------------------

test_that("confint is a Wald interval built from coef and vcov", {
  ci <- confint(fit_bc)
  cf <- coef(fit_bc)
  se <- sqrt(diag(vcov(fit_bc)))
  expect_identical(dim(ci), c(length(cf), 2L))
  expect_identical(rownames(ci), names(cf))
  expect_equal(unname(ci[, 1]), unname(cf - qnorm(0.975) * se))
  expect_equal(unname(ci[, 2]), unname(cf + qnorm(0.975) * se))
})

test_that("level widens the interval and labels the columns", {
  expect_identical(colnames(confint(fit_bc)), c("2.5 %", "97.5 %"))
  expect_identical(colnames(confint(fit_bc, level = 0.9)), c("5.0 %", "95.0 %"))
  wide <- diff(as.numeric(confint(fit_bc, level = 0.99)))
  narrow <- diff(as.numeric(confint(fit_bc, level = 0.90)))
  expect_gt(wide, narrow)
})

test_that("parm subsets by name or index, and rejects unknown names", {
  nm <- names(coef(fit_bc))[1]
  expect_identical(rownames(confint(fit_bc, parm = nm)), nm)
  expect_equal(confint(fit_bc, parm = nm), confint(fit_bc, parm = 1))
  expect_error(confint(fit_bc, parm = "not_a_covariate"), "not found")
})

test_that("intervals are on the beta scale, so exp() gives the hazard ratio", {
  ci <- confint(fit_bc)
  hr <- exp(ci)
  expect_true(all(hr > 0))
  expect_lt(hr[1, 1], exp(coef(fit_bc))[1])
  expect_gt(hr[1, 2], exp(coef(fit_bc))[1])
})

# ---------------------------------------------------------------------------
# Documented degenerate case: with no covariates coxph_mpl() carries a phantom
# column of zeros, marks it as out of the parameter set, and ends up with no
# usable covariance for it.  vcov() must pass that through rather than invent a
# number or fail, and it must agree with the standard error the fit already
# stores.  Recorded here so that changing any of it is a deliberate decision.
# ---------------------------------------------------------------------------

test_that("an intercept-only fit degrades quietly rather than erroring", {
  fit0 <- coxph_mpl(
    Surv(left, right, type = "interval2") ~ 1,
    data   = bcos2,
    basis  = "uniform",
    smooth = 0
  )
  expect_identical(unname(coef(fit0)), 0)

  V <- vcov(fit0)
  expect_identical(dim(V), c(1L, 1L))
  expect_equal(unname(sqrt(diag(V))), unname(fit0$se$Beta$M2QM2))
  expect_identical(dim(confint(fit0)), c(1L, 2L))
})
