library(survivalMPL)

# Smoke tests — must pass before AND after every change.
# Numerical snapshots were taken from installed survivalMPL 0.2-4.
# Re-run devtools::test() to compare; any diff signals a regression.
#
# Coverage:
#   1. Right-censored case (lung, msplines + uniform)
#   2. S3 method structure (class, column names)
#   3. Interval-censored method smoke (melanoma — run-without-error only)

lung_clean <- local({
  d <- survival::lung
  d[complete.cases(d[, c("time", "status", "age", "sex", "ph.karno", "wt.loss")]), ]
})

# ---------------------------------------------------------------------------
# 1. Right-censored numerical stability
# ---------------------------------------------------------------------------

test_that("lung msplines: Beta and ploglik stable", {
  fit <- coxph_mpl(
    Surv(time, status == 2) ~ age + sex + ph.karno + wt.loss,
    data   = lung_clean,
    basis  = "msplines",
    smooth = 0
  )
  expect_snapshot(round(coef(fit, "Beta"), 6))
  expect_snapshot(round(fit$ploglik[1], 4))
  expect_snapshot(fit$dim[c("n", "p", "m")])
})

test_that("lung uniform: Beta and ploglik stable", {
  fit <- coxph_mpl(
    Surv(time, status == 2) ~ age + sex + ph.karno + wt.loss,
    data   = lung_clean,
    basis  = "uniform",
    smooth = 0
  )
  expect_snapshot(round(coef(fit, "Beta"), 6))
  expect_snapshot(round(fit$ploglik[1], 4))
})

# ---------------------------------------------------------------------------
# 2. S3 method output structure (class + column names)
# ---------------------------------------------------------------------------

test_that("S3 methods return correct classes and structure", {
  fit <- coxph_mpl(
    Surv(time, status == 2) ~ age + sex + ph.karno + wt.loss,
    data   = lung_clean,
    basis  = "msplines",
    smooth = 0
  )

  # summary
  s <- summary(fit)
  expect_s3_class(s, "summary.coxph_mpl")
  expect_named(s, c("Beta", "Theta", "inf"))

  # predict
  p <- predict(fit)
  expect_s3_class(p, "predict.coxph_mpl")
  expect_named(p, c("time", "risk", "se", "low", "high"))

  # predict survival
  ps <- predict(fit, type = "survival")
  expect_named(ps, c("time", "survival", "se", "low", "high"))

  # residuals
  r <- residuals(fit)
  expect_s3_class(r, "residuals.coxph_mpl")
  expect_named(r, c("time1", "time2", "censoring", "coxsnell", "martingale"))

  # coef
  b <- coef(fit, "Beta")
  expect_named(b, c("age", "sex", "ph.karno", "wt.loss"))
  expect_named(coef(fit, "Theta"), as.character(seq_len(fit$dim$m)))
})

test_that("print and plot run without error", {
  fit <- coxph_mpl(
    Surv(time, status == 2) ~ age + sex + ph.karno + wt.loss,
    data   = lung_clean,
    basis  = "msplines",
    smooth = 0
  )
  expect_invisible(print(fit))
  expect_invisible(print(summary(fit)))
  withr::with_tempfile("f", fileext = ".png", code = {
    grDevices::png(f)
    expect_no_error(plot(fit, ask = FALSE))
    expect_no_error(plot(residuals(fit), ask = FALSE))
    expect_no_error(plot(predict(fit)))
    grDevices::dev.off()
  })
})

# ---------------------------------------------------------------------------
# 3. Interval-censored smoke (melanoma) — run-without-error only
# ---------------------------------------------------------------------------

test_that("melanoma msplines interval-censored: runs without error", {
  data(melanoma)
  expect_no_error({
    fit <- coxph_mpl(
      Surv(t_L, t_R, type = "interval2") ~
        Arm + Leg + Trunk + mm1to2 + mm2to4 + mm4plus + Female + Age_centred,
      data  = melanoma,
      basis = "msplines",
      smooth = 0
    )
    summary(fit)
  })
})
