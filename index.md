# survivalMPL

Penalised maximum likelihood fits for Cox proportional hazards models with right, left, and interval censoring. Baseline hazards are smoothed via flexible bases (step, Gaussian, M-splines, Epanechnikov) with non-negativity enforced during optimisation.

## Features
- Cox PH models with right, left, and interval censoring
- Smoothed baseline hazard with multiple basis choices
- Penalised likelihood with REML or fixed smoothing
- Inference for coefficients, baseline hazard/survival, predictions, and residuals

## Installation

```r
install.packages("survivalMPL")
```

## Quick start

A right-censored fit on `survival::lung`:

```r
library(survivalMPL)
library(survival)

lung <- na.omit(lung[, c("time", "status", "age", "sex", "ph.karno", "wt.loss")])

fit_lung <- coxph_mpl(
  Surv(time, status == 2) ~ age + sex + ph.karno + wt.loss,
  data = lung,
  control = coxph_mpl.control(
    n.obs = sum(lung$status == 2),
    max.iter = c(40, 2000, 4000),
    smooth = 0
  )
)

summary(fit_lung)
plot(predict(fit_lung, type = "survival"))
```

Because `coxph_mpl()` estimates the baseline hazard rather than profiling it out,
`predict()` returns absolute survival and hazard estimates with standard errors.

See the [Getting Started](articles/getting-started.html) article for
installation notes, the `coxph_mpl()` interface, how each censoring scheme is
encoded, and the control arguments. The remaining articles work through one
censoring scheme at a time and compare the modelling choices.

## Resources
- Reference documentation: `?coxph_mpl`, `?coxph_mpl.control`, `?predict.coxph_mpl`, `?residuals.coxph_mpl`
