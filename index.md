# survivalMPL

Penalised maximum likelihood fits for Cox proportional hazards models with right, left, and interval censoring. Baseline hazards are smoothed via flexible bases (step, Gaussian, M-splines, Epanechnikov) with non-negativity enforced during optimisation.

## Features
- Cox PH models with right, left, and interval censoring
- Smoothed baseline hazard with multiple basis choices
- Penalised likelihood with REML or fixed smoothing
- Inference for coefficients, baseline hazard/survival, predictions, and residuals

## Installation
```
# CRAN (when available)
install.packages("survivalMPL")
```

## Quick start

### Right-censored example (`survival::lung`)
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

### Interval-censored example (`bcos2`)
```r
data(bcos2)

fit_bcos <- coxph_mpl(
  Surv(left, right, type = "interval2") ~ treatment,
  data = bcos2,
  control = coxph_mpl.control(
    basis = "msplines",
    n.obs = nrow(bcos2),
    max.iter = c(40, 2000, 4000),
    smooth = 0
  )
)

summary(fit_bcos)
plot(predict(fit_bcos, type = "survival", i = 1))
```

## Resources
- Vignette: `vignettes/getting-started.Rmd`
- Reference documentation: `?coxph_mpl`, `?coxph_mpl.control`, `?predict.coxph_mpl`, `?residuals.coxph_mpl`
