# survivalMPL

Penalised maximum likelihood fits for Cox proportional hazards models
with right, left, and interval censoring. Baseline hazards are smoothed
via flexible bases (step, Gaussian, M-splines, Epanechnikov) with
non-negativity enforced during optimisation.

## Features

- Cox PH models with right, left, and interval censoring
- Smoothed baseline hazard with multiple basis choices
- Penalised likelihood with REML or fixed smoothing
- Inference for coefficients, baseline hazard/survival, predictions, and
  residuals

## Installation

[`install.packages`](https://rdrr.io/r/utils/install.packages.html)`(``"survivalMPL"``)`

## Quick start

A right-censored fit on
[`survival::lung`](https://rdrr.io/pkg/survival/man/lung.html):

[`library`](https://rdrr.io/r/base/library.html)`(`[`survivalMPL`](https://CRAN.R-project.org/package=survivalMPL)`)`` `[`library`](https://rdrr.io/r/base/library.html)`(`[`survival`](https://github.com/therneau/survival)`)`` `` ``lung`` ``<-`` `[`na.omit`](https://rdrr.io/r/stats/na.fail.html)`(``lung``[``, `[`c`](https://rdrr.io/r/base/c.html)`(``"time"``, ``"status"``, ``"age"``, ``"sex"``, ``"ph.karno"``, ``"wt.loss"``)``]``)`` `` ``fit_lung`` ``<-`` `[`coxph_mpl`](https://CRAN.R-project.org/package=survivalMPL/reference/coxph_mpl.md)`(`` `` `[`Surv`](https://rdrr.io/pkg/survival/man/Surv.html)`(``time``, ``status`` ``==`` ``2``)`` ``~`` ``age`` ``+`` ``sex`` ``+`` ``ph.karno`` ``+`` ``wt.loss``,`` `` data ``=`` ``lung``,`` `` tol ``=`` ``1e-05`` ``)`` `` `[`summary`](https://rdrr.io/r/base/summary.html)`(``fit_lung``)`` `[`plot`](https://rdrr.io/r/graphics/plot.default.html)`(`[`predict`](https://rdrr.io/r/stats/predict.html)`(``fit_lung``, type ``=`` ``"survival"``)``)`

Because
[`coxph_mpl()`](https://CRAN.R-project.org/package=survivalMPL/reference/coxph_mpl.md)
estimates the baseline hazard rather than profiling it out,
[`predict()`](https://rdrr.io/r/stats/predict.html) returns absolute
survival and hazard estimates with standard errors.

See the [Getting
Started](https://CRAN.R-project.org/package=survivalMPL/articles/getting-started.md)
article for installation notes, the
[`coxph_mpl()`](https://CRAN.R-project.org/package=survivalMPL/reference/coxph_mpl.md)
interface, how each censoring scheme is encoded, and the control
arguments. The remaining articles work through one censoring scheme at a
time and compare the modelling choices.

## Resources

- Reference documentation:
  [`?coxph_mpl`](https://CRAN.R-project.org/package=survivalMPL/reference/coxph_mpl.md),
  [`?coxph_mpl.control`](https://CRAN.R-project.org/package=survivalMPL/reference/coxph_mpl.control.md),
  [`?predict.coxph_mpl`](https://CRAN.R-project.org/package=survivalMPL/reference/predict.coxph_mpl.md),
  [`?residuals.coxph_mpl`](https://CRAN.R-project.org/package=survivalMPL/reference/residuals.coxph_mpl.md)
