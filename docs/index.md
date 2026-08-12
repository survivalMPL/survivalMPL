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

    install.packages("survivalMPL")

## Quick start

### Right-censored example (`survival::lung`)

[`library`](https://rdrr.io/r/base/library.html)`(`[`survivalMPL`](https://CRAN.R-project.org/package=survivalMPL)`)`` `[`library`](https://rdrr.io/r/base/library.html)`(`[`survival`](https://github.com/therneau/survival)`)`` `` ``lung`` ``<-`` `[`na.omit`](https://rdrr.io/r/stats/na.fail.html)`(``lung``[``, `[`c`](https://rdrr.io/r/base/c.html)`(``"time"``, ``"status"``, ``"age"``, ``"sex"``, ``"ph.karno"``, ``"wt.loss"``)``]``)`` `` ``fit_lung`` ``<-`` `[`coxph_mpl`](https://CRAN.R-project.org/package=survivalMPL/reference/coxph_mpl.md)`(`` `` `[`Surv`](https://rdrr.io/pkg/survival/man/Surv.html)`(``time``, ``status`` ``==`` ``2``)`` ``~`` ``age`` ``+`` ``sex`` ``+`` ``ph.karno`` ``+`` ``wt.loss``,`` `` data ``=`` ``lung``,`` `` control ``=`` `[`coxph_mpl.control`](https://CRAN.R-project.org/package=survivalMPL/reference/coxph_mpl.control.md)`(`` `` n.obs ``=`` `[`sum`](https://rdrr.io/r/base/sum.html)`(``lung``$``status`` ``==`` ``2``)``,`` `` max.iter ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``40``, ``2000``, ``4000``)``,`` `` smooth ``=`` ``0`` `` ``)`` ``)`` `` `[`summary`](https://rdrr.io/r/base/summary.html)`(``fit_lung``)`` `[`plot`](https://rdrr.io/r/graphics/plot.default.html)`(`[`predict`](https://rdrr.io/r/stats/predict.html)`(``fit_lung``, type ``=`` ``"survival"``)``)`

### Interval-censored example (`bcos2`)

[`data`](https://rdrr.io/r/utils/data.html)`(``bcos2``)`` `` ``fit_bcos`` ``<-`` `[`coxph_mpl`](https://CRAN.R-project.org/package=survivalMPL/reference/coxph_mpl.md)`(`` `` `[`Surv`](https://rdrr.io/pkg/survival/man/Surv.html)`(``left``, ``right``, type ``=`` ``"interval2"``)`` ``~`` ``treatment``,`` `` data ``=`` ``bcos2``,`` `` control ``=`` `[`coxph_mpl.control`](https://CRAN.R-project.org/package=survivalMPL/reference/coxph_mpl.control.md)`(`` `` basis ``=`` ``"msplines"``,`` `` n.obs ``=`` `[`nrow`](https://rdrr.io/r/base/nrow.html)`(``bcos2``)``,`` `` max.iter ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``40``, ``2000``, ``4000``)``,`` `` smooth ``=`` ``0`` `` ``)`` ``)`` `` `[`summary`](https://rdrr.io/r/base/summary.html)`(``fit_bcos``)`` `[`plot`](https://rdrr.io/r/graphics/plot.default.html)`(`[`predict`](https://rdrr.io/r/stats/predict.html)`(``fit_bcos``, type ``=`` ``"survival"``, i ``=`` ``1``)``)`

## Resources

- Vignette: `vignettes/getting-started.Rmd`
- Reference documentation:
  [`?coxph_mpl`](https://CRAN.R-project.org/package=survivalMPL/reference/coxph_mpl.md),
  [`?coxph_mpl.control`](https://CRAN.R-project.org/package=survivalMPL/reference/coxph_mpl.control.md),
  [`?predict.coxph_mpl`](https://CRAN.R-project.org/package=survivalMPL/reference/predict.coxph_mpl.md),
  [`?residuals.coxph_mpl`](https://CRAN.R-project.org/package=survivalMPL/reference/residuals.coxph_mpl.md)
