# Getting Started with survivalMPL

## Overview

`survivalMPL` fits Cox proportional hazards models by **maximum
penalised likelihood** (MPL). Unlike partial likelihood, MPL estimates
the baseline hazard $`h_0`$ jointly with the regression coefficients
$`\boldsymbol{\beta}`$, by expanding it in a basis of non-negative
functions

``` math
h_0(t) = \sum_{u=1}^{m} \theta_u \psi_u(t), \qquad \boldsymbol{\theta} \geq \mathbf{0},
```

and penalising the roughness of the result. Because $`h_0`$ is estimated
rather than profiled out, fits yield *absolute* hazard, cumulative
hazard and survival estimates, and the method extends naturally to
censoring schemes the partial likelihood cannot handle.

Supported data:

- right censoring, including left truncation (delayed entry)
- left censoring
- interval censoring
- any mixture of the above, together with exact event times

## Installation

From CRAN:

[`install.packages`](https://rdrr.io/r/utils/install.packages.html)`(``"survivalMPL"``)`

The development version:

`# install.packages("remotes")`` ``remotes``::``install_github``(``"survivalMPL/survivalMPL"``)`

Then attach it alongside `survival`, which provides the
[`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) response
constructor:

[`library`](https://rdrr.io/r/base/library.html)`(`[`survivalMPL`](https://CRAN.R-project.org/package=survivalMPL)`)`` ``#> Loading required package: survival`` ``#> Loading required package: MASS`` `[`library`](https://rdrr.io/r/base/library.html)`(`[`survival`](https://github.com/therneau/survival)`)`

## The main function

[`coxph_mpl()`](https://CRAN.R-project.org/package=survivalMPL/reference/coxph_mpl.md)
is the entry point. Its interface follows
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html):

[`coxph_mpl`](https://CRAN.R-project.org/package=survivalMPL/reference/coxph_mpl.md)`(``formula``, ``data``, ``subset``, ``na.action``, ``entry``, ``control``, ``...``)`

| Argument | Purpose |
|----|----|
| `formula` | A model formula whose response is a [`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) object |
| `data` | Data frame containing the variables in `formula` |
| `entry` | Optional left-truncation (delayed entry) times |
| `control` | A list from [`coxph_mpl.control()`](https://CRAN.R-project.org/package=survivalMPL/reference/coxph_mpl.control.md) |

The response tells
[`coxph_mpl()`](https://CRAN.R-project.org/package=survivalMPL/reference/coxph_mpl.md)
which censoring scheme applies. Right-censored data can use the familiar
two-argument form; everything else uses `type = "interval2"`, where `NA`
marks an unobserved endpoint:

| Data              | Response                             |
|-------------------|--------------------------------------|
| Right censored    | `Surv(time, status)`                 |
| Left censored     | `Surv(NA, t_R, type = "interval2")`  |
| Interval censored | `Surv(t_L, t_R, type = "interval2")` |
| Exact event       | `Surv(t, t, type = "interval2")`     |

### Tuning the fit

[`coxph_mpl.control()`](https://CRAN.R-project.org/package=survivalMPL/reference/coxph_mpl.control.md)
collects everything governing the baseline hazard and the optimiser. The
arguments most often changed are:

| Argument | Meaning |
|----|----|
| `basis` | Basis for $`h_0`$: `"uniform"`, `"gaussian"`, `"msplines"`, `"epanechikov"` |
| `n.knots` | Length-2 vector: numbers of quantile and equally spaced internal knots |
| `smooth` | Smoothing parameter $`\lambda`$; `NULL` estimates it by REML |
| `n.obs` | Number of exact events, used to set default knot counts |
| `max.iter` | Iteration limits for the three optimisation stages |

The registered bases are:

[`list_bases`](https://CRAN.R-project.org/package=survivalMPL/reference/list_bases.md)`(``)`` ``#> name label`` ``#> 1 bsplines B-Splines`` ``#> 2 epanechikov Epanechnikov`` ``#> 3 gaussian Gaussian`` ``#> 4 msplines M-Splines`` ``#> 5 uniform Uniform`

### Methods

A fitted `coxph_mpl` object supports the usual generics:

| Call | Returns |
|----|----|
| `summary(fit)` | Coefficients, standard errors, and $`\boldsymbol{\theta}`$ |
| `coef(fit)` | Regression coefficients |
| `predict(fit, type = "risk")` | Estimated hazard $`\hat h(t)`$ |
| `predict(fit, type = "survival")` | Estimated survival $`\hat S(t)`$ |
| `plot(fit)` | Baseline hazard and survival with confidence bands |
| `residuals(fit)` | Martingale residuals |

## Where to go next

The tutorials work through one censoring scheme at a time — **right
censoring**, **left censoring**, and **interval censoring** — each with
a complete worked example. Two further articles cover the modelling
choices: **basis functions for the baseline hazard**, and a side-by-side
comparison of **`coxph` and `coxph_mpl`** on the same data.

See
[`?coxph_mpl`](https://CRAN.R-project.org/package=survivalMPL/reference/coxph_mpl.md)
and
[`?coxph_mpl.control`](https://CRAN.R-project.org/package=survivalMPL/reference/coxph_mpl.control.md)
for the full argument lists.
