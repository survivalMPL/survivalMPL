# Fit Cox Proportional Hazards Regression Model Via MPL

Simultaneously estimate the regression coefficients and the baseline
hazard function of proportional hazard Cox models using maximum
penalised likelihood (MPL).

## Usage

``` r
coxph_mpl(formula, data, subset, na.action, control, ...)

# S3 method for class 'coxph_mpl'
print(x, ...)
```

## Arguments

- formula:

  A survival formula with the response on the left-hand side and
  covariates on the right. The response must be built with
  \[survival::Surv()\] using `type = "interval2"` for interval-censored
  data (right-censored responses are converted internally).

- data:

  Optional data frame in which to evaluate `formula`.

- subset:

  Optional expression specifying a subset of observations to use.

- na.action:

  Optional missing-data filter function applied to the model frame;
  defaults to `options()\$na.action`.

- control:

  Optional list returned by \[coxph_mpl.control()\] specifying basis
  choice, smoothing value, iteration limits, and related options. When
  missing, defaults are built from `...`.

- ...:

  Additional arguments passed to \[coxph_mpl.control()\].

- x:

  An object of class `"coxph_mpl"`.

## Value

An object of class `"coxph_mpl"`; see \[coxph_mpl.object\] for
components.

Invisibly returns `x`.

## Details

`coxph_mpl` fits a Cox proportional hazards model allowing right, left,
and interval censoring by maximising a penalised likelihood in which a
penalty term smooths the baseline hazard estimate. Optimisation combines
a Newton step for regression coefficients with a multiplicative step for
the baseline hazard parameters while enforcing non-negativity
constraints (see Ma, Couturier, Heritier and Marschner (2021)). The
covariate matrix is centred during optimisation; baseline estimates and
covariance matrices are corrected afterwards using a delta-method
adjustment.

## See also

\[coxph_mpl.object()\], \[coxph_mpl.control()\],
\[summary.coxph_mpl()\], \[plot.coxph_mpl()\], \[predict.coxph_mpl()\]

## Examples

``` r
if (FALSE) { # \dontrun{
## Right-censored example: survival::lung
data(lung, package = "survival")
fit_mpl <- coxph_mpl(Surv(time, status == 2) ~ age + sex + ph.karno +
 wt.loss,
                     data = lung)
summary(fit_mpl)

## Interval-censored example: bcos2
data(bcos2)
fit_mpl <- coxph_mpl(Surv(left, right, type = "interval2") ~ treatment,
                     data = bcos2, basis = "m")
summary(fit_mpl)
} # }
```
