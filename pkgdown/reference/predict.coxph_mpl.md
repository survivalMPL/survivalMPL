# Predictions for a Cox Model Fit via MPL

Compute predicted instantaneous risk or survival probabilities for a
fitted `coxph_mpl` model.

## Usage

``` r
# S3 method for class 'coxph_mpl'
predict(
  object,
  se = "M2QM2",
  type = "risk",
  i = NULL,
  time = NULL,
  upper.quantile = 0.95,
  ...
)

# S3 method for class 'predict.coxph_mpl'
plot(x, ...)
```

## Arguments

- object:

  A fitted model of class `"coxph_mpl"`.

- se:

  Inference method for confidence intervals. One of `"H"`, `"M2QM2"`, or
  `"M2HM2"`. Default `"M2QM2"`.

- type:

  Prediction type: `"risk"` for instantaneous risk or `"survival"` for
  survival probability. Default `"risk"`.

- i:

  Optional integer index of the observation whose covariates are used.
  Defaults to mean covariates.

- time:

  Optional numeric vector of times at which to predict. Defaults to 1000
  equally spaced times over the outcome range.

- upper.quantile:

  Quantile of the response used to bound the x-axis when plotting.
  Default `0.95`.

- ...:

  Additional plotting parameters.

- x:

  An object of class `"predict.coxph_mpl"`.

## Value

A data frame of class `"predict.coxph_mpl"` with columns `time`, `risk`
or `survival`, `se`, `low`, and `high`.

## Details

Predictions incorporate the baseline hazard or cumulative baseline
hazard, giving absolute (not relative) risk and survival estimates.
Standard errors and confidence intervals are computed via the delta
method and truncated to the parameter range.

## See also

\[coxph_mpl()\], \[coxph_mpl.control()\], \[residuals.coxph_mpl()\],
\[summary.coxph_mpl()\]

## Examples

``` r
if (FALSE) { # \dontrun{
data(lung)
fit_mpl <- coxph_mpl(Surv(time, status == 2) ~ age + sex + ph.karno + wt.loss,
                     data = lung)
plot(predict(fit_mpl))
} # }
```
