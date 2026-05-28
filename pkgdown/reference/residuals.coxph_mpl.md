# Residuals for a Cox Model Fit via MPL

Compute martingale and Cox-Snell residuals for a `coxph_mpl` model. The
returned object has a plot method.

## Usage

``` r
# S3 method for class 'coxph_mpl'
residuals(object, ...)

# S3 method for class 'residuals.coxph_mpl'
plot(x, ask = TRUE, which = 1:2, upper.quantile = 0.95, ...)
```

## Arguments

- object:

  A fitted model of class `"coxph_mpl"`.

- ...:

  Additional plotting parameters.

- x:

  An object of class `"residuals.coxph_mpl"`.

- ask:

  Logical; whether to prompt before each plot. Default `TRUE`.

- which:

  Integer vector selecting residual plots (`1:2`).

- upper.quantile:

  Quantile used to bound the y-axis for Cox-Snell residuals when
  `which == 3`. Default `0.95`.

## Value

A data frame of class `"residuals.coxph_mpl"` with columns `time1`,
`time2`, `censoring`, `coxsnell`, and `martingale`.

## References

Farrington (2000), Collett (2003), Moeschberger (2003).

## See also

\[coxph_mpl()\], \[predict.coxph_mpl()\], \[summary.coxph_mpl()\]

## Examples

``` r
if (FALSE) { # \dontrun{
data(lung, package = "survival")
fit_mpl <- coxph_mpl(Surv(time, status == 2) ~ age + sex + ph.karno + wt.loss,
                     data = lung)
par(mfrow = c(1, 2))
plot(residuals(fit_mpl), which = 1:2, ask = FALSE)
} # }
```
