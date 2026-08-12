# Extract Coefficients from `coxph_mpl` Fits

Extract Coefficients from `coxph_mpl` Fits

## Usage

``` r
# S3 method for class 'summary.coxph_mpl'
coef(object, parameters = "Beta", ...)

# S3 method for class 'coxph_mpl'
coef(object, parameters = "Beta", ...)
```

## Arguments

- object:

  An object of class `"coxph_mpl"` or `"summary.coxph_mpl"`.

- parameters:

  Parameter set of interest: `"Beta"` for regression coefficients or
  `"Theta"` for baseline hazard parameters. Default `"Beta"`.

- ...:

  Additional arguments passed to methods.

## Value

A vector of coefficients or a matrix with estimates, standard errors,
z-statistics, and p-values.

## Details

For `summary.coxph_mpl` inputs with `parameters == "Theta"`, only
baseline hazard estimates exceeding `min.theta` are reported.

## See also

\[coxph_mpl()\], \[summary.coxph_mpl()\]

## Examples

``` r
if (FALSE) { # \dontrun{
data(lung, package = "survival")
fit_mpl <- coxph_mpl(Surv(time, status == 2) ~ age + sex + ph.karno + wt.loss,
                     data = lung)
coef(fit_mpl)
coef(fit_mpl, parameters = "Theta")
coef(summary(fit_mpl))
} # }
```
