# Print Method for `summary.coxph_mpl`

Extracts additional information for a fitted model and returns an object
suitable for printing. Baseline hazard parameters smaller than
`min.theta` are omitted unless `full = TRUE`.

## Usage

``` r
# S3 method for class 'summary.coxph_mpl'
print(x, se = "M2QM2", ...)

# S3 method for class 'coxph_mpl'
summary(object, se = "M2QM2", full = FALSE, ...)
```

## Arguments

- x:

  An object of class `"summary.coxph_mpl"`.

- se:

  Inference method. One of `"H"`, `"M2QM2"`, or `"M2HM2"`. Default is
  `"M2QM2"`.

- ...:

  Additional arguments passed to methods.

- object:

  A fitted model of class `"coxph_mpl"`.

- full:

  Logical; if `TRUE`, include inference for baseline hazard parameters.
  Default `FALSE`.

## Value

Invisibly returns `x`.

An object of class `"summary.coxph_mpl"` with components:

- Beta:

  Matrix of regression estimates, standard errors, z-statistics, and
  p-values.

- Theta:

  Baseline hazard estimates (or, if `full = TRUE`, a matrix of estimates
  with standard errors, z-statistics, and p-values).

- inf:

  List with convergence details, penalised likelihood value, and control
  settings.

## See also

\[summary.coxph_mpl()\], \[coxph_mpl()\], \[coxph_mpl.control()\]

\[coxph_mpl()\], \[coxph_mpl.control()\], \[plot.coxph_mpl()\]

## Examples

``` r
if (FALSE) { # \dontrun{
data(lung)
fit_mpl <- coxph_mpl(Surv(time, status == 2) ~ age + sex + ph.karno + wt.loss,
                     data = lung)
summary(fit_mpl, full = TRUE)
summary(fit_mpl, se = "M2HM2")
} # }
```
