# Plot a `coxph_mpl` Object

Plot the bases used to estimate the baseline hazard along with the
estimated baseline hazard, cumulative baseline hazard, and baseline
survival functions. Each plot can be toggled with `which`.

## Usage

``` r
# S3 method for class 'coxph_mpl'
plot(x, se = "M2QM2", ask = TRUE, which = 1:4, upper.quantile = 0.95, ...)
```

## Arguments

- x:

  A fitted model of class `"coxph_mpl"`.

- se:

  Inference method for confidence intervals. One of `"H"`, `"M2QM2"`, or
  `"M2HM2"`. Default is `"M2QM2"`.

- ask:

  Logical; whether to prompt before each plot. See \[graphics::par()\].
  Default `TRUE`.

- which:

  Integer vector selecting plots to produce (subset of `1:4`).

- upper.quantile:

  Quantile of the response used to set the upper x-axis limit for
  baseline function plots. Default `0.95`.

- ...:

  Additional arguments passed to plotting functions.

## Details

Bases whose estimates are near zero (below `min.theta` in
\[coxph_mpl.control()\]) are drawn with dashed lines. Confidence
intervals are obtained via the delta method.

## See also

\[coxph_mpl()\], \[coxph_mpl.control()\], \[coxph_mpl.object()\],
\[summary.coxph_mpl()\]

## Examples

``` r
if (FALSE) { # \dontrun{
data(lung, package = "survival")
fit_mpl <- coxph_mpl(Surv(time, status == 2) ~ age + sex + ph.karno + wt.loss,
                     data = lung)
par(mfrow = c(2, 2))
plot(fit_mpl, ask = FALSE, cex.main = 0.75)
} # }
```
