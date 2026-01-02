# Breast Cosmesis Data

Interval-censored breast cosmesis data from Finkelstein and Wolfe
(1985), discussed by Moore (2016, example 12.2) and available in the
interval package. Compared to that version, `bcos2` recodes the lower
bound of left-censored data (`NA` instead of `0`) and the upper bound of
right-censored data (`NA` instead of `Inf`) to ease identification via
`Surv(type = "interval2")`.

## Format

A data frame with 94 observations on 3 variables:

- left:

  Numeric lower bound.

- right:

  Numeric upper bound.

- treatment:

  Factor with levels `Rad` and `RadChem`.

## Source

Finkelstein, D.M., and Wolfe, R.A. (1985). A semiparametric model for
regression analysis of interval-censored failure time data. *Biometrics*
41: 731-740.

Moore, D.K. (2016). *Applied Survival Analysis Using R*. Springer.

## Examples

``` r
data(bcos2)
```
