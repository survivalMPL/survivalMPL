# Getting Started with survivalMPL

## Overview

`survivalMPL` fits Cox proportional hazards models with right, left, and
interval censoring using maximum penalised likelihood (MPL). The
baseline hazard is estimated with a smooth basis (step function,
Gaussian, M-splines, or Epanechnikov) and non-negativity constraints are
enforced during optimisation.

Below are two small examples: a right-censored fit on the
[`survival::lung`](https://rdrr.io/pkg/survival/man/lung.html) data and
an interval-censored fit on the bundled `bcos2` data. To keep the
vignette fast to build, we use modest iteration limits in
[`coxph_mpl.control()`](https://CRAN.R-project.org/package=survivalMPL/reference/coxph_mpl.control.md).

## Right-censored example (`lung`)

``` r

library(survivalMPL)
#> Loading required package: survival
#> Loading required package: MASS
library(survival)

lung <- na.omit(lung[, c("time", "status", "age", "sex", "ph.karno", "wt.loss")])

fit_lung <- coxph_mpl(
  Surv(time, status == 2) ~ age + sex + ph.karno + wt.loss,
  data = lung,
  control = coxph_mpl.control(
    n.obs = sum(lung$status == 2),
    max.iter = c(40, 2000, 4000),
    smooth = 0
  )
)

summary(fit_lung)
#> 
#> coxph_mpl(formula = Surv(time, status == 2) ~ age + sex + ph.karno + 
#>     wt.loss, data = lung, control = coxph_mpl.control(n.obs = sum(lung$status == 
#>     2), max.iter = c(40, 2000, 4000), smooth = 0))
#> 
#> -----
#> 
#> Cox Proportional Hazards Model Fit Using MPL 
#> 
#> 
#> Penalized log-likelihood  :  -1052.062
#> Fixed smoothing value     :  0
#> Convergence               :  Yes (8 iter.) 
#> 
#> Data             : lung
#> Number of obs.   : 214
#> Number of events : 152 (71.02804%)
#> Number of cens.  :  62 (28.97196%)
#> 
#> Regression parameters : Surv(time, status == 2) ~ age + sex + ph.karno + wt.loss
#>            Estimate Std. Error z-value Pr(>|z|)   
#> age       0.0151013  0.0101696  1.4849 0.137559   
#> sex      -0.5148282  0.1691462 -3.0437 0.002337 **
#> ph.karno -0.0126768  0.0086042 -1.4733 0.140664   
#> wt.loss  -0.0020624  0.0069833 -0.2953 0.767745   
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Baseline hasard parameters approximated using Uniform :
#>  (11 equal events bins)
#>            1            2            3            4            5            6 
#> 3.091190e-03 3.073782e-03 8.005268e-03 6.378672e-03 5.284837e-03 6.891828e-03 
#>            7            8            9           10           11 
#> 6.688682e-03 6.777236e-03 1.262414e-02 9.950425e-03 2.312947e-10 
#> 
#> -----
```

You can inspect the baseline hazard components or predict survival
curves:

``` r

pred_lung <- predict(fit_lung, type = "survival")
head(pred_lung)
#>        time  survival           se       low      high
#> 1  4.999000 1.0000000 0.0000000000 1.0000000 1.0000000
#> 2  6.017019 0.9986404 0.0002927774 0.9980548 0.9992259
#> 3  7.035038 0.9972826 0.0005847586 0.9961131 0.9984521
#> 4  8.053057 0.9959267 0.0008759453 0.9941748 0.9976785
#> 5  9.071076 0.9945726 0.0011663392 0.9922399 0.9969052
#> 6 10.089095 0.9932203 0.0014559417 0.9903084 0.9961322
```

## Interval-censored example (`bcos2`)

``` r

data(bcos2)

fit_bcos <- coxph_mpl(
  Surv(left, right, type = "interval2") ~ treatment,
  data = bcos2,
  control = coxph_mpl.control(
    basis = "msplines",
    n.obs = nrow(bcos2),
    max.iter = c(40, 2000, 4000),
    smooth = 0
  )
)

summary(fit_bcos)
#> 
#> coxph_mpl(formula = Surv(left, right, type = "interval2") ~ treatment, 
#>     data = bcos2, control = coxph_mpl.control(basis = "msplines", 
#>         n.obs = nrow(bcos2), max.iter = c(40, 2000, 4000), smooth = 0))
#> 
#> -----
#> 
#> Cox Proportional Hazards Model Fit Using MPL 
#> 
#> 
#> Penalized log-likelihood  :  -140.254
#> Fixed smoothing value     :  0
#> Convergence               :  NO 
#> 
#> Data             : bcos2
#> Number of obs.   : 94
#> Number of events :  0 (  0%)
#> Number of cens.  : 94 (100%)
#> 
#> Regression parameters : Surv(left, right, type = "interval2") ~ treatment
#>                  Estimate Std. Error z-value Pr(>|z|)   
#> treatmentRadChem  0.87311    0.33348  2.6181 0.008841 **
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Baseline hasard parameters approximated using M-Splines :
#>  (2 (min/max) + 8 quantile knots + 2 equally spaced knots + 3 (order) - 2 = 13 parameters) 
#>            1            2            3            4            5            6 
#> 4.284005e-02 1.608827e-02 5.030707e-02 3.116483e-02 1.317693e-01 1.213149e-01 
#>            7            8            9           10           11           12 
#> 2.396608e-02 2.067590e-01 2.068148e-01 6.402848e-10 4.310156e+00 6.628441e-01 
#>           13 
#> 6.628441e-01 
#> 
#> -----
```

Predicted survival for the two treatment groups:

``` r

pred_bcos <- predict(fit_bcos, type = "survival", i = 1:2)
#> Warning: only the first observation will be considered
head(pred_bcos)
#>       time  survival          se       low high
#> 1 4.000000 1.0000000 0.000000000 1.0000000    1
#> 2 4.056056 0.9983046 0.001302038 0.9957005    1
#> 3 4.112112 0.9966486 0.002543082 0.9915625    1
#> 4 4.168168 0.9950313 0.003724858 0.9875816    1
#> 5 4.224224 0.9934518 0.004849072 0.9837537    1
#> 6 4.280280 0.9919096 0.005917414 0.9800747    1
```

For richer diagnostics, see `plot(fit_bcos)`, `plot(predict(fit_bcos))`,
and `plot(residuals(fit_lung))`.
