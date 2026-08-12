# Pseudo-Melanoma Survival Data

Simulated interval-censored survival data based on the design in Moore
(2016, Examples 2.6–2.7). The true baseline hazard is Weibull with shape
0.5: \\h_0(t) = t^{-1/2}\\, \\H_0(t) = 2\sqrt{t}\\.

## Format

A data frame with 300 observations on 10 variables:

- t_L:

  Numeric lower bound of the observed interval (0 for left-censored
  observations).

- t_R:

  Numeric upper bound of the observed interval (`Inf` for right-censored
  observations).

- Arm:

  Indicator: tumour located on the arm (reference: head/neck).

- Leg:

  Indicator: tumour located on the leg.

- Trunk:

  Indicator: tumour located on the trunk.

- mm1to2:

  Indicator: tumour thickness 1–2 mm (reference: \< 1 mm).

- mm2to4:

  Indicator: tumour thickness 2–4 mm.

- mm4plus:

  Indicator: tumour thickness \> 4 mm.

- Female:

  Indicator: patient is female.

- Age_centred:

  Age in decades, centred at the sample mean.

## Source

Moore, D.K. (2016). *Applied Survival Analysis Using R*. Springer.

## Details

Observations are classified as exact events (`t_L == t_R`), left-
censored (`t_L == 0`), right-censored (`t_R == Inf`), or
interval-censored otherwise. Pass to `coxph_mpl` via
`Surv(t_L, t_R, type = "interval2")`.

True regression coefficients (`beta_true`) used in the simulation:

|               |       |
|---------------|-------|
| `Arm`         | -0.56 |
| `Leg`         | 0.01  |
| `Trunk`       | -0.22 |
| `mm1to2`      | 0.22  |
| `mm2to4`      | 0.87  |
| `mm4plus`     | 1.13  |
| `Female`      | -0.17 |
| `Age_centred` | 0.14  |

Generated with `set.seed(1)`; see `dev/data-raw/melanoma.R` to
reproduce.

## Examples

``` r
data(melanoma)
if (FALSE) { # \dontrun{
fit <- coxph_mpl(
  Surv(t_L, t_R, type = "interval2") ~
    Arm + Leg + Trunk + mm1to2 + mm2to4 + mm4plus + Female + Age_centred,
  data = melanoma, basis = "m"
)
summary(fit)
} # }
```
