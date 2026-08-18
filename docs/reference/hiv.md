# Pseudo-HIV Seroconversion Data

Simulated left- and right-censored survival data for a cohort of people
who inject drugs, followed for HIV seroconversion. Time is measured in
years from the start of injecting drug use.

## Format

A data frame with 300 observations on 6 variables:

- t_L:

  Numeric lower bound of the observed interval, in years; `NA` for
  left-censored subjects.

- t_R:

  Numeric upper bound of the observed interval, in years; `NA` for
  right-censored subjects. Equal to `t_L` for exact events.

- Sharing:

  Indicator: shared injecting equipment at enrolment.

- Prison:

  Indicator: history of incarceration.

- Female:

  Indicator: subject is female.

- Age_centred:

  Age in decades, centred at the sample mean.

## Source

None: these data are entirely synthetic, with no real cohort and no
published analysis behind them. The design was written for this
package's left-censoring tutorial. The covariates are ones a study of
this kind would plausibly record, but the coefficient values above are
chosen for illustration and are not estimates from any study. Do not
cite these numbers as evidence about HIV seroconversion.

## Details

Each subject has a first HIV test some time after the time origin. A
subject already seropositive at that test seroconverted somewhere in
\\(0, t_R\]\\ and is **left censored**; a subject negative at the first
test is retested frequently, so a later seroconversion is recorded as an
**exact event**; a subject still negative when follow-up ends is **right
censored**. There is no interval censoring by construction, which makes
this dataset a clean illustration of left censoring on its own.

Of the 300 subjects, 96 are left censored, 117 are exact events and 87
are right censored. Pass to `coxph_mpl` via
`Surv(t_L, t_R, type = "interval2")`.

The true baseline hazard is Weibull with shape 1.5 and scale 6, \\H_0(t)
= (t/6)^{1.5}\\, and the true regression coefficients used in the
simulation are:

|               |       |
|---------------|-------|
| `Sharing`     | 0.90  |
| `Prison`      | 0.45  |
| `Female`      | -0.25 |
| `Age_centred` | -0.15 |

Generated with `set.seed(11)`; see `dev/data-raw/hiv.R` to reproduce.

## Examples

``` r
data(hiv)
if (FALSE) { # \dontrun{
fit <- coxph_mpl(
  Surv(t_L, t_R, type = "interval2") ~
    Sharing + Prison + Female + Age_centred,
  data = hiv, basis = "msplines"
)
summary(fit)
} # }
```
