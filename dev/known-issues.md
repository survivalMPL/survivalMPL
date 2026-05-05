# Known Issues & Quirks

## Non-finite Hessian weights — PENDING FIX

**Status: investigated, fix not yet applied.**

When `M_S1_ni1` / `M_S2_ni1` underflow to 0 and `M_S1mS2_ni1` is clamped to epsilon,
the interval-censored weight expression:

```r
M_S1_ni1 * M_S2_ni1 / M_S1mS2_ni1^2 * (M_H2_ni1 - M_H1_ni1)^2 +
  (M_S1_ni1 * M_H1_ni1 - M_S2_ni1 * M_H2_ni1) / M_S1mS2_ni1
```

produces `NaN`/`Inf`. This is passed directly into `crossprod()` at **`R/coxph.r:298`**
(inner Newton loop) and **`R/coxph.r:523`** (outer loop), causing:

```
Error in chol.default(M_hessbeta_p1) :
  the leading minor of order 1 is not positive definite
```

A secondary failure mode: if the log-likelihood itself becomes non-finite (`thetaRtheta`
overflows to `NaN`), the step-halving `if (s_lik < s_lik_OLD)` throws
"missing value where TRUE/FALSE needed".

**Proposed fix** (in `non-infinite_issue.R`): guard the weight before `crossprod()`:
```r
w_ni_guarded <- ifelse(is.finite(w_ni), w_ni, 0)
```
This is a stopgap. A proper fix should address the root cause (extreme cumulative hazards
from near-zero survival values) rather than silently zeroing the contribution.

## Seed sensitivity in melanoma simulation

`simulate_mel()` uses `h0(t) = t^(-1/2)`, which is steep near 0. Seeds 7, 22, 41, 55, 64, 94
produce extreme draws that make plots look bad. Use `set.seed(1)` in the vignette.

## `data(lung)` does not work

`data(lung)` fails silently or errors. Use `survival::lung` directly, or
`data(lung, package = "survival")`.

## `coxph.r` filename

The main fitting engine file has a lowercase `.r` extension — this is intentional, do not rename it.
