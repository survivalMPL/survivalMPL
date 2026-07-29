# Fix: numerical-stability crashes in `coxph_mpl` for interval-censored M-spline fits

> Notes / proposed plan — saved for later. Not yet implemented.
> Companion to `dev/experiments/non-infinite_issue.R` and the failing-seed notes in
> `dev/experiments/book_sim.R` (lines 190–194).

## Context

`book_sim.R` (lines 190–194) documents that `coxph_mpl()` crashes on specific RNG seeds when
fitting interval-censored data with the M-spline basis:

- `set.seed(7)` / `(41)` → `Error in chol.default(M_hessbeta_p1) : the leading minor of order 1 is not positive`
- `set.seed(630)/(22)/(55)/(64)/(94)` → `Error in if (s_lik < s_lik_OLD) { : missing value where TRUE/FALSE needed`

Both are symptoms of the Newton-MI iteration **diverging** (β and/or θ running away) and the
math then overflowing. A trace of `R/coxph.r` (the whole fitting engine lives in this one file)
identifies **three independent root causes**, each mapping to a distinct part of the failure.

### Root cause 1 — line-search predicate mishandles non-finite likelihood → error 1
In both update blocks the full Newton step is taken at ω=1, then `s_lik` is recomputed and compared
via `if (s_lik < s_lik_OLD)` (lines 334, 337, 434, 437). When a step drives `μ = exp(Xβ)` up,
cumulative hazards overflow, `S1,S2 → 0`, `S1-S2` is clamped to `ε=1e-10`, and the interval Hessian
weight `S1*S2/(S1-S2)^2*(H2-H1)^2 + (S1*H1-S2*H2)/(S1-S2)` becomes `NaN/Inf`. The comparison against
`NaN` is `NA`, which both throws here *and* (with the naive `is.finite(...) && <` guard from
`non-infinite_issue.R`) would **skip/exit** step-halving and *accept* the divergent iterate —
surfacing as the chol failure on the next iteration. **The naive guard is backwards: it treats
non-finite as "stop", which locks in the runaway.** The runaway is never rejected.

### Root cause 2 — smoothing-parameter update is unconstrained → error 2
Lines 576–578: `s_df = m - trace(HRinv·Rstar)/sigma2_old`, `s_sigma2 = θ'Rθ / s_df`,
`s_lambda = 1/(2·s_sigma2)`. Since `θ'Rθ ≥ 0`, an ill-conditioned `H` can push `s_df < 0`, making
`sigma2 < 0` and **`s_lambda < 0`** (the `-49955049` in the bug notes). A negative λ *flips the
penalty sign* — the objective now rewards roughness — so θ explodes, `θ'Rθ` overflows to `NaN`,
`s_lik` is `NaN`, and the comparison throws. `df ∈ (0, m]` and `λ ≥ 0` are the *valid statistical
ranges*; the current code never enforces them. This is plausibly the true origin of the error-2
seeds, and the `non-infinite_issue.R` band-aids miss it entirely.

### Root cause 3 — β-Hessian is not guaranteed PD and is never checked before `chol()`
Line 305 calls `chol(M_hessbeta_p1)` directly. The censored contributions (lines 296–298) are not
sign-definite, and overflow can make entries non-finite, so the matrix can be indefinite or `NaN`
even on otherwise-recoverable iterations. The existing `if (p==1 && hess==0)` patch (lines 300–304)
handles only a trivial scalar case.

**Goal (user-confirmed):** primary = never crash, emit a non-convergence warning and return the
last stable iterate on genuinely degenerate data; secondary = converge as many of the failing seeds
as possible **without regressing currently-passing seeds**. All three root causes in scope.

## Approach

All edits in `R/coxph.r` unless noted. Reuse the existing `try(chol2inv(chol(...)))`→fallback idiom
already used at lines 569–575 / 646–665 rather than inventing new machinery.

### Fix 1 — correct line-search acceptance (β block 334–382, θ block 434–488)
- Define an "improved" predicate: `accept <- is.finite(s_lik) && s_lik >= s_lik_OLD`.
  - Outer trigger: `if (!accept) { ... }`.
  - Inner loop: `while (!accept && i <= maxhalve)`, recomputing `accept` each pass.
  - Treat non-finite `s_lik` as **reject** (keep shrinking ω), never as "stop".
- If halving ends without an accepted step, **revert** the parameter to its `_OLD` value
  (β←`M_beta_p1_OLD`, θ←`M_theta_m1_OLD`), recompute dependent `μ/H/S/s_lik`, set
  `converged <- FALSE`, and break the inner `k` loop. Currently the code keeps the last (bad)
  iterate — this is what feeds the next-iteration chol failure.
- Keep the existing `i > 500` ω-floor as `maxhalve`.

### Fix 2 — constrain the smoothing update (lines 576–579)
- Clamp df to the valid range: `s_df <- min(max(s_df, df_floor), m)`, small positive `df_floor`
  (e.g. `1e-3`).
- Guard `s_sigma2` finite and `> 0`; if `θ'Rθ` or `s_df` is degenerate, **retain `s_lambda_old`**
  instead of updating.
- Enforce `s_lambda <- max(s_lambda, 0)` and finiteness before it re-enters `TwoLRtheta`
  (line 579) and the next outer-iteration penalty (line 256). A correct line search (Fix 1) stops
  the *crash*; only this keeps the fit *correct* (penalty sign preserved).

### Fix 3 — make the β-step solve robust (line 305, with guards on 293–298 & 272–276)
- Sanitize inner-loop weights/objects: replace non-finite entries of the interval/left-censored
  weight vectors (and of `M_gradbeta_p1`) with `0` before forming the step —
  `w <- ifelse(is.finite(w), w, 0)`. A `NaN` gradient otherwise makes `step = Hinv %*% NaN`.
- Replace bare `chol2inv(chol(M_hessbeta_p1))` with the file's own pattern:
  `temp <- try(chol2inv(chol(M_hessbeta_p1)), silent=TRUE)`; on `try-error`/non-finite, **damp**
  (Levenberg-Marquardt): add `tau * diag(p)`, increasing `tau` until PD; else fall back to
  `MASS::ginv(M_hessbeta_p1)`.
- **Keep inference clean:** damp only the local step-solve variable `M_hessbeta_p1`. Do **not**
  touch `H[1:p,1:p]` (518–523), `M_2` (610), or the SE Cholesky/`ginv` blocks (646–665) — those are
  recomputed fresh at the final iterate and feed the reported `M2HM2`/`M2QM2`/`H` SEs the simulation
  records. Separate variables ⇒ separation is automatic; just verify no leakage.

### Crash-proofing / reporting
- Thread a `converged` logical through the loops (set `FALSE` on: unimproved step after halving,
  `full.iter > max.iter[3]`, or non-finite final `s_lik`). Emit
  `warning("coxph_mpl did not converge ...")` when `FALSE`.
- Add `fit$converged <- converged` to the return list (near 670–684). Non-breaking.

### Optional secondary (helps "converge where possible", not required to stop crashes)
- Stable interval term: `S1 - S2 = exp(-H1) * (1 - exp(-(H2-H1)))`; cap cumulative hazards
  (`H <= 700`) so legitimate large-hazard points don't spuriously overflow before the line search
  reacts. Implement only if seeds still fail to *converge* after Fixes 1–3.

## Critical files
- `R/coxph.r` — all three fixes + convergence flag (lines 305, 334–382, 434–488, 576–579, 670–684).
- `R/control.R` — only if exposing `df_floor`/`maxhalve` as control args (defer; internal consts first).
- `tests/testthat/test-smoke.R` — extend the interval-censored test (≈ lines 103–115) with a
  regression case on a previously-failing seed.

## Why the `non-infinite_issue.R` band-aids are insufficient
- `w_ni_guarded <- ifelse(is.finite(w_ni), w_ni, 0)` — partial: zeroing weights can still leave the
  Hessian singular/indefinite, and doesn't address the gradient or the divergence source.
- `is.finite(s_lik) && is.finite(s_lik_OLD) && (s_lik < s_lik_OLD)` — **backwards**: makes the
  predicate `FALSE` on `NaN`, which skips/exits step-halving and *accepts* the divergent step,
  converting error 2 into error 1. Correct predicate keeps halving while the new point is *not* a
  finite improvement.
- Neither touches root cause 2 (negative λ), the likely true origin of the error-2 seeds.

## Verification
1. **No crash on documented failing seeds.** Run `simulate_mel(300, 0.37, 0.6, 1.2)` fit
   (`book_sim.R` 200–208) under `devtools::load_all()` for seeds **7, 41** (chol) and
   **630, 22, 55, 64, 94** (NaN), plus the weibull loop. Each returns a fit or a clean
   non-convergence warning — never an error.
2. **No regression on passing seeds.** Re-run a sample of currently-passing seeds; confirm
   `fit$coef$Beta` and `fit$se$Beta$M2HM2` unchanged within tolerance.
3. **λ stays valid.** Instrument a failing seed: `s_lambda ≥ 0` and `s_df ∈ (0, m]` throughout.
4. `devtools::test()` / `R CMD check` green; new regression test passes.
