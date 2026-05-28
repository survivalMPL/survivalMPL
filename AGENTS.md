# AGENTS.md — survivalMPL

AI agent entry point. Read this first, then follow links for detail.

---

## What this package does

`survivalMPL` fits Cox proportional hazards models with right, left, and interval censoring
using maximum penalised likelihood. The baseline hazard is estimated non-parametrically via
smooth basis functions (step, Gaussian, M-splines, Epanechnikov, B-splines).

---

## Orientation

```r
git log --oneline -10   # recent work
git branch              # current branch (dev: test-dev, stable: main)
devtools::load_all()
devtools::test()        # must pass before any commit
devtools::check()       # full CRAN check
```

---

## Testing protocol — change → test → change

**Run smoke tests before and after every change. Only commit if they pass.**

```r
# Fast check (~5 s) — run before starting work and after each change
devtools::test(filter = "smoke")

# Full suite — must pass before any commit
devtools::test()
```

Smoke tests cover (`tests/testthat/test-smoke.R`):
- Right-censored numerical stability: `lung` × msplines + uniform (Beta, ploglik)
- S3 method structure: correct class and column names for `predict`, `residuals`, `coef`, `summary`
- Print and plot run without error and return invisibly
- Interval-censored run-without-error: `melanoma` × msplines

Snapshots in `tests/testthat/_snaps/` are the numerical baseline anchored to
installed v0.2-4. Use `testthat::snapshot_accept()` only for confirmed
platform-level floating-point noise — not to paper over a real value shift.

---

## Experiments

Use `dev/experiments/` for exploratory scripts, temporary investigation folders, and
prototype outputs that should not be part of the R package build. The entire `dev/` tree is
excluded via `.Rbuildignore`; move anything production-ready into the appropriate
package location before documenting or testing it.

---

## Topic guides

| Topic | File |
|-------|------|
| Basis registry architecture | [dev/architecture.md](dev/architecture.md) |
| Development workflow & rules | [dev/workflow.md](dev/workflow.md) |
| Vignette guidelines | [dev/vignettes.md](dev/vignettes.md) |
| pkgdown site configuration | [dev/pkgdown.md](dev/pkgdown.md) |
| Known bugs & quirks | [dev/known-issues.md](dev/known-issues.md) |
| Recent changelog | [dev/DEVLOG.md](dev/DEVLOG.md) |
