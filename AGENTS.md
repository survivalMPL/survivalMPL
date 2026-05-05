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

## Experiments

Use `experiments/` for exploratory scripts, temporary investigation folders, and
prototype outputs that should not be part of the R package build. The folder is
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
| Recent changelog | [old_sources/DEVLOG.md](old_sources/DEVLOG.md) |
| Refactoring plan (phases 0–5) | [old_sources/REFACTORING_PLAN.md](old_sources/REFACTORING_PLAN.md) |
