# AGENTS.md — survivalMPL

AI agent entry point. Read this first, then follow links for detail.

------------------------------------------------------------------------

## What this package does

`survivalMPL` fits Cox proportional hazards models with right, left, and
interval censoring using maximum penalised likelihood. The baseline
hazard is estimated non-parametrically via smooth basis functions (step,
Gaussian, M-splines, Epanechnikov, B-splines).

------------------------------------------------------------------------

## Orientation

``` r
git log --oneline -10   # recent work
git branch              # current branch (dev: test-dev, stable: main)
devtools::load_all()
devtools::test()        # must pass before any commit
devtools::check()       # full CRAN check
```

------------------------------------------------------------------------

## Topic guides

| Topic | File |
|----|----|
| Basis registry architecture | [dev/architecture.md](https://CRAN.R-project.org/package=survivalMPL/dev/architecture.md) |
| Development workflow & rules | [dev/workflow.md](https://CRAN.R-project.org/package=survivalMPL/dev/workflow.md) |
| Vignette guidelines | [dev/vignettes.md](https://CRAN.R-project.org/package=survivalMPL/dev/vignettes.md) |
| pkgdown site configuration | [dev/pkgdown.md](https://CRAN.R-project.org/package=survivalMPL/dev/pkgdown.md) |
| Known bugs & quirks | [dev/known-issues.md](https://CRAN.R-project.org/package=survivalMPL/dev/known-issues.md) |
| Recent changelog | [old_sources/DEVLOG.md](https://CRAN.R-project.org/package=survivalMPL/old_sources/DEVLOG.md) |
| Refactoring plan (phases 0–5) | [old_sources/REFACTORING_PLAN.md](https://CRAN.R-project.org/package=survivalMPL/old_sources/REFACTORING_PLAN.md) |
