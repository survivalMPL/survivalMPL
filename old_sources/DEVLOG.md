# survivalMPL Development Log

This file tracks development work across released and in-progress versions.

## 0.2-4.9001 (development, 2026-04 / 2026-05)

### Added
- Refactoring roadmap and simulation material (`REFACTORING_PLAN.md`, `book_sim.R`).
- Regression safety net for baseline behavior across basis types (`tests/testthat/test-baseline.R`, snapshot file).
- Basis registry infrastructure and startup registration (`R/basis.R`, `R/zzz.R`).
- Dedicated basis implementation files: uniform, Gaussian, Epanechnikov, M-splines, and B-splines.

### Changed
- Refactored monolithic method code into per-concern files (`control.R`, `plot.R`, `summary.R`, `print.R`, `coef.R`, `residuals.R`, `predict.R`) with no intended logic change in phase 1.
- Wired basis dispatching through core call sites in model fitting, prediction, plotting, summaries, and residual workflows.
- Cleaned up dead internal paths and updated basis file contracts to improve extensibility.

### Fixed
- Made pseudo melanoma simulation loop reproducible and consistent by locking the loop seed to `41` (`book_sim.R`).

### Investigated, not yet applied
- Non-finite interval-censored Hessian weights: when `M_S1_ni1`/`M_S2_ni1` underflow to 0
  and `M_S1mS2_ni1` is clamped to epsilon, the weight expression passed to `crossprod()` at
  `coxph.r:298` and `coxph.r:523` becomes `NaN`/`Inf`, causing `chol(M_hessbeta_p1)` to fail
  with "leading minor of order 1 is not positive definite".
  Proposed fix (`ifelse(is.finite(w_ni), w_ni, 0)`) and two related failure modes documented
  in `non-infinite_issue.R`. Fix has **not** been applied to `R/coxph.r` — needs a cleaner
  solution before merging.
- Fixed pkgdown code blocks rendering with leading blank lines caused by badge-stripping; resolved by introducing `index.md` at the package root as the pkgdown homepage source, keeping `README.md` for GitHub only.
- Removed stale `pkgdown/DEVLOG.md` and `pkgdown/REFACTORING_PLAN.md` build artifacts that caused pkgdown to re-build those files on every run.

### Added (2026-05)
- Interval-censoring tutorial vignette (`vignettes/interval-censoring.Rmd`) reproducing Examples 2.6 and 2.7 from Ma, Webb & Hudson (2024) with minimal corrections.
- BibTeX citation entry for the companion book (`vignettes/references.bib`).
- Getting-started tutorial vignette (`vignettes/getting-started.Rmd`) with right-censored (`lung`) and interval-censored (`bcos2`) quick examples.
- pkgdown Tutorials dropdown navbar replacing the default Articles link; tutorials index configured in `_pkgdown.yml`.
- `index.md` at package root as the pkgdown homepage source (badge-free, avoids downlit blank-line artifact).
- AI agent documentation: `AGENTS.md` (index) and topic files under `dev/` (architecture, workflow, vignettes, pkgdown, known-issues).

### Changed (2026-05)
- Moved `DEVLOG.md` and `REFACTORING_PLAN.md` from package root to `old_sources/` to prevent pkgdown from building them as site pages.

## 0.2-4 (2025-03-23)
- Changed maintainer email address.
- Added contributor.

## 0.2-3 (2022-11-21)
- Sorted minor typos in documentation.

## 0.2-2 (2022-11-21)
- Extended the example section of `coxph_mpl`.
- Amended code to allow models without covariates.
- Removed a duplicated item in `coxph_mpl.object`.

## 0.2-1 (2021-11-14)
- Change of maintainer (D-L Couturier).
- Minor fixes (warning messages from `try` and `model.matrix`).

## 0.2 (2017-12-11)
- Interval censoring implementation.

## 0.1.2 (2017-10-13)
- Change of maintainer (Maurizio Manuguerra).

## 0.1.0 (2014-04-03)
- Initial package release (from archived source metadata).
