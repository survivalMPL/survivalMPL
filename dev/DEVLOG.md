---
tags: [survivalMPL, devlog]
created: 2026-05-08
status: active
---

# survivalMPL — Devlog

## Package snapshot

| Field | Value |
|---|---|
| Name | survivalMPL |
| Repo path | `C:\GIT\survivalMPL.package` |
| Current branch | `test-dev` |
| Version | 0.2-4.9001 |
| Collaborators | Dominique-Laurent Couturier (aut, cre), Jun Ma (aut), Stephane Heritier (aut), Maurizio Manuguerra (aut), Serigne Lo (aut); Nikita Mozgunov (dev contributor) |

---

## Feature branches

| Branch | Status | Description |
|---|---|---|
| `main` | 🟢 done | Stable CRAN-tagged baseline (0.2-4) |
| `test-dev` | 🟡 in progress | Full refactor: basis registry, file split, vignettes, pkgdown |
| `origin/dev` | 🔴 blocked | Pre-refactor scratch; superseded by `test-dev` |

---

## Log

### 2026-05-23 — Vignette framing fixes and pkgdown blank-lines resolution

**What was done**

- `README.md`: fixed a double blank line (`\r\n\r\n\r\n`) inside the Quick Start code block
  between the closing `)` and `summary(fit_lung)`. The Edit tool could not match CRLF sequences;
  applied via Python binary patch. This was one of the two root causes of the pkgdown
  blank-lines accumulation bug.
- `vignettes/coxph-comparison.Rmd` — Baseline survival section rewritten for accuracy:
  - Previous framing "three estimators of the baseline survival function, evaluated at mean
    covariate values" was wrong: Kaplan–Meier is marginal and unadjusted, not at mean covariates.
  - New framing: KM = marginal reference, Breslow and M-splines = same quantity
    $\hat{S}(t \mid \bar{\mathbf{x}})$, should track each other closely.
  - y-axis label changed from `hat(S)[0](t)` (baseline only) to `hat(S)(t)` (general survival).
  - Plot title and legend updated to flag KM as marginal.
  - Closing paragraph explains why KM sits above the adjusted curves.
- `vignettes/coxph-comparison.Rmd` — Coefficient table: added a `Covariate` column with
  explicit unit labels; units footnote moved from paragraph text into the `kable` caption
  (`^a^ The forest plot rescales these covariates to per-10-unit changes`).

**Design decisions**
- Breslow and MPL both use mean covariates; this is confirmed in `predict.coxph_mpl` line 55:
  `xTB = if(is.null(i)){apply(object$data$X, 2, mean)}else{...}` — no code change needed.
- KM is kept in the plot as a reference curve, not as a competing estimator — the text now
  makes this distinction explicit.

**pkgdown blank-lines fix summary**
Two root causes, both addressed:
1. `pkgdown/index.html` tracked in git at stale pkgdown format — use `pkgdown::clean_site()`
   before `build_site()` after any pkgdown version change (documented in `dev/pkgdown.md`).
2. Double blank line in `README.md` — fixed.

---

### 2026-05-21 — Design pattern review

**What was discussed**
- Audited the codebase against the refactoring.guru design pattern catalog.
- Identified patterns already in use: Strategy (basis registry with pluggable `knots_fn`/`matrix_fn`/`penalty_fn`), Registry (`.basis_registry` environment with `register_basis`/`get_basis`), Template Method (S3 generics `print`/`summary`/`predict`/`coef`/`residuals`/`plot`).

**Patterns considered for future work**
- Builder for `coxph_mpl.control()`: would make the 17-parameter function more ergonomic and testable. Rejected for now — backward compatibility concern (existing callers pass positional/named args directly). Could be added as an opt-in alternative API without breaking anything, but not a priority.
- Observer for iteration callbacks: selected for implementation (see TODO #8). A `callbacks = NULL` parameter on `coxph_mpl()` is fully additive. Motivating use case: long fits on large datasets have no visibility into convergence; callbacks enable progress bars, custom stopping criteria, and convergence logging.
- Null Object for `smooth = NULL` vs `smooth = 0`: noted but deferred — the `is.null(smooth)` checks are contained and not causing active bugs.

---

### 2026-05-21 — Basis comparison and coxph vs coxph_mpl tutorials

**What was built**
- `vignettes/basis-comparison.Rmd` — tutorial covering all four production basis types (uniform, Gaussian, M-splines, Epanechnikov); includes MathJax equations from the GKS technical note, visual comparison of basis shapes using `melanoma` event times for knot placement, and a four-way fitting comparison on the full melanoma dataset with true Weibull baseline overlay
- `vignettes/coxph-comparison.Rmd` — tutorial comparing `survival::coxph()` (partial likelihood) against `coxph_mpl()` (penalised full likelihood) on the `lung` dataset; includes coefficient/SE table, a scaled forest plot (age and Karnofsky per-10-unit) with right-side numeric annotations, and a three-way baseline survival comparison (KM, Breslow step function, M-spline smooth)
- Added `@article{Jun2013}` entry to `vignettes/references.bib`
- Updated `_pkgdown.yml`: both articles added to `articles:` contents and Tutorials navbar menu

**Design decisions**
- Basis visualisation x-axis cropped to 97th percentile of event times; y-axis capped at 99th percentile per panel — prevents the right-skewed melanoma data from collapsing all basis functions into a flat line against one tall spike
- Plot code chunks use `echo=FALSE`; non-plot chunks remain visible so the tutorial is readable as a teaching document without being buried in plot-setup boilerplate
- Forest plot scales `age` and `ph.karno` to per-10-unit changes so all four CIs are visually comparable; raw-unit values retained in the coefficient table
- `rmarkdown::html_document` with `toc_float: true` used for both vignettes (allows floating TOC); acceptable file-size trade-off for tutorial-style articles
- B-splines mentioned in basis tutorial as a registered extensibility proof-of-concept, not as a production option, per project convention

---

### 2026-05-21 — Roxygen audit and pseudo-melanoma dataset

**What was built**
- Audited all 10 roxygen-documented help pages against old package (`dev/old_packages/survivalMPL_0.2-4/man/`)
- Fixed four concrete bugs found during the audit:
  - `print.coxph_mpl` and `print.summary.coxph_mpl` both claimed `@return Invisibly returns x` but neither called `invisible(x)` — added the call to both
  - `coxph_mpl.Rd` showed `...` as "passed to `base::print()`" (wrong) because `R/print.R` was adding a competing `@param ...` on the shared page; removed the competing tag so `coxph.r`'s description ("passed to `coxph_mpl.control()`") wins
  - `summary.coxph_mpl.Rd` had a duplicate `\seealso` block from two `@seealso` tags on the same page; removed `@seealso` from `print.summary.coxph_mpl`, merged unique links into `summary.coxph_mpl`'s tag
  - `coef.coxph_mpl.Rd` had no examples; added `\dontrun{}` block covering `coef()`, `coef(parameters="Theta")`, and `coef(summary(...))`
- Fixed `data(lung)` in all six example blocks — changed to `data(lung, package = "survival")` so examples work under `devtools::load_all()` (closes the known gap logged 2026-04-01)
- Patched the corresponding stale `.Rd` files directly (they still reference `R/methods.R` in their headers); run `devtools::document()` to fully regenerate from the fixed sources
- Created pseudo-melanoma package dataset:
  - `dev/data-raw/melanoma.R` — self-contained generation script (`set.seed(1)`, documents how to regenerate)
  - `data/melanoma.rda` — 300×10 data frame, bzip2-compressed (117 exact events, 132 left-censored, 16 right-censored, 35 interval-censored)
  - `R/data-melanoma.R` — roxygen docs: all 10 variables, true `beta_true` table, baseline hazard formula, `coxph_mpl` example

**Design decisions**
- Generation script placed in `dev/data-raw/` (not `data-raw/` at root) to stay within the already-`.Rbuildignore`-d `dev/` tree; closes the open question logged below
- `data(lung, package = "survival")` preferred over adding a compatibility shim — explicit is more robust across load contexts
- Per-function `@author` and `\encoding{utf8}` tags absent from all current files (present in old package) are intentional under the new roxygen2 setup — package-level authorship is in `DESCRIPTION`; not restored

**Known gaps**
- [ ] Run `devtools::document()` to regenerate all `.Rd` headers (they still say `R/methods.R` for several pages)

---

### 2026-05-05 — Documentation, experiments scaffold, and interval-censoring vignette

**What was built**
- Added `experiments/` directory (`.Rbuildignore`-excluded) for exploratory scripts that must not enter the package build
- Created `dev/` suite of topic docs: `architecture.md`, `workflow.md`, `vignettes.md`, `pkgdown.md`, `known-issues.md`
- Added `AGENTS.md` (root + `pkgdown/`) as AI-agent entry points linking to the topic guides
- Added `vignettes/interval-censoring.Rmd` covering Examples 2.6 and 2.7 from the book
- Updated `_pkgdown.yml` tutorials nav to include the new vignette
- Tidied `TODO.R` and `README.md`

**Key commits**
- `86a62ce` Add experiments/ directory, old_sources, and AGENTS docs; tidy TODO and README
- `18e1eed` Update pkgdown tutorials/nav and add interval-censoring vignette; add dev/legacy docs and ignore patterns

**Design decisions**
- `experiments/` chosen over extending `dev/` so that the distinction between *dev docs* (permanent reference) and *scratch scripts* (throwaway) is unambiguous
- `AGENTS.md` rather than `CLAUDE.md` at repo root so the index is recognised by multiple AI toolchains (Claude Code, Cursor, Copilot)
- Vignette content sourced exclusively from the book PDF, not from `book_sim.R`, which contains internal modifications not in the published text

**Known gaps**
- [ ] `old_sources/DEVLOG.md` reference in `AGENTS.md` is stale — `old_sources/` is deleted from the working tree on `test-dev` (staged deletions visible in `git status`)
- [ ] `dev/workflow.md` still references `non-infinite_issue.R` and `book_sim.R` as root-level untracked files; confirm whether they exist or should be moved to `experiments/`

---

### 2026-04-13 — Seed fix in melanoma simulation

**What was built**
- Updated seed in the pseudo-melanoma simulation loop in a vignette or dev script

**Key commits**
- `fd835f5` fix: update seed in pseudo melanoma simulation loop for consistency

**Design decisions**
- `h0(t) = t^(-1/2)` is steep near 0; seeds 7, 22, 41, 55, 64, 94 produce extreme draws making plots misleading. `set.seed(1)` chosen as stable default. Documented in `dev/known-issues.md`.

---

### 2026-04-01 — Basis registry refactor (phases 0–5)

**What was built**
- `R/basis.R` — strategy registry: `register_basis()`, `get_basis()`, `list_bases()`, six dispatcher functions
- `R/zzz.R` — `.onLoad()` hook that registers all four production bases on package load
- `R/basis-uniform.R`, `R/basis-gaussian.R`, `R/basis-msplines.R`, `R/basis-epanechnikov.R` — one file per basis, each exporting a `*_spec` list conforming to the registry contract
- `R/basis-bsplines.R` — B-splines basis registered as proof-of-concept/extensibility demo
- `R/control.R`, `R/plot.R`, `R/summary.R`, `R/print.R`, `R/coef.R`, `R/residuals.R`, `R/predict.R` — file split of the original monolithic `coxph.r` (zero logic change)
- `tests/testthat/_snaps/` — golden snapshots for all four basis types established before any refactor touch

**Key commits**
- `c7817fd` phase 0: testthat safety net with golden snapshots for all 4 basis types
- `63a9c45` phase 1: split into per-concern files, zero logic change
- `9ae917d` phase 2: add registry infrastructure (basis.R + zzz.R)
- `377e467` phase 3: add basis-*.R files, register all 4 bases in .onLoad()
- `a0f7064` phase 4: wire registry dispatchers at all 9 call sites
- `82505b3` phase 5: cleanup dead code + B-splines extensibility proof
- `010cefd` phase 5: update basis file contract and cleanup internal functions

**Design decisions**
- Safety-net-first discipline: snapshots committed before any structural change so regressions are caught at the test level, not by eye
- Each phase kept to a single concern (split / infrastructure / wire-up / cleanup) so bisect remains meaningful
- Old helpers (`basis.name_mpl()`, `knots_mpl()`, `basis_mpl()`, `penalty_mpl()`, `penalty.order_mpl()`) deleted rather than aliased — no backwards-compatibility shims
- `coxph.r` kept with lowercase `.r` extension (original CRAN artefact); renaming would break `R CMD build` on case-sensitive systems and is explicitly flagged in `dev/known-issues.md`
- B-splines registered but treated as experimental; not documented in user-facing help

**Known gaps**
- [ ] Non-finite Hessian weights in interval-censored branch: `NaN`/`Inf` from underflowed `M_S1_ni1`/`M_S2_ni1` propagates into `crossprod()` at `R/coxph.r:298` and `:523` — see `dev/known-issues.md` for proposed stopgap (`ifelse(is.finite(w_ni), w_ni, 0)`) and root-cause note
- [ ] `data(lung)` fails silently; must use `survival::lung` — needs a package-level note or compatibility shim decision

---

### 2026-01-02 — Initial dev branch setup

**Key commits**
- `81bcdf7` initial upd commit

---

### 2025-12-22 — CRAN 0.2-4 baseline

**Key commits**
- `4883028` CRAN 0.2-4

> [!note] This is the last stable release. Everything above is dev-branch work not yet merged to `main`.

---

## Open TODOs

Sourced from `TODO.R`. Tick off here and in `TODO.R` when done.

**Code quality**
- [ ] Extract internal helpers for repeated computations in `R/coxph.r`: (a) wrapper that updates hazards/survival from current `beta`/`theta`, (b) wrapper for penalised log-likelihood calculation
- [ ] Fix formatting / style in `R/coxph.r`
- [x] Audit roxygen output against old package documentation — confirm all arguments, return values, and examples are present and correct *(done 2026-05-21)*
- [x] `data(lung)` fails silently — fixed: all examples now use `data(lung, package = "survival")` *(done 2026-05-21)*

**Content and tutorials**
- [x] Create pseudo-melanoma package dataset: fix seed, save as `data/melanoma.rda`, add `R/data-melanoma.R` docs, verify `data(melanoma)` works; keep generation script under `dev/data-raw/` *(done 2026-05-21)*
- [ ] Add left-censoring examples once Jun and Lo provide them — decide between extending `getting-started.Rmd` or a dedicated vignette
- [x] Add basis-function comparison tutorial using equations from `old_sources/GKS-equations.tex`; add article to `vignettes/` and `_pkgdown.yml` *(done 2026-05-21)*
- [x] Add `survival::coxph()` vs `coxph_mpl()` comparison tutorial (coefficients, SEs, ±1.96 SE intervals; consider KM/Nelson-Aalen baseline hazard comparison) *(done 2026-05-21)*
- [ ] Create draft package logo; keep source assets editable, do not commit generated experiments
- [ ] pkgdown `README.md` — needed for the pkgdown homepage

**Collaboration**
- [ ] Prepare minimal numerical-issue examples for Jun: failing seeds, minimal reproduction code, observed error, line numbers in `R/coxph.r`
- [x] Review `old_sources/LTRC_codes/` (left-truncation paper) and implement left truncation via `entry=` in `coxph_mpl()` *(done 2026-07-15)*
- [ ] Extend `residuals.coxph_mpl()`/`predict.coxph_mpl()` to account for `entry` (currently compute from time 0, ignoring truncation)

---

## Open questions

> [!question] Non-finite Hessian fix: stopgap vs root cause
> The proposed `ifelse(is.finite(w_ni), w_ni, 0)` guard suppresses the error but silently drops the contribution of interval-censored observations with extreme cumulative hazards. Should the fix address root cause (clamp cumulative hazard before computing weights, or regularise the penalty) rather than zeroing post-hoc?

> [!question] ~~LTRC (left-truncation + right-censoring) incorporation~~ *(resolved 2026-07-15)*
> Not a new censoring type: `entry=` differences the existing which=2 basis matrices (`Psi(t) -> Psi(t) - Psi(entry)`) via the same `compute_basis_matrix()` entry point already used for `t_i1`/`t_i2` — no registry changes needed. See `R/coxph.r` and `tests/testthat/test-ltrc.R`. `residuals.coxph_mpl()`/`predict.coxph_mpl()` are a documented follow-up, not yet implemented.

> [!question] ~~Pseudo-melanoma dataset as package data~~ *(resolved 2026-05-21)*
> Script placed in `dev/data-raw/melanoma.R`; `data/melanoma.rda` committed.

> [!question] Left-censoring vignette dependency on collaborator examples
> `TODO.R` item 2: left-censoring examples are blocked on Jun and Lo providing their worked examples. Should a placeholder vignette shell be committed now to hold the nav slot, or wait until content exists?

> [!question] ~~`survival::coxph()` comparison tutorial~~ *(resolved 2026-05-21)*
> Includes coefficient/SE table, scaled forest plot, and three-way baseline survival comparison (KM / Breslow / M-splines). See `vignettes/coxph-comparison.Rmd`.

---

## Related projects

| Repo / resource | Relationship |
|---|---|
| `survival` (CRAN) | Dependency; `coxph_mpl` parallels `survival::coxph()` API by design |
| Ma, Webb & Hudson (2024) book | Canonical algorithm reference; vignettes reproduce its examples |
| `old_sources/LTRC_codes/` | Left-truncation R scripts from a related paper; pending incorporation review |

---

## Quick reference

```r
devtools::load_all()          # reload after any R/ change
devtools::test()              # run testthat suite (must pass before commit)
devtools::check()             # full CRAN check
devtools::build_vignettes()   # rebuild vignettes locally
pkgdown::build_site()         # rebuild pkgdown site
devtools::document()          # regenerate NAMESPACE + man/ from roxygen
```
