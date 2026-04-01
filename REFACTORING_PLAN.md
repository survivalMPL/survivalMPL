# survivalMPL — Modular Basis Registry: Refactoring Plan
**Created:** 2026-04-01  
**Author:** survivalMPL-dev  
**Branch:** `refactor/modular-basis-registry`

---

## Context & Goal

The package currently uses a monolithic two-file pattern:
- `R/coxph.r` (~1,058 lines): fitting engine + all numerical internals
- `R/methods.R` (~500 lines): all 10 S3 methods bundled together

**Goal:** Introduce a **strategy/registry pattern** for the basis system so that
adding a new baseline hazard approximation (e.g. B-splines, Bernstein polynomials)
requires creating **one new file only** — zero changes to existing files.

---

## Chosen Architecture: Option B — Strategy Pattern for Basis Types

### Final File Structure

```
R/
│  # --- Registry infrastructure (MUST load first) ---
├── basis.R                   # register_basis(), get_basis(), list_bases(), dispatchers
│
│  # --- Self-registering basis implementations ---
├── basis-uniform.R
├── basis-gaussian.R
├── basis-msplines.R
├── basis-epanechnikov.R
├── basis-bsplines.R          # Phase 5 proof-of-concept
│
│  # --- Core fitting ---
├── control.R                 # coxph_mpl.control()
├── coxph.R                   # coxph_mpl() main fitting engine (slimmed)
├── optim-helpers.R           # compute_loglik(), step_halve(), update_quantities()
│
│  # --- S3 methods (one file per generic) ---
├── plot.R                    # plot.coxph_mpl()
├── summary.R                 # summary.coxph_mpl() + print.summary.coxph_mpl()
├── print.R                   # print.coxph_mpl()
├── coef.R                    # coef.coxph_mpl()
├── predict.R                 # predict.coxph_mpl()
├── residuals.R               # residuals.coxph_mpl()
│
└── zzz.R                     # .onLoad() — registers all basis types
```

> **Why flat structure?** R's build system does NOT support subfolders inside `R/`.
> File-name prefixes (`basis-`, `optim-`) provide logical grouping within the
> flat directory.

---

## Core Design: The Registry

### `basis.R` — Registry Infrastructure

```
.basis_registry (private environment)
    "msplines" ──► { knots_fn, matrix_fn, penalty_fn, label, ... }
    "m"        ──┘  (alias, points to same entry)
    "mspline"  ──┘  (alias, points to same entry)
    "uniform"  ──► { ... }
    "u"        ──┘
    ...
```

**Key functions:**
- `register_basis(name, aliases, knots_fn, matrix_fn, penalty_fn, label, ...)` — adds entry to registry
- `get_basis(name)` — resolves any alias to full entry; informative error if unknown
- `list_bases()` — returns `data.frame(name, label)` of all registered bases

**Thin dispatchers** (replace old `if/else` helpers):

| New dispatcher            | Replaces               |
|---------------------------|------------------------|
| `resolve_basis_name()`    | `basis.name_mpl()`     |
| `compute_knots()`         | `knots_mpl()`          |
| `compute_basis_matrix()`  | `basis_mpl()`          |
| `compute_penalty()`       | `penalty_mpl()`        |
| `compute_penalty_order()` | `penalty.order_mpl()`  |
| `basis_label()`           | `if/else` in plot/summary |

### Each `basis-*.R` File Contract

Every basis file must define exactly three functions and a spec list:

```r
knots_<name>  <- function(control, events)           → list(m, Alpha, Delta, ...)
matrix_<name> <- function(x, knots, order, which)    → n × m matrix
penalty_<name><- function(control, knots)             → m × m matrix

.<name>_spec <- list(
  name             = "...",
  aliases          = c(...),
  knots_fn         = knots_<name>,
  matrix_fn        = matrix_<name>,
  penalty_fn       = penalty_<name>,
  label            = "...",
  default_n_knots  = c(8, 2),
  penalty_order_fn = function(p, order) ...
)
```

### `zzz.R` — Load Hook (avoids `Collate:` ordering issues)

```r
.onLoad <- function(libname, pkgname) {
  register_basis(.msplines_spec)
  register_basis(.uniform_spec)
  register_basis(.gaussian_spec)
  register_basis(.epanechnikov_spec)
  # Future: register_basis(.bsplines_spec)
}
```

> **Why `.onLoad()` instead of `Collate:`?**  
> Alphabetically `basis-*.R` sorts before `basis.R`, so `register_basis()` would
> be called before it exists. `.onLoad()` defers registration until all functions
> exist, eliminating load-order bugs entirely.

---

## What Changes in S3 Methods

Two recurring patterns are replaced across all method files:

```r
# BEFORE — string branching repeated everywhere
if (basis == "uniform")       "Uniform step function"
else if (basis == "gaussian") "Gaussian splines"
else if (basis == "msplines") "M-splines"
else                          "Epanechikov splines"   # note: typo in original

# AFTER — single registry call
basis_label(control$basis)
```

```r
# BEFORE
basis_mpl(x, knots, basis, order, which = 1)

# AFTER
compute_basis_matrix(x, knots, basis, order, which = 1)
```

---

## Mapping: Current Code → New Files

| Current location                             | Content                        | New file              |
|----------------------------------------------|--------------------------------|-----------------------|
| `coxph.r` L1–30: `basis.name_mpl()`         | Alias resolution               | `basis.R`             |
| `coxph.r` L31–50: `penalty.order_mpl()`     | Penalty order dispatch         | `basis.R` dispatcher  |
| `coxph.r` L51–170: `knots_mpl()`            | 4 branches                     | `basis-*.R` each      |
| `coxph.r` L171–350: `basis_mpl()`           | 4 branches × 2 (ψ/Ψ)          | `basis-*.R` each      |
| `coxph.r` L351–420: `penalty_mpl()`         | 4 branches                     | `basis-*.R` each      |
| `coxph.r` L421–490: `coxph_mpl.control()`  | Control                        | `control.R`           |
| `coxph.r` L491–1058: `coxph_mpl()`         | Fitting + repeated loglik      | `coxph.R` + `optim-helpers.R` |
| `methods.R`: `plot.coxph_mpl()`             | Plot method                    | `plot.R`              |
| `methods.R`: `summary.coxph_mpl()`         | Summary method                 | `summary.R`           |
| `methods.R`: `print.summary.coxph_mpl()`   | Print summary method           | `summary.R`           |
| `methods.R`: `print.coxph_mpl()`           | Print method                   | `print.R`             |
| `methods.R`: `coef.coxph_mpl()`            | Coef method                    | `coef.R`              |
| `methods.R`: `predict.coxph_mpl()`         | Predict method                 | `predict.R`           |
| `methods.R`: `residuals.coxph_mpl()`       | Residuals method               | `residuals.R`         |

---

## Implementation Phases

### Phase 0 — Safety Net ✅ START HERE

```bash
git checkout -b refactor/modular-basis-registry
git status   # confirm clean working tree
```

Create `tests/testthat/test-baseline.R` with numerical snapshots:

```r
lung_fit <- coxph_mpl(
  Surv(time, status == 2) ~ age + sex + ph.karno + wt.loss,
  data = survival::lung,
  control = coxph_mpl.control(basis = "msplines")
)
# snapshot: coef(), loglik, dim
```

```r
devtools::test()   # RUN ONCE to CREATE golden snapshot files
```

**Checkpoint 0:** Snapshots exist. Package loads. No errors.

---

### Phase 1 — File Split, Zero Logic Change

**Actions:**
1. Extract `coxph_mpl.control()` → `R/control.R`
2. Extract repeated computation blocks → `R/optim-helpers.R`
3. Split `R/methods.R` into: `plot.R`, `summary.R`, `print.R`, `coef.R`, `predict.R`, `residuals.R`
4. Delete `R/methods.R`
5. Remove moved sections from `R/coxph.r`

**Checkpoint 1:**
```r
devtools::load_all()   # no errors
devtools::test()       # snapshots match
devtools::check()      # no new warnings
```
```bash
git commit -m "phase 1: split into per-concern files, zero logic change"
```

---

### Phase 2 — Registry Infrastructure (not connected yet)

**Actions:**
1. Create `R/basis.R` with `register_basis()`, `get_basis()`, `list_bases()`, all dispatchers
2. Create `R/zzz.R` with empty `.onLoad()`

**Checkpoint 2:**
```r
devtools::load_all()
list_bases()      # returns empty data.frame — expected
get_basis("m")   # ERROR: "Unknown basis 'm'. Available: " — expected
devtools::test()  # snapshots still match (old code paths untouched)
```
```bash
git commit -m "phase 2: registry infrastructure in place, not yet connected"
```

---

### Phase 3 — Migrate Basis Logic (one basis at a time)

**For each basis (msplines → uniform → gaussian → epanechnikov):**

1. Create `R/basis-<name>.R`:
   - Move `knots_mpl()` branch → `knots_<name>()`
   - Move `basis_mpl()` branch → `matrix_<name>()`
   - Move `penalty_mpl()` branch → `penalty_<name>()`
   - Define `.<name>_spec` list
2. Register in `R/zzz.R`
3. Add targeted test in `tests/testthat/test-basis-registry.R`

**Checkpoint 3.x (after each basis):**
```r
devtools::load_all()
get_basis("m")$label   # "M-splines"
list_bases()           # grows by one row
devtools::test()       # snapshots match (old code still used in fitting)
```
```bash
git commit -m "phase 3: all 4 basis types migrated to registry"
```

---

### Phase 4 — Connect Registry to Fitting Engine

Replace call sites one at a time (lowest risk first):

| Step | File         | Change                                              |
|------|--------------|-----------------------------------------------------|
| 4.1  | `coxph.R`    | `basis.name_mpl()` → `resolve_basis_name()`         |
| 4.2  | `coxph.R`    | `penalty.order_mpl()` → `compute_penalty_order()`   |
| 4.3  | `coxph.R`    | `knots_mpl()` → `compute_knots()`                   |
| 4.4  | `coxph.R`    | `penalty_mpl()` → `compute_penalty()`               |
| 4.5  | `coxph.R`    | `basis_mpl()` → `compute_basis_matrix()`            |
| 4.6  | `plot.R`     | `if/else` label → `basis_label()`                   |
| 4.7  | `summary.R`  | `if/else` label → `basis_label()`                   |
| 4.8  | `residuals.R`| `basis_mpl()` → `compute_basis_matrix()`            |
| 4.9  | `predict.R`  | `basis_mpl()` → `compute_basis_matrix()`            |

**Checkpoint 4 (after EACH step):**
```r
devtools::load_all()
devtools::test()   # MUST match after every single substitution
```
```bash
git commit -m "phase 4: all call sites connected to registry dispatchers"
```

---

### Phase 5 — Cleanup & Proof of Concept

**Actions:**
1. Delete now-empty internal functions from `coxph.R`:
   `basis.name_mpl()`, `penalty.order_mpl()`, `knots_mpl()`, `basis_mpl()`, `penalty_mpl()`
2. Create `R/basis-bsplines.R` (placeholder implementations)
3. Add `register_basis(.bsplines_spec)` to `zzz.R`
4. Address items from `TODO.R`:
   - Fix `data(lung)` issue → use `survival::lung` or `LazyData` in DESCRIPTION
   - Fix formatting in `coxph.R`
   - Verify roxygen documentation matches old package docs
   - Do `README.md` for pkgdown

**Final Checkpoint:**
```r
devtools::load_all()
devtools::test()              # all snapshots match
devtools::check()             # zero errors, zero warnings
list_bases()                  # 5 rows including bsplines
coxph_mpl(..., basis = "b")  # runs without error
```
```bash
git commit -m "phase 5: cleanup + bsplines proof of concept"
git tag v-modular-registry-complete
```

---

## How to Add a New Basis (Post-Refactor)

```
1. CREATE  R/basis-newtype.R    ← implement knots_fn, matrix_fn, penalty_fn
                                    define .<newtype>_spec
2. ADD     one line to zzz.R   ← register_basis(.<newtype>_spec)
3. DONE    — zero changes to any other existing file
```

---

## Known Bugs Fixed by This Refactor

- `"epa"` vs `"epanechikov"` alias mismatch in `penalty.order_mpl()` switch
  → structurally impossible after registry: each basis owns its aliases

---

## Quick Reference: Resuming After Chat Break

1. Check current git branch: `git branch`
2. Check what phase you are in: `git log --oneline -10`
3. Run `devtools::test()` to confirm current status
4. Re-read the relevant Phase section above
5. Share this file with the new chat session for full context

**To resume, tell the AI:**
> "I am continuing a refactoring project in survivalMPL. See REFACTORING_PLAN.md
> for full context. I am currently at Phase X, Step Y."
