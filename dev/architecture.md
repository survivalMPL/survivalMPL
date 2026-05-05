# Architecture — Basis Registry

The package uses a strategy/registry pattern so that adding a new basis requires
creating one new file only — zero changes to existing files.

## File layout

```
R/
├── basis.R              # Registry: register_basis(), get_basis(), list_bases(), dispatchers
├── basis-uniform.R      # Step function basis
├── basis-gaussian.R     # Gaussian basis
├── basis-msplines.R     # M-splines (most common; aliases: "m", "mspline")
├── basis-epanechnikov.R # Epanechnikov basis
├── basis-bsplines.R     # B-splines (proof-of-concept, registered)
├── control.R            # coxph_mpl.control()
├── coxph.r              # Main fitting engine  ← lowercase .r, do not rename
├── plot.R / summary.R / print.R / coef.R / predict.R / residuals.R
└── zzz.R                # .onLoad() — registers all bases
```

## Registry dispatchers

Use these everywhere. The old helpers are removed.

| Use this                  | Replaced                   |
|---------------------------|----------------------------|
| `resolve_basis_name()`    | `basis.name_mpl()`         |
| `compute_knots()`         | `knots_mpl()`              |
| `compute_basis_matrix()`  | `basis_mpl()`              |
| `compute_penalty()`       | `penalty_mpl()`            |
| `compute_penalty_order()` | `penalty.order_mpl()`      |
| `basis_label()`           | `if/else` label strings    |

## Adding a new basis

1. Create `R/basis-<name>.R` implementing:
   - `.<name>_knots_fn(control, events)` → list with knots
   - `.<name>_matrix_fn(x, knots, order, which)` → n × m matrix
   - `.<name>_penalty_fn(control, knots)` → m × m matrix
   - `.<name>_spec` list (name, aliases, fns, label, default_n_knots, penalty_order_fn)
2. Add one line to `zzz.R`: `register_basis(.<name>_spec)`
3. Done — no other file changes needed

See `old_sources/REFACTORING_PLAN.md` for the full spec contract and `penalty_order_fn` values per basis.
