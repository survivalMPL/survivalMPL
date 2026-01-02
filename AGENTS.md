# survivalMPL.package Agent Guidelines

These notes are for Copilot/Codex and other automated assistants working on this
package. Keep changes small, reviewable, and consistent with CRAN practices.

## Scope and stability
- Keep the exported API stable: `coxph_mpl`, `coxph_mpl.control`, and S3 methods
  for `print`, `summary`, `plot`, `predict`, `residuals`, `coef`.
- Do not change the structure or naming of fields on the returned
  "coxph_mpl" object unless explicitly requested.
- Prefer behavior-preserving refactors; performance tweaks must be justified and
  verified with before/after comparisons of coefficients and log-likelihood.

## Documentation workflow (roxygen2)
- Edit only `R/*.r` roxygen comments; regenerate docs with
  `devtools::document()`.
- Do not hand-edit `NAMESPACE` or `man/*.Rd` (they are generated).
- Keep examples runnable and small; avoid long-running examples in docs.

## File and data hygiene
- Keep source in `R/`; add data as `.rda` in `data/` with a matching
  `R/data-*.R` roxygen block.
- Avoid committing generated site output (e.g., `pkgdown/`) unless asked.
- Keep `NEWS` updated when user-facing behavior or docs change.

## Quality checks
- If changing numerical code, compare outputs on small datasets (e.g. a
  right-censored and an interval-censored example) before and after.
- Avoid large reformatting; minimize whitespace-only changes.

If anything here conflicts with a direct user request, the user request wins.
