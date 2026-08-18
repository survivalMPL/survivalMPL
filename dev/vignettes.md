# Vignettes

## Current tutorials

| File | Content |
|------|---------|
| `vignettes/getting-started.Rmd` | The `coxph_mpl()` interface, how each censoring scheme is encoded, the bases (listed in markdown, not via `list_bases()`), and pointers to the other articles |
| `vignettes/right-censoring.Rmd` | `lung` fitted with defaults plus `tol = 1e-5` (no `control` block), baseline hazard, predicted survival, and a short residuals section — martingale and Cox-Snell plots and how to read them. Deliberately stops there: the PH diagnostics of the book's Section 5.2.1 (log-log survival, `tt()`, `cox.zph()`) were written and then removed to keep the article simple |
| `vignettes/left-censoring.Rmd` | `hiv` — left- and right-censored only, no interval censoring; fit vs true coefficients, baseline hazard, predicted survival by sharing status. A section quantifying the bias from imputing at `t_R` or dropping left-censored subjects was written and then removed to keep the article simple |
| `vignettes/interval-censoring.Rmd` | Bundled `melanoma` (no inline simulation); the book's own fit call (`basis = "m"`, `tol = 1e-5`, `n.knots = c(7, 0)`, automatic smoothing), hazard-ratio table as in its Table 2.6, baseline hazard as in its Figure 2.5 |
| `vignettes/left-truncation.Rmd` | `hiroshima` with `entry=` (defaults otherwise), baseline hazard, and the non-uniform basis fallback warning. The "What ignoring the truncation costs" comparison against a fit without `entry=` was written and then removed to keep the article simple |
| `vignettes/control-parameters.Rmd` | `coxph_mpl.control()` argument by argument: basis, knots, `smooth` (REML / `0` / fixed), optimiser limits, ties |
| `vignettes/basis-comparison.Rmd` | MathJax definitions of all four basis types; shape visualisation using `melanoma` event times; four-way fit comparison |
| `vignettes/coxph-comparison.Rmd` | `coxph` (PL) vs `coxph_mpl` (MPL) on `lung`: coefficient/SE table, forest plot, baseline survival comparison, binned-Breslow vs M-spline instantaneous hazard |
| `vignettes/articles/book.Rmd` | Website-only page for Ma, Webb & Hudson (2024): cover, DOI, citation, article-to-chapter map |

## Rules

- `book_sim.R` is a **debug harness**, not a source: it wraps the book's examples in
  seed loops, `tictoc` timing and old-vs-new comparisons. Its `coxph_mpl()` calls do
  match the book (checked for the melanoma example against Section 2.11 of the PDF),
  but always confirm against the PDF before copying anything from it
- Citation key `@MaWebbHudson2024` is for the book; `@Ma2014` is for the CSDA journal article; `@BhaskaranLiquetMa` is the left-truncation paper — all defined in `vignettes/references.bib`. Keep the `.bib` free of entries nothing cites
- `basis-comparison.Rmd` setup chunk includes `pkgload::load_all()` — required because `melanoma` is a dev-only dataset not installed from CRAN
- Fits that use `ties = "epsilon"` (the default) jitter duplicated event times, so any vignette whose output depends on the knots needs `set.seed()` in the setup chunk
- `vignettes/articles/` is website-only and is listed in `.Rbuildignore`
- Vignette build cost is dominated by the outer iterations of the automatic smoothing selection. `msplines` with automatic smoothing on `melanoma` takes ~100 s and on `hiv` with default knots ~225 s; `n.knots = c(5, 0)` on `hiv` takes 7 s, and `smooth = 0` on `melanoma` is instant for the same estimates

## Book reference

Ma, Webb & Hudson (2024) *Likelihood Methods in Survival Analysis with R Examples*,
Chapman & Hall/CRC. DOI: 10.1201/9781351109710.
PDF: `old_sources/Penalized_likelihood_based_survival_data_analysis_using_R_05_04_24.pdf`

## Journal reference

Ma, Jun, Héritier, Stephane & Lo, Serigne N. (2014). On the maximum penalised likelihood
approach for proportional hazard models with right censored survival data.
*Computational Statistics & Data Analysis*, 74, 142–156. DOI: 10.1016/j.csda.2014.01.005
