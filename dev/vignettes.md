# Vignettes

## Current tutorials

| File | Content |
|------|---------|
| `vignettes/getting-started.Rmd` | Right-censored (`lung`) and interval-censored (`bcos2`) quick examples |
| `vignettes/interval-censoring.Rmd` | Examples 2.6 and 2.7 from the book |
| `vignettes/basis-comparison.Rmd` | MathJax definitions of all four basis types; shape visualisation using `melanoma` event times; four-way fit comparison with true Weibull baseline overlay |
| `vignettes/coxph-comparison.Rmd` | `coxph` (PL) vs `coxph_mpl` (MPL) on `lung`: coefficient/SE table with unit footnote, scaled forest plot, three-way baseline survival comparison (KM marginal, Breslow step-function, M-splines smooth) |

## Rules

- Do **not** use `book_sim.R` as a source — it contains internal modifications not in the book
- Citation key `@MaWebbHudson2024` is for the book; `@Ma2014` is for the CSDA journal article — both defined in `vignettes/references.bib`
- `basis-comparison.Rmd` setup chunk includes `pkgload::load_all()` — required because `melanoma` is a dev-only dataset not installed from CRAN

## Book reference

Ma, Webb & Hudson (2024) *Likelihood Methods in Survival Analysis with R Examples*,
Chapman & Hall/CRC. DOI: 10.1201/9781351109710.
PDF: `old_sources/Penalized_likelihood_based_survival_data_analysis_using_R_05_04_24.pdf`

## Journal reference

Ma, Jun, Héritier, Stephane & Lo, Serigne N. (2014). On the maximum penalised likelihood
approach for proportional hazard models with right censored survival data.
*Computational Statistics & Data Analysis*, 74, 142–156. DOI: 10.1016/j.csda.2014.01.005
