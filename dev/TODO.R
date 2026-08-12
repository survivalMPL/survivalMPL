#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TODO(survivalMPL-dev):
#// lung dataset is not loading with data(lang) but lung = survival::lung
#* Extract small internal helpers for repeated computations
#    - e.g. wrapper for updating hazards/survival from current beta/theta
#    - e.g. wrapper for penalized log-likelihood calculation
#* Fix formatting in coxph.r
#* check roxygen results and compare them with the old package documentation
#// Do README.md for the pkgdown
#
#* Planned package/data/tutorial work:
# // 1. Create a pseudo-melanoma generated dataset and save it as package data so
# //   examples can use data(melanoma) instead of inline simulation code.
# //   - Keep the data-generation seed and script under data-raw/ or dev/.
# //   - Add R/data-melanoma.R documentation and confirm data(melanoma) works.
# // 2. Add left-censoring examples once Jun and Lo provide their examples.
# //   - Decide whether these belong in getting-started.Rmd or a dedicated
# //     left-censoring vignette.
# //3. Add a website tutorial comparing the supported basis functions.
# //   - Use formulas from old_sources/GKS-equations.tex.
# //  - Add the article to vignettes/ and _pkgdown.yml tutorials dropdown.
# //4. Add a tutorial comparing survival::coxph() and coxph_mpl().
# //  - Compare beta coefficients, standard errors, +/- 1.96 SE intervals, and
# //    other relevant fitted quantities in a table.
# //   - Consider an additional comparison with Kaplan-Meier / Nelson-Aalen style
# //     non-parametric hazard or cumulative hazard estimates, if statistically
# //     coherent for the chosen example.
# 5. Create a draft package logo.
#    - Keep source assets editable and avoid committing generated experiments
#      unless they are selected for the package.
# 6. Prepare examples for Jun that demonstrate numerical issues.
#    - Record failing seeds, minimal reproduction code, observed failure, and
#      the corresponding line numbers in old_sources/survivalMPL_0.2-4.
#// 7. [DONE] Left truncation (entry=) added to coxph_mpl() in R/coxph.r, based
# //   on review of old_sources/LTRC_codes (left-truncation paper). Implements
# //   H*(t) = H(t) - H(entry) by differencing the which=2 basis matrices;
# //   density matrix (which=1) untouched; knot range extended to cover entry
# //   times. Verified against dev/experiments/old_sources/LTRC_codes/ reference
# //   implementation in tests/testthat/test-ltrc.R (loose tolerance — the two
# //   estimators use different knot/bin placement schemes); smoke-tested in
# //   tests/testthat/test-smoke.R.
# //   - Follow-up (not yet done): residuals.coxph_mpl() and predict.coxph_mpl()
# //     do not account for entry (they compute from time 0, not from each
#  //    subject's entry time).
# 8. Add Observer pattern for iteration callbacks in coxph_mpl().
#    - Add a `callbacks = NULL` parameter to coxph_mpl() — fully backward-compatible.
#    - Define a callback contract: on_iter(iter, ploglik, delta) called after each
#      inner iteration; on_converge(iter, ploglik) called on exit.
#    - Implementation lives inside coxph.r; no changes to other files.
#    - Enables: progress bars, custom stopping criteria, convergence logging.



