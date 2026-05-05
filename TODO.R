#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TODO(survivalMPL-dev):
#^ lung dataset is not loading with data(lang) but lung = survival::lung
#* Extract small internal helpers for repeated computations
#    - e.g. wrapper for updating hazards/survival from current beta/theta
#    - e.g. wrapper for penalized log-likelihood calculation
#* Fix formatting in coxph.r
#* check roxygen results and compare them with the old package documentation
#// Do README.md for the pkgdown
#
#* Planned package/data/tutorial work:
# 1. Create a pseudo-melanoma generated dataset and save it as package data so
#    examples can use data(melanoma) instead of inline simulation code.
#    - Keep the data-generation seed and script under data-raw/ or dev/.
#    - Add R/data-melanoma.R documentation and confirm data(melanoma) works.
# 2. Add left-censoring examples once Jun and Lo provide their examples.
#    - Decide whether these belong in getting-started.Rmd or a dedicated
#      left-censoring vignette.
# 3. Add a website tutorial comparing the supported basis functions.
#    - Use formulas from old_sources/GKS-equations.tex.
#    - Add the article to vignettes/ and _pkgdown.yml tutorials dropdown.
# 4. Add a tutorial comparing survival::coxph() and coxph_mpl().
#    - Compare beta coefficients, standard errors, +/- 1.96 SE intervals, and
#      other relevant fitted quantities in a table.
#    - Consider an additional comparison with Kaplan-Meier / Nelson-Aalen style
#      non-parametric hazard or cumulative hazard estimates, if statistically
#      coherent for the chosen example.
# 5. Create a draft package logo.
#    - Keep source assets editable and avoid committing generated experiments
#      unless they are selected for the package.
# 6. Prepare examples for Jun that demonstrate numerical issues.
#    - Record failing seeds, minimal reproduction code, observed failure, and
#      the corresponding line numbers in old_sources/survivalMPL_0.2-4.
# 7. Review old_sources/LTRC_codes from the left-truncation paper / Lifetime
#    Data Analysis paper.
#    - Understand the code structure and algorithmic assumptions.
#    - Produce an incorporation plan only; defer package implementation until
#      the design is agreed.



