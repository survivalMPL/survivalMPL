# survivalMPL 0.2-4.9005 (development)
- **`entry=` no longer errors on a non-uniform basis.** Left truncation is only
  implemented for the piecewise-constant basis, so asking for any other basis
  together with `entry=` now falls back to `"uniform"` with a warning instead of
  stopping. The knot count reverts to the uniform default when it was the
  kernel-basis default, and the penalty order is recomputed for the new basis.
- **B-splines removed.** `R/basis-bsplines.R` was a proof-of-concept never
  intended as a production basis; it is gone, along with its registration in
  `zzz.R` and its mentions in the articles. The registered bases are uniform,
  M-splines, Gaussian and Epanechnikov.
- **New dataset `hiv`.** 300 simulated HIV seroconversion times with left- and
  right-censored observations and exact events, and no interval censoring - a
  clean illustration of left censoring on its own. See `dev/data-raw/hiv.R`.
- **`melanoma` regenerated.** Two thirds of the sample (200 subjects) are now
  exact events and one third is censored: 58 interval, 23 right and 19 left
  censored. Previously it was 117 exact, 35 interval, 16 right and 132 left
  censored, i.e. dominated by left censoring. The simulation design, true
  coefficients and baseline hazard are unchanged; only `pi_E`, `a1` and `a2`
  differ. Fits on this dataset will give different numbers than in 0.2-4.9004.
- New article **Control Parameters**, covering `coxph_mpl.control()`: the basis,
  the knots, the smoothing value and its automatic selection, the optimiser
  limits and tie handling.
- **The default `tol` makes the outer loop hard to satisfy.** The automatic
  smoothing parameter selection stops when the degrees of freedom $\nu$ change
  by less than `10 * tol`, i.e. $10^{-6}$ by default, where Ma, Webb and Hudson
  (2024, Section 2.10) call $\nu$ settled at changes of 0.5 to 1. On `lung` the
  default fit therefore runs out its 150 outer iterations and reports
  `Convergence: NO`, though its coefficients differ from a converged fit by at
  most a twentieth of a standard error. The examples now pass `tol = 1e-05`,
  the value the book uses throughout; the default itself is unchanged.
- New website page **Book**, with the cover, DOI and citation for Ma, Webb and
  Hudson (2024), and a map from each article to the material it follows.
- **Right Censoring** drops the `control = coxph_mpl.control(...)` block from
  its example and gains a residuals section: how to read the martingale and
  Cox-Snell plots, and why the index axis of both reflects the row order of
  `lung` rather than the fit.
- **Left Censoring** rewritten around the new `hiv` dataset, with no interval
  censoring in sight: the fit checked against the coefficients that generated
  the data, the baseline hazard, and predicted survival by sharing status.
- **Interval Censoring** now uses the bundled `melanoma` data instead of
  simulating them inline, and follows the book's own example (Section 2.11): its
  fit call unchanged, including the automatic smoothing parameter selection, a
  hazard ratio table as in its Table 2.6, and the baseline hazard of its
  Figure 2.5.
- **Comparing `coxph()` and `coxph_mpl()`** gains an instantaneous hazard
  comparison: the binned Breslow hazard against the smooth M-spline estimate.
- Baseline hazard plots added to the right-censoring, left-censoring,
  interval-censoring and left-truncation articles.
- **Left Truncation** drops its redundant `control = coxph_mpl.control(basis =
  "uniform", n.obs = ...)` blocks - `"uniform"` is the default and `n.obs` is
  not read by the fitting code - and its "What ignoring the truncation costs"
  comparison, keeping the article to the `entry=` fit, its baseline hazard and
  the basis fallback.
- **Getting Started** lists the bases in the text instead of printing
  `list_bases()`, and points at the new Control Parameters article.

# survivalMPL 0.2-4.9004 (development)
- **`plot()` legends now describe only what is drawn.** `plot.coxph_mpl()`
  listed "Observed", "Censored" and a quantile reference line in its legends,
  none of which the code draws: the points were never plotted and the quantile
  line and shaded region are commented out. The legends now show the knots and
  the line type marking $\hat\theta_u = 0$, and, for the baseline function
  panels, the estimate and its confidence band. Panels gained vertical headroom
  and the legends an opaque background so they no longer sit on the curves.
  `plot.residuals.coxph_mpl()` gained the same treatment.
- **`plot.coxph_mpl()` gains an `xlim` argument** for narrowing the displayed
  time range, which the baseline function panels need for right-skewed data
  where the full knot range compresses everything into the left edge. Only the
  display changes; the fit and the knots do not. `upper.quantile` is documented
  as having no effect, which has been true since the quantile line was disabled.
- New vignette **Left Truncation**, following the profile-likelihood treatment
  of Bhaskaran, Liquet and Ma. It shows on the `hiroshima` data that ignoring
  delayed entry biases the dose and city coefficients upwards by 1.3 and 2.1
  standard errors, and explains why `entry=` is restricted to the uniform basis.
- The tutorials are regrouped into "Censoring and truncation" and "Modelling
  choices". The home page now carries installation and one right-censored
  example, pointing at Getting Started for the rest. All seven articles share
  one output format (`html_vignette` with a table of contents) and one set of
  knitr defaults, use `Example: ...` headings, and end with a Next steps
  section.
- Documentation regenerated: `man/coxph_mpl.Rd` had not been rebuilt after the
  `entry=` basis restriction was documented.

# survivalMPL 0.2-4.9003 (development)
- **`hiroshima` dataset rebuilt, breaking change.** It is now one row per
  subject (86,611 subjects), reconstructed from the RERF LSS14 grouped
  person-time table, rather than one pseudo-subject per grouped-table cell
  (53,782 rows). The pseudo-subject version was not usable for analysis: it
  produced a *negative* colon-dose coefficient for all-cause mortality, an
  artefact of high-dose cells having more person-years and hence a longer mean
  attained age. Changes for existing code:
  - `entry` is now the attained age at the 1950-10-01 start of follow-up, the
    age at which a subject actually enters the risk set, rather than a cell
    mean age at exposure. Age at exposure is available as the new `agex`
    column (exactly `entry - 5.15`, so the two are collinear).
  - `dose` is now in **Gy**, not mGy.
  - `status` is a per-subject death indicator (50,620 deaths), not a
    "any death in this cell" indicator.
  - reproduces the source table's subjects, deaths and colon-cancer deaths
    exactly, and its person-years to within 0.33%.
- **`entry=` now requires `basis = "uniform"`.** Left truncation is implemented
  by differencing the cumulative basis, `H(t) - H(entry)`, which is only
  correct for the piecewise-constant basis. Any other basis combined with
  `entry=` previously returned a silently wrong fit and now raises an error.
- New vignettes **Right Censoring** and **Left Censoring**; left censoring was
  previously undocumented. `Getting Started` is now general information only:
  installation, the `coxph_mpl()` interface, response encodings, control
  arguments and methods, with its worked examples moved into those two articles.
  The pkgdown site groups the tutorials into sections.
- The pkgdown site is now built into `docs/` rather than `pkgdown/`. The old
  destination collided with pkgdown's own source lookup: a generated
  `pkgdown/index.md` outranked the real `index.md` and the home page rendered
  its own previous output.

# survivalMPL 0.2-4 (2025-03-23)
- changed maintainer email address
- added contributor

# survivalMPL 0.2-3 (2022-11-21)
- sorted minor typos in documentation

# survivalMPL 0.2-2 (2022-11-21)
- extended the example section of `coxph_mpl`
- amended code to allow models without covariates
- removed a duplicated item in `coxph_mpl.object`

# survivalMPL 0.2-1 (2021-11-14)
- change of maintainer (D-L Couturier)
- minor fixes (warning messages from `try` and `model.matrix`)

# survivalMPL 0.2 (2017-12-11)
- interval censoring implementation

# survivalMPL 0.1.2 (2017-10-13)
- change of maintainer (Maurizio Manuguerra)
