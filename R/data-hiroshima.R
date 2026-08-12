#' Atomic Bomb Survivor Mortality Data (Life Span Study Report 14)
#'
#' A left-truncation example derived from the Radiation Effects Research
#' Foundation's (RERF) Life Span Study (LSS) Report 14 cancer and noncancer
#' disease mortality data, 1950-2003 (Ozasa et al. 2012). One row per subject,
#' for all 86,611 subjects in the cohort.
#'
#' @details
#' \strong{Left truncation.} LSS follow-up begins on 1950-10-01, five years
#' after the 1945 bombings, so a subject enters the risk set only at the
#' attained age they had reached on that date: anyone exposed in 1945 who died
#' before follow-up began is not in the cohort at all. Analyses on the age
#' scale must therefore account for delayed entry via \code{entry} (or, in
#' \code{\link[survival]{coxph}}, \code{Surv(entry, time, status)}). Entry ages
#' range from 5.2 to 93.8 years, so the truncation is substantial and ignoring
#' it biases the estimates - see the examples.
#'
#' \strong{Reconstruction from grouped data.} The RERF source file is a grouped
#' person-time table stratified by city, sex, ground distance, Adult Health
#' Study participation, age at exposure, attained age, calendar time and colon
#' dose; it does not distribute individual subject records. Individual-level
#' records are reconstructed from it, which is possible because the table's
#' subject counts are recorded once per subject (in each subject's first
#' calendar-time cell) and so sum to the cohort size. Each death in a cell
#' becomes one record, with the attained age at death drawn uniformly within
#' that cell's five-year attained-age band; the remaining subjects in each
#' stratum are censored administratively at their attained age on 2003-12-31.
#' Within-cell variation in age at exposure and in dose is collapsed to the
#' cell mean, and the pairing of a particular death age to a particular subject
#' within a stratum is arbitrary. The reconstruction reproduces the source
#' table's 86,611 subjects, 50,620 deaths and 621 colon-cancer deaths exactly,
#' and its 3,294,282 person-years to within 0.33\%. See
#' \code{dev/data-raw/hiroshima.R} to reproduce and
#' \code{dev/experiments/hiroshima_dataset/} for the original grouped table and
#' its documentation.
#'
#' \strong{Not a reanalysis of LSS Report 14.} Ozasa et al. do not ungroup the
#' table at all: they fit Poisson excess-relative-risk models
#' (\eqn{\lambda_0(\mathrm{strata})(1 + \beta d)}) directly to the person-years,
#' with the baseline stratified on 1200 cells. For analysing that table their
#' approach is better; this reconstruction exists only because a
#' proportional-hazards fit needs one row per subject.
#'
#' How far the coefficients can be trusted has been checked against the
#' grouped-data Poisson fit, which is the correct analysis of the source table
#' (see \code{dev/experiments/hiroshima_reconstruction_checks.R} and
#' \code{dev/experiments/hiroshima_dataset/UNGROUPING.md}). Conclusions about
#' \emph{left truncation} are robust: the difference between the truncated and
#' naive fits is positive for every coefficient under all 12 reconstruction
#' variants tried (RNG seed, within-band placement rule, dose level, entry
#' offset, end of follow-up), and stable to four decimal places across seeds,
#' because it depends only on the entry-age distribution, which the source table
#' reports directly. Dose-response estimates are also recovered: fitted on this
#' dataset with a matched stratification, the all-cause ERR per Gy is 0.145 for
#' males and 0.276 for females, inside Ozasa et al.'s published intervals of
#' 0.150 (0.104-0.199) and 0.298 (0.244-0.354). What the dataset does
#' \emph{not} support is quoting a coefficient to three decimals as an LSS14
#' estimate: the same 12 variants move the absolute coefficients by 1.8 to 4.1
#' standard errors, uncertainty that no reported standard error captures.
#'
#' @format A data frame with 86,611 observations (one per subject) on 9
#' variables:
#' \describe{
#'   \item{city}{Factor: \code{"Hiroshima"} or \code{"Nagasaki"}.}
#'   \item{sex}{Factor: \code{"Male"} or \code{"Female"}.}
#'   \item{gd3}{Factor: ground distance, \code{"<3km"} or \code{"3-10km"}.}
#'   \item{agex}{Numeric age at exposure in 1945, in years. Exactly
#'   \code{entry - 5.15} for every subject (5.15 years separate the bombings
#'   from the start of follow-up), so \code{agex} and \code{entry} are
#'   collinear and must not both enter a model as covariates.}
#'   \item{entry}{Numeric left-truncation (delayed entry) time: attained age on
#'   1950-10-01, in years.}
#'   \item{time}{Numeric event/censoring time: attained age at death, or at the
#'   2003-12-31 end of follow-up for survivors, in years. Strictly greater than
#'   \code{entry} for every subject.}
#'   \item{status}{Integer event indicator: 1 if the subject died of any cause
#'   during follow-up, 0 if censored (50,620 deaths).}
#'   \item{status_colon}{Integer event indicator for the cause-specific colon
#'   cancer hazard: 1 if the subject died of colon cancer, 0 otherwise, with
#'   deaths from other causes treated as censored at their age at death (621
#'   colon-cancer deaths).}
#'   \item{dose}{Numeric DS02-weighted colon dose (gamma + 10*neutron), in
#'   \strong{Gy}. Note the unit: the RERF source column is in mGy.}
#' }
#' @source
#' Radiation Effects Research Foundation (RERF) Life Span Study Report 14:
#' Ozasa K, Shimizu Y, Suyama A, Kasagi F, Soda M, Grant EJ, Sakata R,
#' Sugiyama H, Kodama K. Studies of the mortality of atomic bomb survivors,
#' Report 14, 1950-2003: An overview of cancer and noncancer diseases.
#' \emph{Radiat Res} 2012 [March]; 177(3):229-43.
#'
#' This report makes use of data obtained from the Radiation Effects Research
#' Foundation (RERF), Hiroshima and Nagasaki, Japan. RERF is a private,
#' non-profit foundation funded by the Japanese Ministry of Health, Labour and
#' Welfare (MHLW) and the U.S. Department of Energy (DOE), the latter in part
#' through DOE Award DE-HS0000031 to the National Academy of Sciences. The
#' conclusions in this report are those of the authors and do not necessarily
#' reflect the scientific judgment of RERF or its funding agencies.
#' @examples
#' data(hiroshima)
#' \dontrun{
#' # Delayed entry is handled by entry=, which currently requires the
#' # piecewise-constant baseline hazard, basis = "uniform".
#' fit <- coxph_mpl(
#'   Surv(time, status) ~ dose + sex + city,
#'   data   = hiroshima,
#'   entry  = entry,
#'   basis  = "uniform",
#'   smooth = 0
#' )
#' summary(fit)
#'
#' # Cross-check: coxph()'s partial likelihood handles (start, stop] exactly.
#' coef(coxph(Surv(entry, time, status) ~ dose + sex + city, data = hiroshima))
#'
#' # Ignoring the truncation biases the dose and city coefficients upwards
#' # (here by 1.3 and 2.1 standard errors respectively), because entry age is
#' # associated with both. The sex coefficient is essentially unaffected, entry
#' # age being unrelated to sex.
#' coef(coxph_mpl(
#'   Surv(time, status) ~ dose + sex + city,
#'   data = hiroshima, basis = "uniform", smooth = 0
#' ), "Beta")
#' }
#' @docType data
#' @name hiroshima
#' @keywords dataset
NULL
