#' Pseudo-HIV Seroconversion Data
#'
#' Simulated left- and right-censored survival data for a cohort of people who
#' inject drugs, followed for HIV seroconversion.  Time is measured in years
#' from the start of injecting drug use.
#'
#' Each subject has a first HIV test some time after the time origin.  A
#' subject already seropositive at that test seroconverted somewhere in
#' \eqn{(0, t_R]} and is \strong{left censored}; a subject negative at the
#' first test is retested frequently, so a later seroconversion is recorded as
#' an \strong{exact event}; a subject still negative when follow-up ends is
#' \strong{right censored}.  There is no interval censoring by construction,
#' which makes this dataset a clean illustration of left censoring on its own.
#'
#' Of the 300 subjects, 96 are left censored, 117 are exact events and 87 are
#' right censored.  Pass to \code{coxph_mpl} via
#' \code{Surv(t_L, t_R, type = "interval2")}.
#'
#' @format A data frame with 300 observations on 6 variables:
#' \describe{
#'   \item{t_L}{Numeric lower bound of the observed interval, in years;
#'   \code{NA} for left-censored subjects.}
#'   \item{t_R}{Numeric upper bound of the observed interval, in years;
#'   \code{NA} for right-censored subjects.  Equal to \code{t_L} for exact
#'   events.}
#'   \item{Sharing}{Indicator: shared injecting equipment at enrolment.}
#'   \item{Prison}{Indicator: history of incarceration.}
#'   \item{Female}{Indicator: subject is female.}
#'   \item{Age_centred}{Age in decades, centred at the sample mean.}
#' }
#' @details
#' The true baseline hazard is Weibull with shape 1.5 and scale 6,
#' \eqn{H_0(t) = (t/6)^{1.5}}, and the true regression coefficients used in the
#' simulation are:
#' \tabular{lr}{
#'   \code{Sharing}     \tab  0.90 \cr
#'   \code{Prison}      \tab  0.45 \cr
#'   \code{Female}      \tab -0.25 \cr
#'   \code{Age_centred} \tab -0.15 \cr
#' }
#' Generated with \code{set.seed(11)}; see \code{dev/data-raw/hiv.R} to
#' reproduce.
#' @source
#' None: these data are entirely synthetic, with no real cohort and no
#' published analysis behind them. The design was written for this package's
#' left-censoring tutorial. The covariates are ones a study of this kind would
#' plausibly record, but the coefficient values above are chosen for
#' illustration and are not estimates from any study. Do not cite these numbers
#' as evidence about HIV seroconversion.
#' @examples
#' data(hiv)
#' \dontrun{
#' fit <- coxph_mpl(
#'   Surv(t_L, t_R, type = "interval2") ~
#'     Sharing + Prison + Female + Age_centred,
#'   data = hiv, basis = "msplines"
#' )
#' summary(fit)
#' }
#' @docType data
#' @name hiv
#' @keywords dataset
NULL
