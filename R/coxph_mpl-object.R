#' MPL Proportional Hazards Regression Object
#'
#' Objects returned by [coxph_mpl()] represent proportional hazards models
#' fitted via maximum penalised likelihood. They provide methods for
#' \code{print}, \code{summary}, \code{plot}, \code{residuals}, and
#' \code{predict}.
#'
#' All components listed below must be present in a valid
#' \code{"coxph_mpl"} object.
#'
#' @name coxph_mpl.object
#' @format A list with the following elements:
#' \describe{
#'   \item{coef}{List with regression estimates (\code{Beta}) of length \eqn{p}
#'   and baseline hazard estimates (\code{Theta}) of length \eqn{m}.}
#'   \item{se}{List with standard-error matrices for \code{Beta} and
#'   \code{Theta} across inference methods.}
#'   \item{covar}{List of covariance matrices for the available inference
#'   methods.}
#'   \item{ploglik}{Length-2 vector with penalised log-likelihood components.}
#'   \item{iter}{Length-3 vector with iteration counts for smoothing and
#'   parameter updates.}
#'   \item{knots}{List with basis parameters: \code{m}, \code{Alpha},
#'   \code{Delta}, and (for Gaussian bases) \code{Sigma}.}
#'   \item{control}{Control settings as returned by [coxph_mpl.control()].}
#'   \item{dim}{List with \code{n}, number of events, ties, number of covariates
#'   (\code{p}), and number of baseline parameters (\code{m}).}
#'   \item{call}{The matched call.}
#'   \item{data}{List with outcome times, censoring indicators, and design
#'   matrix \code{X}.}
#' }
#' @seealso [coxph_mpl()], [summary.coxph_mpl()], [coef.coxph_mpl()],
#'   [plot.coxph_mpl()], [residuals.coxph_mpl()], [predict.coxph_mpl()]
#' @keywords survival
NULL
