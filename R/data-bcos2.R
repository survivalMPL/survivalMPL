#' Breast Cosmesis Data
#'
#' Interval-censored breast cosmesis data from Finkelstein and Wolfe (1985),
#' discussed by Moore (2016, example 12.2) and available in the \pkg{interval}
#' package. Compared to that version, \code{bcos2} recodes the lower bound of
#' left-censored data (\code{NA} instead of \code{0}) and the upper bound of
#' right-censored data (\code{NA} instead of \code{Inf}) to ease identification
#' via \code{Surv(type = "interval2")}.
#'
#' @format A data frame with 94 observations on 3 variables:
#' \describe{
#'   \item{left}{Numeric lower bound.}
#'   \item{right}{Numeric upper bound.}
#'   \item{treatment}{Factor with levels \code{Rad} and \code{RadChem}.}
#' }
#' @source
#' Finkelstein, D.M., and Wolfe, R.A. (1985). A semiparametric model for
#' regression analysis of interval-censored failure time data. \emph{Biometrics}
#' 41: 731-740.
#'
#' Moore, D.K. (2016). \emph{Applied Survival Analysis Using R}. Springer.
#' @examples
#' data(bcos2)
#' @docType data
#' @name bcos2
#' @keywords dataset
NULL
