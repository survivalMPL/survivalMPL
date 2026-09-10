#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Variance-Covariance Matrix and Confidence Intervals for \code{coxph_mpl} Fits
#'
#' Extract the estimated variance-covariance matrix of a fitted
#' \code{coxph_mpl} model, and the Wald confidence intervals that follow from
#' it.
#'
#' @param object A fitted model of class \code{"coxph_mpl"}.
#' @param parameters Parameter block of interest: \code{"Beta"} for the
#'   regression coefficients or \code{"Theta"} for the baseline hazard
#'   parameters. Default \code{"Beta"}.
#' @param se Inference method, matching the choices offered by
#'   [summary.coxph_mpl()] and [predict.coxph_mpl()]: one of \code{"M2QM2"},
#'   \code{"M2HM2"} or \code{"H"}. Default \code{"M2QM2"}.
#' @param parm Parameters to report an interval for, given as names or as
#'   indices into the selected block. Defaults to all of them.
#' @param level Confidence level. Default \code{0.95}.
#' @param ... Additional arguments passed to methods.
#'
#' @details \code{coxph_mpl} estimates the regression coefficients
#'   \eqn{\boldsymbol\beta} and the baseline hazard parameters
#'   \eqn{\boldsymbol\theta} jointly, and stores the covariance of the whole
#'   \eqn{(\boldsymbol\beta, \boldsymbol\theta)} vector. \code{vcov} returns one
#'   diagonal block of it, so \code{vcov(fit)} is the \eqn{p \times p} matrix a
#'   reader of \code{coef(fit)} expects. The cross-block covariances between
#'   \eqn{\boldsymbol\beta} and \eqn{\boldsymbol\theta} are not discarded, only
#'   unreported here; they are available as \code{fit$covar[[se]]}.
#'
#'   Rows and columns for baseline parameters that were driven to the boundary
#'   (at or below \code{min.theta}, see [coxph_mpl.control()]) are zero, as are
#'   all entries when the covariance could not be inverted; in the latter case
#'   the matrix is filled with \code{NA} and the intervals follow suit.
#'
#'   \code{confint} forms symmetric Wald intervals on the scale the parameters
#'   were estimated on, so for regression coefficients they are intervals for
#'   \eqn{\beta}, not for the hazard ratio; exponentiate to obtain the latter.
#'
#' @return \code{vcov} returns a symmetric matrix with one row and column per
#'   parameter in the selected block. \code{confint} returns a matrix with two
#'   columns giving the lower and upper limits.
#' @seealso [coxph_mpl()], [coef.coxph_mpl()], [summary.coxph_mpl()]
#' @examples
#' \dontrun{
#' data(bcos2)
#' fit <- coxph_mpl(Surv(left, right, type = "interval2") ~ treatment,
#'                  data = bcos2, basis = "msplines")
#' vcov(fit)
#' confint(fit)
#' exp(confint(fit))          # hazard ratio scale
#' }
#' @export
#' @rdname vcov.coxph_mpl
#' @method vcov coxph_mpl
vcov.coxph_mpl <- function(object, parameters = c("Beta", "Theta"),
                           se = c("M2QM2", "M2HM2", "H"), ...) {
  parameters <- match.arg(parameters)
  se <- match.arg(se)
  V <- object$covar[[se]]
  if (is.null(V)) {
    stop(gettextf("no covariance matrix %s in this fit", dQuote(se)),
         domain = NA, call. = FALSE)
  }
  p <- object$dim$p
  m <- object$dim$m
  idx <- if (parameters == "Beta") seq_len(p) else p + seq_len(m)
  out <- V[idx, idx, drop = FALSE]
  nm <- if (parameters == "Beta") {
    colnames(object$data$X)
  } else {
    as.character(seq_len(m))
  }
  if (!is.null(nm)) dimnames(out) <- list(nm, nm)
  out
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @export
#' @rdname vcov.coxph_mpl
#' @method confint coxph_mpl
#' @importFrom stats qnorm
confint.coxph_mpl <- function(object, parm, level = 0.95,
                              parameters = c("Beta", "Theta"),
                              se = c("M2QM2", "M2HM2", "H"), ...) {
  parameters <- match.arg(parameters)
  se <- match.arg(se)
  cf <- coef(object, parameters = parameters)
  ses <- sqrt(diag(vcov(object, parameters = parameters, se = se)))
  if (missing(parm)) {
    parm <- seq_along(cf)
  } else if (is.character(parm)) {
    parm <- match(parm, names(cf))
    if (anyNA(parm)) stop("'parm' names not found in the fit", call. = FALSE)
  }
  a <- (1 - level) / 2
  z <- qnorm(c(a, 1 - a))
  out <- cf[parm] + ses[parm] %o% z
  dimnames(out) <- list(names(cf)[parm],
                        sprintf("%.1f %%", 100 * c(a, 1 - a)))
  out
}
