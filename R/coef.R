#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Extract Coefficients from \code{coxph_mpl} Fits
#'
#' @param object An object of class \code{"coxph_mpl"} or
#'   \code{"summary.coxph_mpl"}.
#' @param parameters Parameter set of interest: \code{"Beta"} for regression
#'   coefficients or \code{"Theta"} for baseline hazard parameters. Default
#'   \code{"Beta"}.
#' @param ... Additional arguments passed to methods.
#'
#' @details For \code{summary.coxph_mpl} inputs with \code{parameters == "Theta"},
#'   only baseline hazard estimates exceeding \code{min.theta} are reported.
#' @return A vector of coefficients or a matrix with estimates, standard errors,
#'   z-statistics, and p-values.
#' @seealso [coxph_mpl()], [summary.coxph_mpl()]
#' @examples
#' \dontrun{
#' data(lung, package = "survival")
#' fit_mpl <- coxph_mpl(Surv(time, status == 2) ~ age + sex + ph.karno + wt.loss,
#'                      data = lung)
#' coef(fit_mpl)
#' coef(fit_mpl, parameters = "Theta")
#' coef(summary(fit_mpl))
#' }
#' @export
#' @rdname coef.coxph_mpl
#' @method coef summary.coxph_mpl
coef.summary.coxph_mpl=function(object, parameters = "Beta", ...) {
  object[[parameters]]
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @export
#' @rdname coef.coxph_mpl
#' @method coef coxph_mpl
coef.coxph_mpl=function(object, parameters = "Beta", ...) {
  out = object$coef[[parameters]]
  if(parameters == "Beta") names(out) = colnames(object$data$X)
  else names(out) = 1:object$dim$m
  out
}
