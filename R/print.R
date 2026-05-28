#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @rdname coxph_mpl
#' @param x An object of class \code{"coxph_mpl"}.
#' @return Invisibly returns \code{x}.
#' @export
#' @method print coxph_mpl
print.coxph_mpl=function(x, ...) {
  cat("\n")
  print(x$call)
  cat("\nLog-likelihood : ",x$ploglik[1],"\n",sep="")
  cat("\nRegression parameters :\n")
  vect=c(x$coef$Beta)
  names(vect)=dimnames(x$data$X)[[2]]
  print(vect, ...)
  cat("\nBaseline hasard parameters : \n")
  print(x$coef$Theta, ...)
  cat("\n")
  invisible(x)
}
