#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Print Method for \code{summary.coxph_mpl}
#'
#' @rdname summary.coxph_mpl
#' @param x An object of class \code{"summary.coxph_mpl"}.
#' @param se Inference method to display. One of \code{"H"}, \code{"M2QM2"}, or
#'   \code{"M2HM2"}. Default is \code{"M2QM2"}.
#' @param ... Additional arguments passed to [base::print()].
#'
#' @return Invisibly returns \code{x}.
#' @seealso [summary.coxph_mpl()], [coxph_mpl()], [coxph_mpl.control()]
#' @export
#' @method print summary.coxph_mpl
print.summary.coxph_mpl=function(x,se="M2QM2",...) {
  inf = x$inf
  cat("\n")
  print(inf$call)
  cat("\n-----\n\n")
  cat("Cox Proportional Hazards Model Fit Using MPL","\n\n\n")
  cat("Penalized log-likelihood  :  ",inf$ploglik,"\n",sep="")
  cat(ifelse(inf$control$max.iter[1]==1,
             "Fixed smoothing value     :  ",
             "Estimated smoothing value :  "),
      inf$control$smooth,"\n",sep="")
  cat("Convergence               : ",
      if(inf$control$max.iter[1]==1){
        ifelse(inf$iter[2]<inf$control$max.iter[2],
               paste(c("Yes (",inf$iter[2]," iter.)"),collapse=""),"NO")        
      }else{
        ifelse(inf$iter[1]<inf$control$max.iter[1]&inf$iter[2]<inf$control$max.iter[3],
               paste(c("Yes (",inf$iter[1]," + ",inf$iter[2]," iter.)"),collapse=""),"NO")
      },"\n\n")
  cat("Data             : ",inf$data,"\n",sep="")    
  number=format(c(inf$dim$n,inf$dim$n.obs,inf$dim$n-inf$dim$n.obs))
  percentage=format(c(0,inf$dim$n.obs,inf$dim$n-inf$dim$n.obs)/inf$dim$n*100)
  cat("Number of obs.   : ",number[1],"\n",sep="")
  cat("Number of events : ",number[2]," (",percentage[2],"%)\n",sep="")
  cat("Number of cens.  : ",number[3]," (",percentage[3],"%)\n\n",sep="")
  cat("Regression parameters : ",deparse(inf$call[[2]]),"\n",sep="")
  printCoefmat(x$Beta, P.values=TRUE, has.Pvalue=TRUE,...)    
  cat("\nBaseline hasard parameters approximated using",
      basis_label(inf$control$basis),":\n")
  if(inf$control$basis=="uniform"){
    cat(paste(" (",inf$dim$m," equal events bins)\n",sep=""))
  }else{
    cat(paste(" (2 (min/max) + ",inf$control$n.knots[1]," quantile knots + ",inf$control$n.knots[2]," equally spaced knots",
              if(inf$control$basis=="msplines"|inf$control$basis=="epanechikov"){
                paste(" +",inf$control$order,"(order) - 2")}," = ",inf$dim$m," parameters)",sep=""),"\n")}
  if(inf$full){
    printCoefmat(x$Theta, P.values=TRUE, has.Pvalue=TRUE,...)    
  }else{print(x$Theta,...)}
  cat("\n-----\n\n")
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Summarise a \code{coxph_mpl} Object
#'
#' Extracts additional information for a fitted model and returns an object
#' suitable for printing. Baseline hazard parameters smaller than
#' \code{min.theta} are omitted unless \code{full = TRUE}.
#'
#' @param object A fitted model of class \code{"coxph_mpl"}.
#' @param se Inference method. One of \code{"H"}, \code{"M2QM2"}, or
#'   \code{"M2HM2"}. Default is \code{"M2QM2"}.
#' @param full Logical; if \code{TRUE}, include inference for baseline hazard
#'   parameters. Default \code{FALSE}.
#' @param ... Additional arguments passed to methods.
#'
#' @return An object of class \code{"summary.coxph_mpl"} with components:
#'   \item{Beta}{Matrix of regression estimates, standard errors, z-statistics,
#'   and p-values.}
#'   \item{Theta}{Baseline hazard estimates (or, if \code{full = TRUE}, a matrix
#'   of estimates with standard errors, z-statistics, and p-values).}
#'   \item{inf}{List with convergence details, penalised likelihood value, and
#'   control settings.}
#' @seealso [coxph_mpl()], [coxph_mpl.control()], [plot.coxph_mpl()]
#' @importFrom stats printCoefmat
#' @examples
#' \dontrun{
#' data(lung)
#' fit_mpl <- coxph_mpl(Surv(time, status == 2) ~ age + sex + ph.karno + wt.loss,
#'                      data = lung)
#' summary(fit_mpl, full = TRUE)
#' summary(fit_mpl, se = "M2HM2")
#' }
#' @export
#' @method summary coxph_mpl
summary.coxph_mpl=function(object,se="M2QM2",full=FALSE,...) {
  col.names = c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
  seB   = object$se$Beta[[se]]
  matxB = cbind(object$coef$Beta,seB,object$coef$Beta/seB,2*(1-pnorm(abs(object$coef$Beta/seB))))
  dimnames(matxB)=list(colnames(object$data$X),col.names)
  pos = object$coef$Theta>object$control$min.theta
  if(full){
    seT   = object$se$Theta[[se]]
    matxT = cbind(object$coef$Theta,seT,object$coef$Theta/seT,2*(1-pnorm(abs(object$coef$Theta/seT))))
    dimnames(matxT)=list(format(1:object$dim$m,just="right"),col.names)
    matxT = matxT[pos,]
  }else{    matxT = object$coef$Theta[pos]
  names(matxT) = seq(1,object$dim$m)[pos]}
  out = list(Beta = matxB, Theta = matxT, inf = list(iter = object$iter, call=object$call, data = object$data$name,
                                                     control = object$control, dim = object$dim, trunc = trunc, full = full, ploglik = object$ploglik[1]))
  class(out) = "summary.coxph_mpl"
  out
}
