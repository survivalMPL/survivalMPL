#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Residuals for a Cox Model Fit via MPL
#'
#' Compute martingale and Cox-Snell residuals for a \code{coxph_mpl} model. The
#' returned object has a plot method.
#'
#' @param object A fitted model of class \code{"coxph_mpl"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame of class \code{"residuals.coxph_mpl"} with columns
#'   \code{time1}, \code{time2}, \code{censoring}, \code{coxsnell}, and
#'   \code{martingale}.
#' @references Farrington (2000), Collett (2003), Moeschberger (2003).
#' @seealso [coxph_mpl()], [predict.coxph_mpl()], [summary.coxph_mpl()]
#' @examples
#' \dontrun{
#' data(lung)
#' fit_mpl <- coxph_mpl(Surv(time, status == 2) ~ age + sex + ph.karno + wt.loss,
#'                      data = lung)
#' par(mfrow = c(1, 2))
#' plot(residuals(fit_mpl), which = 1:2, ask = FALSE)
#' }
#' @export
#' @method residuals coxph_mpl
residuals.coxph_mpl=function(object,...) {
  control  = object$control
  out = as.data.frame(matrix(NA,object$dim$n,5))
  colnames(out)  = c("time1", "time2","censoring","coxsnell","martingale")
  out$time1      = object$data$time[,1L]
  out$time2      = object$data$time[,2L]
  out$censoring  = object$data$censoring
  i.r   = which(object$data$censoring==0)
  i.obs = which(object$data$censoring==1)
  i.l   = which(object$data$censoring==2)
  i.ic  = which(object$data$censoring==3)
  M_Psi_Xm1 = basis_mpl(out$time1,object$knots,control$basis,control$order,which=2)
  S1 = exp(-exp(object$data$X%*%object$coef$Beta)*M_Psi_Xm1%*%object$coef$Theta)
  r = rep(NA, object$dim$n)
  if(length(i.r)>0) r[i.r]   = 1 - log(S1[i.r])
  if(length(i.obs)>0) r[i.obs] = -log(S1[i.obs])
  if(length(i.l)>0) r[i.l]   = (1-S1[i.l]*(1-log(S1[i.l])))/(1-S1[i.l])
  if(length(i.ic)>0) {
    M_Psi_Xm2 = basis_mpl(out$time2[i.ic],object$knots,control$basis,control$order,which=2)
    S2 = exp(-exp(object$data$X[i.ic,,drop=F]%*%object$coef$Beta)*M_Psi_Xm2%*%object$coef$Theta)
    r[i.ic]  = (S1[i.ic]*(1-log(S1[i.ic])) - S2*(1-log(S2)))/(S1[i.ic] - S2)
  }

  out$coxsnell   = r
  out$martingale = 1 - out$coxsnell
  class(out) =c("residuals.coxph_mpl","data.frame")
  out
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @export
#' @rdname residuals.coxph_mpl
#' @method plot residuals.coxph_mpl
#' @param x An object of class \code{"residuals.coxph_mpl"}.
#' @param ask Logical; whether to prompt before each plot. Default \code{TRUE}.
#' @param which Integer vector selecting residual plots (\code{1:2}).
#' @param upper.quantile Quantile used to bound the y-axis for Cox-Snell
#'   residuals when \code{which == 3}. Default \code{0.95}.
#' @param ... Additional plotting parameters.
plot.residuals.coxph_mpl=function(x,ask=TRUE,which=1:2,upper.quantile=.95,...){
  prob = upper.quantile
  which.plot=rep(TRUE,2)
  if(!is.null(which)){which.plot[-which]=FALSE}
  if(sum(which.plot)==1){ask=FALSE}
  if(ask){oask <- devAskNewPage(TRUE)
  on.exit(devAskNewPage(oask))
  }
  # martingale
  if(which.plot[1]){
    plot(1:nrow(x),x$martingale,col=(!(x$censoring==1))+1,ylab="",
         xlab="Index",main="Martingale Residuals",axes=FALSE,
         ylim=c(min(x$martingale),1+(1-min(x$martingale))*.075),...)
    mtext(expression(delta[i]-plain(e)^(x[i]^T*hat(beta))*widehat(H)[0](t[i])),2,padj=-2,...)
    abline(h=1,col="light gray")
    abline(h=0,col="light gray",lty=2)
    axis(1,...)
    axis(2,...)
    legend("top",ncol=2,legend=c("observed","censored"),col=c(1,2),pch=1,cex=.75)
  }
  # coxsnell
  if(which.plot[2]){
    plot(1:nrow(x),x$coxsnell,col=(!(x$censoring==1))+1,ylab="",
         xlab="Index",main="Adjusted Cox & Snell Residuals",axes=FALSE,
         ylim=c(0,max(x$coxsnell)*1.075),...)
    mtext(expression(plain(e)^(x[i]^T*hat(beta))*widehat(H)[0](t[i])),2,padj=-2,...)
    abline(h=1,col="light gray",lty=2)
    axis(1,pos=0,...)
    axis(2,...)
    legend("top",ncol=2,legend=c("observed","censored"),col=c(1,2),pch=1,cex=.75)
  }
}
