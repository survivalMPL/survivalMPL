#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Predictions for a Cox Model Fit via MPL
#'
#' Compute predicted instantaneous risk or survival probabilities for a fitted
#' \code{coxph_mpl} model.
#'
#' @param object A fitted model of class \code{"coxph_mpl"}.
#' @param se Inference method for confidence intervals. One of \code{"H"},
#'   \code{"M2QM2"}, or \code{"M2HM2"}. Default \code{"M2QM2"}.
#' @param type Prediction type: \code{"risk"} for instantaneous risk or
#'   \code{"survival"} for survival probability. Default \code{"risk"}.
#' @param i Optional integer index of the observation whose covariates are used.
#'   Defaults to mean covariates.
#' @param time Optional numeric vector of times at which to predict. Defaults to
#'   1000 equally spaced times over the outcome range.
#' @param upper.quantile Quantile of the response used to bound the x-axis when
#'   plotting. Default \code{0.95}.
#' @param ... Additional arguments passed to methods.
#'
#' @return A data frame of class \code{"predict.coxph_mpl"} with columns
#'   \code{time}, \code{risk} or \code{survival}, \code{se}, \code{low}, and
#'   \code{high}.
#' @details Predictions incorporate the baseline hazard or cumulative baseline
#'   hazard, giving absolute (not relative) risk and survival estimates.
#'   Standard errors and confidence intervals are computed via the delta method
#'   and truncated to the parameter range.
#' @seealso [coxph_mpl()], [coxph_mpl.control()], [residuals.coxph_mpl()],
#'   [summary.coxph_mpl()]
#' @examples
#' \dontrun{
#' data(lung)
#' fit_mpl <- coxph_mpl(Surv(time, status == 2) ~ age + sex + ph.karno + wt.loss,
#'                      data = lung)
#' plot(predict(fit_mpl))
#' }
#' @export
#' @method predict coxph_mpl
predict.coxph_mpl=function(object,se="M2QM2",type="risk",i=NULL,time=NULL,upper.quantile=.95,...) {
  prob  = upper.quantile
  covar = object$covar[[se]]
  Beta  = object$coef$Beta
  Theta = object$coef$Theta
  p     = object$dim$p
  m     = object$dim$m
  if(length(i)>1){warning("only the first observation will be considered\n",call. = FALSE)}
  # time
  if(is.null(time)){
    n.x   = 1000
    V_x_X = seq(object$knots$Alpha[1],max(object$knots$Alpha),length=n.x)
  }else{	
    n.x   = length(time)
    V_x_X = time
  }
  # x
  xTB = if(is.null(i)){apply(object$data$X,2,mean)}else{object$data$X[i[1],,drop=FALSE]}
  Mu  = c(exp(xTB%*%Beta))
  # risk
  out = data.frame(time = V_x_X, mid = NA, se = NA, low = NA, high = NA)
  if(type=="risk"){
    M_psi_Xm = basis_mpl(V_x_X,object$knots,object$control$basis,object$control$order,which=1)
    out$mid  = Mu*M_psi_Xm%*%Theta
    # correction factor
    M_corr_mpm = matrix(c(rep(M_psi_Xm%*%Theta*Mu,p)*rep(Beta,each=n.x),Mu*M_psi_Xm),ncol=p+m)
    out$se     = sqrt(diag(M_corr_mpm%*%object$covar[[se]]%*%t(M_corr_mpm)))
    out$low    = out$mid - 2*out$se; out$low[out$low<0] = 0
    out$high   = out$mid + 2*out$se
    # survival
  }else{
    M_Psi_Xm = basis_mpl(V_x_X,object$knots,object$control$basis,object$control$order,which=2)
    out$mid  = exp(-Mu*M_Psi_Xm%*%Theta)
    # correction factor
    M_corr_mpm = rep(-out$mid*Mu,p+m)*matrix(c(rep(M_Psi_Xm%*%Theta,p)*rep(Beta,each=n.x),M_Psi_Xm),ncol=p+m)
    out$se     = sqrt(diag(M_corr_mpm%*%object$covar[[se]]%*%t(M_corr_mpm)))
    out$low    = out$mid - 2*out$se; out$low[out$low<0] = 0
    out$high   = out$mid + 2*out$se; out$high[out$high>1] = 1
  }
  # out
  times=c(object$data$time[,1L],object$data$time[which(object$data$censoring==3),2L])
  attributes(out)$inf = list(i=i[1], upper.quantile=prob, upper.value = quantile(times,prob),
                             max = max(times), user.time = !is.null(time), m=m, risk = type=="risk")
  colnames(out)[2]=type
  class(out) =c("predict.coxph_mpl","data.frame")
  out
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @export
#' @rdname predict.coxph_mpl
#' @method plot predict.coxph_mpl
#' @param x An object of class \code{"predict.coxph_mpl"}.
#' @param ... Additional plotting parameters.
plot.predict.coxph_mpl=function(x,...){
  inf     = attr(x,"inf")
  colw    = terrain.colors(3)[1:2]
  pos.var = x$se>0
  main    = paste(ifelse(inf$risk,"Predicted instantaneous risk at time t",
                         "Predicted probability of survival after time t"),"\n",
                  ifelse(is.null(inf$i),"for 'average' covariates",paste("for observation",inf$i)),sep="")
  plot(1,1,pch="",axes=FALSE,xlab="Time",ylim=c(0,max(x$high)),main=main,
       xlim=if(inf$user.time){c(0.5,nrow(x)+.5)}else{c(min(x$time),inf$upper.value)},ylab="",...)
  axis(2,...)
  abline(h=0,col="light gray")
  mtext(if(inf$risk){expression(hat(h)(t[i]))}else{expression(widehat(S)(t[i]))},2,padj=-2,...)
  if(inf$user.time){
    axis(1,at=1:nrow(x),labels=x$time,tick=TRUE,pos=0,...)
    for(tw in 1:nrow(x)){
      if(x$se[tw]>0){
        arrows(tw, y0=x$low[tw], y1=x$high[tw], angle=90, code=3, col=colw[2],...)
      }
    }
    points(1:nrow(x),x[,2],col=colw[1],...)
    legend("topleft",legend=c("Estimate","95% conf. interval"),pch=c(NA,15),lty=c(1,0),
           col=c(colw[1],confcol),ncol=1,bty="n",cex=.75)
  }else{
    axis(1,pos=0,...)
    rect(inf$upper.value,0,inf$max,max(x$high)*1.5,col=gray(.9),border = NA)
    if(inf$upper.quantile<1)abline(v=inf$upper.value,col=gray(.5),lty=2)
    xx = c(x$time[pos.var],x$time[pos.var][length(x$time[pos.var]):1])
    yy = c(x$high[pos.var],x$low[pos.var][length(x$time[pos.var]):1])
    confcol = paste(substr(colw[2],1,7),70,sep="")
    polygon(xx,yy,col=confcol,border = "gray",...)
    lines(x$time,x[,2],lwd=1.1,col=colw[1],...)
    legend(ifelse(inf$risk,"topleft","topright"),legend=c("Estimate","95% conf. interval",
                                                          as.expression(bquote(widehat(F)^-1*(.(inf$upper.quantile))))),pch=c(NA,15,NA),lty=c(1,0,3),
           col=c(colw[1],confcol,gray(.5)),ncol=1,bty="n",cex=.75)
  }
}
