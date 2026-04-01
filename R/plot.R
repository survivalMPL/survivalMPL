#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Plot a \code{coxph_mpl} Object
#'
#' Plot the bases used to estimate the baseline hazard along with the estimated
#' baseline hazard, cumulative baseline hazard, and baseline survival functions.
#' Each plot can be toggled with \code{which}.
#'
#' @param x A fitted model of class \code{"coxph_mpl"}.
#' @param se Inference method for confidence intervals. One of \code{"H"},
#'   \code{"M2QM2"}, or \code{"M2HM2"}. Default is \code{"M2QM2"}.
#' @param ask Logical; whether to prompt before each plot. See
#'   [graphics::par()]. Default \code{TRUE}.
#' @param which Integer vector selecting plots to produce (subset of \code{1:4}).
#' @param upper.quantile Quantile of the response used to set the upper x-axis
#'   limit for baseline function plots. Default \code{0.95}.
#' @param ... Additional arguments passed to plotting functions.
#'
#' @details Bases whose estimates are near zero (below \code{min.theta} in
#'   [coxph_mpl.control()]) are drawn with dashed lines. Confidence intervals
#'   are obtained via the delta method.
#' @seealso [coxph_mpl()], [coxph_mpl.control()], [coxph_mpl.object()],
#'   [summary.coxph_mpl()]
#' @examples
#' \dontrun{
#' data(lung)
#' fit_mpl <- coxph_mpl(Surv(time, status == 2) ~ age + sex + ph.karno + wt.loss,
#'                      data = lung)
#' par(mfrow = c(2, 2))
#' plot(fit_mpl, ask = FALSE, cex.main = 0.75)
#' }
#' @export
#' @method plot coxph_mpl
#' @importFrom grDevices devAskNewPage gray terrain.colors
#' @importFrom graphics abline arrows axis legend lines mtext plot points polygon rect
plot.coxph_mpl=function(x,se="M2QM2",ask=TRUE,which=1:4,upper.quantile=.95,...){
  which.plot=rep(TRUE,4)
  if(!is.null(which)){which.plot[-which]=FALSE}
  if(sum(which.plot)==1){ask=FALSE}    
  if(ask){oask <- devAskNewPage(TRUE)
  on.exit(devAskNewPage(oask))
  }
  control = x$control
  knots   = x$knots
  pos     = x$coef$Theta<x$control$min.theta
  n.x       = 1000
  V_x_X     = seq(x$knots$Alpha[1],max(x$knots$Alpha),length=n.x)    
  colw      = terrain.colors(x$dim$m+1)
  prob      = upper.quantile
  #quant     = quantile(x$data$time,prob=prob)
  M_psi_Xm  = basis_mpl(V_x_X,knots,control$basis,control$order,which=1)    
  if(which.plot[1]){
    plot(1,1,pch="",xlim=range(V_x_X),ylim=max(M_psi_Xm)*c(-.05,1),axes=FALSE,
         xlab="Survival time",ylab=expression(psi[u]^{o}*(t)),
         main=paste(if(control$basis=="uniform"){"Uniform"}else{
           if(control$basis=="gaussian"){"Gaussian"}else{
             if(control$basis=="msplines"){"M-spline"}else{
               if(control$basis=="epanechikov"){"Epanechikov"}}}},
           " bases used to approximate the baseline hazard\n",
           "(",x$dim$m," bases)",sep=""))
    abline(v=knots$Alpha,col=gray(.9))
    #abline(v=quant,col=gray(.5),lty=2)
    for(u in 1:x$dim$m){lines(V_x_X,M_psi_Xm[,u],col=colw[u],lty=2*pos[u]+1)}
    axis(2,las=2)
    axis(1)
    if(any(pos)){
      legend("topright",legend=c("Observed","Censored","Knots",as.expression(bquote(widehat(F)^-1*(.(prob)))),
                                 expression(hat(theta)[u]==0),expression(hat(theta)[u]>0)),
             pch=c(1,1,NA,NA,NA,NA),lty=c(0,0,1,2,3,1),col=c("black","red",gray(.9),gray(.5),colw[1],colw[1]),ncol=1,bty="n",cex=.75)
    }else{
      legend("topright",legend=c("Observed","Censored","Knots",as.expression(bquote(widehat(F)^-1*(.(prob))))),
             pch=c(1,1,NA,NA),lty=c(0,0,1,2),col=c("black","red",gray(.9),gray(.5)),ncol=1,bty="n",cex=.75)            
    }
  }
  if(any(which.plot[2:4])){
    cov_ThetaTheta = x$covar[[se]][-c(1:x$dim$p),-c(1:x$dim$p)]    
    xlim = range(x$knots$Alpha)
    plot_bh <- function(j,V_x_X,Theta,covar,control,knots,pos,prob,xlim,...){
      M_Ppsi_Xm   = basis_mpl(V_x_X,knots,control$basis,control$order,which=as.numeric(j>1)+1)
      V_sd2.Hh0.X = diag(M_Ppsi_Xm[,]%*%covar[,]%*%t(M_Ppsi_Xm[,]))
      pos.var     = V_sd2.Hh0.X>0
      V_Hh0_X     = c(M_Ppsi_Xm%*%matrix(Theta,ncol=1))[pos.var]
      V_x_X       = V_x_X[pos.var]
      V_sd.Hh0.X  = sqrt(V_sd2.Hh0.X[pos.var])
      upper       = V_Hh0_X+2*V_sd.Hh0.X
      lower       = V_Hh0_X-2*V_sd.Hh0.X
      lower[lower<0] = 0
      if(j==3){
        V_Hh0_X = exp(-V_Hh0_X)
        lower   = exp(-lower)
        upper   = exp(-upper)
      }
      plot(1,1,pch="",xlim=xlim,axes=FALSE,
           ylim=c(ifelse(j<3,0,min(upper[V_x_X<xlim[2]])),ifelse(j<3,max(upper[V_x_X<xlim[2]]),max(lower[V_x_X<xlim[2]]))),
           main=paste("Estimate of the",c(" baseline hazard"," cumulative baseline hazard"," baseline survival")[j]," function",sep=""),
           xlab="Survival time",ylab=c(expression(h[0]*(t)),expression(H[0]*(t)),expression(S[0]*(t)))[j])
      #rect(quant,0,max(V_x_X),max(upper)*1.5,col=gray(.9),border = NA)
      #abline(v=quant,col=gray(.5),lty=2)            
      axis(2,pos=0,las=2)
      axis(1,pos=0)
      xx = c(V_x_X,V_x_X[length(V_x_X):1])
      yy = c(upper,lower[length(V_x_X):1])
      confcol = paste(substr(colw[length(Theta)],1,7),70,sep="")
      polygon(xx,yy,col=confcol,border = "gray")
      lines(V_x_X,V_Hh0_X,lwd=1.1,col=colw[1])
      legend(ifelse(j<3,"topleft","topright"),legend=c("Estimate","95% conf. interval",
                                                       as.expression(bquote(widehat(F)^-1*(.(prob))))),pch=c(NA,15,NA),lty=c(1,0,3),
             col=c(colw[1],confcol,gray(.5)),ncol=1,bty="n",cex=.75)        
    }
    if(which.plot[2]){plot_bh(1,V_x_X,x$coef$Theta,cov_ThetaTheta,control,knots,pos,prob,xlim)}
    if(which.plot[3]){plot_bh(2,V_x_X,x$coef$Theta,cov_ThetaTheta,control,knots,pos,prob,xlim,...)}    
    if(which.plot[4]){plot_bh(3,V_x_X,x$coef$Theta,cov_ThetaTheta,control,knots,pos,prob,xlim,...)}                
  }
}
