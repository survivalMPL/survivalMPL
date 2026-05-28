% Part of the survivalMPL-package
% Distributed under GPL 2 or later
% predict.coxph_mpl 
\name{predict.coxph_mpl}
\encoding{utf8}
\alias{predict.coxph_mpl}
\alias{plot.predict.coxph_mpl}
\title{Predictions for a Cox model}
\description{
Compute (and plot) predictions of instantaneous risk and of probability of survival 
for a model fitted by \code{\link{coxph_mpl}}.
}
\usage{
\method{predict}{coxph_mpl}(object,se="M2QM2",type="risk",i=NULL,time=NULL,upper.quantile=.95,\dots)

\method{plot}{predict.coxph_mpl}(x,\dots)
}
\arguments{
\item{object}{
an object inheriting from class \code{\link{coxph_mpl}}, representing
a fitted Cox proportional hazard model.
}
\item{se}{
a character string indicating a method to build confidence intevals for the predictions. 
Possibilites are \code{"H"}, \code{"M2QM2"} and \code{"M2HM2"}. 
Refer to the \bold{Details} Section of \code{\link{coxph_mpl}}.  
Default is \code{se="M2QM2"}. 
}
\item{type}{
character string indicating the type of wished predictions. Possibilies are \code{"risk"},
for predicted instantaneous risk at time \eqn{t} (see argument \code{"time"}); and 
\code{"survival"}, for predicted probability of survival above time \eqn{t}.
Default is \code{type="risk"}.
}
\item{i}{
an integer indicating the covariate vector of interest (i.e., line of the 
\bold{X} matrix). If \code{i==NULL}, the mean of each covariate will be used.  
Default is \code{i=NULL}.   
}
\item{time}{
a double-precision vector indicating at which time the preditions should be computed.
If \code{time==NULL}, predictions are computed for 1000 equally spaced times 
in the range of the outcomes. Default is \code{time=NULL}.
}
\item{upper.quantile}{
quantile of the model response defining the upper limit of the x-axis of the
plot of the predictions. This argument is passed through to 
\code{\link{plot.predict.coxph_mpl}}. Default is \code{upper.quantile=.95}.
}
\item{x}{
an object inheriting from class \code{predict.coxph_mpl}, representing
the prediction of a fitted Cox proportional hazard model.
}
\item{\dots}{
other parameters to be passed through to printing or plotting functions.
}
}
\details{
The available predictions incorporate the baseline hazard  
(instantaneous risk) or cumulated baseline hazard estimate (survival function)
and are thus absolute instead of relative (see \code{\link[survival:predict.coxph]{predict.coxph}}).

Prediction standard errors and confidence intervals are obtained by use of the delta 
method. 

In the plots, the confidence intervals are forced to belong to the parameter 
range, which is \eqn{[0,inf[} for instantaneous risk, and \eqn{[0,1]} 
for survival probabilities.
}
\value{
a data.frame of class \code{predict.coxph_mpl} with following columns: 
\code{'time'}, the prediction time (as defined in argument \code{'time'}); 
\code{'risk'} or \code{'survival'}, the wished predictions; \code{'se'}, 
the standard error of each prediction; \code{'lower'} and \code{'upper'}, 
the lower and upper bound of the prediction confidence interval. 
}
\seealso{
\code{\link{coxph_mpl}}, \code{\link{coxph_mpl.control}},
\code{\link{coxph_mpl.object}}, \code{\link{residuals.coxph_mpl}} and
\code{\link{summary.coxph_mpl}}.
}
\author{Dominique-Laurent Couturier, Maurizio Manuguerra}
\examples{
\dontrun{
data(lung)
fit_mpl <- coxph_mpl(Surv(time, status == 2) ~ age + sex + ph.karno + wt.loss, data = lung)

plot(predict(fit_mpl))
}
}
\keyword{survival}
