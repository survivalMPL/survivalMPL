% Part of the survivalMPL-package
% Distributed under GPL 2 or later
% plot.coxph_mpl 
\name{plot.coxph_mpl}
\encoding{utf8}
\alias{plot.coxph_mpl}
\title{Plot a coxph_mpl Object}
\description{
Plot the bases used to estimate the baseline hazard parameters, as well as 
the estimate and confidence interval of the baseline hazard, cumulative 
baseline hazard and baseline survival functions (plots are selectable by
\code{which}).
}
\usage{
\method{plot}{coxph_mpl}(x, se="M2QM2", ask=TRUE, which=1:4, upper.quantile=.95,\dots)
}
\arguments{
\item{x}{
an object inheriting from class \code{\link{coxph_mpl}}, representing
a fitted Cox proportional hazard model.
}
\item{se}{
an inference method (to build confidence intevals for the baseline hazard, 
cumulative baseline hazard and baseline survival functions). 
Possibilites are \code{"H"}, \code{"M2QM2"} and \code{"M2HM2"}. 
Refer to the \bold{Details} Section of \code{\link{coxph_mpl}}.  
Default is \code{se="M2QM2"}. 
}
\item{ask}{
logical. If \code{TRUE}, the user is asked to hit the \kbd{enter keyboard} 
before each plot. See \code{\link{par}(ask=.)}. Default is \code{ask=TRUE}.   
}
\item{which}{
integer vector indicating the list of wished plots. If a subset of the plots 
is required, specify a subset of the numbers \code{1:4}. By default, all
plots are provided.
}
\item{upper.quantile}{
quantile of the model response defining the upper limit of the x-axis of the
plots of the baseline hazard, cumulative baseline hazard and baseline survival 
functions. Default is \code{upper.quantile=.95}.
}
\item{\dots}{
other parameters to be passed through to plotting functions.
}
}
\details{
In the first plot, the bases corresponding to zero (or close to zero) estimates
appear in dashed line. An estimate is considered as a zero if it is smaller than 
\code{min.Theta} (See \code{\link{coxph_mpl.control}}).

Confidence intervals for the baseline hazard, cumulative baseline hazard and
baseline survival functions are obtained using the delta method.
}
\seealso{
\code{\link{coxph_mpl}},  \code{\link{coxph_mpl.control}},
\code{\link{coxph_mpl.object}}  and  \code{\link{summary.coxph_mpl}}.
}
\author{Dominique-Laurent Couturier, Maurizio Manuguerra}
\examples{
\dontrun{
data(lung)
fit_mpl <- coxph_mpl(Surv(time, status == 2) ~ age + sex + ph.karno + wt.loss, data = lung)

par(mfrow=c(2,2))
plot(fit_mpl, ask=FALSE, cex.main=.75)
}
}
\keyword{survival}
