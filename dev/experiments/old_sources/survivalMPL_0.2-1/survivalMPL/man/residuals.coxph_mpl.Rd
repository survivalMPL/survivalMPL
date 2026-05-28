% Part of the survivalMPL-package
% Distributed under GPL 2 or later
% residuals.coxph_mpl 
\name{residuals.coxph_mpl}
\encoding{utf8}
\alias{residuals.coxph_mpl}
\alias{plot.residuals.coxph_mpl}
\title{Residuals for a Cox model}
\description{
Compute martingale and Cox and Snell residuals for a model fitted 
by \code{\link{coxph_mpl}}. Return objects are of class 
\code{residuals.coxph_mpl} and have methods for \code{plot}.
}
\usage{
\method{residuals}{coxph_mpl}(object, \dots)

\method{plot}{residuals.coxph_mpl}(x, ask=TRUE, which=1:2, upper.quantile=.95, \dots)
}
\arguments{
\item{object}{
an object inheriting from class \code{\link{coxph_mpl}}, representing
a fitted Cox proportional hazard model. 
}
\item{x}{
an object inheriting from class \code{\link{residuals.coxph_mpl}}, representing
the residuals of a Cox proportional hazard model fit with \code{\link{coxph_mpl}}.
}
\item{ask}{
logical. If \code{TRUE}, the user is asked to hit the \kbd{enter keyboard} 
before each plot. See \code{\link{par}(ask=.)}. Default is \code{ask=TRUE}.   
}
\item{which}{
integer vector indicating the list of wished plots. If a subset of the plots 
is required, specify a subset of the numbers \code{1:2}. See \bold{Details}. 
By default, all plots are provided.
}
\item{upper.quantile}{
quantile of the Cox and Snell residuals used when \code{which==3}. 
Default is \code{upper.quantile=.95}.
}
\item{\dots}{
other parameters to be passed through to plotting or printing functions.
}
}
\details{
Refer to Collet (2003, Chapter 4) for a review of model check in the Cox 
regression model, and specifically to Farrington (2000) for an overview 
on residuals with interval-censored survival data.

For object of class \code{residuals.coxph_mpl}, the available residual plots
are, respectively, the martingale residual plot (\code{which==1}) and
the Cox and Snell residual plot (\code{which==2}).

}
\value{
A data.frame of class \code{residuals.coxph_mpl} of \eqn{n} rows 
with following columns: 
\code{'time1'}, the model outcome (with a random noise
added to event ties if \code{ties=='epsilon'} in \code{\link{coxph_mpl.control}}); 
\code{'time2'}, ending time of the interval for interval censored data only (unused otherwise);
\code{'censoring'}, the status indicator as in the \code{Surv()} function, i.e. 0=right censored, 1=event at time, 2=left censored, 3=interval censored; 
\code{'coxsnell'}, the Cox and Snell residuals; 
\code{'martingale'}, the martingale residuals.
}
\seealso{
\code{\link{coxph_mpl}}, \code{\link{coxph_mpl.control}},
\code{\link{coxph_mpl.object}}, \code{\link{predict.coxph_mpl}} and
\code{\link{summary.coxph_mpl}}.
}
\references{
Farrington C.P. (2000),
\emph{Residuals for Proportional Hazard Models with Interval-Censored Data}, 
Biometrics 56, 473-482.  

Collett, D. (2003), and Moeschberger, M. L. (2003),
\emph{Modelling Survival Data in Medical Research}, Chapman and All.  
}
\author{Dominique-Laurent Couturier, Maurizio Manuguerra}
\examples{
\dontrun{
### lung data of the survival package (see ?lung)
data(lung)
fit_mpl <- coxph_mpl(Surv(time, status == 2) ~ age + sex + ph.karno + wt.loss, data = lung)

par(mfrow=c(1,2))
plot(residuals(fit_mpl), which=1:2, ask=FALSE)
}
}
\keyword{survival}
