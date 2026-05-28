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
\method{residuals}{coxph_mpl}(object, klein=FALSE, \dots)

\method{plot}{residuals.coxph_mpl}(x, ask=TRUE, which=1:3, upper.quantile=.95, \dots)
}
\arguments{
\item{object}{
an object inheriting from class \code{\link{coxph_mpl}}, representing
a fitted Cox proportional hazard model. 
}
\item{klein}{
a logical indicating if a cumulated hazard estimate for the Cox and Snell
residuals should be compute or not. May be (slightly) time consuming 
if \code{klein==TRUE}. This is required to perform Klein and Moeschberger 
plot of Cox and Snell residuals (See \bold{Details}).
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
is required, specify a subset of the numbers \code{1:3}. See \bold{Details}. 
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
regression model.

For object of class \code{residuals.coxph_mpl}, the available residual plots
are, respectively, the martingale residual plot (\code{which==1}), 
the Cox and Snell residual plot (\code{which==2}) and Klein and Moeschberger 
plot of Cox and Snell residuals (\code{which==3}, see Klein and Moeschberger 
(2005, p. 356)). 

The residual cumulative hasard function estimate used in the Klein and 
Moeschberger plot of Cox and Snell residuals is computed with the same control 
arguments (see \code{\link{coxph_mpl.control}}) as in \code{object}.
}
\value{
A data.frame of class \code{residuals.coxph_mpl} of \eqn{n} rows 
with following columns: \code{'time'}, the model outcome (with a random noise
added to event ties if \code{ties=='epsilon'} in \code{\link{coxph_mpl.control}}); 
\code{'observed'}, a logical vector indicating if observations are fully observed
(i.e., non-censored); \code{'coxsnell'}, the Cox and Snell residuals; 
\code{'martingale'}, the martingale residuals; \code{'H0coxsnell'}, the 
residual cumulative hasard function estimate (only available if 
\code{klein==TRUE}).
}
\seealso{
\code{\link{coxph_mpl}}, \code{\link{coxph_mpl.control}},
\code{\link{coxph_mpl.object}}, \code{\link{predict.coxph_mpl}} and
\code{\link{summary.coxph_mpl}}.
}
\references{Klein, J. P. and Moeschberger, M. L. (2003),
\emph{Survival Analysis: Techniques for Censored and Truncated Data}, 
2nd edition, Springer.  

Collett, D. (2003), and Moeschberger, M. L. (2003),
\emph{Modelling Survival Data in Medical Research}, Chapman and All.  
}
\author{Dominique-Laurent Couturier}
\examples{
### lung data of the survival package (see ?lung)
data(lung)
fit_mpl <- coxph_mpl(Surv(time, status == 2) ~ age + sex + ph.karno + wt.loss, data = lung)

par(mfrow=c(1,2))
plot(residuals(fit_mpl), which=1:2, ask=FALSE)

}
\keyword{survival}
