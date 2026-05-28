% Part of the survivalMPL-package
% Distributed under GPL 2 or later
% summary.coxph_mpl function and print
\name{summary.coxph_mpl}
\encoding{utf8}
\alias{summary.coxph_mpl}
\alias{print.summary.coxph_mpl}
\title{Summarise a coxph_mpl Object}
\description{
Additional information about the Cox proportional hazard model fit represented
by \code{object} is extracted and included in the returned object, which
is suitable for printing with the generic \code{print} function. The generic
\code{coef} function will extract the matrix of coefficients of interest
with standard errors, \eqn{z}-statistics and \eqn{p}-values. 
See \code{\link{coef.summary.coxph_mpl}}.

Only the baseline hazard parameters larger than \code{min.Theta} (see 
\code{\link{coxph_mpl.control}}) are reported.  
}
\usage{
\method{summary}{coxph_mpl}(object, se="M2QM2", full=FALSE, \dots)

\method{print}{summary.coxph_mpl}(x, se="M2QM2", \dots)
}
\arguments{
\item{object}{
In an object inheriting from class \code{\link{coxph_mpl}}, representing
a fitted Cox proportional hazard model.
}
\item{se}{
an inference method.  Possibilites are \code{"H"}, \code{"M1QM1"}, \code{"M2QM2"}, \code{"M1HM1"} 
and \code{"M2HM2"}. Refer to the \bold{Details} Section of \code{\link{coxph_mpl}}.  
Default is \code{se="M2QM2"}. 
}
\item{full}{
logical. If \code{TRUE} inference for the baseline hasard parameters is 
provided. Default is \code{full=FALSE}.   
}
\item{x}{
an object inheriting from class \code{\link{summary.coxph_mpl}}, representing
the summary of a fitted Cox proportional hazard model.
}
\item{\dots}{
Other arguments  passed through to printing functions.
}
}
\value{
an object of class \code{summary.coxph_mpl} representing the fit and additional
information.
\item{Beta}{
a matrix of \eqn{p} rows indicating the regression parameter estimates, standard
errors, z-statistics values and according p-values.
}
\item{Theta}{
If \code{full == TRUE}, a matrix of \eqn{m} rows indicating the baseline hazard
parameter estimates, standard errors, z-statistics values and according p-values. 
If \code{full == FALSE}, the baseline hazard estimates.
}
\item{inf}{
a list of elements extracted from the object of class \code{\link{coxph_mpl}}
including the number of iterations and the penalised likelihood 
value, for example.
}
}
\seealso{
\code{\link{coxph_mpl}},  \code{\link{coxph_mpl.control}},
\code{\link{coxph_mpl.object}}  and  \code{\link{plot.coxph_mpl}}.
}
\author{Dominique-Laurent Couturier.}
\examples{
data(lung)

fit_mpl <- coxph_mpl(Surv(time, status == 2) ~ age + sex + ph.karno + wt.loss, data = lung)

summary(fit_mpl, full = TRUE)
summary(fit_mpl, se = "M1HM1")
}
\keyword{survival}
