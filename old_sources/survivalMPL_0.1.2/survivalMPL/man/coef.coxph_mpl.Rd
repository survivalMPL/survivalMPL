% Part of the survivalMPL-package
% Distributed under GPL 2 or later
% coef.coxph_mpl function
\name{coef.coxph_mpl}
\encoding{utf8}
\alias{coef.coxph_mpl}
\alias{coef.summary.coxph_mpl}
\title{Extract coefficients of a coxph_mpl Object or of its \code{summary}}
\description{
Extract the coefficients of the model part of interest of a \code{coxph_mpl} 
object, and the matrix of coefficients of the model part of interest and 
corresponding standard errors, \eqn{z}-statistics and \eqn{p}-values of a 
\code{summary.coxph_mpl} object. 
}
\usage{
\method{coef}{coxph_mpl}(object, parameters="Beta",\dots)

\method{coef}{summary.coxph_mpl}(object, parameters="Beta",\dots)
}
\arguments{
\item{object}{
an object inheriting from class \code{\link{coxph_mpl}} or \code{\link{summary.coxph_mpl}}.
}
\item{parameters}{
the set of parameters of interest.  Indicate \code{parameters="Beta"} for the
regression parameters and \code{parameters="Theta"} for the baseline hazard 
parameters.  Default is \code{parameters="Beta"}. 
}
\item{\dots}{
other arguments. 
}
}
\value{
a vector of coefficients or a matrix of coefficients with standard errors, 
z-statistics and corresponding p-values.
}
\details{
When the input is of class \code{summary.coxph_mpl} and \code{parameters=="Theta"},
only the parameter estimates larger than \code{min.Theta} (see 
\code{\link{coxph_mpl.control}}) are reported. 
}
\seealso{
\code{\link{coxph_mpl}} and  \code{\link{summary.coxph_mpl}}.
}
\author{Dominique-Laurent Couturier.}
\examples{
data(lung)

fit_mpl <- coxph_mpl(Surv(time, status == 2) ~ age + sex + ph.karno + wt.loss, data = lung)

coef(fit_mpl)
coef(summary(fit_mpl, se = "M1HM1", full=TRUE))
}
\keyword{survival}
