% Part of the survivalMPL-package
% Distributed under GPL 2 or later
% coxph_mpl function
\name{coxph_mpl}
\encoding{utf8}
\alias{coxph_mpl}
\alias{print.coxph_mpl}
\title{Fit Cox Proportional Hazards Regression Model Via MPL}
\description{
Simultaneously estimate the regression coefficients and the 
baseline hazard function of proportional hazard Cox models 
using maximum penalised likelihood (MPL).
}
\usage{
coxph_mpl(formula, data, subset, na.action, control, \dots)

\method{print}{coxph_mpl}(x,\dots)
}
\arguments{
\item{formula}{
a formula object, with the response on the left of a \code{~} operator, and 
the terms on the right.  The response must be a survival object as 
returned by the \code{\link[survival]{Surv}} function. 
}
\item{data}{
a data.frame in which to interpret the variables named in 
the \code{formula}, or in the \code{subset} argument. If no dataset is indicated,
variables will be taken from the global environment. 
}
\item{subset}{
expression indicating which subset of the rows of data should be used in 
the fit. All observations are included by default. 
}
\item{na.action}{
a missing-data filter function.  This is applied to the model.frame
after any subset argument has been used.  Default is \code{options()\$na.action}.
}
\item{x}{
an object inheriting from class \code{coxph_mpl}, representing
a fitted Cox proportional hazard model.
}
\item{control}{
Object of class \code{\link{coxph_mpl.control}} specifying control options like 
basis choice, smoothing parameter value and maximum number of itereations, 
for example. Refer to \code{\link{coxph_mpl.control}} to see the defaults.
}
\item{\dots}{Other arguments. In \code{coxph_mpl}, these elements, 
will be passed to \code{\link{coxph_mpl.control}}. In \code{print.coxph_mpl}, 
these elements will be passed to the \code{print} function.
}
} 
\value{
an object of class \code{coxph_mpl} representing the fit.
See \code{\link{coxph_mpl.object}} for details.
}
\details{
\code{coxph_mpl} allows to simultaneously estimate the regression
coefficients and baseline hazard function of Cox proportional hazard models,
with right censored data and independent censoring, by maximising a penalised
likelihood, in which a penalty function is used to smooth the baseline hazard 
estimate. 

Optimisation is achieved using a new iterative algorithm, which combines 
Newton's method and the multiplicative iterative algorithm proposed by 
Ma (2010), and respects the non-negativity constraints on the baseline 
hazard estimate (refer to Jun, Heritier and Lo (2014)).  

The centered \bold{X} matrix is used in the optimisation process to get a 
better shaped (penalised) log-likelihood. Baseline hazard parameter estimates 
and covariance matrix are then respectively corrected using a correction factor
and the delta method.

When the chosen basis is not uniform, estimates of zero are possible for 
baseline hazard parameters and will correspond to active constraints as 
defined by Moore and Sadler (2008). Inference, as described by Ma, Heritier and 
Lo (2014), is then corrected accordingly (refer to Moore and Sadler (2008)) by
adequately 'cutting' the corresponding covariance matrix.

There are currently 5 ways to perform inference on model parameters:  

Let \eqn{H} denote the Hessian matrix of the unpenalised likelihood, 
\eqn{M_{1}} denote the product of the first order derivative of the penalised 
likelihood by the transpose of the first order derivative of the unpenalised 
likelihood, \eqn{Q} denote the product of the first order derivative of 
the penalised likelihood by its transpose, and \eqn{M_{2}} denote the second 
order derivative of the penalised likelihood. Then, 
\itemize{
\item \code{'H'} refers to \eqn{H^{-1}}, the inverse of the Hessian matrix, 
\item \code{'M1QM1'}, refers to the sandwich formula \eqn{M_{1}^{-1} Q (M_{1}^{-1})^{T}},
\item \code{'M2QM2'}, refers to the sandwich formula \eqn{M_{2}^{-1} Q M_{2}^{-1}},
\item \code{'M1HM1'}, refers to the sandwich formula \eqn{M_{1}^{-1} H (M_{1}^{-1})^{T}},
\item \code{'M2HM2'}, refers to the sandwich formula \eqn{M_{2}^{-1} H M_{2}^{-1}}.
}
Simulations analysing the coverage levels of confidence intervals for the 
regression parameters seem to indicate that \eqn{M_{2}^{-1} H M_{2}^{-1}} 
performs better when using the uniform basis, and that 
\eqn{M_{2}^{-1} Q M_{2}^{-1}} performs when using other bases.
}
\seealso{
\code{\link{coxph_mpl.object}},  \code{\link{coxph_mpl.control}},
\code{\link{summary.coxph_mpl}}  and  \code{\link{plot.coxph_mpl}}.
}
\references{Ma, J. and Heritier, S. and Lo, S. (2014), On the Maximum Penalised
Likelihood Approach for Proportional Hazard Models with Right Censored 
Survival Data. \emph{Computational Statistics and Data Analysis} \bold{74}, 142-156.

Ma, J. (2010), Positively constrained multiplicative iterative algorithm for 
maximum penalised likelihood tomographic reconstruction. 
\emph{IEEE Transactions On Signal Processing} \bold{57}, 181-192.

Moore, T. J. and Sadler, B. M. and Kozick R. J. (2008), Maximum-Likelihood 
Estimation, the Cramer-Rao Bound, and the Method of Scoring With Parameter 
Constraints, \emph{IEEE Transactions On Signal Processing} \bold{56, 3}, 895-907. 
}
\author{Dominique-Laurent Couturier, Jun Ma, Stephane Heritier. Design inspired by the function 
\code{\link[survival]{coxph}} of the \code{survival} package.}
\examples{
## data lung of the survival package. Refer to ?lung.
data(lung)

fit_mpl <- coxph_mpl(Surv(time, status == 2) ~ age + sex + ph.karno + wt.loss, data = lung)

fit_mpl
}
\keyword{survival}
