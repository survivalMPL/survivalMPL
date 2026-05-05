%%
%% coxph_mpl function
%%
\name{coxph_mpl.object}
\alias{coxph_mpl.object}
\encoding{utf8}
\title{
  MPL Proportional Hazards Regression Object
}
\description{
This class of objects is returned by the \code{\link{coxph_mpl}} class of functions
to represent a proportional hazards model fitted by maximum penalised likelihood.
Objects of this class have methods for the functions \code{print}, 
\code{summary}, \code{plot}, \code{residuals} and \code{predict}.

All components described under \bold{Arguments} must be included in a legitimate
\code{\link{coxph_mpl}} object.
}
\arguments{
\item{coef}{
a list of length 2 containg the parameter estimates of each model part.
The first list, named \code{'Beta'}, contains the vector of regression parameter
estimates of length \eqn{p}. The second list, named \code{'Theta'}, contains the 
vector of the baseline hazard parameter estimates of length \eqn{m}.
}
\item{se}{
a list of length 2 containg the parameter standard errors of each model part.
The first list, named \code{'Beta'}, is a (\eqn{p} x \eqn{5}) matrix indicating 
the standard errors of each regression parameter according to 5 inference methods. 
See \bold{Details} in \code{\link{coxph_mpl}}. The second list, named \code{'Theta'}, 
is a (\eqn{m} x \eqn{5}) matrix indicating the standard errors of each baseline 
hazard parameter according to 5 inference methods. See \bold{Details} in 
\code{\link{coxph_mpl}}.
}
\item{covar}{
a list of length 5 containg the (\eqn{p+m} x \eqn{p+m}) covariance matrix according to 
5 inference methods. See \bold{Details} in \code{\link{coxph_mpl}}.
}
\item{ploglik}{
a vector of length 2. The first element is the penalised log-likelihood 
with the final values of the coefficients. (The second element is a correction 
factor for the baseline hazard parameters due to the use of a centered \bold{X} matrix 
in the estimation process.)
}
\item{iter}{
a vector of length 3 indicating the number of iterations used to estimate the 
smoothing parameter (first value, equal to \code{1} when the user specified a 
chosen value), the Beta and Theta parameters during the entire process (second 
value), and Beta and Theta parameters during the last smoothing parameter iteration
(third value). 
}
\item{knots}{
list of length 3 to 4 containg parameters of the chosen basis: \code{'m'}, the number
of used bases; \code{'Alpha'}, the knot sequence of length \eqn{m+1} for the uniform 
basis, and of length \eqn{m} otherwise; \code{'Delta'}, the value of the integral of
each base over the data support (which equals 1 when \code{basis != 'uniform'}); 
\code{'Sigma'}, only available for the Gaussian basis, corresponds to the standard 
deviation of each truncated Gaussian base.
}
\item{control}{
Object of class \code{\link{coxph_mpl.control}} specifying the basis, smoothing parameter 
value and other options. See \code{\link{coxph_mpl.control}}.
}
\item{dim}{
a list of length 5 with following elements: \code{'n'}, the sample size; \code{'n.events'}, 
the number of events; \code{'n.ties'}, the number of duplicated observations; \code{'p'}, 
the number of regression parameters; and \code{'m'}, the number of baseline hazard 
parameters. 
}
\item{iter}{
number of iterations used. 
}
\item{call}{
the matched call.
}
\item{data}{
a list of length 3 with following elements: \code{'time'}, the outcome vector with 
an added noise applied  to duplicated observation if \code{ties == "epsilon"} 
in \code{\link{coxph_mpl.control}}; \code{'observed'}, a logical vector indicating 
if outcomes are fully observed or censored; \code{'X'}, the \bold{X} matrix corresponding 
to the model formula indicated in \code{\link{coxph_mpl}}.
}
}
\seealso{
\code{\link{coxph_mpl}}, \code{\link{summary.coxph_mpl}}, \code{\link{coef.coxph_mpl}},
\code{\link{plot.coxph_mpl}},\code{\link{residuals.coxph_mpl}}
and \code{\link{predict.coxph_mpl}}.   
}
\author{Dominique-Laurent Couturier.}
\keyword{survival}
