% Part of the survivalMPL-package
% Distributed under GPL 2 or later
% coxph_mpl.control function
\name{coxph_mpl.control}
\alias{coxph_mpl.control}
\title{Ancillary arguments for controling coxph_mpl fits}
\description{
This is used to set various numeric parameters controling a Cox model fit 
using \code{\link{coxph_mpl}}. Typically it would only be used in a call 
to \code{coxph_mpl}. Some basic checks are performed on inputs, such that 
impossible argument values (like a negative number of events per base, 
for example) are avoided.
}
\usage{
coxph_mpl.control(n.obs=NULL, basis = "uniform", smooth = NULL, max.iter=c(50,1e+05),
        tol=1e-7, n.knots = NULL, n.events_basis = NULL, 
        range.quant = c(0.075,.9), cover.sigma.quant = .25, 
        cover.sigma.fixed=.25, min.theta = 1e-10, penalty = 2L,
        order = 3L, kappa = 1/.6, epsilon = 1e-50, ties = "unique", 
        seed = NULL)
}
\arguments{
\item{n.obs}{
the number of fully observed (i.e., non censored) outcomes. This
argument is only required when \code{basis=="uniform"} to define an acceptable 
range of values for \code{n.events_basis}.
}
\item{basis}{
the name of the basis to use to approximate the baseline hazard function. 
Available options are \code{"uniform"}, for a step function approximation); 
\code{"gaussian"}, using truncated Gaussian densities; \code{"msplines"}, 
as defined by Ramsay (1988); and \code{"epanechikov"}. Default is 
\code{basis="uniform"}.
}
\item{smooth}{
the smoothing parameter value. When specified, it should be larger or equal 
to zero. By default, the smoothing value is set to \code{NULL} and its optimal 
value is estimated via REML. Maximum likelihood estimates are obtained by 
specifying \code{smooth=0}. The effect of the smooting parameter over the estimates 
(currently) depends on the response range, such that its value is (currently)
difficult to interpret.    
}
\item{max.iter}{
a vector of 2 integers defining the maximum number of iterations for the smooth 
parameter (first value) and for the Beta and Theta (second value) parameters 
to attempt for convergence. Default is \code{max.iter=c(50,1e+05)}.
}
\item{tol}{
the convergence tolerence value. Convergence is achieved when the maximum 
absolute difference between the parameter estimates at iteration k and iteration 
k-1 is smaller than \code{tol}. Default is \code{tol=1e-7}.  
}
\item{n.knots}{
a vector of 2 integers defining how the internal knot sequence (the minimum and
maximum observations define the external knots) of non-uniform bases 
should be set. The first value specify the number of quantile knots to be set
between the \code{range.quant} quantiles of the fully observed (i.e., non 
censored) outcomes. The second value specify the number of equally spaced knots
to be set outside the range of the quantile knots. The first and last equally 
spaced knots equal the minimum and maximum response values. When the number of 
quantile knots is larger than 0, other equally spaced knots are set between the 
largest quantile knot and the maximum outcome value. The minimal total number 
of knots is 3. Default is \code{n.knots=c(8,2)} when \code{basis=="msplines"} 
and \code{n.knots=c(0,20)} otherwise.
}
\item{n.events_basis}{
an integer specifing the number of fully observed (i.e., non censored) outcome 
per uniform base. The value has to be larger or equal to one and smaller than 
\code{n.obs} divided by 2. Default is \code{round(3.5*log(n.obs)-7.5)} if it 
belongs to the accepted range of values.  
}
\item{range.quant}{
a vector of length 2 defining the range of the quantile knots when a non uniform
basis is chosen. By default, \code{range.quant = c(0.075,.9)} such that 
\code{n.knots[1]} quantile knots are set between the quantiles \code{0.075} 
and \code{0.9} of the fully observed (i.e., non censored) outcomes. 
}
\item{cover.sigma.quant}{
the proportion of fully observed (i.e., non censored) outcomes that should belong
to the interval defined by the quantiles 0.025 and 0.975 of each truncated
Gaussian base corresponding to a quantile knot (see \code{n.kots}). 
This value allows to define the standard deviation of these bases. 
Default is \code{cover.sigma.quant=.25}.
}
\item{cover.sigma.fixed}{
the proportion of the outcome range that should belong to the interval 
defined by the quantiles 0.025 and 0.975 of each \emph{un}truncated 
Gaussian base corresponding to each fixed knot (see \code{n.kots}). Default is 
\code{cover.sigma.fixed=.25}.
}
\item{min.theta}{
a value indicating the minimal baseline hazard parameter value in the output 
(i.e., after the fit). Baseline hazard parameter estimates lower than \code{min.theta}
will be considered as zero. Consequenlty, in the inference, these zero estimates
will correspond to active constraints as defined by Moore and Sadler (2008).
Default is \code{1e-10}.
}
\item{penalty}{
an integer specifying the order of the penalty matrix (see Ma, Heritier and
Lo (2008)). Currently, the first and second order penalty matrices are available 
for the \code{"uniform"} and \code{"gaussian"} bases, the second order penalty 
matrix is available for the \code{"epanechikov"} basis, and the penalty matrix
of the \code{"msplines"} basis is set to \code{order-1} (see \code{order} below).
Default is \code{penalty=2}.
}
\item{order}{
an integer specifying the order of the \code{"msplines"} (as defined by Ramsay (1988))
and \code{"epanechikov"} basis. Default is \code{order=3}. M-splines of order 1
correspond to a uniform base (with density equal to one) and M-splines of order
2 correspond to a triangular base.
}
\item{kappa}{
a value larger than 1 used in the fitting algorithm to decrease the step size 
when the penalised likelihood doesn't increase during the iterative process. 
Default is \code{kappa=1/.6}.
}
\item{epsilon}{
a value indicating the minimum value of baseline parameter estimates \bold{in} the 
fitting algorithm in order to avoid problems with logarithms. Default is 
\code{epsilon=1e-50}
}
\item{ties}{
a character string indicating a method to handle duplicated outcomes when defining the 
knots sequence (see \code{n.events_basis} and \code{n.knots}). Current available 
options are \code{"epsilon"} which add a random noise smaller than 1e-10 to each
duplicate fully observed (i.e., non censored) outcomes, and \code{"unique"} which
delete duplicated fully observed (i.e., non censored) outcomes when defing the 
knot sequence. Default is \code{ties="unique"}.
}
\item{seed}{
\code{NULL} or an integer vector compatible with \code{\link{.Random.seed}}: 
the seed to be used when adding a random noise to duplicate events when 
\code{ties="epsilon"}. The current value of \code{.Random.seed} will be
preserved if \code{seed} is set, i.e. non-\code{NULL};
otherwise, as by default, \code{.Random.seed} will be used and
modified as usual from calls to \code{\link{runif}()} etc.  
}
}
\value{
a list containing the values of each of the above arguments (except \code{n.obs}).
}
\references{
Ma, J. and Heritier, S. and Lo, S. (2014), On the Maximum Penalised
Likelihood Approach for Proportional Hazard Models with Right Censored 
Survival Data. \emph{Computational Statistics and Data Analysis} \bold{74}, 142-156.

Moore, T. J. and Sadler, B. M. and Kozick R. J. (2008), Maximum-Likelihood 
Estimation, the Cramer-Rao Bound, and the Method of Scoring With Parameter 
Constraints, \emph{IEEE Transactions On Signal Processing} \bold{56, 3}, 895-907. 

Ramsay, J. O. (1988), Monotone Regression Splines in Action, 
\emph{Statistical Science} \bold{3, 4}, 425-441. 
}
\author{Dominique-Laurent Couturier.}
\seealso{\code{\link{coxph_mpl}}
}
\keyword{survival}
