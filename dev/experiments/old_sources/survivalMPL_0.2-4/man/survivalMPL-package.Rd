% Part of the survivalMPL-package
% Distributed under GPL 2 or later
\name{survivalMPL-package}
\alias{survivalMPL-package}
\alias{survivalMPL}
\docType{package}
\title{
Penalised Maximum Likelihood for Survival Analysis Models
}
\description{
Simultaneously estimate the regression coefficients and provide a 
'non-parametric' smooth estimate of the baseline hazard function for 
proportional hazard Cox models using maximum penalised likelihood (MPL).
}
\details{
This package allows to perform simultaneous estimation of the regression
coefficients and baseline hazard function in Cox proportional hazard models,
with right, left and interval censored data and independent censoring, by maximising a penalised
likelihood, in which a penalty function is used to smooth the baseline hazard 
estimate. 

Optimisation is achieved using a new iterative algorithm, which combines 
Newton's method and the multiplicative iterative algorithm by Ma (2010), and 
respects the non-negativity constraints on the baseline hazard estimate 
(refer to Ma, Couturier, Heritier and Marschner (2021)).  

Valid inferences for the regression coefficients and the baseline hazard, 
cumulative baseline hazard and survival functions as well as for their 
predictions are available.

This software is accepted by users "as is" and without warranties or guarantees
of any kind.
}
\author{
Dominique-Laurent Couturier, Jun Ma, Stephane Heritier, Maurizio Manuguerra. 

Maintainer: Dominique-Laurent Couturier \email{dominique-laurent.couturier@cruk.cam.ac.uk}.
}
\references{
Ma, J. and Couturier, D.-L., and Heritier, S. and Marschner, I.C. (2021), Penalized likelihood estimation of the proportional hazards model for survival data with interval censoring. \emph{International Journal of Biostatistics},doi:10.1515/ijb-2020-0104.

Ma, J. and Heritier, S. and Lo, S. (2014), On the Maximum Penalised
Likelihood Approach for Proportional Hazard Models with Right Censored 
Survival Data. \emph{Computational Statistics and Data Analysis} \bold{74}, 142-156.

Ma, J. (2010), Positively constrained multiplicative iterative algorithm for 
maximum penalised likelihood tomographic reconstruction. 
\emph{IEEE Transactions On Signal Processing} \bold{57}, 181-192.
}
\keyword{ package }
\keyword{ survival }
