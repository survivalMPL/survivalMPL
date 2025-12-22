\name{bcos2}
\alias{bcos2}
\docType{data}
\title{Breast Cosmesis Data}
\description{
Interval censored data described and given in full in Finkelstein and Wolfe (1985), discussed by More (2016, example 12.2) and available in the R package 'interval' (refer to ?bcos). Compared to the interval package version, bcos2 simply recode the lower value of left-censored data (NA instead of 0) and upper value of right-censored data (NA instead of Inf) to allow an easy identification of left-censored data by means of the function \code{\link[survival]{Surv}} when \code{type=="interval2"}. 
}
\usage{data(bcos2)}
\format{
  A data frame with 94 observations on the following 3 variables.
  \describe{
    \item{\code{left}}{a numeric vector}
    \item{\code{right}}{a numeric vector}
    \item{\code{treatment}}{a factor with levels \code{Rad} and \code{RadChem}}
  }
}
%\details{}
\source{
Finkelstein, D.M., and Wolfe, R.A. (1985). A semiparametric model for regression analysis of interval-censored 
failure time data. Biometrics 41: 731-740.

Moore, D. K. (2016), Applied Survival Analysis Using R, \emph{Springer}. 

}
\examples{
data(bcos2)
}
\keyword{dataset}
