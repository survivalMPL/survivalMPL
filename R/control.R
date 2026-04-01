#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
basis.name_mpl <- function(k){
  if(k == "discr"| k == "discretized" | k == "discretised" | k == "unif" | k == "uniform"){"uniform"
  }else{if(k == "m" | k == "msplines" | k == "mspline"){"msplines"
  }else{if(k == "gauss" | k == "gaussian"){"gaussian"
  }else{if(k == "epa" | k == "epanechikov"){"epanechikov"
  }else{stop("Unknown basis choice", call. = FALSE)}}}}}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
penalty.order_mpl <- function(p,basis,order){
  p = as.integer(p)
  switch(basis,
         'uniform'     = ifelse(p>0 & p<3,p,2),
         'gaussian'    = ifelse(p>0 & p<3,p,2),
         'msplines'    = order-1,
         'epanechikov' = 2)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Ancillary Arguments for Controlling \code{coxph_mpl} Fits
#'
#' Set numeric and algorithmic controls for \code{coxph_mpl} fits. The function
#' validates inputs (e.g., number of events per basis element, iteration limits)
#' to avoid impossible settings.
#'
#' @param n.obs Number of fully observed (non-censored) outcomes. Required when
#'   \code{basis == "uniform"} to derive an acceptable range for
#'   \code{n.events_basis}.
#' @param basis Basis used to approximate the baseline hazard. One of
#'   \code{"uniform"}, \code{"gaussian"}, \code{"msplines"}, or
#'   \code{"epanechikov"}. Defaults to \code{"uniform"}.
#' @param smooth Smoothing parameter value. Defaults to \code{NULL} (estimated
#'   via REML). Set to \code{0} for maximum-likelihood estimates. Must be
#'   non-negative.
#' @param max.iter Integer vector of length 3 giving maximum iterations for (1)
#'   smoothing parameter updates, (2) inner beta/theta updates, and (3) total
#'   inner iterations. Defaults to \code{c(150, 7.5e4, 1e6)}.
#' @param tol Convergence tolerance on parameter change between iterations.
#'   Defaults to \code{1e-7}.
#' @param n.knots Integer vector of length 2 controlling internal knots for
#'   non-uniform bases. The first entry sets quantile knots between
#'   \code{range.quant}; the second sets equally spaced knots outside that
#'   range. Defaults to \code{c(8, 2)} for M-splines and \code{c(0, 20)}
#'   otherwise.
#' @param n.events_basis Integer giving the number of fully observed outcomes
#'   per uniform basis element. Must lie in \code{[1, floor(n.obs/2)]}. Defaults
#'   to \code{round(3.5 * log(n.obs) - 7.5)} when valid.
#' @param range.quant Length-2 numeric vector giving the quantile range used
#'   when setting quantile knots for non-uniform bases. Defaults to
#'   \code{c(0.075, 0.9)}.
#' @param cover.sigma.quant Proportion of fully observed outcomes targeted
#'   within the 0.025-0.975 interval of truncated Gaussian bases tied to
#'   quantile knots. Defaults to \code{0.25}.
#' @param cover.sigma.fixed Proportion of the outcome range targeted within the
#'   0.025-0.975 interval of untruncated Gaussian bases tied to fixed knots.
#'   Defaults to \code{0.25}.
#' @param min.theta Minimum baseline hazard parameter value reported; estimates
#'   below are treated as zero (active constraints). Defaults to \code{1e-10}.
#' @param penalty Integer specifying penalty order. First- and second-order
#'   penalties are available for \code{"uniform"} and \code{"gaussian"} bases;
#'   \code{"epanechikov"} uses second-order; \code{"msplines"} uses
#'   \code{order - 1}. Defaults to \code{2}.
#' @param order Integer order for \code{"msplines"} and \code{"epanechikov"}
#'   bases (default \code{3}). Order 1 M-splines correspond to uniform bases;
#'   order 2 to triangular bases.
#' @param kappa Step-size reduction factor (>1) used when the penalised
#'   likelihood fails to increase. Defaults to \code{1 / 0.6}.
#' @param epsilon Length-2 numeric vector giving safeguards for survival and
#'   baseline hazard values to avoid logarithm issues. Defaults to
#'   \code{c(1e-16, 1e-10)}.
#' @param ties Strategy for handling duplicated fully observed outcomes when
#'   defining knot sequences. Use \code{"epsilon"} to jitter duplicates with
#'   small random noise; use \code{"unique"} to drop duplicates. Defaults to
#'   \code{"epsilon"}.
#' @param seed Optional seed (integer vector compatible with
#'   \code{.Random.seed}) used when \code{ties == "epsilon"}; preserves the
#'   current RNG state when set.
#'
#' @return A list with validated control settings (class
#'   \code{"coxph_mpl.control"}).
#' @seealso [coxph_mpl()]
#' @export
coxph_mpl.control <- function(n.obs=NULL, basis = "uniform", smooth = NULL, max.iter=c(1.5e+2,7.5e+4,1e+6), tol=1e-7,
                              n.knots = NULL, n.events_basis = NULL, range.quant = c(0.075,.9),
                              cover.sigma.quant = .25, cover.sigma.fixed=.25, min.theta = 1e-10,
                              penalty = 2L, order = 3L, kappa = 1/.6, epsilon = c(1e-16,1e-10), ties = "epsilon", seed = NULL){
  basis        = resolve_basis_name(basis)
  max.iter     = c(ifelse(is.null(smooth),ifelse(max.iter[1]>0,as.integer(max.iter[1]),1.5e+2),1L),
                   ifelse(max.iter[2]>0,as.integer(max.iter[2]),7.5e+4),
                   ifelse(length(max.iter)==2,1e+6,
                          ifelse(max.iter[3]>ifelse(max.iter[2]>0,as.integer(max.iter[2]),7.5e+4),
                                 as.integer(max.iter[3]),1e+6)))
  tol          = ifelse(tol>0 & tol<1,tol,1e-7)
  order        = ifelse(order>0 & order<6,as.integer(order),3L)
  min.theta    = ifelse(min.theta>0 & min.theta<1e-3,min.theta,1e-10)
  penalty      = compute_penalty_order(basis, penalty, order)
  kappa        = ifelse(kappa>1, kappa, 1/.6)
  cover.sigma.quant  = ifelse(cover.sigma.quant>0 & cover.sigma.quant<0.4,cover.sigma.quant,.75)
  cover.sigma.fixed  = ifelse(cover.sigma.fixed>0 & cover.sigma.fixed<0.4,cover.sigma.fixed,.75)
  if(all(range.quant<=1) & all(range.quant>=0) & length(range.quant)==2){
    range.quant = range.quant[order(range.quant)]
  }else{range.quant = c(0.075,.9)}
  if(is.null(n.knots)|sum(n.knots)<3|length(n.knots)!=2){
    n.knots    = if(basis!='uniform' & basis!='msplines'){c(0,20)}else{c(8,2)}
  }
  if(!is.null(n.events_basis)){
    n.events_basis = ifelse(n.events_basis<1|n.events_basis>floor(n.obs/2),
                            max(round(3.5*log(n.obs)-7.5),1L),round(n.events_basis))
  }else{n.events_basis = max(round(3.5*log(n.obs)-7.5),1L)}
  if(!is.null(smooth)){
    smooth = ifelse(smooth<0,0,smooth)
  }else{smooth=0}
  out = list(basis = basis, smooth = smooth, max.iter = max.iter, tol = tol,
             order = order, penalty = penalty, n.knots = n.knots, range.quant = range.quant,
             cover.sigma.quant = cover.sigma.quant, cover.sigma.fixed = cover.sigma.fixed,
             n.events_basis = as.integer(n.events_basis), min.theta = min.theta, ties = ties,
             seed = as.integer(seed), kappa = kappa, epsilon = epsilon)
  class(out) = "coxph_mpl.control"
  out
}
