## Fixed-smoothing fits for the "what the penalty does" panel of FIG 2.
##
## lambda is not on a dataset-independent scale.  melanoma selects lambda ~ 1e-5
## and the multiplicative update already fails to converge when that same value
## is imposed from the start, so melanoma has no room to show the penalty doing
## anything.  bcos2 selects lambda ~ 3e9, is small (n = 94), and is the dataset
## already quoted in the call and output panels, so the grid is fitted there.
##
## Every fit is kept, with the diagnostics needed to see whether it converged:
## a run that reaches the total iteration cap comes back with a penalised
## log-likelihood around -1e12, which is a divergence, not a fit.

suppressMessages(pkgload::load_all(quiet = TRUE))
library(survival)

cache_file <- file.path("dev", "poster", "cache", "fits.rds")
cache <- readRDS(cache_file)
data(bcos2)

CAP <- 5e4
n_obs <- sum(!is.na(bcos2$left) & !is.na(bcos2$right))

lambdas <- c(0, 1e2, 1e4, 1e6, 1e7, 1e8, 1e9, 3e9, 1e11)
lam_fits <- list()

for (lam in lambdas) {
  t0 <- Sys.time()
  f <- try(coxph_mpl(Surv(left, right, type = "interval2") ~ treatment,
    data = bcos2,
    control = coxph_mpl.control(n.obs = n_obs, basis = "msplines",
                                n.knots = c(8, 2), smooth = lam,
                                max.iter = c(1, 2e4, CAP))), silent = TRUE)
  if (inherits(f, "try-error")) {
    message(sprintf("lambda=%-9g  ERROR", lam)); next
  }
  ok <- is.finite(f$ploglik) && f$ploglik > -1e6 && f$iter[2] < CAP
  message(sprintf("lambda=%-9g %4.0f s  %s  df=%5.2f  ploglik=%12.3f  beta=%.3f",
                  lam, as.numeric(difftime(Sys.time(), t0, units = "secs")),
                  if (ok) "ok      " else "DIVERGED", f$df, f$ploglik, coef(f)[1]))
  if (ok) lam_fits[[sprintf("bcos2_%g", lam)]] <- f
}

cache$lam_fits <- lam_fits
saveRDS(cache, cache_file)
message("cached ", length(lam_fits), " convergent fixed-lambda bcos2 fits")
