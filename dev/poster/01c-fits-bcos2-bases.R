## The four bases fitted on bcos2, so FIG 2 can use the same dataset as the
## worked-example panel.  melanoma's event times are heavily bunched (median
## 0.07 years against a 4-year maximum), which crushes every curve into the left
## edge of the panel whatever is done with the styling; bcos2 runs 4 to 60
## months evenly and is real trial data.
suppressMessages(pkgload::load_all(quiet = TRUE))
library(survival)

cache_file <- file.path("dev", "poster", "cache", "fits.rds")
cache <- readRDS(cache_file)
data(bcos2)

n_obs <- sum(!is.na(bcos2$left) & !is.na(bcos2$right))
form <- Surv(left, right, type = "interval2") ~ treatment
bases <- c("uniform", "msplines", "gaussian", "epanechikov")

fits <- list()
for (b in bases) {
  t0 <- Sys.time()
  f <- try(coxph_mpl(form, data = bcos2,
    control = coxph_mpl.control(n.obs = n_obs, basis = b, n.knots = c(8, 2),
                                max.iter = c(40, 2e4, 5e4))), silent = TRUE)
  if (inherits(f, "try-error")) { message(b, ": ERROR"); next }
  message(sprintf("%-12s %5.0f s  lambda = %-10.3g df = %5.2f  ploglik = %9.3f",
                  b, as.numeric(difftime(Sys.time(), t0, units = "secs")),
                  f$control$smooth, f$df, f$ploglik))
  fits[[b]] <- f
}

cache$fits_bcos2 <- fits
saveRDS(cache, cache_file)
message("cached ", length(fits), " bcos2 basis fits")
