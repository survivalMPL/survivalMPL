## Fits used by the ICTMC 2026 poster figures.  A few of them take minutes, so
## everything is cached in dev/poster/cache/fits.rds and the plotting scripts
## only read that.
##
##   bcos2    - real trial data (breast cosmesis), every observation censored;
##              used for FIG 1, the "what a trial actually records" figure.
##   melanoma - simulated, 8 covariates, all four observation types; used for
##              FIG 2 and FIG 3, where the basis comparison needs more structure
##              than a single treatment indicator.

suppressMessages(pkgload::load_all(quiet = TRUE))
library(survival)

out_dir <- file.path("dev", "poster", "cache")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------- melanoma --
data(melanoma)

form_mel <- Surv(t_L, t_R, type = "interval2") ~ Arm + Leg + Trunk +
  mm1to2 + mm2to4 + mm4plus + Female + Age_centred
n_obs <- sum(melanoma$t_L == melanoma$t_R & is.finite(melanoma$t_R))

## Same knot specification for every basis, so the comparison is between basis
## shapes and not between different numbers of parameters.  Smoothing is left at
## its default (selected automatically), which is what the poster claims.
base_ctrl <- list(n.obs = n_obs, n.knots = c(8, 2), max.iter = c(40, 2e4, 5e4))

bases <- c("uniform", "msplines", "gaussian", "epanechikov")
fits <- list()
for (b in bases) {
  t0 <- Sys.time()
  fits[[b]] <- coxph_mpl(form_mel, data = melanoma,
    control = do.call(coxph_mpl.control, c(base_ctrl, list(basis = b))))
  message(sprintf("melanoma %-12s %5.0f s   lambda = %.3g", b,
                  as.numeric(difftime(Sys.time(), t0, units = "secs")),
                  fits[[b]]$control$smooth))
}

## ------------------------------------------------------------------- bcos2 --
data(bcos2)
bc <- bcos2
bc$type <- with(bc, ifelse(is.na(left), "left",
                    ifelse(is.na(right), "right", "interval")))

t0 <- Sys.time()
fit_bc <- coxph_mpl(Surv(left, right, type = "interval2") ~ treatment,
                    data = bcos2, basis = "msplines", n.knots = c(8, 2))
message(sprintf("bcos2    %-12s %5.0f s   lambda = %.3g", "msplines",
                as.numeric(difftime(Sys.time(), t0, units = "secs")),
                fit_bc$control$smooth))

## Naive comparator: collapse each subject to one time and hand it to the
## partial likelihood.  (0, r] -> r/2, [l, r] -> midpoint, [l, Inf) -> l.
bc_mid <- with(bcos2, ifelse(is.na(left), right / 2,
                      ifelse(is.na(right), left, (left + right) / 2)))
bc_status <- as.integer(!is.na(bcos2$right))
cox_bc <- coxph(Surv(bc_mid, bc_status) ~ treatment,
                data = data.frame(bc_mid, bc_status, treatment = bcos2$treatment))

## Same collapse for melanoma, used nowhere on the poster but kept so the
## comparison can be repeated there.
mel_mid <- with(melanoma, ifelse(!is.finite(t_R), t_L,
                          ifelse(t_L == t_R, t_L, (t_L + t_R) / 2)))
mel_status <- as.integer(is.finite(melanoma$t_R))

saveRDS(list(
  fits       = fits,          # melanoma, one per basis
  fit_bcos2  = fit_bc,        # bcos2, msplines
  cox_bcos2  = cox_bc,        # bcos2, midpoint-imputed partial likelihood
  bcos2_mid  = bc_mid,
  bcos2_status = bc_status,
  mel_mid    = mel_mid,
  mel_status = mel_status,
  beta_true  = c(Arm = -0.56, Leg = 0.01, Trunk = -0.22, mm1to2 = 0.22,
                 mm2to4 = 0.87, mm4plus = 1.13, Female = -0.17,
                 Age_centred = 0.14)),
  file.path(out_dir, "fits.rds"))
message("saved ", file.path(out_dir, "fits.rds"))
