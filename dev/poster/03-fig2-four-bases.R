## FIG 2 - the same model fitted with all four basis families.
## Right-hand column of the poster, under the basis cards: 16.0 x 4.55 in.
##
## Three panels: the baseline hazard (where the bases visibly disagree), the
## baseline survival (where they do not), and the regression coefficients
## (which is what a trialist actually reports).

suppressMessages(pkgload::load_all(quiet = TRUE))
source(file.path("dev", "poster", "_theme.R"))

cache <- readRDS(file.path("dev", "poster", "cache", "fits.rds"))
fits  <- cache$fits
data(melanoma)

bases <- names(fits)
cols  <- unname(BASIS_COL[bases])

## Confine the picture to where there is data to estimate from: beyond the 95th
## percentile of the event times the tails are noise, not signal.
events <- sort(c(melanoma$t_L[melanoma$t_L > 0],
                 melanoma$t_R[is.finite(melanoma$t_R)]))
t_lo <- as.numeric(quantile(events, 0.05))
t_hi <- as.numeric(quantile(events, 0.95))
tt   <- seq(t_lo, t_hi, length.out = 400)

haz <- lapply(fits, baseline, time = tt, what = "hazard")
srv <- lapply(fits, baseline, time = tt, what = "survival")

out <- file.path("dev", "poster", "figures", "fig2-four-bases.png")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)

poster_png(out, width = 16.04, height = 4.40, pointsize = 17)
layout(matrix(1:3, nrow = 1), widths = c(1, 1, 1.06))
par(mar = c(3.6, 4.4, 3.0, 0.8), bg = PAL$paper, family = SANS,
    lend = 1, ljoin = 1)

## ---- 1. baseline hazard, log scale ----------------------------------------
ylim <- range(unlist(lapply(haz, function(d) d$est)))
plot(NA, xlim = range(tt), ylim = ylim, log = "y", axes = FALSE, xlab = "", ylab = "")
soft_grid()
for (i in seq_along(fits)) lines(tt, haz[[i]]$est, col = cols[i], lwd = 4)
clean_axes(xlab = "t  (years)", ylab = expression(hat(h)[0](t)), cex.axis = 0.92)
panel_title("Baseline hazard", cex = 1.05, line = 1.45)
panel_sub("log scale - where the bases disagree", cex = 0.82, line = 0.35)

## ---- 2. baseline survival --------------------------------------------------
plot(NA, xlim = range(tt), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
soft_grid()
for (i in seq_along(fits)) {
  lines(tt, srv[[i]]$est, col = cols[i], lwd = 4)
}
clean_axes(xlab = "t  (years)", ylab = expression(hat(S)[0](t)), cex.axis = 0.92)
panel_title("Baseline survival", cex = 1.05, line = 1.45)
panel_sub("where they do not", cex = 0.82, line = 0.35)
legend("topright", legend = unname(BASIS_LAB[bases]), col = cols, lwd = 4,
       bty = "n", cex = 0.92, seg.len = 1.6, text.col = PAL$ink2, y.intersp = 1.25)

## ---- 3. regression coefficients -------------------------------------------
cf <- sapply(fits, coef)
nm <- rownames(cf)
np <- length(nm)
xr <- range(cf) + c(-0.18, 0.18) * diff(range(cf))
par(mar = c(3.6, 7.6, 3.0, 0.8))
plot(NA, xlim = xr, ylim = c(0.4, np + 0.6), axes = FALSE, xlab = "", ylab = "")
abline(v = axTicks(1), col = PAL$grid, lwd = 1.4)
abline(v = 0, col = PAL$rule, lwd = 2)
off <- seq(-0.21, 0.21, length.out = length(bases))
for (j in seq_len(np)) {
  y <- np - j + 1
  segments(min(cf[j, ]), y, max(cf[j, ]), y, col = PAL$rule, lwd = 2)
  for (i in seq_along(bases))
    points(cf[j, i], y + off[i], pch = 19, cex = 1.15, col = cols[i])
}
axis(1, col = PAL$rule, col.axis = PAL$mute, lwd = 1.6, cex.axis = 0.92,
     family = SANS, tck = -0.02, mgp = c(3, 0.6, 0))
axis(2, at = np:1, labels = nm, tick = FALSE, las = 1, cex.axis = 0.88,
     col.axis = PAL$ink2, family = MONO, mgp = c(3, 0.4, 0))
title(xlab = expression(hat(beta)), col.lab = PAL$ink2, family = SANS,
      line = 2.0, cex.lab = 1.05)
panel_title("Regression coefficients", cex = 1.05, line = 1.45)
panel_sub("spread across bases \u2264 1 standard error", cex = 0.82, line = 0.35)

invisible(dev.off())
message("wrote ", out)
