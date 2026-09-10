## FIG 2, production version: the original two panels - baseline hazard and
## baseline survival, four bases each - made legible at poster distance.
##
## Presentation changes from the version on the deck, and why:
##
##   log time axis   melanoma's event times sit at median 0.07 of a four-year
##                   range, so on a linear axis every curve is crushed into the
##                   left tenth of the panel.  This is the change that matters;
##                   nothing else helps while the axis is linear.
##   x from 0.005    rather than the 5th percentile, 0.0009 years, which reads
##                   as nine hours.  Every summary number is unchanged by the
##                   trim: 2.9x on the hazard, 6.3 pp on survival.
##   round ticks     0.005 / 0.02 / 0.1 / 0.5.  An earlier attempt included
##                   0.01, only a factor of two from 0.005, and R silently
##                   suppressed the label.
##   no subtitles    the panel titles carry it; the agreement figures (2.9x on
##                   the hazard, 6.3 pp on survival) belong in the caption on
##                   the deck rather than inside the plot.
##   thin lines      lwd 3, not 5.5.  Heavy strokes on four overlapping curves
##                   read as a smear rather than as four estimates.
##
## SIZE.  The earlier draft read as squashed because more than half its height
## went on chrome: at 4.40 in tall, 3.9 lines of bottom margin, 3.1 of top and
## 1.9 of outer margin for the legend left only about 2.05 in of plot, so each
## panel had a 3.1:1 plot region.  Fixed three ways - trim the margins, move the
## legend inside the survival panel (its lower-left corner is empty, since
## survival is still near 1 at the earliest times), and take the figure to
## 5.10 in tall.  Each panel plot region is now about 6.4 x 3.9 in, or 1.6:1.
##
## The picture box on the deck has to match, or PowerPoint stretches the image:
## x 15.82, y 34.80, w 15.56, h 5.10.  That is 0.83 in taller than the box in
## v6; the height comes from dead space at the foot of the SMOOTHING panel and
## a 0.22 in nudge to the closing block - see 08-update-v6.R.
##
## melanoma, not bcos2: with 94 subjects and no exact event times the four bases
## on bcos2 disagree by 31 pp on survival, which is the opposite of the point.

suppressMessages(pkgload::load_all(quiet = TRUE))
source(file.path("dev", "poster", "_theme.R"))

FIG_W <- 15.56
FIG_H <- 5.10
LWD   <- 3.0

cache <- readRDS(file.path("dev", "poster", "cache", "fits.rds"))
fits <- cache$fits
data(melanoma)

bases <- names(fits)
cols <- unname(BASIS_COL[bases])
labs <- unname(BASIS_LAB[bases])

ev <- sort(c(melanoma$t_L[melanoma$t_L > 0],
             melanoma$t_R[is.finite(melanoma$t_R)]))
lo <- 0.005
hi <- as.numeric(quantile(ev, 0.95))
tt <- exp(seq(log(lo), log(hi), length.out = 600))

H <- sapply(fits, function(f) baseline(f, tt, "hazard")$est)
S <- sapply(fits, function(f) baseline(f, tt, "survival")$est)

haz_ratio <- max(apply(H, 1, function(x) max(x) / min(x)))
spread_pp <- 100 * max(apply(S, 1, function(x) diff(range(x))))
message(sprintf("hazard differs by up to %.1fx; survival agrees to %.1f pp",
                haz_ratio, spread_pp))

XAT <- c(0.005, 0.02, 0.1, 0.5)
XLB <- c("0.005", "0.02", "0.1", "0.5")
XLAB <- "t  (years, log scale)"

out <- file.path("dev", "poster", "figures", "fig2-readable.png")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)

poster_png(out, width = FIG_W, height = FIG_H, pointsize = 19)
par(mfrow = c(1, 2), mar = c(3.0, 4.6, 1.6, 0.9), oma = c(0, 0, 0, 0),
    bg = PAL$paper, family = SANS, lend = 1, ljoin = 1)

## ---- baseline hazard: where the bases visibly differ ----------------------
## Headroom at the top: these curves are highest at the earliest times, and
## without it they run into the panel title.
plot(NA, xlim = range(tt), ylim = range(H) * c(1, 1.18), log = "xy",
     axes = FALSE, xlab = "", ylab = "")
soft_grid()
for (i in seq_along(bases)) lines(tt, H[, i], col = cols[i], lwd = LWD)
clean_axes(xlab = XLAB, ylab = expression(hat(h)[0](t)), cex.axis = 0.95,
           xat = XAT, xlabels = XLB)
panel_title("Baseline hazard", cex = 1.1, line = 0.45)

## ---- baseline survival: where they do not --------------------------------
plot(NA, xlim = range(tt), ylim = c(0, 1), log = "x", axes = FALSE,
     xlab = "", ylab = "")
soft_grid()
for (i in seq_along(bases)) lines(tt, S[, i], col = cols[i], lwd = LWD)
clean_axes(xlab = XLAB, ylab = expression(hat(S)[0](t)), cex.axis = 0.95,
           xat = XAT, xlabels = XLB)
panel_title("Baseline survival", cex = 1.1, line = 0.45)
## One legend for both panels, in the empty lower-left corner of this one: it
## costs no plot height, and the colours mean the same thing in each panel.
legend("bottomleft", legend = labs, col = cols, lwd = LWD, bty = "n",
       cex = 0.92, seg.len = 1.5, text.col = PAL$ink2, ncol = 2,
       x.intersp = 0.6, y.intersp = 1.15)

invisible(dev.off())
message(sprintf("wrote %s  (%.2f x %.2f in)", out, FIG_W, FIG_H))
