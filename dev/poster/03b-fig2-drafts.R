## FIG 2 candidates, for review before anything goes near the poster.
##
## Two problems with the version currently on the deck:
##
##   1. The survival panel draws four curves on top of one another, which reads
##      as "identical" - either a bug or a wasted panel.  They are not
##      identical; the spread reaches 6.3 percentage points.  Four overlapping
##      lines simply cannot show a 6 pp difference on a 0-1 axis.
##   2. melanoma's event times are bunched into the first weeks of a four-year
##      range (median 0.07, 75th percentile 0.28, max 4.0), so on a linear axis
##      every curve piles into the left edge.  A log time axis fixes that; no
##      amount of styling does.
##
## It has to stay melanoma.  bcos2 has the nicer time axis and is the dataset
## the worked-example panel uses, but with 94 subjects and no exact event times
## its four bases disagree by 31 percentage points on survival and by a factor
## of 5e6 on the hazard - the opposite of what this panel exists to show.  Set
## POSTER_FIG2_DATA=bcos2 to see that for yourself.
##
## All candidates are drawn at 16.04 x 4.40 in, the size the picture occupies on
## the deck, into dev/poster/figures/drafts/.

suppressMessages(pkgload::load_all(quiet = TRUE))
source(file.path("dev", "poster", "_theme.R"))

cache <- readRDS(file.path("dev", "poster", "cache", "fits.rds"))

DATA <- Sys.getenv("POSTER_FIG2_DATA", "melanoma")
LOGX <- !identical(Sys.getenv("POSTER_FIG2_LOGX", "yes"), "no")

if (DATA == "bcos2") {
  stopifnot(!is.null(cache$fits_bcos2))
  fits <- cache$fits_bcos2
  xlab <- "months since treatment"
  tt <- seq(min(fits[[1]]$knots$Alpha), max(fits[[1]]$knots$Alpha),
            length.out = 500)
} else {
  fits <- cache$fits
  xlab <- "t  (years, log scale)"
  data(melanoma)
  ev <- sort(c(melanoma$t_L[melanoma$t_L > 0],
               melanoma$t_R[is.finite(melanoma$t_R)]))
  ## Start at 0.005 years rather than the 5th percentile (0.0009 = about nine
  ## hours), which reads as nonsense on a poster.  Every headline number is
  ## unchanged by the trim: 3.5 pp, 6.3 pp, 2.87x either way.
  lo <- 0.005
  hi <- as.numeric(quantile(ev, 0.95))
  tt <- if (LOGX) exp(seq(log(lo), log(hi), length.out = 500)) else
    seq(lo, hi, length.out = 500)
}
if (!LOGX) xlab <- sub(", log scale", "", xlab, fixed = TRUE)
tag <- paste0(DATA, if (LOGX) "-logx" else "-linx")

bases <- names(fits)
cols <- unname(BASIS_COL[bases])
labs <- unname(BASIS_LAB[bases])

H <- sapply(fits, function(f) baseline(f, tt, "hazard")$est)
S <- sapply(fits, function(f) baseline(f, tt, "survival")$est)

spread_pp <- 100 * max(apply(S, 1, function(x) diff(range(x))))
haz_ratio <- max(apply(H, 1, function(x) max(x) / min(x)))
message(sprintf("%s: max survival spread %.1f pp, max hazard ratio %.3g x",
                tag, spread_pp, haz_ratio))

out_dir <- file.path("dev", "poster", "figures", "drafts")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(a, b) if (is.null(a)) b else a
XLOG <- if (LOGX && DATA != "bcos2") "x" else ""

## --------------------------------------------------------------------------
## Shared hazard panel: the half that already works, on a log time axis.
## --------------------------------------------------------------------------
panel_hazard <- function(legend = TRUE, sub = NULL, pos = "bottomleft") {
  plot(NA, xlim = range(tt), ylim = range(H), log = paste0(XLOG, "y"),
       axes = FALSE, xlab = "", ylab = "")
  soft_grid()
  for (i in seq_along(bases)) lines(tt, H[, i], col = cols[i], lwd = 5)
  clean_axes(xlab = xlab, ylab = expression(hat(h)[0](t)), cex.axis = 0.95)
  panel_title("Baseline hazard", cex = 1.1, line = 1.45)
  panel_sub(sub %||% sprintf("the four differ by up to %.1fx", haz_ratio),
            cex = 0.86, line = 0.35)
  if (legend)
    legend(pos, legend = labs, col = cols, lwd = 5, bty = "n",
           cex = 0.95, seg.len = 1.5, text.col = PAL$ink2, y.intersp = 1.22)
}

## --------------------------------------------------------------------------
## A: one survival curve, plus a ribbon whose thickness IS the disagreement
## --------------------------------------------------------------------------
f <- file.path(out_dir, sprintf("fig2-A-ribbon-%s.png", tag))
poster_png(f, width = 16.04, height = 4.40, pointsize = 19)
par(mfrow = c(1, 2), mar = c(3.9, 5.0, 3.1, 1.0), bg = PAL$paper,
    family = SANS, lend = 1, ljoin = 1)
panel_hazard()
band_lo <- apply(S, 1, min); band_hi <- apply(S, 1, max)
plot(NA, xlim = range(tt), ylim = c(0, 1), log = XLOG, axes = FALSE,
     xlab = "", ylab = "")
soft_grid()
## Drawn at true width, not exaggerated; the caption carries the number.
ribbon(tt, band_lo, band_hi, PAL$blue, alpha = 0.32)
lines(tt, S[, 1], col = PAL$blue, lwd = 5)
clean_axes(xlab = xlab, ylab = expression(hat(S)[0](t)), cex.axis = 0.95)
panel_title("Baseline survival", cex = 1.1, line = 1.45)
panel_sub(sprintf("band spans all four - never wider than %.1f pp", spread_pp),
          cex = 0.86, line = 0.35)
invisible(dev.off())
message("wrote ", f)

## --------------------------------------------------------------------------
## B: the same four curves, magnified - signed difference from their mean
## --------------------------------------------------------------------------
f <- file.path(out_dir, sprintf("fig2-B-difference-%s.png", tag))
poster_png(f, width = 16.04, height = 4.40, pointsize = 19)
par(mfrow = c(1, 2), mar = c(3.9, 5.0, 3.1, 1.0), bg = PAL$paper,
    family = SANS, lend = 1, ljoin = 1)
panel_hazard()
D <- 100 * (S - rowMeans(S))
plot(NA, xlim = range(tt), ylim = range(D) * 1.18, log = XLOG, axes = FALSE,
     xlab = "", ylab = "")
soft_grid()
abline(h = 0, col = PAL$rule, lwd = 3)
for (i in seq_along(bases)) lines(tt, D[, i], col = cols[i], lwd = 5)
clean_axes(xlab = xlab, ylab = "difference from their mean  (pp)",
           cex.axis = 0.95)
panel_title("Baseline survival", cex = 1.1, line = 1.45)
panel_sub(sprintf("the same four curves, magnified - all within %.1f pp",
                  max(abs(D))), cex = 0.86, line = 0.35)
invisible(dev.off())
message("wrote ", f)

## --------------------------------------------------------------------------
## C: both panels on the hazard - level, then relative to the mean
## --------------------------------------------------------------------------
f <- file.path(out_dir, sprintf("fig2-C-hazard-ratio-%s.png", tag))
poster_png(f, width = 16.04, height = 4.40, pointsize = 19)
par(mfrow = c(1, 2), mar = c(3.9, 5.0, 3.1, 1.0), bg = PAL$paper,
    family = SANS, lend = 1, ljoin = 1)
panel_hazard(sub = "log scale, same knots for every basis")
R <- H / exp(rowMeans(log(H)))
plot(NA, xlim = range(tt), ylim = range(R), log = paste0(XLOG, "y"),
     axes = FALSE, xlab = "", ylab = "")
soft_grid()
abline(h = 1, col = PAL$rule, lwd = 3)
for (i in seq_along(bases)) lines(tt, R[, i], col = cols[i], lwd = 5)
clean_axes(xlab = xlab, ylab = "relative to their mean", cex.axis = 0.95)
panel_title("Where they disagree", cex = 1.1, line = 1.45)
panel_sub("hazard, as a ratio to the four-basis mean", cex = 0.86, line = 0.35)
invisible(dev.off())
message("wrote ", f)

## --------------------------------------------------------------------------
## D: one wide hazard panel; the agreement stated rather than plotted
## --------------------------------------------------------------------------
f <- file.path(out_dir, sprintf("fig2-D-single-%s.png", tag))
poster_png(f, width = 16.04, height = 4.40, pointsize = 19)
layout(matrix(1:2, nrow = 1), widths = c(2.35, 1))
par(mar = c(3.9, 5.4, 3.1, 1.0), bg = PAL$paper, family = SANS,
    lend = 1, ljoin = 1)
panel_hazard(sub = "log scale, same knots for every basis", pos = "topright")
par(mar = c(3.9, 0.6, 3.1, 0.4))
plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
text(0, 0.98, "And what it costs", adj = c(0, 1), font = 2, cex = 1.15,
     col = PAL$ink, family = SANS)
text(0, 0.80, sprintf("%.1f×", haz_ratio), adj = c(0, 1), cex = 2.8,
     col = PAL$blue, family = SANS)
text(0, 0.63, "widest gap between two\nbases on the hazard at\nany one time",
     adj = c(0, 1), cex = 0.92, col = PAL$ink2, family = SANS)
text(0, 0.36, sprintf("%.1f pp", spread_pp), adj = c(0, 1), cex = 2.8,
     col = PAL$blue, family = SANS)
text(0, 0.19, "widest gap on baseline\nsurvival, anywhere in\nthe follow-up",
     adj = c(0, 1), cex = 0.92, col = PAL$ink2, family = SANS)
invisible(dev.off())
message("wrote ", f)

## --------------------------------------------------------------------------
## FINAL (B, refined): one legend serving both panels, placed where there is
## actually room, and explicit log ticks instead of R's irregular defaults.
## Not written into the deck by this script - see 08-update-v6.R for that.
## --------------------------------------------------------------------------
if (DATA == "melanoma" && LOGX) {
  XAT <- c(0.005, 0.01, 0.05, 0.1, 0.5)
  XLB <- c("0.005", "0.01", "0.05", "0.1", "0.5")

  f <- file.path(out_dir, "fig2-final.png")
  poster_png(f, width = 16.04, height = 4.40, pointsize = 19)
  ## The legend goes in an outer margin spanning both panels: four labels in a
  ## row need more width than either panel has on its own, and one legend for
  ## two panels that share a colour scheme is one legend too few, not too many.
  par(mfrow = c(1, 2), mar = c(3.9, 5.0, 3.1, 1.0), oma = c(1.9, 0, 0, 0),
      bg = PAL$paper, family = SANS, lend = 1, ljoin = 1)

  plot(NA, xlim = range(tt), ylim = range(H), log = "xy", axes = FALSE,
       xlab = "", ylab = "")
  soft_grid()
  for (i in seq_along(bases)) lines(tt, H[, i], col = cols[i], lwd = 5)
  clean_axes(xlab = xlab, ylab = expression(hat(h)[0](t)), cex.axis = 0.95,
             xat = XAT, xlabels = XLB)
  panel_title("Baseline hazard", cex = 1.1, line = 1.45)
  panel_sub(sprintf("the four bases differ by up to %.1fx", haz_ratio),
            cex = 0.86, line = 0.35)

  D <- 100 * (S - rowMeans(S))
  plot(NA, xlim = range(tt), ylim = range(D) * 1.2, log = "x",
       axes = FALSE, xlab = "", ylab = "")
  soft_grid()
  abline(h = 0, col = PAL$rule, lwd = 3)
  for (i in seq_along(bases)) lines(tt, D[, i], col = cols[i], lwd = 5)
  ## Short ylab: the panel title and subtitle already say what is plotted, and
  ## a longer rotated label does not fit the panel height.
  clean_axes(xlab = xlab, ylab = "difference  (pp)",
             cex.axis = 0.95, xat = XAT, xlabels = XLB,
             yat = pretty(D, 5))
  panel_title("Baseline survival", cex = 1.1, line = 1.45)
  panel_sub(sprintf("and yet they agree to within %.1f pp", max(abs(D))),
            cex = 0.86, line = 0.35)
  ## One legend for both panels, drawn across the whole figure.
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0),
      new = TRUE)
  plot.new()
  legend("bottom", legend = labs, col = cols, lwd = 5, bty = "n",
         cex = 0.95, seg.len = 1.6, text.col = PAL$ink2, horiz = TRUE,
         x.intersp = 0.6, text.width = NA)
  invisible(dev.off())
  message("wrote ", f)
}
