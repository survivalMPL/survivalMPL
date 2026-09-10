## Shared look-and-feel for the ICTMC 2026 poster figures.
## Colours are taken straight from the poster palette so the plots sit on the
## page as part of the design rather than as pasted-in screenshots.

PAL <- list(
  paper   = "#FDFCF9",  # poster background
  ink     = "#121417",  # body text
  ink2    = "#3B3E44",
  mute    = "#606369",  # captions, axis labels
  rule    = "#D3D8E0",  # hairlines
  grid    = "#E7E9ED",
  blue    = "#025894",  # primary
  blue_lt = "#CFE0F2",
  navy    = "#003A68",
  orange  = "#EB883B",  # accent
  orange_d= "#903300",
  orange_lt = "#FCE6D9"
)

## One colour per basis family, used identically in every figure.
BASIS_COL <- c(
  uniform      = "#025894",
  msplines     = "#B4530A",
  gaussian     = "#12776F",
  epanechikov  = "#6A4C93"
)
BASIS_LAB <- c(uniform = "uniform", msplines = "msplines",
               gaussian = "gaussian", epanechikov = "epanechnikov")

FONT_DIR <- Sys.getenv("POSTER_FONT_DIR", unset = "")

register_poster_fonts <- function() {
  if (!nzchar(FONT_DIR) || !dir.exists(FONT_DIR)) return(invisible(FALSE))
  f <- function(x) file.path(FONT_DIR, x)
  ok <- all(file.exists(f(c("IBMPlexSans-Regular.ttf", "IBMPlexSans-SemiBold.ttf",
                            "IBMPlexSans-Italic.ttf", "IBMPlexMono-Regular.ttf"))))
  if (!ok) return(invisible(FALSE))
  systemfonts::register_font("PlexSans",
    plain = f("IBMPlexSans-Regular.ttf"), bold = f("IBMPlexSans-SemiBold.ttf"),
    italic = f("IBMPlexSans-Italic.ttf"), bolditalic = f("IBMPlexSans-SemiBold.ttf"))
  systemfonts::register_font("PlexMono",
    plain = f("IBMPlexMono-Regular.ttf"), bold = f("IBMPlexMono-SemiBold.ttf"))
  invisible(TRUE)
}

SANS <- if (register_poster_fonts()) "PlexSans" else "sans"
MONO <- if (SANS == "PlexSans") "PlexMono" else "mono"

## agg_png at poster scale: sizes are in inches at the size the figure occupies
## on the A0 sheet, so text set in points here is text in points on the poster.
poster_png <- function(file, width, height, res = 300, pointsize = 18) {
  ragg::agg_png(file, width = width, height = height, units = "in", res = res,
                background = PAL$paper, scaling = 1, pointsize = pointsize)
}

## Panel title in the poster's voice: small, semibold, left-aligned above the box.
panel_title <- function(txt, cex = 1.15, line = 1.1) {
  mtext(txt, side = 3, line = line, adj = 0, font = 2, cex = cex,
        col = PAL$ink, family = SANS)
}
panel_sub <- function(txt, cex = 0.95, line = 0.05) {
  mtext(txt, side = 3, line = line, adj = 0, cex = cex,
        col = PAL$mute, family = SANS)
}

## Light dotted grid drawn behind the data, never over it.
soft_grid <- function(h = TRUE, v = TRUE) {
  u <- par("usr")
  if (v) abline(v = axTicks(1), col = PAL$grid, lwd = 1.4)
  if (h) abline(h = axTicks(2), col = PAL$grid, lwd = 1.4)
}

## Axis box reduced to two hairlines (left + bottom), the poster's line weight.
clean_axes <- function(xlab = "", ylab = "", cex.axis = 1.0, las = 1,
                       xat = NULL, yat = NULL, ylabels = TRUE,
                       xlabels = TRUE) {
  axis(1, at = xat, labels = xlabels, col = PAL$rule, col.axis = PAL$mute,
       lwd = 1.6, cex.axis = cex.axis, family = SANS, tck = -0.02,
       mgp = c(3, 0.6, 0))
  axis(2, at = yat, labels = ylabels, col = PAL$rule, col.axis = PAL$mute,
       lwd = 1.6, las = las, cex.axis = cex.axis, family = SANS,
       tck = -0.02, mgp = c(3, 0.6, 0))
  title(xlab = xlab, col.lab = PAL$ink2, family = SANS, line = 2.0, cex.lab = 1.05)
  title(ylab = ylab, col.lab = PAL$ink2, family = SANS, line = 2.6, cex.lab = 1.05)
}

## Mix a colour towards white (amount 0 = unchanged, 1 = white).
lighten <- function(col, amount = 0.6) {
  rgb_ <- grDevices::col2rgb(col) / 255
  grDevices::rgb(t(rgb_ + (1 - rgb_) * amount))
}
## Mix a colour towards black.
darken <- function(col, amount = 0.25) {
  rgb_ <- grDevices::col2rgb(col) / 255
  grDevices::rgb(t(rgb_ * (1 - amount)))
}

## --- baseline quantities at x = 0 -------------------------------------------
## coxph_mpl() returns Theta already corrected back to the uncentred covariate
## scale, and fit$covar holds the covariance of that corrected (Beta, Theta), so
## the theta block gives baseline standard errors directly.  predict() would
## instead evaluate everything at the *mean* covariate vector.
baseline <- function(fit, time, what = c("hazard", "survival"), se = "M2QM2") {
  what  <- match.arg(what)
  which <- if (what == "hazard") 1 else 2
  B <- survivalMPL:::compute_basis_matrix(time, fit$knots, fit$control$basis,
                                          fit$control$order, which = which)
  theta <- fit$coef$Theta
  p <- fit$dim$p
  V <- fit$covar[[se]][(p + 1):(p + fit$dim$m), (p + 1):(p + fit$dim$m), drop = FALSE]
  lin   <- drop(B %*% theta)
  varl  <- pmax(rowSums((B %*% V) * B), 0)
  sel   <- sqrt(varl)
  if (what == "hazard") {
    data.frame(time = time, est = lin,
               lo = pmax(lin - 1.96 * sel, 0), hi = lin + 1.96 * sel)
  } else {
    ## S0 = exp(-H0); delta method on the cumulative hazard
    s <- exp(-lin)
    data.frame(time = time, est = s,
               lo = pmax(exp(-(lin + 1.96 * sel)), 0),
               hi = pmin(exp(-pmax(lin - 1.96 * sel, 0)), 1))
  }
}

## Confidence band as a soft ribbon rather than a pair of lines.
ribbon <- function(x, lo, hi, col, alpha = 0.16) {
  ok <- is.finite(x) & is.finite(lo) & is.finite(hi)
  polygon(c(x[ok], rev(x[ok])), c(lo[ok], rev(hi[ok])),
          col = adjustcolor(col, alpha.f = alpha), border = NA)
}
