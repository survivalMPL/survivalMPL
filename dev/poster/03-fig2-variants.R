## FIG 2, in two versions, so the poster can be seen either way before
## committing.  Both are drawn at 16.04 x 4.40 in, the size the picture already
## occupies on the v3 deck, so swapping them is a byte swap.
##
##   fig2a-two-panels.png  baseline hazard | baseline survival, both large
##   fig2b-with-lambda.png the same two, plus what the roughness penalty does
##
## Melanoma throughout: 300 subjects, eight covariates, all four observation
## schemes, and a knot specification shared by every basis so the comparison is
## between basis shapes and not between parameter counts.

suppressMessages(pkgload::load_all(quiet = TRUE))
source(file.path("dev", "poster", "_theme.R"))

cache <- readRDS(file.path("dev", "poster", "cache", "fits.rds"))
fits  <- cache$fits
data(melanoma)

bases <- names(fits)
cols  <- unname(BASIS_COL[bases])

## Beyond the 95th percentile of the event times there is almost nothing left to
## estimate from, so the tails there are noise rather than signal.
events <- sort(c(melanoma$t_L[melanoma$t_L > 0],
                 melanoma$t_R[is.finite(melanoma$t_R)]))
tt <- seq(as.numeric(quantile(events, 0.05)),
          as.numeric(quantile(events, 0.95)), length.out = 400)

haz <- lapply(fits, baseline, time = tt, what = "hazard")
srv <- lapply(fits, baseline, time = tt, what = "survival")

out_dir <- file.path("dev", "poster", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## --------------------------------------------------------------- panels ----
panel_hazard <- function(cex_title = 1.05, legend = FALSE) {
  ylim <- range(unlist(lapply(haz, function(d) d$est)))
  plot(NA, xlim = range(tt), ylim = ylim, log = "y", axes = FALSE,
       xlab = "", ylab = "")
  soft_grid()
  for (i in seq_along(fits)) lines(tt, haz[[i]]$est, col = cols[i], lwd = 4.4)
  clean_axes(xlab = "t  (years)", ylab = expression(hat(h)[0](t)), cex.axis = 0.95)
  panel_title("Baseline hazard", cex = cex_title, line = 1.45)
  panel_sub("log scale - where the bases disagree", cex = 0.84, line = 0.35)
  if (legend)
    legend("topright", legend = unname(BASIS_LAB[bases]), col = cols, lwd = 4.4,
           bty = "n", cex = 0.95, seg.len = 1.6, text.col = PAL$ink2,
           y.intersp = 1.25)
}

panel_survival <- function(cex_title = 1.05, legend = TRUE) {
  plot(NA, xlim = range(tt), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  soft_grid()
  for (i in seq_along(fits)) lines(tt, srv[[i]]$est, col = cols[i], lwd = 4.4)
  clean_axes(xlab = "t  (years)", ylab = expression(hat(S)[0](t)), cex.axis = 0.95)
  panel_title("Baseline survival", cex = cex_title, line = 1.45)
  panel_sub("where they do not", cex = 0.84, line = 0.35)
  if (legend)
    legend("topright", legend = unname(BASIS_LAB[bases]), col = cols, lwd = 4.4,
           bty = "n", cex = 0.95, seg.len = 1.6, text.col = PAL$ink2,
           y.intersp = 1.25)
}

## The penalty panel holds the dataset and the basis fixed and varies only
## lambda, so the roughness penalty is the only thing moving.  It is drawn on
## bcos2: melanoma selects lambda ~ 1e-5 and stops converging once that value is
## imposed from the start, so it has no room to show the penalty working, while
## bcos2 selects lambda ~ 3e9 and is the dataset the call and output panels
## already quote.
##
## The largest lambda shown is not a cold-start fixed-lambda fit: forcing
## lambda ~ 1e9 from the start lands the multiplicative update in a much worse
## optimum (penalised log-likelihood -173 against -144 for the warm-started
## automatic fit).  The automatic fit itself is used for that curve instead, so
## every curve in the panel is one the package would actually return.
LAM_SHOW <- as.numeric(strsplit(Sys.getenv("POSTER_LAM_SHOW", "0,10000"),
                                ",")[[1]])

panel_lambda <- function(cex_title = 1.05) {
  lf <- cache$lam_fits
  key <- function(lam) sprintf("bcos2_%g", lam)
  keep <- LAM_SHOW[vapply(LAM_SHOW, function(l) !is.null(lf[[key(l)]]), logical(1))]
  stopifnot(length(keep) >= 1)
  fs <- lapply(keep, function(l) lf[[key(l)]])

  auto <- cache$fit_bcos2
  keep <- c(keep, auto$control$smooth)
  fs   <- c(fs, list(auto))
  LAM_AUTO <- auto$control$smooth

  ## bcos2 is in months and its own knot range, not melanoma's
  tb <- seq(min(fs[[1]]$knots$Alpha), max(fs[[1]]$knots$Alpha), length.out = 400)
  hs <- lapply(fs, baseline, time = tb, what = "hazard")

  ## light -> dark as lambda grows: the eye reads "more penalty" as "more settled"
  lcol <- colorRampPalette(c(PAL$orange, PAL$blue, darken(PAL$navy, 0.2)))(length(keep))

  ## Log scale, as in the first panel: unpenalised, the tail spike sits two
  ## orders of magnitude above the rest and flattens everything else.
  ylim <- range(unlist(lapply(hs, function(d) pmax(d$est, 1e-4))))
  plot(NA, xlim = range(tb), ylim = ylim, log = "y", axes = FALSE,
       xlab = "", ylab = "")
  soft_grid()
  for (i in seq_along(hs)) lines(tb, hs[[i]]$est, col = lcol[i], lwd = 4.4)
  clean_axes(xlab = "months", ylab = expression(hat(h)[0](t)), cex.axis = 0.95)
  panel_title("What the penalty does", cex = cex_title, line = 1.45)
  panel_sub("bcos2, msplines - only lambda changes", cex = 0.84, line = 0.35)

  fmt <- function(l) {
    if (l == 0) return("0")
    e <- floor(log10(l)); m <- l / 10^e
    if (isTRUE(all.equal(m, 1, tolerance = 5e-2))) sprintf("1e%d", e)
    else sprintf("%.1fe%d", m, e)
  }
  labs <- vapply(seq_along(keep), function(i)
    sprintf("lambda = %s   df %.1f%s", fmt(keep[i]), fs[[i]]$df,
            if (isTRUE(all.equal(keep[i], LAM_AUTO))) "  (chosen)" else ""),
    character(1))
  legend("topleft", legend = labs, col = lcol, lwd = 4.4, bty = "n",
         cex = 0.88, seg.len = 1.6, text.col = PAL$ink2, y.intersp = 1.25)
}

## ------------------------------------------------------------- variant a ----
f <- file.path(out_dir, "fig2a-two-panels.png")
poster_png(f, width = 16.04, height = 4.40, pointsize = 19)
par(mfrow = c(1, 2), mar = c(3.6, 4.8, 3.0, 1.0), bg = PAL$paper,
    family = SANS, lend = 1, ljoin = 1)
panel_hazard(cex_title = 1.05, legend = FALSE)
panel_survival(cex_title = 1.05, legend = TRUE)
invisible(dev.off())
message("wrote ", f)

## ------------------------------------------------------------- variant b ----
f <- file.path(out_dir, "fig2b-with-lambda.png")
poster_png(f, width = 16.04, height = 4.40, pointsize = 17)
par(mfrow = c(1, 3), mar = c(3.6, 4.6, 3.0, 0.9), bg = PAL$paper,
    family = SANS, lend = 1, ljoin = 1)
panel_hazard(cex_title = 1.05, legend = FALSE)
panel_survival(cex_title = 1.05, legend = TRUE)
panel_lambda(cex_title = 1.05)
invisible(dev.off())
message("wrote ", f)
