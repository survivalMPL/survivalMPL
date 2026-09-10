## FIG 3 - the four basis families, drawn rather than described.
## Occupies the full-width strip at the foot of the poster: 24.4 x 2.75 in.
## Chrome is stripped to the minimum: the point is the shape of psi_u(t), not
## its level, so the y axis carries no numbers.

suppressMessages(pkgload::load_all(quiet = TRUE))
source(file.path("dev", "poster", "_theme.R"))

data(melanoma)
events <- sort(c(melanoma$t_L[melanoma$t_L > 0],
                 melanoma$t_R[is.finite(melanoma$t_R)]))

## Deliberately few knots: at the package defaults the curves overlap into an
## unreadable thicket at poster viewing distance.  Fewer knots change how many
## basis functions there are, not the shape of any one of them.
n_knots_vis <- c(0, 6)
bases <- c("uniform", "msplines", "gaussian", "epanechikov")

out <- file.path("dev", "poster", "figures", "fig3-basis-shapes.png")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)

poster_png(out, width = 24.80, height = 2.35, pointsize = 15)
par(mfrow = c(1, 4), mar = c(2.7, 0.9, 2.2, 0.9), bg = PAL$paper,
    family = SANS, xaxs = "i", yaxs = "i", lend = 1, ljoin = 1)

for (b in bases) {
  ctrl  <- coxph_mpl.control(basis = b, n.obs = 200, n.knots = n_knots_vis)
  knots <- survivalMPL:::compute_knots(ctrl, events)
  tv    <- seq(min(events[events > 0]), max(events), length.out = 1200)
  psi   <- survivalMPL:::compute_basis_matrix(tv, knots, b, ctrl$order, which = 1)

  m    <- ncol(psi)
  ymax <- max(psi[is.finite(psi)]) * 1.15
  ## One hue per basis family, shaded light -> dark with the basis index, so a
  ## panel still says which family it is even without reading the label.
  base_col <- BASIS_COL[[b]]
  cols <- colorRampPalette(c(lighten(base_col, 0.66), base_col,
                             darken(base_col, 0.30)))(m)

  plot(NA, xlim = range(tv), ylim = c(0, ymax), axes = FALSE, xlab = "", ylab = "")

  ## knots as faint uprights, so the reader can see what the basis is built on
  a <- knots$Alpha[knots$Alpha >= min(tv) & knots$Alpha <= max(tv)]
  abline(v = a, col = PAL$rule, lty = 3, lwd = 1.4)

  ## The uniform basis is a set of indicators; drawn with lines() the risers of
  ## neighbouring functions land on top of each other, so fill them instead.
  if (b == "uniform") {
    for (u in seq_len(m)) {
      inside <- which(psi[, u] > 0)
      if (!length(inside)) next
      rect(tv[min(inside)], 0, tv[max(inside)], psi[inside[1], u],
           col = adjustcolor(cols[u], alpha.f = 0.45), border = cols[u], lwd = 2.6)
    }
  } else {
    for (u in seq_len(m)) lines(tv, psi[, u], col = cols[u], lwd = 3.4)
  }

  axis(1, col = PAL$rule, col.axis = PAL$mute, lwd = 1.8, cex.axis = 0.9,
       family = SANS, tck = -0.035, mgp = c(3, 0.55, 0))
  mtext("t  (years)", side = 1, line = 1.7, adj = 1, cex = 0.78,
        col = PAL$mute, family = SANS)
  mtext(BASIS_LAB[[b]], side = 3, line = 0.75, adj = 0, font = 2, cex = 1.02,
        col = base_col, family = MONO)
  mtext(sprintf("m = %d", m), side = 3, line = 0.85, adj = 1, cex = 0.8,
        col = PAL$mute, family = SANS)
}
invisible(dev.off())
message("wrote ", out)
