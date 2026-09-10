## FIG 1 - the problem, and what the package returns, on real trial data.
## Left-hand column of the poster, under THE PROBLEM: 15.4 x 2.73 in.
##
## bcos2 is the breast-cosmesis trial: 94 women, radiotherapy alone against
## radiotherapy plus chemotherapy, cosmetic deterioration assessed at clinic
## visits.  Not one event time in it is observed exactly, which is the point.
##
## A: what the trial records.
## B: the survival curves that come out, by arm, with 95% bands.
## C: the same data collapsed to one time per subject and handed to coxph().

suppressMessages(pkgload::load_all(quiet = TRUE))
library(survival)
source(file.path("dev", "poster", "_theme.R"))

cache <- readRDS(file.path("dev", "poster", "cache", "fits.rds"))
fit   <- cache$fit_bcos2
data(bcos2)

COL <- list(
  atrisk = PAL$blue,       # under observation, known event-free
  window = PAL$orange,     # the event is known to lie in here
  rad    = PAL$blue,
  chem   = PAL$orange_d
)

out <- file.path("dev", "poster", "figures", "fig1-what-is-observed.png")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)

poster_png(out, width = 15.4, height = 2.55, pointsize = 16)
layout(matrix(1:3, nrow = 1), widths = c(1.42, 1, 1))
par(mar = c(2.9, 8.6, 3.0, 0.5), bg = PAL$paper, family = SANS,
    lend = 1, ljoin = 1)

## ---- A: what the trial records ---------------------------------------------
bc <- transform(bcos2,
  type = ifelse(is.na(left), "left", ifelse(is.na(right), "right", "interval")))

set.seed(4)
pick <- function(ty, n) {
  ok <- which(bc$type == ty)
  ok <- sample(ok, min(n, length(ok)))
  key <- ifelse(bc$type[ok] == "left", bc$right[ok], bc$left[ok])
  ok[order(key)]
}
groups <- list(interval = pick("interval", 11), right = pick("right", 9),
               left = pick("left", 5))
glab <- c(interval = "interval censored", right = "right censored",
          left = "left censored")

rows <- list(); y <- 0; ticks <- c()
for (g in names(groups)) {
  ys <- y + seq_along(groups[[g]])
  rows[[g]] <- data.frame(i = groups[[g]], y = ys)
  ticks[g] <- mean(ys)
  y <- max(ys) + 2.0
}

x_max <- 62
plot(NA, xlim = c(0, x_max), ylim = c(y - 1.2, -0.2), axes = FALSE,
     xlab = "", ylab = "", yaxs = "i")
abline(v = axTicks(1), col = PAL$grid, lwd = 1.4)

bar <- 0.40
for (g in names(rows)) {
  for (k in seq_len(nrow(rows[[g]]))) {
    i <- rows[[g]]$i[k]; yy <- rows[[g]]$y[k]
    tL <- bc$left[i]; tR <- bc$right[i]
    if (g == "interval") {
      segments(0, yy, tL, yy, col = COL$atrisk, lwd = 2.6)
      rect(tL, yy - bar, tR, yy + bar, col = adjustcolor(COL$window, 0.45),
           border = COL$window, lwd = 1.6)
    } else if (g == "right") {
      segments(0, yy, tL, yy, col = COL$atrisk, lwd = 2.6)
      rect(tL, yy - bar, x_max, yy + bar, col = adjustcolor(COL$window, 0.18),
           border = NA)
      arrows(tL, yy, x_max * 0.99, yy, length = 0.05, angle = 22,
             col = COL$window, lwd = 2)
    } else {
      rect(0, yy - bar, tR, yy + bar, col = adjustcolor(COL$window, 0.45),
           border = COL$window, lwd = 1.6)
    }
  }
}
axis(1, col = PAL$rule, col.axis = PAL$mute, lwd = 1.6, cex.axis = 0.88,
     family = SANS, tck = -0.03, mgp = c(3, 0.45, 0))
axis(2, at = ticks, labels = glab[names(ticks)], tick = FALSE, las = 1,
     cex.axis = 0.88, col.axis = PAL$ink2, family = MONO, mgp = c(3, 0.3, 0))
mtext("months since treatment", side = 1, line = 1.6, adj = 1, cex = 0.7,
      col = PAL$mute, family = SANS)
panel_title("What the trial records", cex = 0.98, line = 1.75)
panel_sub("bcos2, 25 of 94 women - not one event time is observed exactly",
          cex = 0.74, line = 0.5)

## ---- B and C: what comes back ----------------------------------------------
tt <- seq(min(fit$knots$Alpha), max(fit$knots$Alpha), length.out = 400)

## predict() at a chosen subject's covariates: one woman from each arm.
i_rad  <- which(bcos2$treatment == "Rad")[1]
i_chem <- which(bcos2$treatment == "RadChem")[1]
s_rad  <- predict(fit, type = "survival", i = i_rad,  time = tt)
s_chem <- predict(fit, type = "survival", i = i_chem, time = tt)

par(mar = c(2.9, 4.4, 3.0, 0.6))
plot(NA, xlim = range(tt), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
soft_grid()
ribbon(tt, s_rad$low,  s_rad$high,  COL$rad,  alpha = 0.16)
ribbon(tt, s_chem$low, s_chem$high, COL$chem, alpha = 0.16)
lines(tt, s_rad$survival,  col = COL$rad,  lwd = 4)
lines(tt, s_chem$survival, col = COL$chem, lwd = 4)
clean_axes(xlab = "months", ylab = expression(hat(S)(t)), cex.axis = 0.88)
panel_title("Time to cosmetic deterioration", cex = 0.98, line = 1.75)
panel_sub("interval-censored fit, 95% bands", cex = 0.74, line = 0.5)
legend("bottomleft", legend = c("radiotherapy", "radiotherapy + chemo"),
       col = c(COL$rad, COL$chem), lwd = 4, bty = "n", cex = 0.82,
       seg.len = 1.4, text.col = PAL$ink2, y.intersp = 1.2)

## ---- C: against collapsing every interval to one time ----------------------
s0 <- baseline(fit, tt, "survival")
bc_mid <- data.frame(t = cache$bcos2_mid, d = cache$bcos2_status,
                     treatment = bcos2$treatment)
cox_mid <- coxph(Surv(t, d) ~ treatment, data = bc_mid)
bh <- basehaz(cox_mid, centered = FALSE)
bh <- rbind(data.frame(hazard = 0, time = 0), bh[bh$time <= max(tt), ])

plot(NA, xlim = range(tt), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
soft_grid()
ribbon(tt, s0$lo, s0$hi, PAL$blue, alpha = 0.16)
lines(bh$time, exp(-bh$hazard), type = "s", col = PAL$orange_d, lwd = 3)
lines(tt, s0$est, col = PAL$blue, lwd = 4)
clean_axes(xlab = "months", ylab = expression(hat(S)[0](t)), cex.axis = 0.88)
panel_title("Same answer, no imputation", cex = 0.98, line = 1.75)
panel_sub("MPL on the intervals vs midpoint + coxph()", cex = 0.74, line = 0.5)
legend("bottomleft", legend = c("coxph_mpl()", "midpoint + coxph()"),
       col = c(PAL$blue, PAL$orange_d), lwd = c(4, 3), bty = "n", cex = 0.82,
       seg.len = 1.4, text.col = PAL$ink2, y.intersp = 1.2)

invisible(dev.off())
message("wrote ", out)
