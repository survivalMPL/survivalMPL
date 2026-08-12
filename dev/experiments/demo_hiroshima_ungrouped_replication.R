# Follow-up to demo_hiroshima_left_truncation.R.
#
# Two questions a single fit cannot answer:
#   1. Is the naive-vs-truncated difference SYSTEMATIC, or just the noise of
#      one particular sample? -> paired replication over subsamples: the same
#      rows are fitted with and without entry=, so sampling variability
#      cancels within each pair.
#   2. Does coxph_mpl(entry=, basis="uniform") stay accurate and affordable as
#      n grows to the full 86,611-subject cohort? -> accuracy/timing sweep
#      against survival::coxph() on (start, stop] data.
#
# Also fits the colon-cancer-specific outcome (621 events), where the LSS dose
# response is meant to live and where the rarity of events makes getting the
# risk sets right matter more.

devtools::load_all()
library(survival)

data(hiroshima)
COVS <- c("dose", "sex", "city")

subsample <- function(n, seed) {
  set.seed(seed)
  hiroshima[sort(sample.int(nrow(hiroshima), min(n, nrow(hiroshima)))), , drop = FALSE]
}

fit_mpl <- function(data_name, status_var, truncated) {
  cl <- call("coxph_mpl",
    formula = reformulate(COVS, response = bquote(Surv(time, .(as.name(status_var))))),
    data = as.name(data_name), basis = "uniform", smooth = 0
  )
  if (truncated) cl$entry <- quote(entry)
  coef(eval(cl, envir = globalenv()), "Beta")
}

fit_ref <- function(data, status_var, truncated) {
  lhs <- if (truncated) {
    bquote(Surv(entry, time, .(as.name(status_var))))
  } else {
    bquote(Surv(time, .(as.name(status_var))))
  }
  coef(coxph(reformulate(COVS, response = lhs), data = data))
}

## ============================================================ 1. paired ====
cat("\n== paired replication: 20 subsamples of n = 8000 subjects ==\n")
R <- 20
res <- vector("list", R)
for (r in seq_len(R)) {
  samp <- subsample(8000, seed = 100 + r)
  assign("samp", samp, envir = globalenv())
  res[[r]] <- c(
    trunc = fit_mpl("samp", "status", TRUE),
    naive = fit_mpl("samp", "status", FALSE)
  )
  cat(".")
}
cat("\n")
res <- do.call(rbind, res)

pars <- c("dose", "sexFemale", "cityNagasaki")
summ <- rbind(
  `MPL entry=` = colMeans(res[, paste0("trunc.", pars)]),
  `MPL naive` = colMeans(res[, paste0("naive.", pars)])
)
colnames(summ) <- pars
cat("\nmean over replications:\n")
print(round(summ, 4))

cat("\nwithin-sample difference (naive - entry=): mean [min, max], against the\n")
cat("Monte Carlo SD of the entry= estimate for scale.\n")
for (p in pars) {
  d <- res[, paste0("naive.", p)] - res[, paste0("trunc.", p)]
  cat(sprintf(
    "  %-13s %+.4f  [%+.4f, %+.4f]   MC SD %.4f   sign %s\n",
    p, mean(d), min(d), max(d), sd(res[, paste0("trunc.", p)]),
    if (all(d > 0)) "always +" else if (all(d < 0)) "always -" else "mixed"
  ))
}

## ============================================================== 2. sweep ===
cat("\n== sample-size sweep (all-cause), entry= vs coxph LTRC ==\n")
cat(sprintf(
  "%8s %6s  %8s %9s %9s  %7s  %6s\n",
  "n", "events", "dose", "sexFemale", "cityNag", "maxdiff", "secs"
))
for (n in c(2000, 8000, 20000, 50000, nrow(hiroshima))) {
  samp <- subsample(n, seed = 1)
  assign("samp", samp, envir = globalenv())
  t0 <- proc.time()[3]
  b <- fit_mpl("samp", "status", TRUE)
  el <- proc.time()[3] - t0
  bref <- fit_ref(samp, "status", TRUE)
  cat(sprintf(
    "%8d %6d  %8.4f %9.4f %9.4f  %7.4f  %6.1f\n",
    n, sum(samp$status), b[1], b[2], b[3], max(abs(b - bref)), el
  ))
}

## ======================================================== 3. colon cancer ==
cat("\n== colon-cancer-specific outcome, full cohort ==\n")
assign("samp", hiroshima, envir = globalenv())
t0 <- proc.time()[3]
tab <- rbind(
  `MPL entry=` = fit_mpl("samp", "status_colon", TRUE),
  `coxph LTRC` = fit_ref(hiroshima, "status_colon", TRUE),
  `MPL naive` = fit_mpl("samp", "status_colon", FALSE),
  `coxph naive` = fit_ref(hiroshima, "status_colon", FALSE)
)
print(round(tab, 4))
cat(sprintf("(%.0fs)\n", proc.time()[3] - t0))
