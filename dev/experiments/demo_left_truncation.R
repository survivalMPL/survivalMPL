# Demonstration: left truncation (entry=) in coxph_mpl()
#
# A single simulation draw is not evidence either way (the naive estimator
# can land closer to the truth than the corrected one by chance), and if
# entry time is independent of the covariates, ignoring truncation barely
# biases beta at all (it mostly distorts the baseline hazard shape instead).
# So this demo:
#   1. Makes truncation INFORMATIVE: entry time depends on x1 (like a cohort
#      where more severe cases, x1=1, are enrolled later) - this is the
#      scenario that genuinely biases beta if truncation is ignored.
#   2. Runs many replications and reports average bias/RMSE, not one draw.
#
# Three fits per replication:
#   1. survival::coxph(Surv(entry, time, status) ~ ...) - trusted reference.
#   2. coxph_mpl(..., entry = entry) - our new argument.
#   3. coxph_mpl(...) with entry omitted - naive, ignores truncation.

devtools::load_all()
library(survival)

true_beta <- c(x1 = 0.7, x2 = -0.5)
n_per_rep <- 800
n_reps <- 100

simulate_one <- function(seed) {
  set.seed(seed)
  n <- n_per_rep
  x1 <- rbinom(n, 1, 0.5)
  x2 <- rnorm(n)
  lp <- true_beta["x1"] * x1 + true_beta["x2"] * x2

  event_time <- rexp(n, rate = 0.5 * exp(lp))

  # Informative truncation: x1 = 1 subjects enter later on average.
  entry_time <- runif(n, 0, 0.5 + 1.5 * x1)

  keep <- event_time > entry_time
  while (any(!keep)) {
    idx <- which(!keep)
    event_time[idx] <- rexp(length(idx), rate = 0.5 * exp(lp[idx]))
    keep <- event_time > entry_time
  }

  cens_time <- entry_time + rexp(n, rate = 0.3)
  time <- pmin(event_time, cens_time)
  status <- as.numeric(event_time <= cens_time)

  data.frame(time = time, status = status, entry = entry_time, x1 = x1, x2 = x2)
}

fit_all_three <- function(df) {
  fit_ref <- coxph(Surv(entry, time, status) ~ x1 + x2, data = df)
  fit_entry <- coxph_mpl(Surv(time, status) ~ x1 + x2,
    data = df, entry = entry, basis = "msplines", smooth = 0
  )
  fit_naive <- coxph_mpl(Surv(time, status) ~ x1 + x2,
    data = df, basis = "msplines", smooth = 0
  )
  list(
    reference = coef(fit_ref),
    entry     = coef(fit_entry, "Beta"),
    naive     = coef(fit_naive, "Beta")
  )
}

# Occasional replications hit a numerical failure in the naive (severely
# mis-specified) fit; skip those rather than letting one bad draw abort the
# whole Monte Carlo run.
fit_all_three_safe <- function(df) {
  tryCatch(fit_all_three(df), error = function(e) NULL)
}

# ---------------------------------------------------------------------------
# 1. One illustrative run, printed in full
# ---------------------------------------------------------------------------
one <- fit_all_three(simulate_one(2026))
cat("One replication (n =", n_per_rep, "), true beta:", true_beta, "\n\n")
cat("survival::coxph reference        :", one$reference, "\n")
cat("coxph_mpl WITH entry=            :", one$entry, "\n")
cat("coxph_mpl WITHOUT entry= (naive) :", one$naive, "\n\n")

# ---------------------------------------------------------------------------
# 2. Monte Carlo: average bias and RMSE over many replications
# ---------------------------------------------------------------------------
results <- lapply(seq_len(n_reps), function(i) fit_all_three_safe(simulate_one(1000 + i)))
n_failed <- sum(vapply(results, is.null, logical(1)))
results <- Filter(Negate(is.null), results)
if (n_failed > 0) {
  cat("(", n_failed, "of", n_reps, "replications skipped due to a numerical fit failure)\n\n")
}

entry_mat <- do.call(rbind, lapply(results, `[[`, "entry"))
naive_mat <- do.call(rbind, lapply(results, `[[`, "naive"))

bias_rmse <- function(mat, truth) {
  bias <- colMeans(mat) - truth
  rmse <- sqrt(colMeans((sweep(mat, 2, truth))^2))
  rbind(bias = bias, rmse = rmse)
}

cat("Monte Carlo over", n_reps, "replications (n =", n_per_rep, "each):\n\n")
cat("coxph_mpl WITH entry= vs true beta:\n")
print(bias_rmse(entry_mat, true_beta))
cat("\ncoxph_mpl WITHOUT entry= (naive) vs true beta:\n")
print(bias_rmse(naive_mat, true_beta))
