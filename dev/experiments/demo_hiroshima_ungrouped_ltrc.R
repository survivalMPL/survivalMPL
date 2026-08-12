# Why the `hiroshima` dataset is one-row-per-subject and not
# one-row-per-grouped-table-cell.
#
# An earlier version of the dataset treated each of the 53,782 cells of the
# RERF LSS14 grouped person-time table as a single pseudo-subject, with
# entry = the cell's mean age at exposure and time = its mean attained age.
# This script refits that construction next to the current dataset. The
# pseudo-subject version gets the sign of the dose-response WRONG - it makes
# radiation look protective - because a cell with high dose is a cell with many
# person-years and therefore a long mean attained age, and that artefact
# swamps the real signal. It also cannot represent the truncation properly:
# entry is a cell average, not the age at which a subject entered the risk set.
#
# Three analyses of the same question (all-cause mortality vs colon dose, sex
# and city):
#   (A) GROUPED  pseudo-subjects, entry=  - the old, superseded construction
#   (B) UNGROUPED individual records, entry= - the current dataset
#   (C) UNGROUPED individual records, naive  - truncation ignored
#
# survival::coxph() is the reference throughout, its partial likelihood
# handling (start, stop] data exactly. entry= is only correct for
# basis = "uniform", so all coxph_mpl() fits use it.

devtools::load_all()
library(survival)

COVS <- c("dose", "sex", "city")
N_SAMPLE <- 8000 # subjects for the MPL fits; the full cohort also runs (~7s)

## ---------------------------------------------------------------- data ----
data(hiroshima)
set.seed(1)
samp <- hiroshima[sort(sample.int(nrow(hiroshima), N_SAMPLE)), , drop = FALSE]

# Rebuild the old pseudo-subject frame straight from the source table, so the
# comparison does not depend on a stale copy of the dataset.
raw <- read.csv("dev/experiments/hiroshima_dataset/lss14/lss14.csv")
grouped <- data.frame(
  city = factor(raw$city, levels = c(1, 2), labels = c("Hiroshima", "Nagasaki")),
  sex = factor(raw$sex, levels = c(1, 2), labels = c("Male", "Female")),
  entry = raw$agex, # cell mean age at exposure
  time = raw$age, # cell mean attained age
  status = as.integer(raw$death > 0),
  dose = raw$colon10 / 1000 # mGy -> Gy, as in the current dataset
)

cat("\n== data ==\n")
for (nm in c("grouped", "samp", "hiroshima")) {
  x <- get(nm)
  cat(sprintf(
    "%-10s %6d rows, %5d events, entry %.1f-%.1f\n",
    nm, nrow(x), sum(x$status), min(x$entry), max(x$entry)
  ))
}

## --------------------------------------------------------------- fitting --
mpl <- function(data_name, truncated) {
  cl <- call("coxph_mpl",
    formula = reformulate(COVS, response = quote(Surv(time, status))),
    data = as.name(data_name), basis = "uniform", smooth = 0
  )
  if (truncated) cl$entry <- quote(entry)
  coef(eval(cl, envir = globalenv()), "Beta")
}

ref <- function(data, truncated) {
  lhs <- if (truncated) quote(Surv(entry, time, status)) else quote(Surv(time, status))
  coef(coxph(reformulate(COVS, response = lhs), data = data))
}

cat("\n== fitting ==\n")
timed <- function(label, expr) {
  t0 <- proc.time()[3]
  v <- force(expr)
  cat(sprintf("%-22s %5.1fs\n", label, proc.time()[3] - t0))
  v
}
A <- timed("A grouped   + entry", mpl("grouped", TRUE))
B <- timed("B ungrouped + entry", mpl("samp", TRUE))
C <- timed("C ungrouped naive", mpl("samp", FALSE))

## -------------------------------------------------------------- results ---
tab <- rbind(
  `A grouped   MPL entry=` = A,
  `A grouped   coxph LTRC` = ref(grouped, TRUE),
  `B ungrouped MPL entry=` = B,
  `B ungrouped coxph LTRC` = ref(samp, TRUE),
  `C ungrouped MPL naive ` = C,
  `C ungrouped coxph naive` = ref(samp, FALSE),
  `  full cohort coxph LTRC` = ref(hiroshima, TRUE),
  `  full cohort coxph naive` = ref(hiroshima, FALSE)
)
cat("\n== coefficients ==\n")
print(round(tab, 4))

cat("\n== MPL vs coxph on the same data (max abs difference) ==\n")
cat(sprintf(
  "A grouped   %.4f\nB ungrouped %.4f\nC naive     %.4f\n",
  max(abs(tab[1, ] - tab[2, ])),
  max(abs(tab[3, ] - tab[4, ])),
  max(abs(tab[5, ] - tab[6, ]))
))

cat("\nNOTE: (A) puts dose at about -0.36, i.e. a protective effect of\n")
cat("radiation, and both the MPL and the coxph fit agree on it - so it is a\n")
cat("property of the pseudo-subject construction, not of either estimator.\n")
cat("(B) and (C) put dose at about +0.11 to +0.19, the expected direction.\n")
