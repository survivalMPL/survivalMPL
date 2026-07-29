# Constructs the `hiroshima` package dataset from RERF's Life Span Study
# Report 14 grouped person-time table (dev/experiments/hiroshima_dataset/).
#
# The source file is a stratified table (city x sex x ground-distance x AHS x
# age-at-exposure-category x attained-age-category x calendar-time-category x
# dose-category), NOT individual-level records: each row gives counts/means
# for a stratum of subjects, not one subject. To use it as a left-truncation
# example for coxph_mpl() (which needs one row per subject), each stratum row
# is treated as a single pseudo-subject:
#   entry        = agex  (person-year weighted mean age at exposure in cell)
#   time         = age   (person-year weighted mean attained age in cell)
#   status       = 1 if any death (any cause) occurred in the cell, else 0
#   status_colon = 1 if any colon-cancer death occurred in the cell, else 0
#   dose         = colon10 (DS02 weighted colon dose, mGy)
# This discards the person-years/subject-count weighting of the original
# grouped table - it is a simplification for demonstration purposes, not a
# substitute for a proper grouped-data (Poisson) reanalysis of LSS Report 14.
#
# Verified before this simplification: age > agex and pyr > 0 for all 53,782
# rows in the source file, so every row satisfies coxph_mpl()'s requirement
# that entry precede the event/censoring time. status_colon has 584 rows with
# an event (rows with colon > 0), matching the count reported for this cohort
# elsewhere (621 is the sum of colon death counts across cells, which
# over-counts cells with more than one colon death).

raw <- read.csv("dev/experiments/hiroshima_dataset/lss14/lss14.csv")

hiroshima <- data.frame(
  city         = factor(raw$city, levels = c(1, 2), labels = c("Hiroshima", "Nagasaki")),
  sex          = factor(raw$sex, levels = c(1, 2), labels = c("Male", "Female")),
  gd3          = factor(raw$gd3, levels = c(1, 2), labels = c("<3km", "3-10km")),
  entry        = raw$agex,
  time         = raw$age,
  status       = as.integer(raw$death > 0),
  status_colon = as.integer(raw$colon > 0),
  dose         = raw$colon10
)

stopifnot(all(hiroshima$time > hiroshima$entry))

save(hiroshima, file = "data/hiroshima.rda", compress = "bzip2")
message(
  "Saved data/hiroshima.rda  (n = ", nrow(hiroshima), " rows, ",
  ncol(hiroshima), " columns)"
)
