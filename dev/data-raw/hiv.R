## Pseudo-HIV seroconversion dataset for survivalMPL
##
## Simulates 300 left- and right-censored survival observations for the
## left-censoring tutorial.  The design mimics a cohort of people who inject
## drugs, enrolled in a needle and syringe programme and followed for HIV
## seroconversion:
##
##   * time origin  = start of injecting drug use (reported at enrolment);
##   * first HIV test happens S years later.  A subject already seropositive
##     at that test seroconverted somewhere in (0, S]: LEFT CENSORED;
##   * a subject negative at the first test is retested frequently, so a
##     seroconversion during follow-up is recorded as an EXACT event;
##   * a subject still negative when follow-up ends is RIGHT CENSORED.
##
## There is no interval censoring by construction, which is what makes this
## dataset a clean illustration of left censoring on its own.
##
## True baseline hazard is Weibull with shape 1.5 and scale 6:
##   H0(t) = (t/6)^1.5,  h0(t) = (1.5/6) * (t/6)^0.5.
##
## Run this script from the package root to regenerate data/hiv.rda:
##   Rscript dev/data-raw/hiv.R

simulate_hiv <- function(n, shape = 1.5, scale = 6,
                         screen.rate = 0.35, dropout.rate = 0.18,
                         max.followup = 8) {
  Sharing     <- rbinom(n, 1, 0.45)
  Prison      <- rbinom(n, 1, 0.30)
  Female      <- rbinom(n, 1, 0.35)
  age         <- rnorm(n, 32, 8)
  Age_centred <- (age - mean(age)) / 10

  X         <- cbind(Sharing, Prison, Female, Age_centred)
  beta_true <- c(0.90, 0.45, -0.25, -0.15)
  eXtB      <- as.numeric(exp(X %*% beta_true))

  # Weibull true seroconversion time under the Cox model
  neglogU <- -log(runif(n))
  y       <- scale * (neglogU / eXtB)^(1 / shape)

  # time of the first HIV test, and end of follow-up after that test
  S <- rexp(n, screen.rate)
  E <- pmin(S + rexp(n, dropout.rate), max.followup)

  left  <- y <= S
  exact <- !left & y <= E

  t_L <- ifelse(left, NA, ifelse(exact, y, E))
  t_R <- ifelse(left, S,  ifelse(exact, y, NA))

  data.frame(t_L, t_R, Sharing, Prison, Female, Age_centred)
}

set.seed(11)
hiv <- simulate_hiv(300)

save(hiv, file = "data/hiv.rda", compress = "bzip2")
message("Saved data/hiv.rda  (n = ", nrow(hiv), " rows, ",
        ncol(hiv), " columns)")
message("  left censored: ", sum(is.na(hiv$t_L)),
        " | exact: ", sum(hiv$t_L == hiv$t_R, na.rm = TRUE),
        " | right censored: ", sum(is.na(hiv$t_R)))
