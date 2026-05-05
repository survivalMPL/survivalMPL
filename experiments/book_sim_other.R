# Example 2.6
# The book prints the midpoint assignment split across two lines:
# BOOK (broken — R would parse this as two statements, syntax error):
#* midpoint[which(is.infinite(midpoint))]
#* <- t_L[which(is.infinite(midpoint))]


# CORRECT (one line):
midpoint[which(is.infinite(midpoint))] <- t_L[which(is.infinite(midpoint))]


# Example 2.7: Missing type = in Surv() call

# Typos in the book

#* mela_fit = coxph_mpl(formula = Surv(t_L, t_R, "interval2") ~ ...

# CORRECT (named argument):
#* mela_fit <- coxph_mpl(formula = Surv(t_L, t_R, type = "interval2") ~ ...

library(survivalMPL)

simulate_mel <- function(n, pi_E, a1, a2) {
  neglogU <- -log(runif(n))

  location <- sample(c(1, 2, 3, 4), n, replace = TRUE,
    prob = c(0.2, 0.15, 0.3, 0.35))
  Arm   <- as.numeric(location == 2)
  Leg   <- as.numeric(location == 3)
  Trunk <- as.numeric(location == 4)

  thickness <- sample(c(1, 2, 3, 4), n, replace = TRUE,
    prob = c(0.15, 0.4, 0.3, 0.15))
  mm1to2  <- as.numeric(thickness == 2)
  mm2to4  <- as.numeric(thickness == 3)
  mm4plus <- as.numeric(thickness == 4)

  Female <- rbinom(n, 1, 0.4)

  x8 <- rnorm(n, 56, 12)
  x8 <- x8 - mean(x8)
  Age_centred <- x8/10

  X <- cbind(Arm, Leg, Trunk, mm1to2, mm2to4, mm4plus, Female, Age_centred)
  beta_true = c(-0.56, 0.01, -0.22, 0.22, 0.87, 1.13, -0.17, 0.14)
  eXtB = exp(X %*% beta_true)

  y <- as.numeric((neglogU/(2 * eXtB))^(2))

  # uniform variables
  U_E <- runif(n)
  U_L <- runif(n, 0, 1)
  U_R <- runif(n, U_L, 1)

  t_L <- y^((U_E < pi_E)) *
    (a1 * U_L)^((pi_E <= U_E & a1 * U_L <= y & y <= a2 * U_R)) *
    (a2 * U_R)^((pi_E <= U_E & a2 * U_R < y)) *
    (0)^((pi_E <= U_E & y < a1 * U_L))

  t_R <- y^((U_E < pi_E)) *
    (a1 * U_L)^((pi_E <= U_E & y < a1 * U_L)) *
    (a2 * U_R)^((pi_E <= U_E & a1 * U_L <= y & y <= a2 * U_R)) *
    Inf^((pi_E <= U_E & a2 * U_R < y))

  df = data.frame(t_L, t_R, X)
  return(df)
}

set.seed(1)
ch2_mel_dat <- simulate_mel(300, 0.37, 0.6, 1.2)
max(ch2_mel_dat$t_L)

mela_fit = coxph_mpl(
  formula = Surv(t_L, t_R, type = "interval2") ~
    Arm + Leg + Trunk + mm1to2 + mm2to4 + mm4plus + Female + Age_centred,
  data = ch2_mel_dat,
  basis = "m", tol = 1e-05, n.knots = c(7, 0),
  max.iter = c(1000, 5000, 1e+05)
)

summary(mela_fit)
plot(mela_fit, which = 2:3)