################################################################################
##### SIM 2: A pseudo melanoma study example ###################################
################################################################################
library(survival)
# library(survivalMPL)
devtools::load_all()
library(splines2)
library(tictoc)


simulate_mel <- function(n, pi_E, a1, a2) {
    neglogU <- -log(runif(n))
    location <- sample(c(1, 2, 3, 4), n,
        replace = TRUE,
        prob = c(0.2, 0.15, 0.3, 0.35)
    )
    Arm <- as.numeric(location == 2)
    Leg <- as.numeric(location == 3)
    Trunk <- as.numeric(location == 4)
    thickness <- sample(c(1, 2, 3, 4), n,
        replace = TRUE,
        prob = c(0.15, 0.4, 0.3, 0.15)
    )
    mm1to2 <- as.numeric(thickness == 2)
    mm2to4 <- as.numeric(thickness == 3)
    mm4plus <- as.numeric(thickness == 4)
    Female <- rbinom(n, 1, 0.4)
    x8 <- rnorm(n, 56, 12)
    x8 <- x8 - mean(x8)
    Age_centred <- x8 / 10
    X <- cbind(
        Arm, Leg, Trunk, mm1to2, mm2to4,
        mm4plus, Female, Age_centred
    )
    beta_true <- c(
        -0.56, 0.01, -0.22, 0.22, 0.87,
        1.13, -0.17, 0.14
    )
    eXtB <- exp(X %*% beta_true)
    y <- as.numeric((neglogU / (2 * eXtB))^(2))
    # uniform variables
    U_E <- runif(n)
    U_L <- runif(n, 0, 1)
    U_R <- runif(n, U_L, 1)
    t_L <- y^((U_E < pi_E)) * (a1 * U_L)^((pi_E <=
        U_E & a1 * U_L <= y & y <= a2 * U_R)) *
        (a2 * U_R)^((pi_E <= U_E & a2 * U_R < y)) *
        (0)^((pi_E <= U_E & y < a1 * U_L))
    t_R <- y^((U_E < pi_E)) * (a1 * U_L)^((pi_E <=
        U_E & y < a1 * U_L)) * (a2 * U_R)^((pi_E <=
        U_E & a1 * U_L <= y & y <= a2 * U_R)) *
        Inf^((pi_E <= U_E & a2 * U_R < y))
    df <- data.frame(t_L, t_R, X)
    return(df)
}
# ^breaks with set.seed(7) (41) 
# [1] 7 Error in chol.default(M_hessbeta_p1) : 
#   the leading minor of order 1 is not positive
# ^breaks with set.seed(630), (22), (55), (64), (94)
# Error in if (s_lik < s_lik_OLD) { : missing value where TRUE/FALSE needed
options(error = recover)

set.seed(41)
# print(i)
ch2_mel_dat <- simulate_mel(300, 0.37, 0.6, 1.2)
max(ch2_mel_dat$t_L)
mela_fit <- coxph_mpl(
    formula = Surv(t_L, t_R, type = "interval2") ~
        Arm + Leg + Trunk + mm1to2 + mm2to4 + mm4plus +
        Female + Age_centred, data = ch2_mel_dat,
    basis = "m", tol = 1e-05, n.knots = c(7, 0),
    max.iter = c(1000, 5000, 1e+05)
)
# summary(mela_fit)
# plot(mela_fit)



# unloadNamespace("survivalMPL")

library(survivalMPL)
fit <- coxph_mpl(Surv(left, right, type = "interval2") ~ treatment,
                 data = bcos2, basis = "msplines")
summary(fit)   # coefficients, lambda, degrees of freedom
plot(fit)      # bases, hazard, cumulative hazard, survival






# on.exit(options(old), add = TRUE)

# for (i in c(1:6, 8:21, 23:40, 42:54, 56:63, 65:93, 95:106)) {
#   message("i = ", i)
#   set.seed(41)

#   ch2_mel_dat <- simulate_mel(300, 0.37, 0.6, 1.2)

#   mela_fit <- coxph_mpl(
#     formula = Surv(t_L, t_R, type = "interval2") ~
#       Arm + Leg + Trunk + mm1to2 + mm2to4 + mm4plus +
#       Female + Age_centred,
#     data = ch2_mel_dat,
#     basis = "m", tol = 1e-05, n.knots = c(7, 0),
#     max.iter = c(1000, 5000, 1e+05)
#   )
# }

# which(!is.finite(M_hessbeta_p1), arr.ind = TRUE)[1:10, , drop = FALSE]

# w_no <- c(M_H_no1)
# w_nr <- c(M_H_nr1)
# w_nl <- c(M_S_nl1/(1 - M_S_nl1)^2 * M_H_nl1^2 -
#           M_S_nl1/(1 - M_S_nl1) * M_H_nl1)
# w_ni <- c(M_S1_ni1 * M_S2_ni1 / M_S1mS2_ni1^2 * (M_H2_ni1 - M_H1_ni1)^2 +
#           (M_S1_ni1 * M_H1_ni1 - M_S2_ni1 * M_H2_ni1) / M_S1mS2_ni1)

# chk <- function(z) c(
#   nonfinite = sum(!is.finite(z)),
#   finite = sum(is.finite(z)),
#   min = if (any(is.finite(z))) min(z[is.finite(z)]) else NA_real_,
#   max = if (any(is.finite(z))) max(z[is.finite(z)]) else NA_real_
# )

# sapply(list(w_no=w_no, w_nr=w_nr, w_nl=w_nl, w_ni=w_ni), chk)


# idx <- which(!is.finite(w_ni))
# head(cbind(
#   H1 = M_H1_ni1[idx], H2 = M_H2_ni1[idx],
#   S1 = M_S1_ni1[idx], S2 = M_S2_ni1[idx], D = M_S1mS2_ni1[idx]
# ), 20)


iter; k
s_lik
s_lik_OLD
is.finite(s_lik)
is.finite(s_lik_OLD)



lik_parts <- c(
  part1 = sum(log(M_mu_no1) + log(M_h0_no1) - M_H_no1),
  part2 = sum(M_H_nr1),
  part3 = sum(log(1 - M_S_nl1)),
  part4 = sum(log(M_S1mS2_ni1)),
  part5 = as.numeric(s_lambda * thetaRtheta)
)
lik_parts
is.finite(lik_parts)


c(
  s_lambda = as.numeric(s_lambda),
  thetaRtheta = as.numeric(thetaRtheta),
  penalty = as.numeric(s_lambda * thetaRtheta)
)
is.finite(s_lambda)
is.finite(thetaRtheta)



# 1) Is theta already contaminated?
sum(!is.finite(M_theta_m1))
range(M_theta_m1[is.finite(M_theta_m1)])

# 2) Is Rtheta contaminated?
sum(!is.finite(Rtheta))
range(Rtheta[is.finite(Rtheta)])

# 3) Recompute thetaRtheta explicitly
thetaRtheta_check <- t(M_theta_m1) %*% Rtheta
thetaRtheta_check
is.finite(thetaRtheta_check)

# 4) Check raw quadratic form against penalty matrix
quad_check <- t(M_theta_m1) %*% M_R_mm %*% M_theta_m1
quad_check
is.finite(quad_check)

# 5) Check smoothing update state
c(
  s_lambda = as.numeric(s_lambda),
  s_df = as.numeric(s_df)
)

# Size/spread of theta
summary(log10(M_theta_m1[M_theta_m1 > 0]))
sum(M_theta_m1 == control$epsilon[2])

# Which component is exploding
j <- which.max(M_theta_m1)
c(j = j, theta = M_theta_m1[j],
  theta_old = M_theta_m1_OLD[j],
  step = M_steptheta_p1[j],
  s_m = M_s_m1[j],
  gradA = M_gradthetaA_m1[j],
  gradB = M_gradthetaB_m1[j],
  grad = M_gradtheta_m1[j])

# Is division by tiny gradB creating huge multiplicative step?
range(M_gradthetaB_m1)
quantile(abs(M_gradthetaB_m1), probs = c(0, .01, .05, .5, .95, .99, 1))
sum(abs(M_gradthetaB_m1) < 1e-12)

# Smoothing/df state
c(s_lambda = as.numeric(s_lambda),
  s_lambda_old = as.numeric(s_lambda_old),
  s_sigma2_old = as.numeric(s_sigma2_old),
  s_df = as.numeric(s_df))



j <- which.max(M_theta_m1)

c(
  term1 = as.numeric(M_tpsi_nom[j, ] %*% (1 / M_h0_no1)),
  term2 = as.numeric(M_tPsi_nlm[j, ] %*% (M_S_nl1 * M_mu_nl1 / (1 - M_S_nl1))),
  term3 = as.numeric(M_tPsi2_nim[j, ] %*% (M_S2_ni1 * M_mu_ni1 / M_S1mS2_ni1)),
  term4 = as.numeric(-TwoLRtheta[j] * (TwoLRtheta[j] < 0) + 0.3),
  gradA = as.numeric(M_gradthetaA_m1[j])
)

range(M_h0_no1)
sum(M_h0_no1 <= 0)
summary(log10(1 / M_h0_no1))


trace_raw <- sum(diag(HRinv %*% M_Rstar_ll))
c(
  m = m,
  trace_raw = trace_raw,
  s_sigma2_old = s_sigma2_old,
  trace_scaled = trace_raw / s_sigma2_old,
  s_df = s_df
)


A <- H[pos, pos] + (1 / s_sigma2_old) * M_Rstar_ll[pos, pos]
A_sym <- (A + t(A)) / 2

ev <- eigen(A_sym, symmetric = TRUE, only.values = TRUE)$values
c(min_ev = min(ev), max_ev = max(ev), kappa = kappa(A_sym))

d <- diag(HRinv %*% M_Rstar_ll)
summary(d)
c(nonfinite = sum(!is.finite(d)), negative = sum(d < 0), trace_raw = sum(d))


exists("temp")
class(temp)
if (exists("temp")) c(nonfinite = sum(!is.finite(temp)))
