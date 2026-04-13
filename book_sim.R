library(survival)
library(survivalMPL)
library(splines2)
library(tictoc)

# Write simulate() function to repeatedly generate a sample
simulate <- function(n, pi_E, g1, g2) {
    # Simulate covariates
    x1 <- runif(n, 0, 1)
    x2 <- 5 * runif(n, 0, 1)
    X <- cbind(x1, x2)
    beta_true <- c(0.75, -0.5)
    eXtB <- exp(X %*% beta_true)
    # Simulate true event times
    neg.log.u <- -log(runif(n, 0, 1))
    y_i <- (neg.log.u / (eXtB))^(1 / 3)

    # Simulate censoring types & times
    u_E <- runif(n, 0, 1)
    u_L <- runif(n, 0, 1)
    u_R <- runif(n, u_L, 1)
    t_L <- (y_i^as.numeric(u_E < pi_E)) *
        (g1 * u_L)^as.numeric(pi_E < u_E &
            g1 * u_L <= y_i & y_i <= g2 * u_R) *
        (g2 * u_R)^as.numeric(pi_E < u_E & g2 * u_R < y_i) *
        0^as.numeric(pi_E < u_E & y_i < g1 * u_L)
    t_R <- (y_i^as.numeric(u_E < pi_E)) *
        (g1 * u_L)^as.numeric(pi_E < u_E & y_i < g1 * u_L) *
        (g2 * u_R)^as.numeric(pi_E < u_E &
            g1 * u_L <= y_i & y_i <= g2 * u_R) *
        Inf^(pi_E < u_E & g2 * u_R < y_i)

    # Compute midpoint and censoring type for midpoint imputation
    midpoint <- t_L + (t_R - t_L) / 2
    event <- rep(1, n)
    event[which(is.infinite(midpoint))] <- 0
    midpoint[which(is.infinite(midpoint))] <- t_L[which(is.infinite(midpoint))]
    df <- data.frame(t_L, t_R, x1, x2, midpoint, event)
    return(df)
}

# Run a simulation study
tic("R example: a simulation study (Old version)")
    ch2_ex26_save <- matrix(0, nrow = 100, ncol = 8)
    baseline_h0 <- NULL
    for (s in 1:100) {
        set.seed(s)
        # Generate data set
        dat <- simulate(n = 500, pi_E = 0.5, g1 = 0.9, g2 = 1.3)
        # Fit MPL model
        fit.mpl <- coxph_mpl(
            Surv(t_L, t_R, type = "interval2") ~
                x1 + x2,
                data = dat,
                basis = "m", 
                n.knots = c(5, 0)
        )
        # Save coefficient estimates & SE
        ch2_ex26_save[s, 1] <- fit.mpl$coef$Beta[1]
        ch2_ex26_save[s, 2] <- fit.mpl$coef$Beta[2]
        ch2_ex26_save[s, 3] <- fit.mpl$se$Beta$M2HM2[1]
        ch2_ex26_save[s, 4] <- fit.mpl$se$Beta$M2HM2[2]
        # # Fit PH model
        # fit.ph <- coxph(Surv(midpoint, event) ~ x1 + x2, data = dat)
        # # Save coefficient estimates & SE
        # ch2_ex26_save[s, 5] <- fit.ph$coefficients[1]
        # ch2_ex26_save[s, 6] <- fit.ph$coefficients[2]
        # ch2_ex26_save[s, 7] <- sqrt(diag(fit.ph$var))[1]
        # ch2_ex26_save[s, 8] <- sqrt(diag(fit.ph$var))[2]

        # psi <- mSpline(survfit(fit.ph)$time,
        #     knots = fit.mpl$knots$Alpha[2:6],
        #     Boundary.knots = c(
        #         fit.mpl$knots$Alpha[1],
        #         fit.mpl$knots$Alpha[7]
        #     )
        # )
        # h0t_MPL <- psi %*% fit.mpl$coef$Theta
        # baseline_h0 <- cbind(
        #     baseline_h0,
        #     basehaz(fit.ph)$time, basehaz(fit.ph)$hazard, h0t_MPL
        # )
    }
toc()



unloadNamespace("survivalMPL")
devtools::load_all()

tic("R example: a simulation study (New version)")
    ch2_ex26_save_new <- matrix(0, nrow = 100, ncol = 8)
    baseline_h0 <- NULL
    for (s in 1:100) {
        # Generate data set
        set.seed(s)
        dat <- simulate(n = 500, pi_E = 0.5, g1 = 0.9, g2 = 1.3)
        # Fit MPL model
        fit.mpl <- coxph_mpl(
            Surv(t_L, t_R, type = "interval2") ~
                x1 + x2,
                data = dat,
                basis = "m", 
                n.knots = c(5, 0)
        )
        # Save coefficient estimates & SE
        ch2_ex26_save_new[s, 1] <- fit.mpl$coef$Beta[1]
        ch2_ex26_save_new[s, 2] <- fit.mpl$coef$Beta[2]
        ch2_ex26_save_new[s, 3] <- fit.mpl$se$Beta$M2HM2[1]
        ch2_ex26_save_new[s, 4] <- fit.mpl$se$Beta$M2HM2[2]
        # # Fit PH model
        # fit.ph <- coxph(Surv(midpoint, event) ~ x1 + x2, data = dat)
        # # Save coefficient estimates & SE
        # ch2_ex26_save_new[s, 5] <- fit.ph$coefficients[1]
        # ch2_ex26_save_new[s, 6] <- fit.ph$coefficients[2]
        # ch2_ex26_save_new[s, 7] <- sqrt(diag(fit.ph$var))[1]
        # ch2_ex26_save_new[s, 8] <- sqrt(diag(fit.ph$var))[2]
        # psi <- mSpline(survfit(fit.ph)$time,
        #     knots = fit.mpl$knots$Alpha[2:6],
        #     Boundary.knots = c(
        #         fit.mpl$knots$Alpha[1],
        #         fit.mpl$knots$Alpha[7]
        #     )
        # )
        # h0t_MPL <- psi %*% fit.mpl$coef$Theta
        # baseline_h0 <- cbind(
        #     baseline_h0,
        #     basehaz(fit.ph)$time, basehaz(fit.ph)$hazard, h0t_MPL
        # )
    }
toc()




################################################################################
##### SIM 2: A pseudo melanoma study example ###################################
################################################################################
library(survival)
library(survivalMPL)
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
tic("A pseudo melanoma study example (Old version)")
for(i in c(1:6, 8:21, 23:40, 42:54, 56:63, 65:93, 95:106)){
set.seed(i)
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
}
toc()

unloadNamespace("survivalMPL")
devtools::load_all()
tic("A pseudo melanoma study example (New version)")
for(i in c(1:6, 8:21, 23:40, 42:54, 56:63, 65:93, 95:106)){
set.seed(i)
# print(i)
ch2_mel_dat <- simulate_mel(300, 0.37, 0.6, 1.2)
max(ch2_mel_dat$t_L)
mela_fit_new <- coxph_mpl(
    formula = Surv(t_L, t_R, type = "interval2") ~
        Arm + Leg + Trunk + mm1to2 + mm2to4 + mm4plus +
        Female + Age_centred, data = ch2_mel_dat,
    basis = "m", tol = 1e-05, n.knots = c(7, 0),
    max.iter = c(1000, 5000, 1e+05)
)
# summary(mela_fit_new)
# plot(mela_fit_new)
}
toc()