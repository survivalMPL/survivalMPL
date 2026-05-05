##  non-finite interval Hessian term in coxph_mpl
#^ [1] 7 Error in chol.default(M_hessbeta_p1) : 
# ^  the leading minor of order 1 is not positive



## Representative values observed in debug:
## M_H1_ni1/M_H2_ni1 were huge, M_S1_ni1/M_S2_ni1 underflowed to 0,
## and M_S1mS2_ni1 was clamped to epsilon.
M_H1_ni1 <- c(3.802983e+153, 1.257154e+154, 7.032951e+153)
M_H2_ni1 <- c(4.094034e+154, 1.353367e+155, 7.571197e+154)
M_S1_ni1 <- c(0, 0, 0)
M_S2_ni1 <- c(0, 0, 0)
M_S1mS2_ni1 <- c(1e-10, 1e-10, 1e-10)

w_ni <- M_S1_ni1 * M_S2_ni1 / M_S1mS2_ni1^2 * (M_H2_ni1 - M_H1_ni1)^2 +
  (M_S1_ni1 * M_H1_ni1 - M_S2_ni1 * M_H2_ni1) / M_S1mS2_ni1
print(w_ni)



#FIX
w_ni_guarded <- ifelse(is.finite(w_ni), w_ni, 0)
print(w_ni_guarded)


###############################################################################
##  if (s_lik < s_lik_OLD)` fails
###############################################################################
#^ Error in if (s_lik < s_lik_OLD) { : missing value where TRUE/FALSE needed

## Values patterned after the debug session.
## s_lik = event_term - right_censored_term + left_censored_term +
##         interval_term - penalty_term
loglik_event_term <- -2.020082e+231
loglik_right_censored_term <- 4.377023e+230
loglik_left_censored_term <- 0
loglik_interval_term <- -4236.757
thetaRtheta <- NaN
s_lambda <- -49955049

penalty_term <- s_lambda * thetaRtheta

s_lik <- loglik_event_term -
  loglik_right_censored_term +
  loglik_left_censored_term +
  loglik_interval_term -
  penalty_term
s_lik_OLD <- 1.232505e+231




 is.finite(s_lik)

## comparison with NaN returns NA (not TRUE/FALSE).
(s_lik < s_lik_OLD)
#FIX
is.finite(s_lik) && is.finite(s_lik_OLD) && (s_lik < s_lik_OLD)



##  why a `thetaRtheta` form can become NaN
M_theta_m1 <- matrix(c(1e229, 1e229), ncol = 1)
M_R_mm <- matrix(c(
  1e245, -1e245,
  -1e245, 1e245
), nrow = 2, byrow = TRUE)
Rtheta <- M_R_mm %*% M_theta_m1
thetaRtheta <- t(M_theta_m1) %*% Rtheta

Rtheta
thetaRtheta
is.finite(thetaRtheta)