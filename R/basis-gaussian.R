#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# basis-gaussian.R — Gaussian kernel basis
#
# Spec: .gaussian_spec
# Registered in zzz.R via register_basis(.gaussian_spec)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.gaussian_knots_fn <- function(control, events) {
  Alpha    <- .compute_alpha_knots(control, events)
  n.Alpha  <- length(Alpha)
  n.events <- length(events)
  rng      <- range(events)

  Sigma <- Delta <- rep(0, n.Alpha)
  for (aw in seq_len(n.Alpha)) {
    if (aw > 1L && aw < (n.Alpha - control$n.knots[2])) {
      while (
        sum(events > (Alpha[aw] - 2 * Sigma[aw]) &
            events < (Alpha[aw] + 2 * Sigma[aw])) <
        n.events * control$cover.sigma.quant
      ) {
        Sigma[aw] <- Sigma[aw] + 0.001
      }
    } else {
      Sigma[aw] <- control$cover.sigma.fixed *
                   (Alpha[n.Alpha] - Alpha[1L]) / 3
    }
    Delta[aw] <- pnorm((rng[2L] - Alpha[aw]) / Sigma[aw]) -
                 pnorm((rng[1L] - control$epsilon[1L] - Alpha[aw]) / Sigma[aw])
  }
  list(m = n.Alpha, Alpha = Alpha, Sigma = Sigma, Delta = Delta)
}

.gaussian_matrix_fn <- function(x, knots, order, which = 1) {
  Alpha        <- knots$Alpha
  Sigma        <- knots$Sigma
  m            <- knots$m
  n            <- length(x)
  which.matrix <- c(1 %in% which, 2 %in% which)
  M_psi_nm     <- matrix(0, n, m)
  M_Psi_nm     <- matrix(0, n, m)

  for (u in seq_len(m)) {
    if (which.matrix[1]) {
      M_psi_nm[, u] <- dnorm((x - Alpha[u]) / Sigma[u]) /
                       (Sigma[u] * knots$Delta[u])
    }
    if (which.matrix[2]) {
      M_Psi_nm[, u] <- (pnorm((x - Alpha[u]) / Sigma[u]) -
                         pnorm((Alpha[1L] - Alpha[u]) / Sigma[u])) /
                       knots$Delta[u]
    }
  }

  if (all(which.matrix)) list(psi = M_psi_nm, Psi = M_Psi_nm)
  else if (which.matrix[1]) M_psi_nm
  else M_Psi_nm
}

# ── Gaussian penalty: analytical Gauss–Hermite integrals ─────────────────────

.gauss_int_rij_1 <- function(x, mu_i, mu_j, sig_i, sig_j, t1, tn) {
  K  <- 4 * (pnorm((t1  - mu_i) / sig_i) - pnorm((tn - mu_i) / sig_i)) *
            (pnorm((t1  - mu_j) / sig_j) - pnorm((tn - mu_j) / sig_j))
  q2 <- dnorm((mu_i - mu_j) / sqrt(sig_i^2 + sig_j^2)) * 2 * pi *
        sig_i * sig_j * (sig_i^2 + sig_j^2 - (mu_i - mu_j)^2) *
        (2 * pnorm(((x - mu_j) * sig_i^2 + (x - mu_i) * sig_j^2) /
                   (sig_i * sig_j * sqrt(sig_i^2 + sig_j^2))) - 1)
  q1 <- 4 * pi * dnorm((x - mu_i) / sig_i) * dnorm((x - mu_j) / sig_j) *
        sqrt(sig_i^2 + sig_j^2) *
        ((mu_i - x) * sig_i^2 + (mu_j - x) * sig_j^2)
  q3 <- pi * sig_i * sig_j * (sig_i^2 + sig_j^2)^(5 / 2) * K
  (q1 + q2) / q3
}

.gauss_int_rij_2 <- function(x, mu_i, mu_j, sig_i, sig_j, t1, tn) {
  K     <- 4 * (pnorm((t1  - mu_i) / sig_i) - pnorm((tn - mu_i) / sig_i)) *
               (pnorm((t1  - mu_j) / sig_j) - pnorm((tn - mu_j) / sig_j))
  q1q2a <- 4 * dnorm(x, mu_i, sig_i) * dnorm(x, mu_j, sig_j) *
           sig_i * sig_j * 2 * pi * sqrt(sig_i^2 + sig_j^2) *
           ((mu_j - x) * sig_i^6 * ((x - mu_i)^2 - sig_i^2) +
            sig_i^4 * sig_j^2 * ((x - 4*mu_j + 3*mu_i) * sig_i^2 -
                                   (x - mu_i)^2 * (3*x - 4*mu_j + mu_i)) -
            sig_i^2 * sig_j^4 * (x - mu_j)^2 * (3*x - 4*mu_i + mu_j) +
            sig_j^6 * ((x + 3*mu_j - 4*mu_i) * sig_i^2 -
                        (x - mu_j)^2 * (x - mu_i)) +
            sig_j^8 * (x - mu_i))
  q1q2b <- 2 * dnorm(mu_j, mu_i, sqrt(sig_i^2 + sig_j^2)) *
           pi * sqrt(sig_i^2 + sig_j^2) * sig_i^3 * sig_j^3 *
           (mu_j^4 - 4*mu_j^3*mu_i + mu_i^4 +
            6*mu_j^2 * (mu_i^2 - sig_i^2 - sig_j^2) -
            6*mu_i^2 * (sig_i^2 + sig_j^2) + 3*(sig_i^2 + sig_j^2)^2 -
            4*mu_i*mu_j * (mu_i^2 - 3*(sig_i^2 + sig_j^2))) *
           (2 * pnorm(((x - mu_j) * sig_i^2 + (x - mu_i) * sig_j^2) /
                      (sig_i * sig_j * sqrt(sig_i^2 + sig_j^2))) - 1)
  q3 <- pi * sig_i^3 * sig_j^3 * (sig_i^2 + sig_j^2)^(9 / 2) * K
  (q1q2a + q1q2b) / q3
}

.gaussian_penalty_fn <- function(control, knots) {
  m       <- knots$m
  penalty <- control$penalty
  t1      <- knots$Alpha[1L]
  tn      <- knots$Alpha[m]
  int_fn  <- if (penalty == 2L) .gauss_int_rij_2 else .gauss_int_rij_1

  M_R_mm <- matrix(0, m, m)
  for (i in seq_len(m)) {
    for (j in i:m) {
      M_R_mm[i, j] <- M_R_mm[j, i] <-
        int_fn(tn, knots$Alpha[i], knots$Alpha[j],
               knots$Sigma[i], knots$Sigma[j], t1, tn) -
        int_fn(t1, knots$Alpha[i], knots$Alpha[j],
               knots$Sigma[i], knots$Sigma[j], t1, tn)
    }
  }
  M_R_mm
}

.gaussian_spec <- list(
  name             = "gaussian",
  aliases          = c("g", "gauss"),
  label            = "Gaussian",
  default_n_knots  = c(8L, 2L),
  penalty_order_fn = function(p, order) { p <- as.integer(p); ifelse(p > 0L & p < 3L, p, 2L) },
  knots_fn         = .gaussian_knots_fn,
  matrix_fn        = .gaussian_matrix_fn,
  penalty_fn       = .gaussian_penalty_fn
)
