#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# basis-uniform.R — Uniform (piecewise-constant) basis
#
# Spec: .uniform_spec
# Registered in zzz.R via register_basis(.uniform_spec)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.uniform_knots_fn <- function(control, events) {
  Alpha <- .compute_alpha_knots(control, events)
  m     <- length(Alpha) - 1L
  Delta <- Alpha[2L:(m + 1L)] - Alpha[1L:m]
  list(m = m, Alpha = Alpha, Delta = Delta)
}

.uniform_matrix_fn <- function(x, knots, order, which = 1) {
  Alpha        <- knots$Alpha
  Delta        <- knots$Delta
  m            <- knots$m
  n            <- length(x)
  M_psi_nm     <- matrix(0, n, m)
  M_Psi_nm     <- matrix(0, n, m)
  which.matrix <- c(1 %in% which, 2 %in% which)

  u_i <- sapply(x, function(y, lim = Alpha[-1L]) sum(lim < y) + 1L)
  for (i in seq_len(n)) {
    if (which.matrix[1]) M_psi_nm[i, u_i[i]] <- 1
    if (which.matrix[2]) {
      M_Psi_nm[i, 1:u_i[i]] <- c(
        if (u_i[i] > 1L) Delta[1L:(u_i[i] - 1L)],
        x[i] - Alpha[u_i[i]]
      )
    }
  }

  if (all(which.matrix)) list(psi = M_psi_nm, Psi = M_Psi_nm)
  else if (which.matrix[1]) M_psi_nm
  else M_Psi_nm
}

.uniform_penalty_fn <- function(control, knots) {
  m       <- knots$m
  penalty <- control$penalty
  D       <- diag(m) * c(1, -2)[penalty]
  E       <- diag(m + 1L)[-(m + 1L), -1L] * c(-1, 1)[penalty]
  B       <- D + E + list(0, t(E))[[penalty]]
  B[c(m + 1L, (m - 1L) * m)[seq_len(penalty)]] <- c(-1, 2)[penalty]
  t(B) %*% B
}

.uniform_spec <- list(
  name             = "uniform",
  aliases          = c("u", "uni"),
  label            = "Uniform",
  default_n_knots  = c(8L, 2L),
  penalty_order_fn = function(p, order) { p <- as.integer(p); ifelse(p > 0L & p < 3L, p, 2L) },
  knots_fn         = .uniform_knots_fn,
  matrix_fn        = .uniform_matrix_fn,
  penalty_fn       = .uniform_penalty_fn
)
