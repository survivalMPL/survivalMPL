#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# basis-bsplines.R — B-Splines basis (proof of extensibility)
#
# This is a skeleton demonstrating that adding a new basis type requires only:
#   1. Creating this file with the three required functions + a spec list.
#   2. Calling register_basis(.bsplines_spec) in zzz.R .onLoad().
#
# The implementations below delegate to splines::splineDesign() for the density
# basis and use a numerical integral for the penalty.  They are functional but
# not tuned for speed — the purpose is to show the registration mechanism works
# end-to-end.
#
# Augmented knot convention: order-1 boundary copies, giving
#   m = n.Alpha + order - 2 B-splines (same m as M-splines).
#
# Spec: .bsplines_spec
# Registered in zzz.R via register_basis(.bsplines_spec)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.bsplines_knots_fn <- function(control, events) {
  Alpha   <- .compute_alpha_knots(control, events)
  n.Alpha <- length(Alpha)
  m       <- n.Alpha + control$order - 2L
  list(m = m, Alpha = Alpha, Delta = rep(1, m))
}

.bsplines_matrix_fn <- function(x, knots, order, which = 1) {
  Alpha        <- knots$Alpha
  n.Alpha      <- length(Alpha)
  m            <- knots$m
  n            <- length(x)
  which.matrix <- c(1 %in% which, 2 %in% which)

  # order-1 boundary copies → length n.Alpha + 2*(order-1) → n.Alpha+order-2 = m B-splines
  knots_aug <- c(rep(Alpha[1L], order - 1L), Alpha,
                 rep(Alpha[n.Alpha], order - 1L))

  M_psi_nm <- matrix(0, n, m)
  M_Psi_nm <- matrix(0, n, m)

  if (which.matrix[1]) {
    M_psi_nm <- splines::splineDesign(knots_aug, x, ord = order,
                                      outer.ok = TRUE)
  }

  if (which.matrix[2]) {
    x_g <- seq(Alpha[1L], max(Alpha), length.out = 500L)
    B_g <- splines::splineDesign(knots_aug, x_g, ord = order,
                                 outer.ok = TRUE)
    for (j in seq_len(n)) {
      idx_g <- x_g <= x[j]
      if (any(idx_g)) {
        M_Psi_nm[j, ] <- apply(B_g[idx_g, , drop = FALSE], 2,
                                function(col) {
                                  xx <- x_g[idx_g]
                                  sum(diff(xx) * (col[-length(col)] + col[-1L]) / 2)
                                })
      }
    }
  }

  if (all(which.matrix)) list(psi = M_psi_nm, Psi = M_Psi_nm)
  else if (which.matrix[1]) M_psi_nm
  else M_Psi_nm
}

.bsplines_penalty_fn <- function(control, knots) {
  Alpha     <- knots$Alpha
  n.Alpha   <- length(Alpha)
  m         <- knots$m
  order     <- control$order

  # Same augmentation as matrix_fn
  knots_aug <- c(rep(Alpha[1L], order - 1L), Alpha,
                 rep(Alpha[n.Alpha], order - 1L))

  x_g <- seq(Alpha[1L], max(Alpha), length.out = 300L)
  h   <- diff(x_g)[1L]
  B2  <- splines::splineDesign(knots_aug, x_g, ord = order,
                               derivs = 2L, outer.ok = TRUE)
  # Gramian: M_R[i,j] = integral of b''_i(t) * b''_j(t) dt
  t(B2) %*% B2 * h
}

.bsplines_spec <- list(
  name             = "bsplines",
  aliases          = c("b", "bs", "bspline"),
  label            = "B-Splines",
  default_n_knots  = c(8L, 2L),
  penalty_order_fn = function(p, order) order - 1L,
  knots_fn         = .bsplines_knots_fn,
  matrix_fn        = .bsplines_matrix_fn,
  penalty_fn       = .bsplines_penalty_fn
)
