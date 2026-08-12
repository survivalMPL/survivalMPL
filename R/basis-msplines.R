#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# basis-msplines.R — M-Splines (I-Splines for cumulative) basis
#
# Spec: .msplines_spec
# Registered in zzz.R via register_basis(.msplines_spec)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.msplines_knots_fn <- function(control, events) {
  Alpha   <- .compute_alpha_knots(control, events)
  n.Alpha <- length(Alpha)
  m       <- n.Alpha + control$order - 2L
  list(m = m, Alpha = Alpha, Delta = rep(1, m))
}

# Private helper — used recursively for the cumulative (I-spline) computation.
# Computes only the density basis matrix (which=1 internally).
.msplines_density_matrix <- function(x, knots, order) {
  Alpha        <- knots$Alpha
  n.Alpha      <- length(Alpha)
  m            <- n.Alpha + order - 2L
  n            <- length(x)
  seq1n        <- seq_len(n)
  Alpha_star   <- c(rep(Alpha[1L], order - 1L), Alpha,
                    rep(Alpha[n.Alpha], order - 1L))
  M_psi_nm     <- matrix(0, n, m + 1L)   # +1 trimmed at end

  Alpha_star_x <- sapply(x, function(y, lim = Alpha[-1L]) sum(lim < y) + 1L) +
                  order - 1L

  # order-1 case: piecewise constant
  M_psi_nm[(Alpha_star_x - 1L) * n + seq1n] <-
    1 / (Alpha_star[Alpha_star_x + 1L] - Alpha_star[Alpha_star_x])

  if (order > 1L) {
    for (ow in 2L:order) {
      uw_x <- Alpha_star_x - ow + 1L
      for (pw in 0L:(ow - 1L)) {
        pos_x <- (uw_x + pw - 1L) * n + seq1n
        M_psi_nm[pos_x] <-
          (ow / ((ow - 1) * (Alpha_star[1:m + ow] - Alpha_star[1:m])))[uw_x + pw] *
          ((x - Alpha_star[uw_x + pw])     * M_psi_nm[pos_x] +
           (Alpha_star[uw_x + pw + ow] - x) * M_psi_nm[pos_x + n])
      }
    }
  }
  M_psi_nm[, 1:m, drop = FALSE]
}

.msplines_matrix_fn <- function(x, knots, order, which = 1) {
  Alpha        <- knots$Alpha
  n.Alpha      <- length(Alpha)
  m            <- knots$m
  n            <- length(x)
  seq1n        <- seq_len(n)
  which.matrix <- c(1 %in% which, 2 %in% which)
  M_psi_nm     <- matrix(0, n, m)
  M_Psi_nm     <- matrix(0, n, m)

  Alpha_star <- c(rep(Alpha[1L], order - 1L), Alpha,
                  rep(Alpha[n.Alpha], order - 1L))

  # ── Density (M-spline) ───────────────────────────────────────────────────
  if (which.matrix[1]) {
    M_psi_nm <- .msplines_density_matrix(x, knots, order)
  }

  # ── Cumulative (I-spline) ────────────────────────────────────────────────
  if (which.matrix[2]) {
    rank.x   <- rank(x)
    x_sorted <- x[order(x)]
    Alpha_x  <- sapply(x_sorted,
                       function(y, lim = Alpha[-1L]) sum(lim < y) + 1L)
    up_u     <- cumsum(tabulate(Alpha_x, n.Alpha - 1L))

    for (uw in 1L:(m - order + 1L)) {
      M_Psi_nm[min(n, up_u[uw] + 1L):n, uw] <- 1
    }

    Alpha_star2 <- c(rep(Alpha[1L], order), Alpha, rep(Alpha[n.Alpha], order))
    factor_v    <- c(
      (Alpha_star2[(order + 2L):length(Alpha_star2)] -
       Alpha_star2[1L:(length(Alpha_star2) - order - 1L)]) / (order + 1L),
      rep(0, order - 1L)
    )

    # Recursive: density at order+1 gives the I-spline increment
    M_psi2_nm <- cbind(
      .msplines_density_matrix(x_sorted, knots, order + 1L),
      matrix(0, n, order - 1L)
    )

    pos_xo  <- rep((Alpha_x - 1L) * n, 1L) + seq1n
    pos_xo1 <- rep(pos_xo, order) + rep(1L:order, each = n) * n

    for (ow in 0L:(order - 1L)) {
      M_Psi_nm[pos_xo + ow * n] <-
        apply(
          matrix(
            M_psi2_nm[pos_xo1 + ow * n] *
              factor_v[rep(Alpha_x, order) + rep((1L:order) + ow, each = n)],
            ncol = order
          ),
          1, sum
        )
    }
    M_Psi_nm <- M_Psi_nm[rank.x, , drop = FALSE]
  }

  if (all(which.matrix)) list(psi = M_psi_nm, Psi = M_Psi_nm)
  else if (which.matrix[1]) M_psi_nm
  else M_Psi_nm
}

.msplines_penalty_fn <- function(control, knots) {
  Alpha        <- knots$Alpha
  n.Alpha      <- length(Alpha)
  m            <- knots$m
  order        <- control$order
  Alpha_star   <- c(rep(Alpha[1L], order - 1L), Alpha,
                    rep(Alpha[n.Alpha], order - 1L))

  seq1n        <- 1L:(n.Alpha - 1L)
  Alpha_star_x <- sapply(Alpha[-1L],
                         function(y, lim = Alpha[-1L]) sum(lim < y) + 1L) +
                  order - 1L
  M_d2f_mm     <- matrix(0, n.Alpha - 1L, n.Alpha + order - 1L)

  M_d2f_mm[(Alpha_star_x - 1L) * (n.Alpha - 1L) + seq1n] <-
    1 / (Alpha_star[Alpha_star_x + 1L] - Alpha_star[Alpha_star_x])

  for (ow in 2L:order) {
    uw_x <- Alpha_star_x - ow + 1L
    for (pw in 0L:(ow - 1L)) {
      M_d2f_mm[(uw_x + pw - 1L) * (n.Alpha - 1L) + seq1n] <-
        (ow / (Alpha_star[1L:(n.Alpha + ow) + ow] -
               Alpha_star[1L:(n.Alpha + ow)]))[uw_x + pw] *
        (M_d2f_mm[(uw_x + pw - 1L) * (n.Alpha - 1L) + seq1n] -
         M_d2f_mm[(uw_x + pw) * (n.Alpha - 1L) + seq1n])
    }
  }

  M_d2f_mm  <- M_d2f_mm[, 1:m, drop = FALSE]
  M_R_mm    <- matrix(0, m, m)
  half_width <- Alpha[-1L] - Alpha[-n.Alpha]
  for (uw in 1:m) {
    for (vw in uw:m) {
      M_R_mm[uw, vw] <- M_R_mm[vw, uw] <-
        sum(M_d2f_mm[, uw] * M_d2f_mm[, vw] * half_width)
    }
  }
  M_R_mm
}

.msplines_spec <- list(
  name             = "msplines",
  aliases          = c("m", "ms"),
  label            = "M-Splines",
  default_n_knots  = c(8L, 2L),
  penalty_order_fn = function(p, order) order - 1L,
  knots_fn         = .msplines_knots_fn,
  matrix_fn        = .msplines_matrix_fn,
  penalty_fn       = .msplines_penalty_fn
)
