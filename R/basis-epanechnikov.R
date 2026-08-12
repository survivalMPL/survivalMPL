#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# basis-epanechnikov.R — Epanechnikov kernel basis
#
# Spec: .epanechnikov_spec
# Registered in zzz.R via register_basis(.epanechnikov_spec)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Shares the same knot structure as M-Splines: m = n.Alpha + order - 2
.epanechnikov_knots_fn <- function(control, events) {
  Alpha   <- .compute_alpha_knots(control, events)
  n.Alpha <- length(Alpha)
  m       <- n.Alpha + control$order - 2L
  list(m = m, Alpha = Alpha, Delta = rep(1, m))
}

.epanechnikov_matrix_fn <- function(x, knots, order, which = 1) {
  Alpha        <- knots$Alpha
  n.Alpha      <- length(Alpha)
  m            <- knots$m
  n            <- length(x)
  seq1n        <- seq_len(n)
  which.matrix <- c(1 %in% which, 2 %in% which)
  M_psi_nm     <- matrix(0, n, m + 1L)  # +1 trimmed at end
  M_Psi_nm     <- matrix(0, n, m + 1L)

  Alpha_star   <- c(rep(Alpha[1L], order - 1L), Alpha,
                    rep(Alpha[n.Alpha], order - 1L))

  # ── Density (Epanechnikov kernel) ─────────────────────────────────────────
  if (which.matrix[1]) {
    Alpha_star_x <- sapply(x,
                           function(y, lim = Alpha[-1L]) sum(lim < y) + 1L) +
                    order - 1L
    uw_x <- Alpha_star_x - order + 1L

    for (pw in 0L:(order - 1L)) {
      pos_x   <- (uw_x + pw - 1L) * n + seq1n
      pos_1   <- (uw_x + pw) == 1L
      pos_m   <- (uw_x + pw) == m
      pos_mid <- !pos_1 & !pos_m

      M_psi_nm[pos_x[pos_mid]] <-
        (6 * (x - Alpha_star[uw_x + pw]) *
             (x - Alpha_star[uw_x + pw + order]) /
         (Alpha_star[uw_x + pw] - Alpha_star[uw_x + pw + order])^3)[pos_mid]

      M_psi_nm[pos_x[pos_1]] <-
        (12 * (x - Alpha_star[uw_x + pw + order]) *
              (x - 2 * Alpha_star[uw_x + pw] + Alpha_star[uw_x + pw + order]) /
         (2 * Alpha_star[uw_x + pw] - 2 * Alpha_star[uw_x + pw + order])^3)[pos_1]

      M_psi_nm[pos_x[pos_m]] <-
        (12 * (x - Alpha_star[uw_x + pw]) *
              (x + Alpha_star[uw_x + pw] - 2 * Alpha_star[uw_x + pw + order]) /
         (2 * Alpha_star[uw_x + pw] - 2 * Alpha_star[uw_x + pw + order])^3)[pos_m]
    }
    M_psi_nm <- M_psi_nm[, 1:m, drop = FALSE]
  }

  # ── Cumulative ────────────────────────────────────────────────────────────
  if (which.matrix[2]) {
    rank.x   <- rank(x)
    x_sorted <- x[order(x)]
    Alpha_x  <- sapply(x_sorted,
                       function(y, lim = Alpha[-1L]) sum(lim < y) + 1L)
    up_u     <- cumsum(tabulate(Alpha_x, n.Alpha - 1L))

    for (uw in 1L:(m - order + 1L)) {
      M_Psi_nm[min(n, up_u[uw] + 1L):n, uw] <- 1
    }

    Alpha_star_x <- sapply(x_sorted,
                           function(y, lim = Alpha[-1L]) sum(lim < y) + 1L) +
                    order - 1L
    uw_x <- Alpha_star_x - order + 1L

    for (pw in 0L:(order - 1L)) {
      pos_x   <- (uw_x + pw - 1L) * n + seq1n
      pos_1   <- (uw_x + pw) == 1L
      pos_m   <- (uw_x + pw) == m
      pos_mid <- !pos_1 & !pos_m

      M_Psi_nm[pos_x[pos_mid]] <-
        ((x_sorted - Alpha_star[uw_x + pw])^2 *
         (2 * x_sorted + Alpha_star[uw_x + pw] - 3 * Alpha_star[uw_x + pw + order]) /
         (Alpha_star[uw_x + pw] - Alpha_star[uw_x + pw + order])^3)[pos_mid]

      M_Psi_nm[pos_x[pos_1]] <-
        ((x_sorted - Alpha_star[uw_x + pw]) *
         (x_sorted^2 - 2 * x_sorted * Alpha_star[uw_x + pw] -
          2 * Alpha_star[uw_x + pw]^2 +
          6 * Alpha_star[uw_x + pw] * Alpha_star[uw_x + pw + order] -
          3 * Alpha_star[uw_x + pw + order]^2) /
         (2 * (Alpha_star[uw_x + pw] - Alpha_star[uw_x + pw + order])^3))[pos_1]

      M_Psi_nm[pos_x[pos_m]] <-
        ((x_sorted - Alpha_star[uw_x + pw])^2 *
         (x_sorted + 2 * Alpha_star[uw_x + pw] -
          3 * Alpha_star[uw_x + pw + order]) /
         (2 * (Alpha_star[uw_x + pw] - Alpha_star[uw_x + pw + order])^3))[pos_m]
    }
    M_Psi_nm <- M_Psi_nm[rank.x, 1:m, drop = FALSE]
  }

  M_psi_nm <- M_psi_nm[, 1:m, drop = FALSE]

  if (all(which.matrix)) list(psi = M_psi_nm, Psi = M_Psi_nm)
  else if (which.matrix[1]) M_psi_nm
  else M_Psi_nm
}

.epanechnikov_penalty_fn <- function(control, knots) {
  Alpha      <- knots$Alpha
  n.Alpha    <- length(Alpha)
  m          <- knots$m
  order      <- control$order
  Alpha_star <- c(rep(Alpha[1L], order - 1L), Alpha,
                  rep(Alpha[n.Alpha], order - 1L))
  M_R_mm     <- matrix(0, m, m)

  for (uw in seq_len(m)) {
    f_u <- ifelse(uw == 1L | uw == m, 4, 1)
    for (vw in uw:m) {
      if (Alpha_star[vw] < Alpha_star[uw + order]) {
        f_v <- ifelse(vw == 1L | vw == m, 4, 1)
        M_R_mm[uw, vw] <- M_R_mm[vw, uw] <-
          (144 * (Alpha_star[uw + order] - Alpha_star[vw])) /
          (f_v * f_u *
           (Alpha_star[uw + order] - Alpha_star[uw])^3 *
           (Alpha_star[vw + order] - Alpha_star[vw])^3)
      }
    }
  }
  M_R_mm
}

.epanechnikov_spec <- list(
  name             = "epanechikov",      # canonical (original spelling kept)
  aliases          = c("e", "epa", "epanechnikov"),
  label            = "Epanechnikov",
  default_n_knots  = c(8L, 2L),
  penalty_order_fn = function(p, order) 2L,
  knots_fn         = .epanechnikov_knots_fn,
  matrix_fn        = .epanechnikov_matrix_fn,
  penalty_fn       = .epanechnikov_penalty_fn
)
