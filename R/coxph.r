#' Fit Cox Proportional Hazards Regression Model Via MPL
#'
#' Simultaneously estimate the regression coefficients and the baseline hazard
#' function of proportional hazard Cox models using maximum penalised
#' likelihood (MPL).
#'
#' \code{coxph_mpl} fits a Cox proportional hazards model allowing right, left,
#' and interval censoring by maximising a penalised likelihood in which a
#' penalty term smooths the baseline hazard estimate. Optimisation combines a
#' Newton step for regression coefficients with a multiplicative step for the
#' baseline hazard parameters while enforcing non-negativity constraints (see
#' Ma, Couturier, Heritier and Marschner (2021)). The covariate matrix is
#' centred during optimisation; baseline estimates and covariance matrices are
#' corrected afterwards using a delta-method adjustment.
#'
#' @param formula A survival formula with the response on the left-hand side and
#'   covariates on the right. The response must be built with
#'   [survival::Surv()] using \code{type = "interval2"} for interval-censored
#'   data (right-censored responses are converted internally).
#' @param data Optional data frame in which to evaluate \code{formula}.
#' @param subset Optional expression specifying a subset of observations to use.
#' @param na.action Optional missing-data filter function applied to the model
#'   frame; defaults to \code{options()\$na.action}.
#' @param control Optional list returned by [coxph_mpl.control()] specifying
#'   basis choice, smoothing value, iteration limits, and related options. When
#'   missing, defaults are built from \code{...}.
#' @param entry Optional numeric vector (or a column reference evaluated in
#'   \code{data}) giving each subject's left-truncation (delayed entry) time.
#'   Every value must be strictly less than that subject's event/interval
#'   lower bound; violations raise an error. Left-truncated and
#'   non-left-truncated subjects may be mixed in the same call. Requires
#'   \code{basis = "uniform"}; any other basis is replaced by \code{"uniform"}
#'   with a warning. Defaults to
#'   \code{NULL} (no truncation), which reproduces prior behaviour exactly.
#' @param ... Additional arguments passed to [coxph_mpl.control()].
#'
#' @return An object of class \code{"coxph_mpl"}; see [coxph_mpl.object] for
#'   components.
#'
#' @section Limitations: \code{entry} is only supported for
#'   \code{basis = "uniform"}; another basis is downgraded to it with a
#'   warning. Also, [residuals.coxph_mpl()] and
#'   [predict.coxph_mpl()] do not yet account for \code{entry} — they compute
#'   cumulative hazard and survival from time 0 rather than from each subject's
#'   entry time. Both are known follow-ups, not yet implemented.
#' @seealso [coxph_mpl.object()], [coxph_mpl.control()], [summary.coxph_mpl()],
#'   [plot.coxph_mpl()], [predict.coxph_mpl()]
#' @examples
#' \dontrun{
#' ## Right-censored example: survival::lung
#' data(lung, package = "survival")
#' fit_mpl <- coxph_mpl(Surv(time, status == 2) ~ age + sex + ph.karno +
#'  wt.loss, data = lung, tol = 1e-05)
#' summary(fit_mpl)
#'
#' ## Interval-censored example: bcos2
#' data(bcos2)
#' fit_mpl <- coxph_mpl(Surv(left, right, type = "interval2") ~ treatment,
#'                      data = bcos2, basis = "m")
#' summary(fit_mpl)
#' }
#' @import survival
#' @importFrom MASS ginv
#' @importFrom stats contrasts dnorm model.extract model.matrix pnorm quantile runif terms
#' @export
coxph_mpl <- function(formula, data, subset, na.action, control, entry, ...) {
  # --- build model frame and response (Surv) ---
  mc <- match.call(expand.dots = FALSE)
  # `entry` is resolved from `data` via the same generic mechanism
  # model.frame.default() uses for `weights`/`offset` in stats::lm(): any extra
  # named argument not among its own formals is evaluated against `data`.
  m <- match(c("formula", "data", "subset", "na.action", "entry"), names(mc), 0)
  mc <- mc[c(1, m)]
  if (m[1] == 0) stop("A formula argument is required")
  data.name <- if (m[2] != 0) {
    deparse(match.call()[[3]])
  } else {
    "-"
  }
  mc[[1]] <- as.name("model.frame")
  mc$formula <- if (missing(data)) {
    terms(formula)
  } else {
    terms(formula, data = data)
  }
  mf <- eval(mc, parent.frame())
  if (any(is.na(mf))) stop("Missing observations in the model variables")
  if (nrow(mf) == 0) stop("No (non-missing) observations")
  mt <- attr(mf, "terms")
  # Y
  y <- model.extract(mf, "response")
  type <- attr(y, "type")
  if (!inherits(y, "Surv")) {
    stop("Response must be a survival object")
  }
  ##
  if (attr(y, which = "type") == "right") {
    left <- y[, 1]
    right <- rep(NA, nrow(y))
    icase <- which(y[, 2] == 1)
    right[icase] <- y[icase, 1]
    y <- Surv(left, right, type = "interval2")
  } else if (type != "interval") {
    stop("\nPlease create the survival object using the option type='interval2' in the Surv function.\n")
  }
  ##
  # --- classify observations by censoring type (right/event/left/interval) ---
  t_i1 <- y[, 1L]
  t_i2 <- y[, 2L]
  n <- length(t_i1)
  ctype <- matrix(NA, nrow = n, ncol = 4)
  colnames(ctype) <- c("r", "e", "l", "i")
  for (tw in 1:4) {
    ctype[, tw] <- y[, 3L] == (tw - 1)
  }
  n.ctype <- apply(ctype, 2, sum)
  ctypeTF <- n.ctype > 0
  observed <- y[, 3L] == 1L
  n.obs <- sum(y[, 3L] != 0)
  # `entry` was already filtered by subset/na.action together with everything
  # else inside model.frame(), so no separate row-alignment is needed here.
  # model.frame() stores extras under a parenthesized name (same convention as
  # "(weights)" in stats::lm()), not the bare name.
  entry <- mf[["(entry)"]]
  if (!is.null(entry)) {
    entry <- as.numeric(entry)
    if (any(entry >= t_i1)) {
      stop("entry time must precede the observed interval/event time")
    }
  }
  # --- control arguments and tie handling ---
  extraArgs <- list(...)
  if (length(extraArgs)) {
    controlargs <- names(formals(coxph_mpl.control))
    m <- pmatch(names(extraArgs), controlargs, nomatch = 0L)
    if (any(m == 0L)) {
      stop(gettextf("Argument(s) %s not matched", names(extraArgs)[m == 0L]),
        domain = NA, call. = FALSE
      )
    }
  }
  if (missing(control)) control <- coxph_mpl.control(n.obs, ...)

  # Left truncation is implemented by differencing the cumulative basis,
  # H(t) - H(entry), which the package supports for the piecewise-constant
  # "uniform" basis only.  Rather than refusing the fit, fall back to
  # "uniform" and say so.
  if (!is.null(entry) && control$basis != "uniform") {
    dropped.basis <- control$basis
    warning(gettextf(
      'entry= (left truncation) is only implemented for basis = "uniform"; basis "%s" was replaced by "uniform".',
      dropped.basis
    ), domain = NA, call. = FALSE)
    # c(0, 20) is the coxph_mpl.control() default for the kernel bases, never a
    # knot count chosen for a step function: revert it to the uniform default.
    if (identical(as.numeric(control$n.knots), c(0, 20))) {
      control$n.knots <- c(8, 2)
    }
    control$basis   <- "uniform"
    control$penalty <- compute_penalty_order("uniform", control$penalty,
                                             control$order)
  }

  # ties
  t_i1.obs <- t_i1[observed]
  ties <- duplicated(t_i1.obs)
  if (any(ties)) {
    if (control$ties == "epsilon") {
      if (length(control$seed) > 0) {
        old <- .Random.seed
        on.exit(.Random.seed <<- old)
        set.seed(control$seed)
      }
      t_i1.obs[ties] <- t_i1.obs[ties] + runif(sum(ties), -1e-11, 1e-11)
      t_i1[observed] <- t_i1.obs
    } else {
      t_i1.obs <- t_i1.obs[!ties]
      n.obs <- length(t_i1.obs)
    }
  }
  # --- design matrix, centering, knots and basis matrices ---
  # X
  X <- model.matrix(mt, mf) # , contrasts)
  X <- X[, !apply(X, 2, function(x) all(x == x[1])), drop = FALSE]
  if (ncol(X) == 0) {
    X <- matrix(0, n, 1)
    noX <- TRUE
  } else {
    noX <- FALSE
  }
  p <- ncol(X)
  mean_j <- apply(X, 2, mean)
  XC <- X - rep(mean_j, each = n)
  # knot sequence and psi matrices
  knots <- compute_knots(
    control,
    c(
      t_i1[ctype[, "i"]],
      t_i2[ctype[, "i"]],
      t_i1[ctype[, "e"]] - 1e-3,
      t_i1[ctype[, "e"]] + 1e-3,
      t_i1[ctype[, "r"]],
      t_i1[ctype[, "l"]],
      # extend knot support down to min(entry): otherwise Psi(entry) can fall
      # outside basis support and truncation silently becomes a no-op.
      if (!is.null(entry)) entry
    )
  )

  ###
  ### Estimation
  ###
  m <- knots$m
  K <- control$max.iter
  s_lambda <- control$smooth
  s_kappa <- control$kappa
  s_t1 <- knots$Alpha[1]
  s_tn <- max(knots$Alpha)
  M_R_mm <- compute_penalty(control, knots)
  M_Rstar_ll <- rbind(matrix(0, p, p + m), cbind(matrix(0, m, p), M_R_mm))
  s_convlimit <- control$tol
  M_X_nop <- XC[ctype[, 2], , drop = FALSE]
  # M_tX_nop    = t(M_X_nop)  # OLD:  replaced with crossprod()
  M_psi_nom <- compute_basis_matrix(t_i1, knots, control$basis,
    control$order,
    which = 1
  )[ctype[, 2], , drop = FALSE]
  # density matrix (which=1): never differenced against entry — truncation
  # changes the cumulative hazard, not the hazard density h0(t).
  M_tpsi_nom <- t(M_psi_nom)
  M_Psi_nom <- compute_basis_matrix(t_i1, knots, control$basis,
    control$order,
    which = 2
  )[ctype[, 2], , drop = FALSE]
  if (!is.null(entry)) {
    # H*(t) = H(t) - H(entry) = (Psi(t) - Psi(entry)) %*% theta
    M_Psi_nom <- M_Psi_nom - compute_basis_matrix(entry, knots, control$basis,
      control$order,
      which = 2
    )[ctype[, 2], , drop = FALSE]
  }
  M_tPsi_nom <- t(M_Psi_nom)
  M_X_nrp <- XC[ctype[, 1], , drop = FALSE]
  # M_tX_nrp    = t(M_X_nrp)  # OLD:  replaced with crossprod()
  M_Psi_nrm <- compute_basis_matrix(t_i1, knots, control$basis,
    control$order,
    which = 2
  )[ctype[, 1], , drop = FALSE]
  if (!is.null(entry)) {
    M_Psi_nrm <- M_Psi_nrm - compute_basis_matrix(entry, knots, control$basis,
      control$order,
      which = 2
    )[ctype[, 1], , drop = FALSE]
  }
  M_tPsi_nrm <- t(M_Psi_nrm)
  M_X_nlp <- XC[ctype[, 3], , drop = FALSE]
  # M_tX_nlp    = t(M_X_nlp)  # OLD:  replaced with crossprod()
  M_Psi_nlm <- compute_basis_matrix(t_i1, knots, control$basis,
    control$order,
    which = 2
  )[ctype[, 3], , drop = FALSE]
  if (!is.null(entry)) {
    M_Psi_nlm <- M_Psi_nlm - compute_basis_matrix(entry, knots, control$basis,
      control$order,
      which = 2
    )[ctype[, 3], , drop = FALSE]
  }
  M_tPsi_nlm <- t(M_Psi_nlm)
  M_X_nip <- XC[ctype[, 4], , drop = FALSE]
  # M_tX_nip    = t(M_X_nip)  # OLD:  replaced with crossprod()
  M_Psi1_nim <- compute_basis_matrix(t_i1, knots, control$basis,
    control$order,
    which = 2
  )[ctype[, 4], , drop = FALSE]
  M_Psi2_nim <- compute_basis_matrix(t_i2, knots, control$basis,
    control$order,
    which = 2
  )[ctype[, 4], , drop = FALSE]
  if (!is.null(entry)) {
    # same entry-time basis reused for both interval bounds
    M_Psi_entry_nim <- compute_basis_matrix(entry, knots, control$basis,
      control$order,
      which = 2
    )[ctype[, 4], , drop = FALSE]
    M_Psi1_nim <- M_Psi1_nim - M_Psi_entry_nim
    M_Psi2_nim <- M_Psi2_nim - M_Psi_entry_nim
  }
  M_tPsi1_nim <- t(M_Psi1_nim)
  M_tPsi2_nim <- t(M_Psi2_nim)

  # --- initialize parameters and derived hazard/survival quantities ---
  M_beta_p1 <- matrix(0, nrow = p, ncol = 1)
  M_theta_m1 <- matrix(1, nrow = m, ncol = 1)
  s_df <- -1

  ## shortcuts
  M_mu_no1 <- exp(M_X_nop %*% M_beta_p1)
  M_mu_nr1 <- exp(M_X_nrp %*% M_beta_p1)
  M_mu_nl1 <- exp(M_X_nlp %*% M_beta_p1)
  M_mu_ni1 <- exp(M_X_nip %*% M_beta_p1)
  M_h0_no1 <- M_psi_nom %*% M_theta_m1
  M_H0_no1 <- M_Psi_nom %*% M_theta_m1
  M_H0_nr1 <- M_Psi_nrm %*% M_theta_m1
  M_H0_nl1 <- M_Psi_nlm %*% M_theta_m1
  M_H01_ni1 <- M_Psi1_nim %*% M_theta_m1
  M_H02_ni1 <- M_Psi2_nim %*% M_theta_m1
  Rtheta <- M_R_mm %*% M_theta_m1
  thetaRtheta <- t(M_theta_m1) %*% Rtheta
  TwoLRtheta <- s_lambda * 2 * Rtheta
  M_H_no1 <- M_H0_no1 * M_mu_no1
  M_H_nr1 <- M_H0_nr1 * M_mu_nr1
  M_H_nl1 <- M_H0_nl1 * M_mu_nl1
  M_H1_ni1 <- M_H01_ni1 * M_mu_ni1
  M_H2_ni1 <- M_H02_ni1 * M_mu_ni1
  M_S_no1 <- exp(-M_H_no1)
  M_S_nl1 <- exp(-M_H_nl1)
  M_S1_ni1 <- exp(-M_H1_ni1)
  M_S2_ni1 <- exp(-M_H2_ni1)
  # avoid division by 0
  M_S_nl1[M_S_nl1 == 1] <- 1 - control$epsilon[1]
  M_S1_ni1[M_S1_ni1 == 1] <- 1 - control$epsilon[1]
  M_S2_ni1[M_S2_ni1 == 1] <- 1 - control$epsilon[1]
  M_S1mS2_ni1 <- M_S1_ni1 - M_S2_ni1
  M_S1mS2_ni1[M_S1mS2_ni1 < control$epsilon[2]] <- control$epsilon[2]


  # --- outer loop over smoothing/df iterations ---

  full.iter <- 0
  K <- ifelse(control$max.iter[1] > 1, control$max.iter[2], control$max.iter[3])
  for (iter in 1:control$max.iter[1]) {
    # loglik
    s_lik <-
      sum(log(M_mu_no1) + log(M_h0_no1) - M_H_no1) -
      sum(M_H_nr1) +
      sum(log(1 - M_S_nl1)) +
      sum(log(M_S1mS2_ni1)) -
      s_lambda * thetaRtheta

    # --- inner loop: alternating beta/theta updates ---

    for (k in 1:K) {
      ## update betas
      s_omega <- 1
      M_beta_p1_OLD <- M_beta_p1
      s_lik_OLD <- s_lik
      ## OLD: using pre-computed transpose M_tX
      # M_gradbeta_p1     =
      #   (M_tX_nop%*%(1-M_H_no1))-
      #   M_tX_nrp%*%M_H_nr1+
      #   M_tX_nlp%*%(M_S_nl1*M_H_nl1/(1-M_S_nl1))+
      #   M_tX_nip%*%((M_H2_ni1*M_S2_ni1-M_H1_ni1*M_S1_ni1)/M_S1mS2_ni1)
      ## NEW: using crossprod()
      M_gradbeta_p1 <-
        crossprod(M_X_nop, (1 - M_H_no1)) -
        crossprod(M_X_nrp, M_H_nr1) +
        crossprod(M_X_nlp, (M_S_nl1 * M_H_nl1 / (1 - M_S_nl1))) +
        crossprod(M_X_nip, ((M_H2_ni1 * M_S2_ni1 - M_H1_ni1 * M_S1_ni1) / M_S1mS2_ni1))
      ## OLD: diag() creates full n x n matrices
      # M_hessbeta_p1 =
      #   M_tX_nop%*%diag(c(M_H_no1),n.ctype[2],n.ctype[2])%*%M_X_nop+
      #   M_tX_nrp%*%diag(c(M_H_nr1),n.ctype[1],n.ctype[1])%*%M_X_nrp+
      #   M_tX_nlp%*%diag(c((M_S_nl1/(1-M_S_nl1)^2*M_H_nl1^2-M_S_nl1/(1-M_S_nl1)*M_H_nl1)),n.ctype[3],n.ctype[3])%*%M_X_nlp+
      #   M_tX_nip%*%diag(c((M_S1_ni1*M_S2_ni1/M_S1mS2_ni1^2*(M_H2_ni1-M_H1_ni1)^2+
      #                        (M_S1_ni1*M_H1_ni1-M_S2_ni1*M_H2_ni1)/M_S1mS2_ni1
      #   )),n.ctype[4],n.ctype[4])%*%M_X_nip
      ## PREV: use element-wise multiplication with pre-computed M_tX
      # M_hessbeta_p1 =
      #   M_tX_nop %*% (c(M_H_no1) * M_X_nop)+
      #   M_tX_nrp %*% (c(M_H_nr1) * M_X_nrp)+
      #   M_tX_nlp %*% (c(M_S_nl1/(1-M_S_nl1)^2*M_H_nl1^2-M_S_nl1/(1-M_S_nl1)*M_H_nl1) * M_X_nlp)+
      #   M_tX_nip %*% (c(M_S1_ni1*M_S2_ni1/M_S1mS2_ni1^2*(M_H2_ni1-M_H1_ni1)^2+
      #                     (M_S1_ni1*M_H1_ni1-M_S2_ni1*M_H2_ni1)/M_S1mS2_ni1) * M_X_nip)
      ## NEW: use crossprod()
      M_hessbeta_p1 <-
        crossprod(M_X_nop, c(M_H_no1) * M_X_nop) +
        crossprod(M_X_nrp, c(M_H_nr1) * M_X_nrp) +
        crossprod(M_X_nlp, c(M_S_nl1 / (1 - M_S_nl1)^2 * M_H_nl1^2 - M_S_nl1 / (1 - M_S_nl1) * M_H_nl1) * M_X_nlp) +
        crossprod(M_X_nip, c(M_S1_ni1 * M_S2_ni1 / M_S1mS2_ni1^2 * (M_H2_ni1 - M_H1_ni1)^2 +
          (M_S1_ni1 * M_H1_ni1 - M_S2_ni1 * M_H2_ni1) / M_S1mS2_ni1) * M_X_nip)
      # avoid division by 0 (leading to the issue spotted by Kenneth Beath [email of 20220112])
      if (p == 1) {
        if (M_hessbeta_p1[1, 1] == 0) {
          M_hessbeta_p1[1, 1] <- control$epsilon[1]
        }
      }

  #     # !Test numerical issues
  # if (!is.finite(s_lik) || s_lambda < 0 || !all(is.finite(M_hessbeta_p1)) ||
  #     M_hessbeta_p1[1, 1] <= 0) browser()


      M_stepbeta_p1 <- chol2inv(chol(M_hessbeta_p1)) %*% M_gradbeta_p1
      M_beta_p1 <- M_beta_p1_OLD + s_omega * M_stepbeta_p1
      M_mu_no1 <- exp(M_X_nop %*% M_beta_p1)
      M_mu_nr1 <- exp(M_X_nrp %*% M_beta_p1)
      M_mu_nl1 <- exp(M_X_nlp %*% M_beta_p1)
      M_mu_ni1 <- exp(M_X_nip %*% M_beta_p1)
      M_H_no1 <- M_H0_no1 * M_mu_no1
      M_H_nr1 <- M_H0_nr1 * M_mu_nr1
      M_H_nl1 <- M_H0_nl1 * M_mu_nl1
      M_H1_ni1 <- M_H01_ni1 * M_mu_ni1
      M_H2_ni1 <- M_H02_ni1 * M_mu_ni1
      M_S_no1 <- exp(-M_H_no1)
      M_S_nl1 <- exp(-M_H_nl1)
      M_S1_ni1 <- exp(-M_H1_ni1)
      M_S2_ni1 <- exp(-M_H2_ni1)
      # avoid division by 0
      M_S_nl1[M_S_nl1 == 1] <- 1 - control$epsilon[1]
      M_S1_ni1[M_S1_ni1 == 1] <- 1 - control$epsilon[1]
      M_S2_ni1[M_S2_ni1 == 1] <- 1 - control$epsilon[1]
      M_S1mS2_ni1 <- M_S1_ni1 - M_S2_ni1
      M_S1mS2_ni1[M_S1mS2_ni1 < control$epsilon[2]] <- control$epsilon[2]
      # loglik
      s_lik <-
        sum(log(M_mu_no1) + log(M_h0_no1) - M_H_no1) -
        sum(M_H_nr1) +
        sum(log(1 - M_S_nl1)) +
        sum(log(M_S1mS2_ni1)) -
        s_lambda * thetaRtheta
      ## if likelihood decreases
      if (s_lik < s_lik_OLD) {
        i <- 0
        s_omega <- 1 / s_kappa
        while (s_lik < s_lik_OLD) {
          M_beta_p1 <- M_beta_p1_OLD + s_omega * M_stepbeta_p1
          M_mu_no1 <- exp(M_X_nop %*% M_beta_p1)
          M_mu_nr1 <- exp(M_X_nrp %*% M_beta_p1)
          M_mu_nl1 <- exp(M_X_nlp %*% M_beta_p1)
          M_mu_ni1 <- exp(M_X_nip %*% M_beta_p1)
          M_H_no1 <- M_H0_no1 * M_mu_no1
          M_H_nr1 <- M_H0_nr1 * M_mu_nr1
          M_H_nl1 <- M_H0_nl1 * M_mu_nl1
          M_H1_ni1 <- M_H01_ni1 * M_mu_ni1
          M_H2_ni1 <- M_H02_ni1 * M_mu_ni1
          M_S_no1 <- exp(-M_H_no1)
          M_S_nl1 <- exp(-M_H_nl1)
          M_S1_ni1 <- exp(-M_H1_ni1)
          M_S2_ni1 <- exp(-M_H2_ni1)
          # avoid division by 0
          M_S_nl1[M_S_nl1 == 1] <- 1 - control$epsilon[1]
          M_S1_ni1[M_S1_ni1 == 1] <- 1 - control$epsilon[1]
          M_S2_ni1[M_S2_ni1 == 1] <- 1 - control$epsilon[1]
          M_S1mS2_ni1 <- M_S1_ni1 - M_S2_ni1
          M_S1mS2_ni1[M_S1mS2_ni1 < control$epsilon[2]] <- control$epsilon[2]
          # loglik
          s_lik <-
            sum(log(M_mu_no1) + log(M_h0_no1) - M_H_no1) -
            sum(M_H_nr1) +
            sum(log(1 - M_S_nl1)) +
            sum(log(M_S1mS2_ni1)) -
            s_lambda * thetaRtheta
          # update value of omega
          if (s_omega >= 1e-2) {
            s_omega <- s_omega / s_kappa
          } else {
            if (s_omega < 1e-2 & s_omega >= 1e-5) {
              s_omega <- s_omega * 5e-2
            } else {
              if (s_omega < 1e-5) {
                s_omega <- s_omega * 1e-5
              }
            }
          }
          i <- i + 1
          if (i > 500) {
            break
          }
        }
      }
      ## update thetas
      s_nu <- 1
      M_theta_m1_OLD <- M_theta_m1
      s_lik_OLD <- s_lik
      M_gradthetaA_m1 <-
        M_tpsi_nom %*% (1 / M_h0_no1) +
        M_tPsi_nlm %*% (M_S_nl1 * M_mu_nl1 / (1 - M_S_nl1)) +
        M_tPsi2_nim %*% (M_S2_ni1 * M_mu_ni1 / (M_S1mS2_ni1)) -
        TwoLRtheta * (TwoLRtheta < 0) + 0.3
      M_gradthetaB_m1 <-
        M_tPsi_nom %*% M_mu_no1 +
        M_tPsi_nrm %*% M_mu_nr1 +
        M_tPsi1_nim %*% (M_S1_ni1 * M_mu_ni1 / (M_S1mS2_ni1)) +
        TwoLRtheta * (TwoLRtheta > 0) + 0.3
      M_gradtheta_m1 <- M_gradthetaA_m1 - M_gradthetaB_m1
      M_s_m1 <- M_theta_m1 / M_gradthetaB_m1
      M_steptheta_p1 <- M_s_m1 * M_gradtheta_m1
      M_theta_m1 <- M_theta_m1_OLD + s_nu * M_steptheta_p1
      M_theta_m1[M_theta_m1 < control$epsilon[2]] <- control$epsilon[2]
      M_h0_no1 <- M_psi_nom %*% M_theta_m1
      M_H0_no1 <- M_Psi_nom %*% M_theta_m1
      M_H0_nr1 <- M_Psi_nrm %*% M_theta_m1
      M_H0_nl1 <- M_Psi_nlm %*% M_theta_m1
      M_H01_ni1 <- M_Psi1_nim %*% M_theta_m1
      M_H02_ni1 <- M_Psi2_nim %*% M_theta_m1
      Rtheta <- M_R_mm %*% M_theta_m1
      thetaRtheta <- t(M_theta_m1) %*% Rtheta
      TwoLRtheta <- s_lambda * 2 * Rtheta
      M_H_no1 <- M_H0_no1 * M_mu_no1
      M_H_nr1 <- M_H0_nr1 * M_mu_nr1
      M_H_nl1 <- M_H0_nl1 * M_mu_nl1
      M_H1_ni1 <- M_H01_ni1 * M_mu_ni1
      M_H2_ni1 <- M_H02_ni1 * M_mu_ni1
      M_S_no1 <- exp(-M_H_no1)
      M_S_nl1 <- exp(-M_H_nl1)
      M_S1_ni1 <- exp(-M_H1_ni1)
      M_S2_ni1 <- exp(-M_H2_ni1)
      # avoid division by 0
      M_S_nl1[M_S_nl1 == 1] <- 1 - control$epsilon[1]
      M_S1_ni1[M_S1_ni1 == 1] <- 1 - control$epsilon[1]
      M_S2_ni1[M_S2_ni1 == 1] <- 1 - control$epsilon[1]
      M_S1mS2_ni1 <- M_S1_ni1 - M_S2_ni1
      M_S1mS2_ni1[M_S1mS2_ni1 < control$epsilon[2]] <- control$epsilon[2]
      # loglik
      s_lik <-
        sum(log(M_mu_no1) + log(M_h0_no1) - M_H_no1) -
        sum(M_H_nr1) +
        sum(log(1 - M_S_nl1)) +
        sum(log(M_S1mS2_ni1)) -
        s_lambda * thetaRtheta
      ## if likelihood decreases
      if (s_lik < s_lik_OLD) {
        i <- 0
        s_omega <- 1 / s_kappa
        while (s_lik < s_lik_OLD) {
          M_theta_m1 <- M_theta_m1_OLD + s_omega * M_steptheta_p1
          M_theta_m1[M_theta_m1 < control$epsilon[2]] <- control$epsilon[2]
          M_h0_no1 <- M_psi_nom %*% M_theta_m1
          M_H0_no1 <- M_Psi_nom %*% M_theta_m1
          M_H0_nr1 <- M_Psi_nrm %*% M_theta_m1
          M_H0_nl1 <- M_Psi_nlm %*% M_theta_m1
          M_H01_ni1 <- M_Psi1_nim %*% M_theta_m1
          M_H02_ni1 <- M_Psi2_nim %*% M_theta_m1
          Rtheta <- M_R_mm %*% M_theta_m1
          thetaRtheta <- t(M_theta_m1) %*% Rtheta
          TwoLRtheta <- s_lambda * 2 * Rtheta
          M_H_no1 <- M_H0_no1 * M_mu_no1
          M_H_nr1 <- M_H0_nr1 * M_mu_nr1
          M_H_nl1 <- M_H0_nl1 * M_mu_nl1
          M_H1_ni1 <- M_H01_ni1 * M_mu_ni1
          M_H2_ni1 <- M_H02_ni1 * M_mu_ni1
          M_S_no1 <- exp(-M_H_no1)
          M_S_nl1 <- exp(-M_H_nl1)
          M_S1_ni1 <- exp(-M_H1_ni1)
          M_S2_ni1 <- exp(-M_H2_ni1)
          # avoid division by 0
          M_S_nl1[M_S_nl1 == 1] <- 1 - control$epsilon[1]
          M_S1_ni1[M_S1_ni1 == 1] <- 1 - control$epsilon[1]
          M_S2_ni1[M_S2_ni1 == 1] <- 1 - control$epsilon[1]
          M_S1mS2_ni1 <- M_S1_ni1 - M_S2_ni1
          M_S1mS2_ni1[M_S1mS2_ni1 < control$epsilon[2]] <- control$epsilon[2]
          # loglik
          s_lik <-
            sum(log(M_mu_no1) + log(M_h0_no1) - M_H_no1) -
            sum(M_H_nr1) +
            sum(log(1 - M_S_nl1)) +
            sum(log(M_S1mS2_ni1)) -
            s_lambda * thetaRtheta
          # update omega
          if (s_omega >= 1e-2) {
            s_omega <- s_omega / s_kappa
          } else {
            if (s_omega < 1e-2 & s_omega >= 1e-5) {
              s_omega <- s_omega * 5e-2
            } else {
              if (s_omega < 1e-5) {
                s_omega <- s_omega * 1e-5
              }
            }
          }
          i <- i + 1
          if (i > 500) {
            break
          }
        }
      }
      if (all(c(abs(M_beta_p1 - M_beta_p1_OLD), abs(M_theta_m1 - M_theta_m1_OLD)) < s_convlimit)) {
        break
      }
      if (control$max.iter[1] == 1) {
        if (any(k == seq(0, control$max.iter[3], control$max.iter[2]))) {
          control$epsilon[2] <- control$epsilon[2] * 10
        }
      }
    }

    # --- Hessian and smoothing parameter update (df-based) ---
    # H matrix
    H <- HRinv <- matrix(0, p + m, p + m)
    ## OLD: diag() creates full n x n matrices
    # H[1:p,1:p] =
    #   M_tX_nop%*%diag(c(M_H_no1),n.ctype[2],n.ctype[2])%*%M_X_nop+
    #   M_tX_nrp%*%diag(c(M_H_nr1),n.ctype[1],n.ctype[1])%*%M_X_nrp+
    #   M_tX_nlp%*%diag(c((M_S_nl1/(1-M_S_nl1)^2*M_H_nl1^2-M_S_nl1/(1-M_S_nl1)*M_H_nl1)),n.ctype[3],n.ctype[3])%*%M_X_nlp+
    #   M_tX_nip%*%diag(c((M_S1_ni1*M_S2_ni1/M_S1mS2_ni1^2*(M_H2_ni1-M_H1_ni1)^2+
    #                        (M_S1_ni1*M_H1_ni1-M_S2_ni1*M_H2_ni1)/M_S1mS2_ni1
    #   )),n.ctype[4],n.ctype[4])%*%M_X_nip
    ## PREV: use element-wise multiplication with pre-computed M_tX
    # H[1:p,1:p] =
    #   M_tX_nop %*% (c(M_H_no1) * M_X_nop)+
    #   M_tX_nrp %*% (c(M_H_nr1) * M_X_nrp)+
    #   M_tX_nlp %*% (c(M_S_nl1/(1-M_S_nl1)^2*M_H_nl1^2-M_S_nl1/(1-M_S_nl1)*M_H_nl1) * M_X_nlp)+
    #   M_tX_nip %*% (c(M_S1_ni1*M_S2_ni1/M_S1mS2_ni1^2*(M_H2_ni1-M_H1_ni1)^2+
    #                     (M_S1_ni1*M_H1_ni1-M_S2_ni1*M_H2_ni1)/M_S1mS2_ni1) * M_X_nip)
    ## NEW: use crossprod()
    H[1:p, 1:p] <-
      crossprod(M_X_nop, c(M_H_no1) * M_X_nop) +
      crossprod(M_X_nrp, c(M_H_nr1) * M_X_nrp) +
      crossprod(M_X_nlp, c(M_S_nl1 / (1 - M_S_nl1)^2 * M_H_nl1^2 - M_S_nl1 / (1 - M_S_nl1) * M_H_nl1) * M_X_nlp) +
      crossprod(M_X_nip, c(M_S1_ni1 * M_S2_ni1 / M_S1mS2_ni1^2 * (M_H2_ni1 - M_H1_ni1)^2 +
        (M_S1_ni1 * M_H1_ni1 - M_S2_ni1 * M_H2_ni1) / M_S1mS2_ni1) * M_X_nip)
    ## OLD: diag() creates full n x n matrices
    # H[1:p,(p+1):(p+m)] =
    #   M_tX_nop%*%diag(c(M_mu_no1),n.ctype[2],n.ctype[2])%*%M_Psi_nom+
    #   M_tX_nrp%*%diag(c(M_mu_nr1),n.ctype[1],n.ctype[1])%*%M_Psi_nrm+
    #   M_tX_nlp%*%diag(c((M_S_nl1/(1-M_S_nl1)^2*M_H_nl1-M_S_nl1/(1-M_S_nl1))*M_mu_nl1),n.ctype[3],n.ctype[3])%*%M_Psi_nlm+
    #   M_tX_nip%*%diag(c(M_S1_ni1*M_S2_ni1/M_S1mS2_ni1^2*(M_H2_ni1-M_H1_ni1)*M_mu_ni1),n.ctype[4],n.ctype[4])%*%(M_Psi2_nim-M_Psi1_nim)+
    #   M_tX_nip%*%diag(c(M_S1_ni1/M_S1mS2_ni1*M_mu_ni1),n.ctype[4],n.ctype[4])%*%M_Psi1_nim-
    #   M_tX_nip%*%diag(c(M_S2_ni1/M_S1mS2_ni1*M_mu_ni1),n.ctype[4],n.ctype[4])%*%M_Psi2_nim
    ## PREV: use element-wise multiplication with pre-computed M_tX
    # H[1:p,(p+1):(p+m)] =
    #   M_tX_nop %*% (c(M_mu_no1) * M_Psi_nom)+
    #   M_tX_nrp %*% (c(M_mu_nr1) * M_Psi_nrm)+
    #   M_tX_nlp %*% (c((M_S_nl1/(1-M_S_nl1)^2*M_H_nl1-M_S_nl1/(1-M_S_nl1))*M_mu_nl1) * M_Psi_nlm)+
    #   M_tX_nip %*% (c(M_S1_ni1*M_S2_ni1/M_S1mS2_ni1^2*(M_H2_ni1-M_H1_ni1)*M_mu_ni1) * (M_Psi2_nim-M_Psi1_nim))+
    #   M_tX_nip %*% (c(M_S1_ni1/M_S1mS2_ni1*M_mu_ni1) * M_Psi1_nim)-
    #   M_tX_nip %*% (c(M_S2_ni1/M_S1mS2_ni1*M_mu_ni1) * M_Psi2_nim)
    ## NEW: use crossprod()
    H[1:p, (p + 1):(p + m)] <-
      crossprod(M_X_nop, c(M_mu_no1) * M_Psi_nom) +
      crossprod(M_X_nrp, c(M_mu_nr1) * M_Psi_nrm) +
      crossprod(M_X_nlp, c((M_S_nl1 / (1 - M_S_nl1)^2 * M_H_nl1 - M_S_nl1 / (1 - M_S_nl1)) * M_mu_nl1) * M_Psi_nlm) +
      crossprod(M_X_nip, c(M_S1_ni1 * M_S2_ni1 / M_S1mS2_ni1^2 * (M_H2_ni1 - M_H1_ni1) * M_mu_ni1) * (M_Psi2_nim - M_Psi1_nim)) +
      crossprod(M_X_nip, c(M_S1_ni1 / M_S1mS2_ni1 * M_mu_ni1) * M_Psi1_nim) -
      crossprod(M_X_nip, c(M_S2_ni1 / M_S1mS2_ni1 * M_mu_ni1) * M_Psi2_nim)
    H[(p + 1):(p + m), 1:p] <- t(H[1:p, (p + 1):(p + m)])
    ## OLD: diag() creates full n x n matrices
    # H[(p+1):(p+m),(p+1):(p+m)] =
    #   M_tpsi_nom%*%diag(c(1/M_h0_no1^2),n.ctype[2],n.ctype[2])%*%M_psi_nom+
    #   M_tPsi_nlm%*%diag(c(M_S_nl1/(1-M_S_nl1)^2*M_mu_nl1^2),n.ctype[3],n.ctype[3])%*%M_Psi_nlm+
    #   (M_tPsi2_nim-M_tPsi1_nim)%*%diag(c(M_S1_ni1*M_S2_ni1/M_S1mS2_ni1^2*M_mu_ni1^2),n.ctype[4],n.ctype[4])%*%(M_Psi2_nim-M_Psi1_nim)
    ## NEW: use element-wise multiplication
    H[(p + 1):(p + m), (p + 1):(p + m)] <-
      M_tpsi_nom %*% (c(1 / M_h0_no1^2) * M_psi_nom) +
      M_tPsi_nlm %*% (c(M_S_nl1 / (1 - M_S_nl1)^2 * M_mu_nl1^2) * M_Psi_nlm) +
      (M_tPsi2_nim - M_tPsi1_nim) %*% (c(M_S1_ni1 * M_S2_ni1 / M_S1mS2_ni1^2 * M_mu_ni1^2) * (M_Psi2_nim - M_Psi1_nim))
    s_lambda_old <- s_lambda
    s_df_old <- s_df
    s_sigma2_old <- 1 / (2 * s_lambda_old)
    # pos            = c(if(noX){FALSE}else{rep(TRUE,p)},M_theta_m1>control$min.theta)
    pos <- c(if (noX) {
      FALSE
    } else {
      rep(TRUE, p)
    }, (M_theta_m1 > control$min.theta & apply(H[(p + 1):(p + m), (p + 1):(p + m)], 2, sum) > 0))
    # MM: (G+Q)^-1
    temp <- try(chol2inv(chol(H[pos, pos] + (1 / s_sigma2_old) * M_Rstar_ll[pos, pos])), silent = T)
    if (class(temp)[1] != "try-error" & !any(is.infinite(temp))) {
      HRinv[pos, pos] <- temp
      ## MM: is ginv(H[pos,pos]) right? Forgot Q?
    } else {
      HRinv[pos, pos] <- MASS::ginv(H[pos, pos])
    }
    s_df <- m - sum(diag(HRinv %*% M_Rstar_ll)) / s_sigma2_old
    s_sigma2 <- c(t(M_theta_m1) %*% M_R_mm %*% M_theta_m1 / s_df)
    s_lambda <- 1 / (2 * s_sigma2)
    TwoLRtheta <- s_lambda * 2 * Rtheta
    full.iter <- full.iter + k
    if ((full.iter / iter) > (control$max.iter[2] * .975)) {
      control$epsilon[2] <- control$epsilon[2] * 10
    }
    if (full.iter > control$max.iter[3]) {
      break
    }
    if ((k < control$max.iter[2]) &
      (abs(s_df - s_df_old) < (control$tol * 10))
    ) {
      break
    }
  }
  s_lambda <- control$smooth <- s_lambda_old
  s_correction <- c(exp(-mean_j %*% M_beta_p1))
  M_thetatilde_m1 <- M_theta_m1
  M_theta_m1 <- M_theta_m1 * s_correction

  ###
  ### Inference
  ###
  # M_corr_ll = cbind(rbind(diag(rep(1,p)),s_correction*matrix(rep(M_thetatilde_m1,p)*rep(-mean_j,each=m),ncol=p)),
  # rbind(matrix(0,ncol=m,nrow=p),diag(rep(s_correction,m))))
  M_corr_ll <- cbind(
    rbind(diag(rep(1, p)), s_correction * matrix(rep(M_thetatilde_m1, p) * rep(-M_beta_p1, each = m), ncol = p)),
    rbind(matrix(0, ncol = m, nrow = p), diag(rep(s_correction, m)))
  )
  M_corr_ll[!pos, ] <- 0


  M_2 <- H + 2 * s_lambda * M_Rstar_ll
  Q <- matrix(NA, n, p + m)
  if (ctypeTF[1]) {
    Q[ctype[, 1], 1:p] <- rep(-M_H_nr1, p) * M_X_nrp
  }
  if (ctypeTF[2]) {
    Q[ctype[, 2], 1:p] <- rep((1 - M_H_no1), p) * M_X_nop
  }
  if (ctypeTF[3]) {
    Q[ctype[, 3], 1:p] <- rep(M_S_nl1 * M_H_nl1 / (1 - M_S_nl1), p) * M_X_nlp
  }
  if (ctypeTF[4]) {
    Q[ctype[, 4], 1:p] <- rep((M_H2_ni1 * M_S2_ni1 - M_H1_ni1 * M_S1_ni1) / M_S1mS2_ni1, p) * M_X_nip
  }
  if (ctypeTF[1]) {
    Q[ctype[, 1], -c(1:p)] <- rep(-M_mu_nr1, m) * M_Psi_nrm
  }
  if (ctypeTF[2]) {
    Q[ctype[, 2], -c(1:p)] <- rep(1 / M_h0_no1, m) * M_psi_nom - rep(M_mu_no1, m) * M_Psi_nom
  }
  if (ctypeTF[3]) {
    Q[ctype[, 3], -c(1:p)] <- rep(M_S_nl1 * M_mu_nl1 / (1 - M_S_nl1), m) * M_Psi_nlm
  }
  if (ctypeTF[4]) {
    Q[ctype[, 4], -c(1:p)] <- rep(M_S2_ni1 * M_mu_ni1 / (M_S1mS2_ni1), m) * M_Psi2_nim -
      rep(M_S1_ni1 * M_mu_ni1 / (M_S1mS2_ni1), m) * M_Psi1_nim
  }
  Sp <- Q - matrix(rep(c(rep(0, p), TwoLRtheta), n), n, byrow = T) / n
  Q <- t(Sp) %*% Sp
  # pos   = c(if(noX){FALSE}else{rep(TRUE,p)},M_theta_m1>control$min.theta)
  pos <- c(if (noX) {
    FALSE
  } else {
    rep(TRUE, p)
  }, (M_theta_m1 > control$min.theta & apply(H[(p + 1):(p + m), (p + 1):(p + m)], 2, sum) > 0))
  Minv_1 <- Minv_2 <- Hinv <- matrix(0, p + m, p + m)
  temp <- try(chol2inv(chol(M_2[pos, pos])), silent = T)
  if (class(temp)[1] != "try-error") {
    Minv_2[pos, pos] <- temp
    cov_NuNu_M2QM2 <- M_corr_ll %*% (Minv_2 %*% Q %*% Minv_2) %*% t(M_corr_ll)
    cov_NuNu_M2HM2 <- M_corr_ll %*% (Minv_2 %*% H %*% Minv_2) %*% t(M_corr_ll)
    se.Eta_M2QM2 <- sqrt(diag(cov_NuNu_M2QM2))
    se.Eta_M2HM2 <- sqrt(diag(cov_NuNu_M2HM2))
  } else {
    cov_NuNu_M2QM2 <- cov_NuNu_M2HM2 <- matrix(NA, p + m, p + m)
    se.Eta_M2QM2 <- se.Eta_M2HM2 <- rep(NA, p + m)
  }
  temp <- try(chol2inv(chol(H[pos, pos])), silent = T)
  if (class(temp)[1] != "try-error") {
    Hinv[pos, pos] <- temp
    cov_NuNu_H <- M_corr_ll %*% Hinv %*% t(M_corr_ll)
    se.Eta_H <- sqrt(diag(cov_NuNu_H))
  } else {
    cov_NuNu_H <- matrix(NA, p + m, p + m)
    se.Eta_H <- rep(NA, p + m)
  }
  mx.seNu.l5 <- as.data.frame(cbind(se.Eta_M2QM2, se.Eta_M2HM2, se.Eta_H))
  colnames(mx.seNu.l5) <- c("M2QM2", "M2HM2", "H")
  rownames(mx.seNu.l5)[(p + 1):(p + m)] <- paste("Theta", 1:m, sep = "")
  rownames(mx.seNu.l5)[(1:p)] <- paste("Beta", 1:p, sep = "")
  fit <- list(
    coef = list(Beta = c(M_beta_p1), Theta = c(M_theta_m1) * (M_theta_m1 > control$min.theta)),
    iter = c(iter, full.iter)
  )
  fit$se <- list(Beta = mx.seNu.l5[1:p, ], Theta = mx.seNu.l5[(p + 1):(p + m), ])
  fit$covar <- list(M2QM2 = cov_NuNu_M2QM2, M2HM2 = cov_NuNu_M2HM2, H = cov_NuNu_H)
  fit$knots <- knots
  fit$control <- control
  fit$call <- match.call()
  fit$dim <- list(n = n, n.obs = sum(observed), n.ties = sum(ties), p = p, m = knots$m)
  fit$data <- list(time = y, censoring = y[, 3L], X = X, name = data.name, entry = entry) # list(name = data.name)#
  fit$df <- s_df
  fit$ploglik <- s_lik
  fit$loglik <- s_lik + s_lambda * thetaRtheta
  class(fit) <- "coxph_mpl"
  fit
}
