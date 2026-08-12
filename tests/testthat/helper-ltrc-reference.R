# Reference LTRC implementation used only to verify coxph_mpl(entry=...)
# against an independent piecewise-hazard profile-likelihood estimator
# (see test-ltrc.R). Copied from dev/experiments/old_sources/LTRC_codes/
# so the comparison test doesn't depend on files outside this branch.

datagen_LT_RC <- function(n, beta, dist, c_max, alpha, psi) {
  nsample <- 0
  nfailed <- 0

  dataMat <- NULL

  while (nsample < n) {
    # Generate covariates:
    x1 <- rbinom(1, 1, 0.5)
    x2 <- rnorm(1, mean = 0, sd = 1)

    x_i <- cbind(x1, x2)

    # Generate standard uniform r.v. for event times:
    U <- runif(1, 0, 1)

    # Generate event times:
    if (dist == "weibull") {
      t_i <- (-log(U) / (psi * exp(x_i %*% beta)))^(1 / alpha)
    } else if (dist == "log-logistic") {
      t_i <- ((exp(-log(U) / exp(x_i %*% beta)) - 1) / psi)^(1 / alpha)
    }

    # Generate left truncation times:
    l_i <- runif(1, 0, 1)

    # Generate right censoring times:
    c_i <- runif(1, 0, c_max)

    # Process samples:
    if (l_i < t_i & t_i < c_i) { # Event time is observed.
      y_i <- t_i
      del_i <- 1

      dataMat <- rbind(dataMat, c(y_i, x_i, l_i, c_i, del_i))
      nsample <- nsample + 1
    } else if (l_i < c_i & c_i < t_i) { # Observation is right-censored.
      y_i <- c_i
      del_i <- 0

      dataMat <- rbind(dataMat, c(y_i, x_i, l_i, c_i, del_i))
      nsample <- nsample + 1
    } else {
      nfailed <- nfailed + 1
    }
  }

  colnames(dataMat) <- c("y_i", "x_i1", "x_i2", "l_i", "c_i", "censor_type")

  censor_prop <- (n - sum(dataMat[, 6])) / n

  list(
    Y = dataMat[, 1], X = dataMat[, 2:3], L = dataMat[, 4], C = dataMat[, 5],
    del = dataMat[, 6],
    censor_prop = censor_prop
  )
}

Cox_LT_RC <- function(Y, X, L, delVec, maxIter, tol, bin_num, bin_option, bin_edges) {
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  L <- as.matrix(L)

  n <- nrow(X)
  p <- ncol(X)

  ### Set starting values for beta: ###
  bVal <- matrix(0, p, 1)

  profLogLikVec <- NULL

  ### Set number of bins if unspecified: ###
  if (missing(bin_num)) bin_num <- floor(sqrt(sum(delVec)))

  num_boundary <- bin_num + 1

  ### Set option for choosing bins if unspecified: ###
  if (missing(bin_option)) {
    bin_option <- "equal-space"
  }

  ### Calculate bin_edges if unspecified: ###
  if (missing(bin_edges)) {
    if (bin_option == "quantile") {
      K <- bin_num
      bin_edges <- c(
        min(L),
        quantile(unique(Y[delVec == 1]), probs = seq(1 / K, 1 - 1 / K, length.out = K - 1)),
        max(Y)
      )
    } else {
      print("Method to place bins is required.")
    }
  }

  upper_b <- bin_edges[2:num_boundary]
  lower_b <- bin_edges[1:(num_boundary - 1)]
  bin_width <- diff(bin_edges)

  h <- hist(Y[delVec == 1], breaks = bin_edges, right = FALSE, plot = FALSE)
  binCount <- matrix(h$counts, nrow = length(h$counts), ncol = 1)

  binIndOne <- t(sapply(1:n, function(i) {
    lower_b <= L[i] & L[i] < Y[i] & Y[i] <= upper_b
  }))
  binIndTwo <- t(sapply(1:n, function(i) {
    lower_b <= L[i] & L[i] < upper_b & upper_b < Y[i]
  }))
  binIndThree <- t(sapply(1:n, function(i) {
    L[i] < lower_b & Y[i] > upper_b
  }))
  binIndFour <- t(sapply(1:n, function(i) {
    L[i] < lower_b & lower_b < Y[i] & Y[i] <= upper_b
  }))

  binLeft <- t(sapply(1:n, function(i) {
    L[i] >= lower_b & L[i] < upper_b
  }))
  binLeftNum <- apply(t(t(binLeft) * 1:bin_num), 1, sum)

  binObs <- t(sapply(1:n, function(i) {
    Y[i] >= lower_b & Y[i] < upper_b
  }))
  binObsNum <- apply(t(t(binObs) * 1:bin_num), 1, sum)

  ############### Check if everything corresponds: ###############
  a <- binObsNum - binLeftNum + 1
  b <- (binIndOne + binIndTwo + binIndThree + binIndFour)
  c <- apply(b, 1, sum)

  if (all(which(a != c) == which(binLeftNum == 0)) != TRUE) stop("bin index mismatch")
  ################################################################

  inputVal <- (as.vector(Y - L) * binIndOne +
    (t(t(binIndTwo) * upper_b) - as.vector(L) * binIndTwo) +
    t(t(binIndThree) * bin_width) +
    (as.vector(Y) * binIndFour - t(t(binIndFour) * lower_b)))

  Xbeta <- X %*% bVal
  A_mat <- as.vector(exp(Xbeta)) * (inputVal)
  sumA_mat <- apply(A_mat, 2, sum)
  W <- t(t(A_mat) / sumA_mat)

  for (k in 1:maxIter) {
    beta_old <- bVal
    beta_new <- NULL

    ############################### Newton Step: ###############################

    gradient <- t(X) %*% (delVec - W %*% binCount)

    Xbar_mat <- t(X) %*% W

    hessPartOne <- Xbar_mat %*% diag(as.vector(binCount)) %*% t(Xbar_mat)
    hessPartTwo <- t(X) %*% diag(as.vector(W %*% binCount)) %*% X

    profHess <- hessPartOne - hessPartTwo

    beta_new <- beta_old - solve(profHess) %*% gradient

    ############################################################################
    bVal <- beta_new

    Xbeta <- X %*% bVal
    A_mat <- as.vector(exp(Xbeta)) * (inputVal)
    sumA_mat <- apply(A_mat, 2, sum)

    sumA_mat_profLik <- sumA_mat

    W <- t(t(A_mat) / sumA_mat)

    proflik <- delVec %*% (Xbeta) - log(sumA_mat_profLik) %*% binCount
    profLogLikVec <- c(profLogLikVec, proflik)

    if (all(abs(beta_old - beta_new) < tol)) {
      cvg <- paste("Simulation converged with", k, "iterations.")
      break
    }
  }

  Xbeta <- X %*% bVal
  A_mat <- as.vector(exp(Xbeta)) * (inputVal)
  sumA_mat <- apply(A_mat, 2, sum)
  W <- t(t(A_mat) / sumA_mat)

  theta <- as.numeric(binCount) / sumA_mat

  ### To find the Hessian of Beta using Profile Likelihood: ###

  Xbar_mat <- t(X) %*% W
  profHessFinal <- Xbar_mat %*% diag(as.vector(binCount)) %*% t(Xbar_mat) - t(X) %*% diag(as.vector(W %*% binCount)) %*% X
  betaVar <- solve(-profHessFinal)

  ##### ASYMPTOTIC CALCULATIONS #####

  F_11 <- t(X) %*% diag(as.vector(inputVal %*% as.vector(theta)) * as.vector(exp(Xbeta))) %*% X

  F_12 <- t(X) %*% diag(as.vector(exp(Xbeta))) %*% inputVal

  F_22 <- diag(as.vector(binCount) / (as.vector(theta)^2))

  blockinv1 <- function(A, B, C, D) {
    topLeft <- solve(A - B %*% solve(D) %*% C)

    topRight <- -topLeft %*% B %*% solve(D)

    bottomLeft <- -solve(D) %*% C %*% topLeft

    bottomRight <- solve(D) - bottomLeft %*% B %*% solve(D)

    rbind(cbind(topLeft, topRight), cbind(bottomLeft, bottomRight))
  }

  asyVarMat <- blockinv1(F_11, F_12, t(F_12), F_22)

  ###################################

  list(
    bVal = bVal, tVal = theta, bin_edges = bin_edges,
    cvg = cvg, grad = gradient, profLogLik = profLogLikVec,
    betaVar = betaVar, asyVarMat = asyVarMat
  )
}
