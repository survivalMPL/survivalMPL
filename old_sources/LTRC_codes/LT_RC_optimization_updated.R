Cox_LT_RC <- function(Y, X, L, delVec, maxIter, tol, bin_num, bin_option, bin_edges){
  
  Y = as.matrix(Y)
  X = as.matrix(X)
  L = as.matrix(L)
  
  n = nrow(X)
  p = ncol(X)
  
  ### Set starting values for beta: ###
  bVal = matrix(0, p, 1)
  
  profLogLikVec = {}
  
  ### Set number of bins if unspecified: ###
  if(missing(bin_num)){bin_num = floor(sqrt(sum(delVec)))}
  
  num_boundary = bin_num + 1
  
  ### Set option for choosing bins if unspecified: ###
  if(missing(bin_option)){
    bin_option = "equal-space"
  }
  
  ### Calculate bin_edges if unspecified: ###
  if(missing(bin_edges)){
    if(bin_option == "quantile"){
      K = bin_num
      bin_edges <- c(min(L), 
                     quantile(unique(Y[delVec == 1]), probs = seq(1/K, 1 - 1/K, length.out = K - 1)), 
                     max(Y))
    }else{
      print("Method to place bins is required.")
    }
  }

  upper_b = bin_edges[2:num_boundary]
  lower_b = bin_edges[1:(num_boundary - 1)]
  bin_width = diff(bin_edges)
  
  h = hist(Y[delVec == 1], breaks = bin_edges, right = FALSE, plot = FALSE)
  binCount = matrix(h$counts, nrow = length(h$counts), ncol = 1)
  
  binIndOne = t(sapply(1:n, function(i){lower_b <= L[i] & L[i] < Y[i] & Y[i] <= upper_b}))
  binIndTwo = t(sapply(1:n, function(i){lower_b <= L[i] & L[i] < upper_b & upper_b < Y[i]}))
  binIndThree = t(sapply(1:n, function(i){L[i] < lower_b & Y[i] > upper_b}))
  binIndFour = t(sapply(1:n, function(i){L[i] < lower_b & lower_b < Y[i] & Y[i] <= upper_b}))
  
  binLeft = t(sapply(1:n, function(i){L[i] >= lower_b & L[i] < upper_b}))
  binLeftNum = apply(t(t(binLeft)*1:bin_num), 1, sum)
  
  binObs = t(sapply(1:n, function(i){Y[i] >= lower_b & Y[i] < upper_b}))
  binObsNum = apply(t(t(binObs)*1:bin_num), 1, sum)
  
  ############### Check if everything corresponds: ###############
  a = binObsNum - binLeftNum + 1
  b = (binIndOne + binIndTwo + binIndThree + binIndFour)
  c = apply(b, 1, sum)
  
  if(all(which(a != c) == which(binLeftNum == 0)) != TRUE){break}
  ################################################################
  
  inputVal = (as.vector(Y-L)*binIndOne + 
                (t(t(binIndTwo)*upper_b) - as.vector(L)*binIndTwo) +
                t(t(binIndThree)*bin_width)+
                (as.vector(Y)*binIndFour - t(t(binIndFour)*lower_b)))
  
  Xbeta = X%*%bVal 
  A_mat = as.vector(exp(Xbeta))*(inputVal)
  sumA_mat = apply(A_mat, 2, sum)
  #sumA_mat[sumA_mat == 0] = 1e-10
  W = t(t(A_mat)/sumA_mat)
  
  for(k in 1:maxIter){
    
    beta_old = bVal
    beta_new = {}
    
    ############################### Newton Step: ###############################
    
    gradient = t(X)%*%(delVec - W%*%binCount)
    
    Xbar_mat = t(X)%*%W
    
    hessPartOne = Xbar_mat%*%diag(as.vector(binCount))%*%t(Xbar_mat)
    hessPartTwo = t(X)%*%diag(as.vector(W%*%binCount))%*%X
    
    profHess = hessPartOne - hessPartTwo
  
    beta_new = beta_old - solve(profHess)%*%gradient
    
    ############################################################################
    bVal = beta_new
    
    Xbeta = X%*%bVal
    A_mat = as.vector(exp(Xbeta))*(inputVal)
    sumA_mat = apply(A_mat, 2, sum)
    
    sumA_mat_profLik = sumA_mat
    #sumA_mat_profLik[sumA_mat_profLik == 0] = 1 ## Check!
    
    #sumA_mat[sumA_mat == 0] = 1e-10
    W = t(t(A_mat)/sumA_mat)
    
    proflik = delVec%*%(Xbeta) - log(sumA_mat_profLik)%*%binCount
    profLogLikVec = c(profLogLikVec, proflik)
    
    all(abs(beta_old - beta_new) < tol)
    
    if(all(abs(beta_old - beta_new) < tol)){
      cvg = paste("Simulation converged with", k, "iterations.")
      break
    }
  }
  
  Xbeta = X%*%bVal
  A_mat = as.vector(exp(Xbeta))*(inputVal)
  sumA_mat = apply(A_mat, 2, sum)
  #sumA_mat[sumA_mat == 0] = 1e-10
  W = t(t(A_mat)/sumA_mat)
  
  theta = as.numeric(binCount)/sumA_mat 
  
  ### To find the Hessian of Beta using Profile Likelihood: ###
  
  Xbar_mat = t(X)%*%W
  profHessFinal = Xbar_mat%*%diag(as.vector(binCount))%*%t(Xbar_mat) - t(X)%*%diag(as.vector(W%*%binCount))%*%X
  betaVar = solve(-profHessFinal)
  
  ##### ASYMPTOTIC CALCULATIONS #####
  
  F_11 = t(X)%*%diag(as.vector(inputVal %*% as.vector(theta))*as.vector(exp(Xbeta)))%*%X
  
  F_12 = t(X)%*%diag(as.vector(exp(Xbeta)))%*%inputVal
  
  F_22 =  diag(as.vector(binCount)/(as.vector(theta)^2))
  
  blockinv1 <- function(A, B, C, D){
    
    topLeft = solve(A - B%*%solve(D)%*%C)
    
    topRight = -topLeft%*%B%*%solve(D)
    
    bottomLeft = -solve(D)%*%C%*%topLeft
    
    bottomRight = solve(D) - bottomLeft%*%B%*%solve(D)
    
    return(rbind(cbind(topLeft,topRight), cbind(bottomLeft, bottomRight)))
  }
  
  asyVarMat = blockinv1(F_11, F_12, t(F_12), F_22)
  
  ###################################

  return(list(bVal = bVal, tVal = theta,  bin_edges = bin_edges, 
              cvg = cvg, grad = gradient, profLogLik = profLogLikVec, 
              betaVar = betaVar, asyVarMat = asyVarMat))
}