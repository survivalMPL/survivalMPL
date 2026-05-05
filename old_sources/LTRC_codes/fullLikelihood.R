library(fastDummies)

####################################################
#### Function to calculate full log-likelihood: ####
####################################################

fullLogLikFun <- function(beta, theta, Y, X, L, delVec, bin_edges){
  
  Y = as.matrix(Y)
  X = as.matrix(X)
  L = as.matrix(L)
  
  n = length(Y)
  bin_num = length(bin_edges) - 1
  
  Xbeta = X%*%beta
  
  bin_width = diff(bin_edges)
  
  num_boundary = bin_num + 1
  upper_b = bin_edges[2:num_boundary]
  lower_b = bin_edges[1:(num_boundary - 1)]
  h = hist(Y[delVec == 1], breaks = bin_edges, right = TRUE, plot = F)
  binCount = matrix(h$counts, bin_num, 1)
  
  binIndOne = t(sapply(1:n, function(i){L[i] >= lower_b & L[i] < upper_b & Y[i] < upper_b}))
  binIndTwo = t(sapply(1:n, function(i){L[i] >= lower_b & L[i] < upper_b & Y[i] > upper_b}))
  binIndThree = t(sapply(1:n, function(i){L[i] < lower_b & Y[i] > upper_b}))
  binIndFour = t(sapply(1:n, function(i){Y[i] >= lower_b & Y[i] < upper_b & lower_b >= L[i]}))
  
  inputVal = (as.vector(Y-L)*binIndOne + 
                (t(t(binIndTwo)*upper_b) - as.vector(L)*binIndTwo) +
                t(t(binIndThree)*bin_width) + t(t(binIndOne)*bin_width) +
                (as.vector(Y)*binIndFour - t(t(binIndFour)*lower_b)))
  
  cumulHaz = apply(t(t(inputVal)*theta), 1, sum) ###CHECK!!!
  
  theta_alt = theta
  theta_alt[theta_alt == 0] = 1
  
  loglik = delVec%*%(Xbeta) + t(binCount)%*%log(theta_alt) - t(exp(Xbeta))%*%as.matrix(cumulHaz)
  
  return(loglik)
}

partialLogLikFun <- function(beta, Y, X, delVec){
  Y = as.matrix(Y)
  X = as.matrix(X)
  L = as.matrix(L)
  
  n = length(Y)
  
  partialLogLik = 0
  
  for(i in 1:n){
    addLogLik = delVec[i]*(X[i]*beta - log(t(exp(X%*%beta))%*%as.vector(Y>=Y[i])))
    partialLogLik = partialLogLik + addLogLik
  }
   
  return(partialLogLik)
}
                  