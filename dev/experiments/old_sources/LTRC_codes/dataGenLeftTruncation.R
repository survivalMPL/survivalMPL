datagen_LT_RC <- function(n, beta, dist, c_max, alpha, psi){
  
  nsample = 0
  nfailed = 0
  
  dataMat = {}
  
  while(nsample < n){
    
    # Generate covariates:
    x1 = rbinom(1, 1, 0.5)
    x2 = rnorm(1, mean = 0 , sd = 1)
    
    x_i = cbind(x1, x2)
    
    # Generate standard uniform r.v. for event times:
    U = runif(1, 0, 1)
    
    # Generate event times:
    if(dist == "weibull"){
      t_i = (-log(U)/(psi*exp(x_i%*%beta)))^(1/alpha)
    }else if(dist == "log-logistic"){
      t_i = ((exp(-log(U)/exp(x_i%*%beta)) - 1)/psi)^(1/alpha)
    }
    
    # Generate left truncation times:
    l_i = runif(1, 0, 1)
    
    # Generate right censoring times:
    c_i = runif(1, 0, c_max)
    
    # Process samples:
    if(l_i < t_i & t_i < c_i){ # Event time is observed.
      y_i = t_i
      del_i = 1
      
      dataMat = rbind(dataMat, c(y_i, x_i, l_i, c_i, del_i))
      nsample = nsample + 1
    }
    else if(l_i <c_i & c_i < t_i){ # Observation is right-censored.
      y_i = c_i 
      del_i = 0
      
      dataMat = rbind(dataMat, c(y_i, x_i, l_i, c_i, del_i))
      nsample = nsample + 1
    }
    else{nfailed = nfailed + 1}
  }
  
  colnames(dataMat) = c("y_i", "x_i1", "x_i2", "l_i", "c_i", "censor_type" )
  
  censor_prop = (n - sum(dataMat[,6]))/n
  
  return(list(Y = dataMat[,1], X = dataMat[,2:3], L = dataMat[,4], C = dataMat[,5], 
              del = dataMat[,6],
              censor_prop = censor_prop))
  
}
