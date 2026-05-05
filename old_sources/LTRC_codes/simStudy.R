simStudy <- function(repeats, n, beta, nbinother, bin_option, 
                     c_max, dist, alpha, psi, maxIter, tol, bin_edge_option){
  
  ### Values to be stored: ###
  i = 1
  repID = {}
  
  overall_prop = {}
  
  allIter = {}
  
  bStore = {}
  tStore = {}
  
  # bStore_eha = {}
  # tStore_eha = {}
  
  bStore_cox = {}
  
  asyVarStore = list()
  asyVarStoreProfBeta = list()
  # asyVarStore_eha = list()
  asyVarStore_cox = list()
  
  allGrad = {}
  profLogLikStore = list()
  
  fullLogLikStore = {}
  # fullLogLikStore_eha = {}
  
  binEdgeStore = {}
  
  allTime = {}
  # allTime_eha = {}
  allTime_cox = {}
  
  maxObs = {}
  
  startTime = Sys.time()
  
  while(i <= repeats){
    
    set.seed(i)
    
    ######################
    ### Simulate data: ###
    ######################
    
    simData =  datagen_LT_RC(n = n, beta = beta, dist = dist, c_max = c_max, 
                             alpha = alpha, psi = psi)
    
    Y = simData$Y
    X = simData$X 
    L = simData$L
    delVec = simData$del
    
    ###############################
    ### Optimization algorithm: ###
    ###############################
    
    cat("Starting replication", i, "now: \n")
    
    
    
    if(missing(bin_edge_option)){
      
      binEdgeSet = FALSE
      
      opt_proflik_start = Sys.time()
      
      optimSet = Cox_LT_RC(Y, X, L, delVec, maxIter, tol, 
                           bin_num = nbinother, bin_option = bin_option)
      
      opt_proflik_end = Sys.time()
      
    }else{
      
      binEdgeSet = TRUE
      
      opt_proflik_start = Sys.time()
      
      optimSet = Cox_LT_RC(Y, X, L, delVec, maxIter, tol, 
                           bin_num = nbinother, bin_option = bin_option,
                           bin_edges = bin_edge_option)
      
      opt_proflik_end = Sys.time()
      
    }
    
    
    bin_edges = optimSet$bin_edges
    
    
    
    #######################################################
    ### Competing optimization algorithm (eha package): ###
    #######################################################
    
    opt_eha_start = Sys.time()
    
    # EHA Package:
    fit <- eha::pchreg(Surv(L, Y, delVec) ~ X[,1] + X[,2], data = simData, cuts = bin_edges)
    
    # pch Package:
    fit2 <- pch::pchreg(Surv(L, Y, delVec) ~ X[,1] + X[,2], data = simData, breaks = bin_edges)
      
    opt_eha_end = Sys.time()
    
    #######################################################
    ### Competing optimization algorithm (coxph): ###
    #######################################################
    
    opt_cox_start = Sys.time()
    
    coxfit <- coxph(Surv(L, Y, delVec) ~ X[,1] + X[,2], data = simData)
    
    opt_cox_end = Sys.time()
     
    #######################
    ### Record outputs: ###
    #######################
    
    repID = rbind(repID, i)
    
    overall_prop = rbind(overall_prop, simData$censor_prop)
    
    allIter = rbind(allIter, optimSet$cvg)
    allGrad = rbind(allGrad, as.vector(optimSet$grad))
    profLogLikStore = c(profLogLikStore, list(optimSet$profLogLik))
    binEdgeStore = rbind(binEdgeStore, bin_edges)
    
    bStore = rbind(bStore, as.vector(optimSet$bVal))
    tStore = rbind(tStore, as.vector(optimSet$tVal))
    fullLogLikStore = rbind(fullLogLikStore, fullLogLikFun(optimSet$bVal, optimSet$tVal, Y, X, L, delVec, bin_edges))
    asyVarStoreProfBeta = c(asyVarStoreProfBeta, list(optimSet$betaVar))
    asyVarStore = c(asyVarStore, list(optimSet$asyVarMat))
    allTime = c(allTime, opt_proflik_end - opt_proflik_start)
    
    # bStore_eha = rbind(bStore_eha, fit$coefficients)
    # tStore_eha = rbind(tStore_eha, as.vector(fit$hazards))
    # fullLogLikStore_eha = rbind(fullLogLikStore_eha, fullLogLikFun(fit$coefficients, as.vector(fit$hazards), Y, X, L, delVec, bin_edges))
    # asyVarStore_eha = c(asyVarStore_eha, list(fit$var))
    # allTime_eha = c(allTime_eha, opt_eha_end - opt_eha_start)
    
    bStore_cox = rbind(bStore_cox, coxfit$coef)
    asyVarStore_cox = c(asyVarStore_cox, list(coxfit$var))
    allTime_cox = c(allTime_cox, opt_cox_end - opt_cox_start)
    
    maxObs = c(maxObs, max(Y[delVec == 1]))

    cat("Replication", i, "is complete: \n")
    i = i + 1
  }
  
  row.names(repID) = NULL
  endTime = Sys.time()
  
  totalTime = endTime - startTime 
  print(totalTime)
  
  
  ###########################################
  ### Save outputs from simulation study: ###
  ###########################################
  
  save(n, dist, nbinother, repeats, repID, allIter, allGrad, 
       profLogLikStore, bStore, tStore, 
       fullLogLikStore, asyVarStoreProfBeta, asyVarStore, allTime, binEdgeStore,
       #bStore_eha, tStore_eha, fullLogLikStore_eha, asyVarStore_eha, allTime_eha,
       bStore_cox, asyVarStore_cox, allTime_cox, overall_prop, maxObs,
       file = paste0(dist,"_n",n,"_nbins_",nbinother,
                     "_cmax",c_max, "_",bin_option,".RData"))
}