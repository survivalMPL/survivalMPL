library(eha)
library(survival)

###############################
### Import necessary files: ###
###############################

source("dataGenLeftTruncation.R")
source("LT_RC_optimization_updated.R")
source("fullLikelihood.R")
source("simStudy.R")

simStudy(repeats = 300, n = 1000, beta = c(1, -1), nbinother = 30, 
         bin_option = "quantile", c_max = 0.95, dist = "weibull", alpha = 3,
         psi = 1, maxIter = 1000, tol = 1e-4)

########## Weibull: ###########
# 30% censoring: c_max = 2.2  #
# 70% censoring: c_max = 0.95 #
###############################

########## Log-logistic: ######
# 30% censoring: c_max = 3.2  #
# 70% censoring: c_max = 1    #
###############################