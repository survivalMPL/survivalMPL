#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
coxph_mpl=function(formula,data,subset,na.action,control,...){
  #
  mc = match.call(expand.dots = FALSE)
  m  = match(c("formula","data","subset","na.action") ,names(mc),0)
  mc = mc[c(1,m)]    
  if (m[1]==0){stop("A formula argument is required")}
  data.name = if(m[2]!=0){deparse(match.call()[[3]])}else{"-"}
  mc[[1]] = as.name("model.frame")
  mc$formula = if(missing(data)) terms(formula)
  else              terms(formula, data=data)
  mf = eval(mc,parent.frame())
  if (any(is.na(mf))) stop("Missing observations in the model variables")
  if (nrow(mf) ==0) stop("No (non-missing) observations")
  mt = attr(mf,"terms")
  # Y
  y    = model.extract(mf, "response")
  type = attr(y, "type")
  if(!inherits(y, "Surv")){stop("Response must be a survival object")}
  ##
  if (attr(y,which = "type")=="right"){
    left=y[,1]
    right=rep(NA, nrow(y))
    icase = which(y[,2]==1)
    right[icase] = y[icase,1]
    y = Surv(left, right, type="interval2")
  } else if (type!="interval"){
    stop("\nPlease create the survival object using the option type='interval2' in the Surv function.\n")
  }
  ##
  t_i1        = y[,1L]
  t_i2        = y[,2L]
  n       = length(t_i1)
  ctype   = matrix(NA, nrow=n, ncol=4)
  colnames(ctype) = c("r","e","l","i")
  for(tw in 1:4){ctype[,tw] = y[,3L]==(tw-1)}
  n.ctype     = apply(ctype,2,sum)
  ctypeTF     = n.ctype>0    
  observed    = y[,3L]==1L
  n.obs       = sum(y[,3L]!=0)    
  # control arguments
  extraArgs <- list(...)
  if (length(extraArgs)) {
    controlargs <- names(formals(coxph_mpl.control)) 
    m <- pmatch(names(extraArgs), controlargs, nomatch=0L)
    if (any(m==0L))
      stop(gettextf("Argument(s) %s not matched", names(extraArgs)[m==0L]),
           domain = NA, call. = F)
  }    
  if (missing(control)) control <- coxph_mpl.control(n.obs, ...)
  
  # ties 
  t_i1.obs  = t_i1[observed]    
  ties     = duplicated(t_i1.obs)
  if(any(ties)){
    if(control$ties=="epsilon"){
      if(length(control$seed)>0){
        old <- .Random.seed
        on.exit({.Random.seed <<- old})
        set.seed(control$seed)
      }
      t_i1.obs[ties] = t_i1.obs[ties]+runif(sum(ties),-1e-11,1e-11)
      t_i1[observed] = t_i1.obs
    }else{    
      t_i1.obs = t_i1.obs[!ties]
      n.obs   = length(t_i1.obs)
    }
  }
  # X
  X           = model.matrix(mt, mf)#, contrasts)
  X           = X[,!apply(X, 2, function(x) all(x==x[1])), drop=FALSE]
  if(ncol(X)==0){
    X   = matrix(0,n,1)
    noX = TRUE
  }else{  noX = FALSE}
  p           = ncol(X)    
  mean_j      = apply(X, 2, mean)    
  XC          = X - rep(mean_j, each=n)    
  # knot sequence and psi matrices
  knots  = knots_mpl(control,
                     c(t_i1[ctype[,"i"]],t_i2[ctype[,"i"]],t_i1[ctype[,"e"]]-1e-3,t_i1[ctype[,"e"]]+1e-3,
                       t_i1[ctype[,"r"]],t_i1[ctype[,"l"]]))
  
  ###                                    
  ### Estimation                         
  ###                                    
  
  m            = knots$m
  K            = control$max.iter
  s_lambda     = control$smooth
  s_kappa      = control$kappa
  s_t1         = knots$Alpha[1]
  s_tn         = max(knots$Alpha)
  M_R_mm       = penalty_mpl(control,knots)
  M_Rstar_ll   = rbind(matrix(0,p,p+m),cbind(matrix(0,m,p),M_R_mm))    
  s_convlimit  = control$tol    
  M_X_nop     = XC[ctype[,2],,drop=F]
  M_tX_nop    = t(M_X_nop)
  M_psi_nom   = basis_mpl(t_i1,knots,control$basis,control$order,which=1)[ctype[,2],,drop=F] 
  M_tpsi_nom  = t(M_psi_nom)
  M_Psi_nom   = basis_mpl(t_i1,knots,control$basis,control$order,which=2)[ctype[,2],,drop=F] 
  M_tPsi_nom  = t(M_Psi_nom)                                                                 
  M_X_nrp     = XC[ctype[,1],,drop=F]
  M_tX_nrp    = t(M_X_nrp)
  M_Psi_nrm   = basis_mpl(t_i1,knots,control$basis,control$order,which=2)[ctype[,1],,drop=F]     
  M_tPsi_nrm  = t(M_Psi_nrm)                                                     
  M_X_nlp     = XC[ctype[,3],,drop=F]
  M_tX_nlp    = t(M_X_nlp)
  M_Psi_nlm   = basis_mpl(t_i1,knots,control$basis,control$order,which=2)[ctype[,3],,drop=F]                                     
  M_tPsi_nlm  = t(M_Psi_nlm)                                         
  M_X_nip     = XC[ctype[,4],,drop=F]
  M_tX_nip    = t(M_X_nip)
  M_Psi1_nim  = basis_mpl(t_i1,knots,control$basis,control$order,which=2)[ctype[,4],,drop=F]                                     
  M_Psi2_nim  = basis_mpl(t_i2,knots,control$basis,control$order,which=2)[ctype[,4],,drop=F]                                     
  M_tPsi1_nim = t(M_Psi1_nim)                                     
  M_tPsi2_nim = t(M_Psi2_nim)                                     
  
  ###
  ### initialise
  ###
  M_beta_p1  = matrix(0,nrow=p,ncol=1)
  M_theta_m1 = matrix(1,nrow=m,ncol=1)
  s_df       = -1    
  
  ## shortcuts
  M_mu_no1   = exp(M_X_nop%*%M_beta_p1)
  M_mu_nr1   = exp(M_X_nrp%*%M_beta_p1)
  M_mu_nl1   = exp(M_X_nlp%*%M_beta_p1)
  M_mu_ni1   = exp(M_X_nip%*%M_beta_p1)
  M_h0_no1    = M_psi_nom%*%M_theta_m1
  M_H0_no1    = M_Psi_nom%*%M_theta_m1
  M_H0_nr1    = M_Psi_nrm%*%M_theta_m1
  M_H0_nl1    = M_Psi_nlm%*%M_theta_m1
  M_H01_ni1   = M_Psi1_nim%*%M_theta_m1
  M_H02_ni1   = M_Psi2_nim%*%M_theta_m1
  Rtheta      = M_R_mm%*%M_theta_m1
  thetaRtheta = t(M_theta_m1)%*%Rtheta
  TwoLRtheta  = s_lambda*2*Rtheta            
  M_H_no1  = M_H0_no1 * M_mu_no1
  M_H_nr1  = M_H0_nr1 * M_mu_nr1
  M_H_nl1  = M_H0_nl1 * M_mu_nl1
  M_H1_ni1 = M_H01_ni1* M_mu_ni1
  M_H2_ni1 = M_H02_ni1* M_mu_ni1
  M_S_no1  = exp(-M_H_no1)
  M_S_nl1  = exp(-M_H_nl1)
  M_S1_ni1 = exp(-M_H1_ni1)
  M_S2_ni1 = exp(-M_H2_ni1)
  # avoid division by 0
  M_S_nl1[M_S_nl1==1]   = 1-control$epsilon[1]
  M_S1_ni1[M_S1_ni1==1] = 1-control$epsilon[1]
  M_S2_ni1[M_S2_ni1==1] = 1-control$epsilon[1]
  M_S1mS2_ni1           = M_S1_ni1-M_S2_ni1  
  M_S1mS2_ni1[M_S1mS2_ni1<control$epsilon[2]] = control$epsilon[2]
  
  
  ### outer loop
  
  full.iter = 0
  K         = ifelse(control$max.iter[1]>1,control$max.iter[2],control$max.iter[3])
  for(iter in 1:control$max.iter[1]){
    # loglik
    s_lik = 
      sum(log(M_mu_no1)+log(M_h0_no1)-M_H_no1)-
      sum(M_H_nr1)+
      sum(log(1-M_S_nl1))+
      sum(log(M_S1mS2_ni1))-
      s_lambda*thetaRtheta
    
    ### inner loop
    
    for(k in 1:K){
      ## update betas
      s_omega       = 1
      M_beta_p1_OLD = M_beta_p1 
      s_lik_OLD     = s_lik 
      M_gradbeta_p1     = 
        (M_tX_nop%*%(1-M_H_no1))-
        M_tX_nrp%*%M_H_nr1+
        M_tX_nlp%*%(M_S_nl1*M_H_nl1/(1-M_S_nl1))+
        M_tX_nip%*%((M_H2_ni1*M_S2_ni1-M_H1_ni1*M_S1_ni1)/M_S1mS2_ni1)
      M_hessbeta_p1 = 
        M_tX_nop%*%diag(c(M_H_no1),n.ctype[2],n.ctype[2])%*%M_X_nop+
        M_tX_nrp%*%diag(c(M_H_nr1),n.ctype[1],n.ctype[1])%*%M_X_nrp+
        M_tX_nlp%*%diag(c((M_S_nl1/(1-M_S_nl1)^2*M_H_nl1^2-M_S_nl1/(1-M_S_nl1)*M_H_nl1)),n.ctype[3],n.ctype[3])%*%M_X_nlp+            
        M_tX_nip%*%diag(c((M_S1_ni1*M_S2_ni1/M_S1mS2_ni1^2*(M_H2_ni1-M_H1_ni1)^2+
                             (M_S1_ni1*M_H1_ni1-M_S2_ni1*M_H2_ni1)/M_S1mS2_ni1
        )),n.ctype[4],n.ctype[4])%*%M_X_nip   
        # avoid division by 0 (leading to the issue spotted by Kenneth Beath [email of 20220112])
        if(p==1){
            if(M_hessbeta_p1[1,1]==0){M_hessbeta_p1[1,1]=control$epsilon[1]}
        }        
      M_stepbeta_p1 = chol2inv(chol(M_hessbeta_p1))%*%M_gradbeta_p1
      M_beta_p1     = M_beta_p1_OLD+s_omega*M_stepbeta_p1
      M_mu_no1   = exp(M_X_nop%*%M_beta_p1)
      M_mu_nr1   = exp(M_X_nrp%*%M_beta_p1)
      M_mu_nl1   = exp(M_X_nlp%*%M_beta_p1)
      M_mu_ni1   = exp(M_X_nip%*%M_beta_p1)
      M_H_no1  = M_H0_no1 * M_mu_no1
      M_H_nr1  = M_H0_nr1 * M_mu_nr1
      M_H_nl1  = M_H0_nl1 * M_mu_nl1
      M_H1_ni1 = M_H01_ni1* M_mu_ni1
      M_H2_ni1 = M_H02_ni1* M_mu_ni1
      M_S_no1  = exp(-M_H_no1)
      M_S_nl1  = exp(-M_H_nl1)
      M_S1_ni1 = exp(-M_H1_ni1)
      M_S2_ni1 = exp(-M_H2_ni1)
      # avoid division by 0
      M_S_nl1[M_S_nl1==1]   = 1-control$epsilon[1]
      M_S1_ni1[M_S1_ni1==1] = 1-control$epsilon[1]
      M_S2_ni1[M_S2_ni1==1] = 1-control$epsilon[1]
      M_S1mS2_ni1           = M_S1_ni1-M_S2_ni1  
      M_S1mS2_ni1[M_S1mS2_ni1<control$epsilon[2]] = control$epsilon[2]
      # loglik
      s_lik = 
        sum(log(M_mu_no1)+log(M_h0_no1)-M_H_no1)-
        sum(M_H_nr1)+
        sum(log(1-M_S_nl1))+
        sum(log(M_S1mS2_ni1))-
        s_lambda*thetaRtheta
      ## if likelihood decreases
      if(s_lik<s_lik_OLD){
        i       = 0            
        s_omega = 1/s_kappa
        while(s_lik<s_lik_OLD){
          M_beta_p1 = M_beta_p1_OLD+s_omega*M_stepbeta_p1
          M_mu_no1   = exp(M_X_nop%*%M_beta_p1)
          M_mu_nr1   = exp(M_X_nrp%*%M_beta_p1)
          M_mu_nl1   = exp(M_X_nlp%*%M_beta_p1)
          M_mu_ni1   = exp(M_X_nip%*%M_beta_p1)
          M_H_no1  = M_H0_no1 * M_mu_no1
          M_H_nr1  = M_H0_nr1 * M_mu_nr1
          M_H_nl1  = M_H0_nl1 * M_mu_nl1
          M_H1_ni1 = M_H01_ni1* M_mu_ni1
          M_H2_ni1 = M_H02_ni1* M_mu_ni1
          M_S_no1  = exp(-M_H_no1)
          M_S_nl1  = exp(-M_H_nl1)
          M_S1_ni1 = exp(-M_H1_ni1)
          M_S2_ni1 = exp(-M_H2_ni1)
          # avoid division by 0
          M_S_nl1[M_S_nl1==1]   = 1-control$epsilon[1]
          M_S1_ni1[M_S1_ni1==1] = 1-control$epsilon[1]
          M_S2_ni1[M_S2_ni1==1] = 1-control$epsilon[1]
          M_S1mS2_ni1           = M_S1_ni1-M_S2_ni1  
          M_S1mS2_ni1[M_S1mS2_ni1<control$epsilon[2]] = control$epsilon[2]
          # loglik
          s_lik = 
            sum(log(M_mu_no1)+log(M_h0_no1)-M_H_no1)-
            sum(M_H_nr1)+
            sum(log(1-M_S_nl1))+
            sum(log(M_S1mS2_ni1))-
            s_lambda*thetaRtheta
          # update value of omega
          if(s_omega>=1e-2){
            s_omega = s_omega/s_kappa
          }else{if(s_omega<1e-2&s_omega>=1e-5){
            s_omega = s_omega*5e-2
          }else{if(s_omega<1e-5){
            s_omega = s_omega*1e-5    
          }}}
          i = i+1
          if(i>500){break}                
        }
      }
      ## update thetas
      s_nu           = 1
      M_theta_m1_OLD = M_theta_m1 
      s_lik_OLD      = s_lik
      M_gradthetaA_m1     = 
        M_tpsi_nom%*%(1/M_h0_no1)+
        M_tPsi_nlm%*%(M_S_nl1*M_mu_nl1/(1-M_S_nl1))+
        M_tPsi2_nim%*%(M_S2_ni1*M_mu_ni1/(M_S1mS2_ni1))-
        TwoLRtheta*(TwoLRtheta<0)+0.3        
      M_gradthetaB_m1     = 
        M_tPsi_nom%*%M_mu_no1+
        M_tPsi_nrm%*%M_mu_nr1+
        M_tPsi1_nim%*%(M_S1_ni1*M_mu_ni1/(M_S1mS2_ni1))+
        TwoLRtheta*(TwoLRtheta>0)+0.3    
      M_gradtheta_m1 = M_gradthetaA_m1-M_gradthetaB_m1
      M_s_m1         = M_theta_m1/M_gradthetaB_m1
      M_steptheta_p1 = M_s_m1 * M_gradtheta_m1
      M_theta_m1     = M_theta_m1_OLD+s_nu*M_steptheta_p1
      M_theta_m1[M_theta_m1<control$epsilon[2]] = control$epsilon[2]
      M_h0_no1    = M_psi_nom%*%M_theta_m1
      M_H0_no1    = M_Psi_nom%*%M_theta_m1
      M_H0_nr1    = M_Psi_nrm%*%M_theta_m1
      M_H0_nl1    = M_Psi_nlm%*%M_theta_m1
      M_H01_ni1   = M_Psi1_nim%*%M_theta_m1
      M_H02_ni1   = M_Psi2_nim%*%M_theta_m1
      Rtheta      = M_R_mm%*%M_theta_m1
      thetaRtheta = t(M_theta_m1)%*%Rtheta
      TwoLRtheta  = s_lambda*2*Rtheta            
      M_H_no1  = M_H0_no1 * M_mu_no1
      M_H_nr1  = M_H0_nr1 * M_mu_nr1
      M_H_nl1  = M_H0_nl1 * M_mu_nl1
      M_H1_ni1 = M_H01_ni1* M_mu_ni1
      M_H2_ni1 = M_H02_ni1* M_mu_ni1
      M_S_no1  = exp(-M_H_no1)
      M_S_nl1  = exp(-M_H_nl1)
      M_S1_ni1 = exp(-M_H1_ni1)
      M_S2_ni1 = exp(-M_H2_ni1)
      # avoid division by 0
      M_S_nl1[M_S_nl1==1]   = 1-control$epsilon[1]
      M_S1_ni1[M_S1_ni1==1] = 1-control$epsilon[1]
      M_S2_ni1[M_S2_ni1==1] = 1-control$epsilon[1]
      M_S1mS2_ni1           = M_S1_ni1-M_S2_ni1  
      M_S1mS2_ni1[M_S1mS2_ni1<control$epsilon[2]] = control$epsilon[2]
      # loglik
      s_lik = 
        sum(log(M_mu_no1)+log(M_h0_no1)-M_H_no1)-
        sum(M_H_nr1)+
        sum(log(1-M_S_nl1))+
        sum(log(M_S1mS2_ni1))-
        s_lambda*thetaRtheta
      ## if likelihood decreases
      if(s_lik<s_lik_OLD){
        i       = 0
        s_omega = 1/s_kappa
        while(s_lik<s_lik_OLD){
          M_theta_m1= M_theta_m1_OLD+s_omega*M_steptheta_p1
          M_theta_m1[M_theta_m1<control$epsilon[2]] = control$epsilon[2]
          M_h0_no1    = M_psi_nom%*%M_theta_m1
          M_H0_no1    = M_Psi_nom%*%M_theta_m1
          M_H0_nr1    = M_Psi_nrm%*%M_theta_m1
          M_H0_nl1    = M_Psi_nlm%*%M_theta_m1
          M_H01_ni1   = M_Psi1_nim%*%M_theta_m1
          M_H02_ni1   = M_Psi2_nim%*%M_theta_m1
          Rtheta      = M_R_mm%*%M_theta_m1
          thetaRtheta = t(M_theta_m1)%*%Rtheta            
          TwoLRtheta  = s_lambda*2*Rtheta            
          M_H_no1  = M_H0_no1 * M_mu_no1
          M_H_nr1  = M_H0_nr1 * M_mu_nr1
          M_H_nl1  = M_H0_nl1 * M_mu_nl1
          M_H1_ni1 = M_H01_ni1* M_mu_ni1
          M_H2_ni1 = M_H02_ni1* M_mu_ni1
          M_S_no1  = exp(-M_H_no1)
          M_S_nl1  = exp(-M_H_nl1)
          M_S1_ni1 = exp(-M_H1_ni1)
          M_S2_ni1 = exp(-M_H2_ni1)
          # avoid division by 0
          M_S_nl1[M_S_nl1==1]   = 1-control$epsilon[1]
          M_S1_ni1[M_S1_ni1==1] = 1-control$epsilon[1]
          M_S2_ni1[M_S2_ni1==1] = 1-control$epsilon[1]
          M_S1mS2_ni1           = M_S1_ni1-M_S2_ni1  
          M_S1mS2_ni1[M_S1mS2_ni1<control$epsilon[2]] = control$epsilon[2]
          # loglik
          s_lik = 
            sum(log(M_mu_no1)+log(M_h0_no1)-M_H_no1)-
            sum(M_H_nr1)+
            sum(log(1-M_S_nl1))+
            sum(log(M_S1mS2_ni1))-
            s_lambda*thetaRtheta
          # update omega
          if(s_omega>=1e-2){
            s_omega = s_omega/s_kappa
          }else{if(s_omega<1e-2&s_omega>=1e-5){
            s_omega = s_omega*5e-2
          }else{if(s_omega<1e-5){
            s_omega = s_omega*1e-5    
          }}}
          i = i+1
          if(i>500){break}                
        }    
      }
      if(all(c(abs(M_beta_p1-M_beta_p1_OLD),abs(M_theta_m1-M_theta_m1_OLD))<s_convlimit)){break}
      if(control$max.iter[1]==1){if(any(k==seq(0,control$max.iter[3],control$max.iter[2]))){
        control$epsilon[2] = control$epsilon[2]*10}}
    }
    
    # H matrix 
    H = HRinv = matrix(0,p+m,p+m)    
    H[1:p,1:p] = 
      M_tX_nop%*%diag(c(M_H_no1),n.ctype[2],n.ctype[2])%*%M_X_nop+
      M_tX_nrp%*%diag(c(M_H_nr1),n.ctype[1],n.ctype[1])%*%M_X_nrp+
      M_tX_nlp%*%diag(c((M_S_nl1/(1-M_S_nl1)^2*M_H_nl1^2-M_S_nl1/(1-M_S_nl1)*M_H_nl1)),n.ctype[3],n.ctype[3])%*%M_X_nlp+            
      M_tX_nip%*%diag(c((M_S1_ni1*M_S2_ni1/M_S1mS2_ni1^2*(M_H2_ni1-M_H1_ni1)^2+
                           (M_S1_ni1*M_H1_ni1-M_S2_ni1*M_H2_ni1)/M_S1mS2_ni1
      )),n.ctype[4],n.ctype[4])%*%M_X_nip                
    H[1:p,(p+1):(p+m)] =  
      M_tX_nop%*%diag(c(M_mu_no1),n.ctype[2],n.ctype[2])%*%M_Psi_nom+
      M_tX_nrp%*%diag(c(M_mu_nr1),n.ctype[1],n.ctype[1])%*%M_Psi_nrm+
      M_tX_nlp%*%diag(c((M_S_nl1/(1-M_S_nl1)^2*M_H_nl1-M_S_nl1/(1-M_S_nl1))*M_mu_nl1),n.ctype[3],n.ctype[3])%*%M_Psi_nlm+
      M_tX_nip%*%diag(c(M_S1_ni1*M_S2_ni1/M_S1mS2_ni1^2*(M_H2_ni1-M_H1_ni1)*M_mu_ni1),n.ctype[4],n.ctype[4])%*%(M_Psi2_nim-M_Psi1_nim)+
      M_tX_nip%*%diag(c(M_S1_ni1/M_S1mS2_ni1*M_mu_ni1),n.ctype[4],n.ctype[4])%*%M_Psi1_nim-
      M_tX_nip%*%diag(c(M_S2_ni1/M_S1mS2_ni1*M_mu_ni1),n.ctype[4],n.ctype[4])%*%M_Psi2_nim
    H[(p+1):(p+m),1:p] = t(H[1:p,(p+1):(p+m)])               
    H[(p+1):(p+m),(p+1):(p+m)] =  
      M_tpsi_nom%*%diag(c(1/M_h0_no1^2),n.ctype[2],n.ctype[2])%*%M_psi_nom+
      M_tPsi_nlm%*%diag(c(M_S_nl1/(1-M_S_nl1)^2*M_mu_nl1^2),n.ctype[3],n.ctype[3])%*%M_Psi_nlm+
      (M_tPsi2_nim-M_tPsi1_nim)%*%diag(c(M_S1_ni1*M_S2_ni1/M_S1mS2_ni1^2*M_mu_ni1^2),n.ctype[4],n.ctype[4])%*%(M_Psi2_nim-M_Psi1_nim)
    s_lambda_old   = s_lambda
    s_df_old       = s_df
    s_sigma2_old   = 1/(2*s_lambda_old)
    #pos            = c(if(noX){FALSE}else{rep(TRUE,p)},M_theta_m1>control$min.theta)
    pos            = c(if(noX){FALSE}else{rep(TRUE,p)},(M_theta_m1>control$min.theta & apply(H[(p+1):(p+m),(p+1):(p+m)],2,sum)>0))
    #MM: (G+Q)^-1
    temp           = try(chol2inv(chol(H[pos,pos]+(1/s_sigma2_old)*M_Rstar_ll[pos,pos])),silent=T)  
    if(class(temp)[1]!="try-error"&!any(is.infinite(temp))){
      HRinv[pos,pos]=temp
    ##MM: is ginv(H[pos,pos]) right? Forgot Q?
    }else{HRinv[pos,pos]=MASS::ginv(H[pos,pos])}
    s_df           = m-sum(diag(HRinv%*%M_Rstar_ll))/s_sigma2_old
    s_sigma2       = c(t(M_theta_m1)%*%M_R_mm%*%M_theta_m1/s_df)    
    s_lambda       = 1/(2*s_sigma2)
    TwoLRtheta     = s_lambda*2*Rtheta        
    full.iter      = full.iter+k
    if((full.iter/iter)>(control$max.iter[2]*.975)){control$epsilon[2] = control$epsilon[2]*10}        
    if(full.iter>control$max.iter[3]){break}                
    if((k<control$max.iter[2])&
       (abs(s_df-s_df_old)<(control$tol*10))
    ){break}
  }
  s_lambda        = control$smooth = s_lambda_old
  s_correction    = c(exp(-mean_j%*%M_beta_p1)) 
  M_thetatilde_m1 = M_theta_m1
  M_theta_m1      = M_theta_m1*s_correction
  
  ###                                    
  ### Inference                          
  ###                                    
  # M_corr_ll = cbind(rbind(diag(rep(1,p)),s_correction*matrix(rep(M_thetatilde_m1,p)*rep(-mean_j,each=m),ncol=p)),
                    # rbind(matrix(0,ncol=m,nrow=p),diag(rep(s_correction,m))))
  M_corr_ll = cbind(rbind(diag(rep(1,p)),s_correction*matrix(rep(M_thetatilde_m1,p)*rep(-M_beta_p1,each=m),ncol=p)),
                    rbind(matrix(0,ncol=m,nrow=p),diag(rep(s_correction,m))))
  M_corr_ll[!pos,] = 0

  
  M_2    = H+2*s_lambda*M_Rstar_ll
  Q = matrix(NA,n,p+m)
  if(ctypeTF[1]){Q[ctype[,1],1:p] = rep(-M_H_nr1,p)*M_X_nrp}
  if(ctypeTF[2]){Q[ctype[,2],1:p] = rep((1-M_H_no1),p)*M_X_nop}
  if(ctypeTF[3]){Q[ctype[,3],1:p] = rep(M_S_nl1*M_H_nl1/(1-M_S_nl1),p)*M_X_nlp}
  if(ctypeTF[4]){Q[ctype[,4],1:p] = rep((M_H2_ni1*M_S2_ni1-M_H1_ni1*M_S1_ni1)/M_S1mS2_ni1,p)*M_X_nip}
  if(ctypeTF[1]){Q[ctype[,1],-c(1:p)] = rep(-M_mu_nr1,m)*M_Psi_nrm}
  if(ctypeTF[2]){Q[ctype[,2],-c(1:p)] = rep(1/M_h0_no1,m)*M_psi_nom-rep(M_mu_no1,m)*M_Psi_nom}
  if(ctypeTF[3]){Q[ctype[,3],-c(1:p)] = rep(M_S_nl1*M_mu_nl1/(1-M_S_nl1),m)*M_Psi_nlm}
  if(ctypeTF[4]){Q[ctype[,4],-c(1:p)] = rep(M_S2_ni1*M_mu_ni1/(M_S1mS2_ni1),m)*M_Psi2_nim-
    rep(M_S1_ni1*M_mu_ni1/(M_S1mS2_ni1),m)*M_Psi1_nim}
  Sp = Q-matrix(rep(c(rep(0,p),TwoLRtheta),n),n,byrow=T)/n
  Q = t(Sp)%*%Sp
  #pos   = c(if(noX){FALSE}else{rep(TRUE,p)},M_theta_m1>control$min.theta)
  pos            = c(if(noX){FALSE}else{rep(TRUE,p)},(M_theta_m1>control$min.theta & apply(H[(p+1):(p+m),(p+1):(p+m)],2,sum)>0))
  Minv_1 = Minv_2 = Hinv = matrix(0,p+m,p+m)                        
  temp = try(chol2inv(chol(M_2[pos,pos])),silent=T) 
  if(class(temp)[1]!="try-error"){
    Minv_2[pos,pos] = temp        
    cov_NuNu_M2QM2  = M_corr_ll%*%(Minv_2%*%Q%*%Minv_2)%*%t(M_corr_ll)
    cov_NuNu_M2HM2  = M_corr_ll%*%(Minv_2%*%H%*%Minv_2)%*%t(M_corr_ll)
    se.Eta_M2QM2    = sqrt(diag(cov_NuNu_M2QM2))
    se.Eta_M2HM2    = sqrt(diag(cov_NuNu_M2HM2))            
  }else{
    cov_NuNu_M2QM2  = cov_NuNu_M2HM2 = matrix(NA,p+m,p+m)
    se.Eta_M2QM2    = se.Eta_M2HM2   = rep(NA,p+m)
  }
  temp = try(chol2inv(chol(H[pos,pos])),silent=T) 
  if(class(temp)[1]!="try-error"){
    Hinv[pos,pos] = temp
    cov_NuNu_H    = M_corr_ll%*%Hinv%*%t(M_corr_ll) 
    se.Eta_H      = sqrt(diag(cov_NuNu_H))
  }else{
    cov_NuNu_H  = matrix(NA,p+m,p+m)
    se.Eta_H    = rep(NA,p+m)
  }
  mx.seNu.l5=as.data.frame(cbind(se.Eta_M2QM2,se.Eta_M2HM2,se.Eta_H))
  colnames(mx.seNu.l5)=c("M2QM2","M2HM2","H")
  rownames(mx.seNu.l5)[(p+1):(p+m)] = paste("Theta",1:m,sep="")        
  rownames(mx.seNu.l5)[(1:p)] = paste("Beta",1:p,sep="")
  fit         = list(coef=list(Beta=c(M_beta_p1),Theta=c(M_theta_m1)*(M_theta_m1>control$min.theta)),
                     iter = c(iter,full.iter))
  fit$se      = list(Beta=mx.seNu.l5[1:p,],Theta=mx.seNu.l5[(p+1):(p+m),])
  fit$covar   = list(M2QM2=cov_NuNu_M2QM2,M2HM2=cov_NuNu_M2HM2,H=cov_NuNu_H)      
  fit$knots   = knots
  fit$control = control
  fit$call    = match.call()
  fit$dim     = list(n = n, n.obs = sum(observed), n.ties = sum(ties), p = p, m = knots$m)
  fit$data    = list(time = y, censoring=y[,3L], X = X, name = data.name)# list(name = data.name)#
  fit$df      = s_df
  fit$ploglik = s_lik
  fit$loglik  = s_lik+s_lambda*thetaRtheta
  class(fit)  = "coxph_mpl"
  fit
}







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
coxph_mpl.control <- function(n.obs=NULL, basis = "uniform", smooth = NULL, max.iter=c(1.5e+2,7.5e+4,1e+6), tol=1e-7, 
                              n.knots = NULL, n.events_basis = NULL, range.quant = c(0.075,.9),
                              cover.sigma.quant = .25, cover.sigma.fixed=.25, min.theta = 1e-10,
                              penalty = 2L, order = 3L, kappa = 1/.6, epsilon = c(1e-16,1e-10), ties = "epsilon", seed = NULL){
  basis        = basis.name_mpl(basis)
  max.iter     = c(ifelse(is.null(smooth),ifelse(max.iter[1]>0,as.integer(max.iter[1]),1.5e+2),1L),
                   ifelse(max.iter[2]>0,as.integer(max.iter[2]),7.5e+4),
                   ifelse(length(max.iter)==2,1e+6,
                          ifelse(max.iter[3]>ifelse(max.iter[2]>0,as.integer(max.iter[2]),7.5e+4),
                                 as.integer(max.iter[3]),1e+6)))    
  tol          = ifelse(tol>0 & tol<1,tol,1e-7)    
  order        = ifelse(order>0 & order<6,as.integer(order),3L)    
  min.theta    = ifelse(min.theta>0 & min.theta<1e-3,min.theta,1e-10)    
  penalty      = penalty.order_mpl(penalty,basis,order)
  kappa        = ifelse(kappa>1, kappa, 1/.6)
  cover.sigma.quant  = ifelse(cover.sigma.quant>0 & cover.sigma.quant<0.4,cover.sigma.quant,.75)
  cover.sigma.fixed  = ifelse(cover.sigma.fixed>0 & cover.sigma.fixed<0.4,cover.sigma.fixed,.75)    
  if(all(range.quant<=1) & all(range.quant>=0) & length(range.quant)==2){
    range.quant = range.quant[order(range.quant)]
  }else{range.quant = c(0.075,.9)}
  if(is.null(n.knots)|sum(n.knots)<3|length(n.knots)!=2){
    n.knots    = if(basis!='uniform' & basis!='msplines'){c(0,20)}else{c(8,2)}
  }
  if(!is.null(n.events_basis)){
    n.events_basis = ifelse(n.events_basis<1|n.events_basis>floor(n.obs/2),
                            max(round(3.5*log(n.obs)-7.5),1L),round(n.events_basis))
  }else{n.events_basis = max(round(3.5*log(n.obs)-7.5),1L)}    
  if(!is.null(smooth)){
    smooth = ifelse(smooth<0,0,smooth)
  }else{smooth=0}
  out = list(basis = basis, smooth = smooth, max.iter = max.iter, tol = tol,
             order = order, penalty = penalty, n.knots = n.knots, range.quant = range.quant,
             cover.sigma.quant = cover.sigma.quant, cover.sigma.fixed = cover.sigma.fixed,
             n.events_basis = as.integer(n.events_basis), min.theta = min.theta, ties = ties,
             seed = as.integer(seed), kappa = kappa, epsilon = epsilon)
  class(out) = "coxph_mpl.control"
  out
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
basis.name_mpl <- function(k){
  if(k == "discr"| k == "discretized" | k == "discretised" | k == "unif" | k == "uniform"){"uniform"
  }else{if(k == "m" | k == "msplines" | k == "mspline"){"msplines"
  }else{if(k == "gauss" | k == "gaussian"){"gaussian"
  }else{if(k == "epa" | k == "epanechikov"){"epanechikov"
  }else{stop("Unkown basis choice", call. = FALSE)}}}}
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
penalty.order_mpl <- function(p,basis,order){
  p = as.integer(p)
  switch(basis,
         'uniform'  = ifelse(p>0 & p<3,p,2),
         'gaussian' = ifelse(p>0 & p<3,p,2),
         'msplines' = order-1,
         'epa'      = 2)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
knots_mpl=function(control,events){
  n.events  = length(events)
  range = range(events)  
  if(control$n.knots[2]==0){
    Alpha    = quantile(events,seq(0,1,length.out=(control$n.knots[1]+2)))
  }else{
    Alpha1   = quantile(events,seq(0,control$range.quant[2],length.out=(control$n.knots[1]+1)))    
    Alpha2   = seq(quantile(events,control$range.quant[2]),range(events)[2],length=control$n.knots[2]+2)
    Alpha    = c(Alpha1,Alpha2[-1])
  }
  n.Alpha  = length(Alpha)
  if(control$basis=="gaussian"){
    Sigma = Delta = rep(0,n.Alpha)
    for(aw in 1:n.Alpha){
      if(aw>1 & aw<(n.Alpha-control$n.knots[2])){
        while(sum(events>(Alpha[aw]-2*Sigma[aw])&events<(Alpha[aw]+2*Sigma[aw]))<(n.events*control$cover.sigma.quant)){
          Sigma[aw] = Sigma[aw] + 0.001}
      }else{Sigma[aw] = control$cover.sigma.fixed*(Alpha[n.Alpha]-Alpha[1])/3}
      Delta[aw]= pnorm((range[2]-Alpha[aw])/Sigma[aw])-
        pnorm((range[1]-control$epsilon[1]-Alpha[aw])/Sigma[aw])
    }
    list(m=n.Alpha, Alpha=Alpha, Sigma=Sigma, Delta=Delta)
  }else{if(control$basis=="msplines"|control$basis=="epanechikov"){
    m = n.Alpha+control$order-2
    list(m=m, Alpha=Alpha, Delta=rep(1,m))
  }else{if(control$basis=="uniform"){
    m          = length(Alpha)-1
    Delta      = Alpha[2L:(m+1L)]-Alpha[1L:m]
    list(m=m,Alpha=Alpha,Delta=Delta)
  }else{stop("Unkown basis choice")}}}
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
basis_mpl = function(x,knots,basis,order,which=c(1,2)){
  which.matrix = rep(T,2)
  which.matrix[-which]=FALSE    
  n        = length(x)
  Alpha    = knots$Alpha
  Delta    = knots$Delta
  n.Alpha  = length(Alpha)
  m        = ifelse(basis=="msplines"|basis=="epanechikov",n.Alpha+order-2,knots$m)
  M_Psi_nm = M_psi_nm = matrix(0,n,m)
  ##
  if(basis=="uniform"){
    u_i = sapply(x,function(y,lim=Alpha[-1L])sum(lim<y)+1L)
    for(i in 1:n){
      M_psi_nm[i,u_i[i]]   = 1
      M_Psi_nm[i,1:u_i[i]] = c(if(u_i[i]>1){Delta[1:(u_i[i]-1)]},
                               x[i]-Alpha[u_i[i]])
    }
    ##
  }else{
    if(basis=="gaussian"){
      Sigma = knots$Sigma
      for(u in 1:m){
        M_psi_nm[,u] =  dnorm((x-Alpha[u])/Sigma[u])/(Sigma[u]*Delta[u])
        M_Psi_nm[,u] = (pnorm((x-Alpha[u])/Sigma[u])-
                          pnorm((Alpha[1]-Alpha[u])/Sigma[u]))/Delta[u]
      }
      ##
    }else{
      seq1n = 1:n
      Alpha_star   = as.numeric(c(rep(Alpha[1],order-1L),Alpha,rep(Alpha[n.Alpha],order-1L)))
      M_psi_nm     = M_Psi_nm = cbind(M_psi_nm,0)        
      if(which.matrix[1]){
        Alpha_star_x = sapply(x,function(y,lim=Alpha[-1L])sum(lim<y)+1L)+order-1L
        if(basis=="msplines"){
          M_psi_nm[(Alpha_star_x-1L)*n+seq1n]=1/(Alpha_star[Alpha_star_x+1]-Alpha_star[Alpha_star_x])
          if(order>1){
            for(ow in 2L:order){
              uw_x = Alpha_star_x-ow+1L
              for(pw in 0:(ow-1L)){
                pos_x = (uw_x+pw-1L)*n+seq1n
                M_psi_nm[pos_x]=(ow/((ow-1)*(Alpha_star[1:m+ow]-Alpha_star[1:m])))[uw_x+pw]*
                  ((x-Alpha_star[uw_x+pw])*M_psi_nm[pos_x]+
                     (Alpha_star[uw_x+pw+ow]-x)*M_psi_nm[pos_x+n])    
              }
            }
          }
        }else{
          uw_x = Alpha_star_x-order+1L
          for(pw in 0:(order-1L)){
            pos_x = (uw_x+pw-1L)*n+seq1n
            pos_1 = (uw_x+pw)==1
            pos_m = (uw_x+pw)==m
            pos_other = pos_1==FALSE & pos_m==FALSE    
            M_psi_nm[pos_x[pos_other]]=(6*(x-Alpha_star[uw_x+pw])*(x-Alpha_star[uw_x+pw+order])/ 
                                          ((Alpha_star[uw_x+pw]-Alpha_star[uw_x+pw+order])^3))[pos_other]
            M_psi_nm[pos_x[pos_1]]=(12*(x-Alpha_star[uw_x+pw+order])*(x-2*Alpha_star[uw_x+pw]+Alpha_star[uw_x+pw+order])/ 
                                      ((2*Alpha_star[uw_x+pw]-2*Alpha_star[uw_x+pw+order])^3))[pos_1]
            M_psi_nm[pos_x[pos_m]]=(12*(x-Alpha_star[uw_x+pw])*(x+Alpha_star[uw_x+pw]-2*Alpha_star[uw_x+pw+order])/ 
                                      ((2*Alpha_star[uw_x+pw]-2*Alpha_star[uw_x+pw+order])^3))[pos_m]
          }
        }
        M_psi_nm = M_psi_nm[,1:m,drop=FALSE]
      }    
      if(which.matrix[2]){
        rank.x   = rank(x)
        x        = x[order(x)]    
        Alpha_x  = sapply(x,function(y,lim=Alpha[-1L])sum(lim<y)+1L)    
        up_u     = cumsum(tabulate(Alpha_x,n.Alpha-1))
        for(uw in 1:(m-order+1)){M_Psi_nm[min(n,up_u[uw]+1):n,uw] = 1}    
        if(basis=="msplines"){
          Alpha_star2 = c(rep(Alpha[1],order),Alpha,rep(Alpha[n.Alpha],order))    
          factor_v    = c((Alpha_star2[(order+2):length(Alpha_star2)]-Alpha_star2[1:(length(Alpha_star2)-order-1)])/
                            (order+1),rep(0,order-1))
          M_psi2_nm   = cbind(basis_mpl(x,knots,basis=basis,order=order+1,which=1),matrix(0,n,order-1))
          pos_xo  = rep((Alpha_x-1L)*n,1)+seq1n
          pos_xo1 = rep(pos_xo,order)+rep(1:order,each=n)*n
          for(ow in 0:(order-1)){
            M_Psi_nm[pos_xo+ow*n] = apply(matrix(M_psi2_nm[pos_xo1+ow*n]*
                                                   factor_v[rep(Alpha_x,order)+rep((1:order)+ow,each=n)],ncol=order),1,sum)
          }
        }else{
          Alpha_star_x = sapply(x,function(y,lim=Alpha[-1L])sum(lim<y)+1L)+order-1L
          uw_x = Alpha_star_x-order+1L
          for(pw in 0:(order-1L)){
            pos_x = (uw_x+pw-1L)*n+seq1n
            pos_1 = (uw_x+pw)==1
            pos_m = (uw_x+pw)==m
            pos_other = pos_1==FALSE & pos_m==FALSE    
            M_Psi_nm[pos_x[pos_other]]=((x-Alpha_star[uw_x+pw])^2*(2*x+Alpha_star[uw_x+pw]-3*Alpha_star[uw_x+pw+order])/ 
                                          ((Alpha_star[uw_x+pw]-Alpha_star[uw_x+pw+order])^3))[pos_other]
            M_Psi_nm[pos_x[pos_1]]=((x-Alpha_star[uw_x+pw])*(x^2-2*x*Alpha_star[uw_x+pw]-2*Alpha_star[uw_x+pw]^2+
                                                               6*Alpha_star[uw_x+pw]*Alpha_star[uw_x+pw+order]-3*Alpha_star[uw_x+pw+order]^2)/ 
                                      (2*(Alpha_star[uw_x+pw]-Alpha_star[uw_x+pw+order])^3))[pos_1]
            M_Psi_nm[pos_x[pos_m]]=((x-Alpha_star[uw_x+pw])^2*(x+2*Alpha_star[uw_x+pw]-3*Alpha_star[uw_x+pw+order])/ 
                                      (2*(Alpha_star[uw_x+pw]-Alpha_star[uw_x+pw+order])^3))[pos_m]
          }
        }
        M_Psi_nm = M_Psi_nm[rank.x,1:m,drop=FALSE]
      }}}
  if(all(which.matrix)){list(psi=M_psi_nm,Psi=M_Psi_nm)
  }else{if(which.matrix[1]){M_psi_nm}else{M_Psi_nm}}
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
penalty_mpl=function(control,knots){
  #
  penalty = control$penalty
  m       = knots$m
  if(control$basis=="uniform"){
    D = diag(m)*c(1,-2)[penalty]
    E = diag(m+1)[-(m+1),-1]*c(-1,1)[penalty]
    B = D+E+list(0,t(E))[[penalty]]
    B[c(m+1,(m-1)*m)[1:penalty]] = c(-1,2)[penalty]
    M_R_mm = t(B)%*%B        
  }else{
    M_R_mm = matrix(0,m,m)
    if(control$basis=="gaussian"){
      int_rij_2.fun=function(x,mu_i,mu_j,sig_i,sig_j,t1,tn){
        K     = 4*(pnorm((t1-mu_i)/sig_i)-pnorm((tn-mu_i)/sig_i))*(pnorm((t1-mu_j)/sig_j)-pnorm((tn-mu_j)/sig_j))
        q1q2a = 4*dnorm(x,mu_i,sig_i)*dnorm(x,mu_j,sig_j)*sig_i*sig_j*2*pi*sqrt(sig_i^2+sig_j^2)*
          ((mu_j-x)*sig_i^6*((x-mu_i)^2-sig_i^2)+
             sig_i^4*sig_j^2*((x-4*mu_j+3*mu_i)*sig_i^2-(x-mu_i)^2*(3*x-4*mu_j+mu_i))-
             sig_i^2*sig_j^4*(x-mu_j)^2*(3*x-4*mu_i+mu_j)+
             sig_j^6*((x+3*mu_j-4*mu_i)*sig_i^2-(x-mu_j)^2*(x-mu_i))+
             sig_j^8*(x-mu_i)
          )
        q1q2b = 2*dnorm(mu_j,mu_i,sqrt(sig_i^2+sig_j^2))*pi*sqrt(sig_i^2+sig_j^2)*sig_i^3*sig_j^3*
          (mu_j^4-4*mu_j^3*mu_i+mu_i^4+6*mu_j^2*(mu_i^2-sig_i^2-sig_j^2)-
             6*mu_i^2*(sig_i^2+sig_j^2)+3*(sig_i^2+sig_j^2)^2-4*mu_i*mu_j*(mu_i^2-3*(sig_i^2+sig_j^2))
          )*
          (2*pnorm(((x-mu_j)*sig_i^2+(x-mu_i)*sig_j^2)/(sig_i*sig_j*sqrt(sig_i^2+sig_j^2)))-1)
        q3    = pi*sig_i^3*sig_j^3*(sig_i^2+sig_j^2)^(9/2)*K
        (q1q2a+q1q2b)/q3
      }
      int_rij_1.fun=function(x,mu_i,mu_j,sig_i,sig_j,t1,tn){
        K  = 4*(pnorm((t1-mu_i)/sig_i)-pnorm((tn-mu_i)/sig_i))*(pnorm((t1-mu_j)/sig_j)-pnorm((tn-mu_j)/sig_j))
        q2 =  dnorm((mu_i-mu_j)/sqrt(sig_i^2+sig_j^2))*2*pi*
          sig_i*sig_j*(sig_i^2+sig_j^2-(mu_i-mu_j)^2)*
          (2*pnorm(((x-mu_j)*sig_i^2+(x-mu_i)*sig_j^2)/(sig_i*sig_j*sqrt(sig_i^2+sig_j^2)))-1)
        q1 =   4*pi*dnorm((x-mu_i)/sig_i)*dnorm((x-mu_j)/sig_j)*
          sqrt(sig_i^2+sig_j^2)*((mu_i-x)*sig_i^2+(mu_j-x)*sig_j^2)
        q3 =  pi*sig_i*sig_j*(sig_i^2+sig_j^2)^(5/2)*K
        (q1+q2)/q3    
      }
      for(i in 1:m){
        for(j in i:m){
          if(penalty==2){
            M_R_mm[i,j] = M_R_mm[j,i] = 
              int_rij_2.fun(knots$Alpha[m],knots$Alpha[i],knots$Alpha[j],knots$Sigma[i],knots$Sigma[j],knots$Alpha[1],knots$Alpha[m])-
              int_rij_2.fun(knots$Alpha[1],knots$Alpha[i],knots$Alpha[j],knots$Sigma[i],knots$Sigma[j],knots$Alpha[1],knots$Alpha[m])
          }else{  M_R_mm[i,j] = M_R_mm[j,i] = 
            int_rij_1.fun(knots$Alpha[m],knots$Alpha[i],knots$Alpha[j],knots$Sigma[i],knots$Sigma[j],knots$Alpha[1],knots$Alpha[m])-
            int_rij_1.fun(knots$Alpha[1],knots$Alpha[i],knots$Alpha[j],knots$Sigma[i],knots$Sigma[j],knots$Alpha[1],knots$Alpha[m])}
        }
      }
    }else{  Alpha      = knots$Alpha
    n.Alpha    = length(Alpha)
    order      = control$order
    Alpha_star = c(rep(Alpha[1],order-1L),Alpha,rep(Alpha[n.Alpha],order-1L))
    if(control$basis=="msplines"){
      seq1n        = 1L:(n.Alpha-1)
      n.Alpha_star = length(Alpha_star)
      Alpha_star_x = sapply(Alpha[-1],function(y,lim=Alpha[-1L])sum(lim<y)+1L)+order-1L                
      M_d2f_mm = matrix(0,n.Alpha-1,n.Alpha+order-1L)
      M_d2f_mm[(Alpha_star_x-1L)*(n.Alpha-1)+seq1n]=1/(Alpha_star[Alpha_star_x+1]-Alpha_star[Alpha_star_x])
      for(ow in 2L:order){
        pw   = 1L:ow 
        uw_x = Alpha_star_x-ow+1L
        for(pw in 0:(ow-1L)){
          M_d2f_mm[(uw_x+pw-1L)*(n.Alpha-1)+seq1n]=
            (ow/(Alpha_star[1:(n.Alpha+ow)+ow]-Alpha_star[1:(n.Alpha+ow)]))[uw_x+pw]*
            (M_d2f_mm[(uw_x+pw-1L)*(n.Alpha-1)+seq1n]-M_d2f_mm[(uw_x+pw)*(n.Alpha-1)+seq1n])
        }
      }
      M_d2f_mm = M_d2f_mm[,1:m,drop=FALSE]
      for(uw in 1:m){
        for(vw in uw:m){
          M_R_mm[uw,vw] = M_R_mm[vw,uw] = 
            sum((M_d2f_mm[,uw]*M_d2f_mm[,vw])*(Alpha[-1]-Alpha[-n.Alpha]))
        }
      }
    }else{    for(uw in 1:m){
      f_u = ifelse(uw==1|uw==m,4,1)
      for(vw in uw:m){
        if(Alpha_star[vw]<Alpha_star[uw+order]){
          f_v = ifelse(vw==1|vw==m,4,1)
          M_R_mm[uw,vw] = M_R_mm[vw,uw] = (144*(Alpha_star[uw+order]-Alpha_star[vw]))/
            (f_v*f_u*(Alpha_star[uw+order]-Alpha_star[uw])^3*(Alpha_star[vw+order]-Alpha_star[vw])^3)
        }
      }
    }
    }
    }}
  M_R_mm
}








