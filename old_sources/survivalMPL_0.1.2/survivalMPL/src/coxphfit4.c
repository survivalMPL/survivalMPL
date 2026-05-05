#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h> 
#include <Rmath.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

void R_init_survivalMPL(DllInfo* info) {
	R_registerRoutines(info, NULL, NULL, NULL, NULL);
	R_useDynamicSymbols(info, TRUE);
}

SEXP coxph_mpl(SEXP status2, SEXP X2, SEXP meanX2, SEXP R2, SEXP psi2, SEXP PSI2, 
	         SEXP Beta0, SEXP Theta0, SEXP gamma2, SEXP kappa2,
	         SEXP convVal2, SEXP minTheta2, SEXP maxiter2)
	{
	/* define objects */
	int i, ii, j, jj, m, n, p, l, u, r, nprotect=0, one_int=1;
	SEXP list_out, list_out_names, list_coef, list_coef_names;
	SEXP list_matrices, list_matrices_names;
	double sum, pen, loglik, varepsilon, correction; 
	double one_double = 1.0;
	double zero_double = 0.0;
	double gamma, minTheta, convVal, kappa;
	int maxiter;
	gamma    = asReal(gamma2);
	minTheta = asReal(minTheta2);
	convVal  = asReal(convVal2);
	kappa    = asReal(kappa2);
	maxiter  = asInteger(maxiter2);
	p = ncols(X2);
	n = nrows(X2);
	m = ncols(R2);
        l = p+m;
        
    	SEXP Beta_old2, StepBeta2, GradBeta2;
     	SEXP Theta_old2, StepTheta2, GradTheta2, sTheta2, PsiTMu2, RTheta2;
    	SEXP Mu2, PsiTheta2, psiTheta2, PsiThetaMu2, StatusMinPsiThetaMu2;
    	SEXP Z2, W2, status_double2;
 	double *Beta_old, *StepBeta, *GradBeta;
	double *Theta_old, *StepTheta, *GradTheta, *sTheta, *PsiTMu, *RTheta;
	double *Mu, *PsiTheta, *psiTheta, *PsiThetaMu, *StatusMinPsiThetaMu;
	double *Z, *W, *status_double;
	PROTECT(Beta_old2 = allocVector(REALSXP, p)); nprotect++;
	Beta_old = REAL(Beta_old2);
	PROTECT(StepBeta2 = allocVector(REALSXP, p)); nprotect++;
	StepBeta = REAL(StepBeta2);
	PROTECT(GradBeta2 = allocVector(REALSXP, p)); nprotect++;
	GradBeta = REAL(GradBeta2);
	PROTECT(Theta_old2 = allocVector(REALSXP, m)); nprotect++;
	Theta_old = REAL(Theta_old2);
	PROTECT(StepTheta2 = allocVector(REALSXP, m)); nprotect++;
	StepTheta = REAL(StepTheta2);
	PROTECT(GradTheta2 = allocVector(REALSXP, m)); nprotect++;
	GradTheta = REAL(GradTheta2);
	PROTECT(sTheta2 = allocVector(REALSXP, m)); nprotect++;	
	sTheta = REAL(sTheta2);
	PROTECT(PsiTMu2 = allocVector(REALSXP, m)); nprotect++;	
	PsiTMu = REAL(PsiTMu2);
 	PROTECT(RTheta2 = allocVector(REALSXP, m)); nprotect++;
	RTheta = REAL(RTheta2);	
	PROTECT(Mu2 = allocVector(REALSXP, n)); nprotect++;
	Mu = REAL(Mu2);
	PROTECT(PsiTheta2 = allocVector(REALSXP, n)); nprotect++;
	PsiTheta = REAL(PsiTheta2);
	PROTECT(psiTheta2 = allocVector(REALSXP, n)); nprotect++;
	psiTheta = REAL(psiTheta2);
	PROTECT(PsiThetaMu2 = allocVector(REALSXP, n)); nprotect++;
	PsiThetaMu = REAL(PsiThetaMu2);
 	PROTECT(StatusMinPsiThetaMu2 = allocVector(REALSXP, n)); nprotect++;
	StatusMinPsiThetaMu = REAL(StatusMinPsiThetaMu2);
 	PROTECT(Z2 = allocVector(REALSXP, n*p)); nprotect++;
	Z = REAL(Z2);
 	PROTECT(W2 = allocVector(REALSXP, n*m)); nprotect++;
	W = REAL(W2);
  	PROTECT(status_double2 = coerceVector(status2, REALSXP)); nprotect++;
	status_double = REAL(status_double2);
      
	double *X, *R, *PSI, *psi, *meanX;
	int *status;
	X = REAL(X2);	
	meanX = REAL(meanX2);
	R = REAL(R2);
	PSI = REAL(PSI2);
	psi = REAL(psi2);
        status = INTEGER(status2);

	SEXP HessianBeta2, inv_inf2, chol_inf2, iter2, ploglik2, Beta2, Theta2, flag2;
	double *ploglik, *Beta, *Theta, *HessianBeta;
	int *chol_inf, *inv_inf, *flag, *iter;
	PROTECT(Beta2 = duplicate(Beta0));nprotect++;
	Beta = REAL(Beta2);
	PROTECT(Theta2 = duplicate(Theta0));nprotect++;
	Theta = REAL(Theta2);	
	PROTECT(HessianBeta2 = allocVector(REALSXP,p*p)); nprotect++;
	HessianBeta = REAL(HessianBeta2);	
	PROTECT(inv_inf2 = allocVector(INTSXP, 1)); nprotect++;
	inv_inf = INTEGER(inv_inf2);
	PROTECT(chol_inf2 = allocVector(INTSXP, 1)); nprotect++;
	chol_inf = INTEGER(chol_inf2);	
	PROTECT(iter2 = allocVector(INTSXP, 1)); nprotect++;
	iter = INTEGER(iter2);	
	PROTECT(ploglik2 = allocVector(REALSXP, 2)); nprotect++;
	ploglik = REAL(ploglik2);
	PROTECT(flag2 = allocVector(INTSXP, 3)); nprotect++;
	flag = INTEGER(flag2);	
	
	/*             */
	/* INITIALISE  */
	/*             */
	
	F77_CALL(dgemm)("N","N",&n,&one_int,&m,&one_double,PSI,&n,Theta,&m,&zero_double,PsiTheta,&n); 
	F77_CALL(dgemm)("N","N",&n,&one_int,&m,&one_double,psi,&n,Theta,&m,&zero_double,psiTheta,&n); 
	F77_CALL(dgemm)("N","N",&n,&one_int,&p,&one_double,X  ,&n,Beta ,&p,&zero_double,Mu      ,&n); 
	F77_CALL(dgemm)("N","N",&m,&one_int,&m,&one_double,R  ,&m,Theta,&m,&zero_double,RTheta  ,&m); 
	
	loglik = 0;
	for(i = 0; i<n; i++){
		Mu[i] = exp(Mu[i]);
		PsiThetaMu[i] = PsiTheta[i]*Mu[i];	
		if(status[i] == 1) loglik += log(psiTheta[i]*Mu[i]);
		loglik -= PsiThetaMu[i];
	}	
	pen = 0; 
	for(u = 0; u<m; u++) pen += Theta[u] * RTheta[u];
	ploglik[0] = (1-gamma)*loglik - gamma*pen;	
	for(i=0; i<3; i++) flag[i] = 0;
	
	/*             		       			*/
	/* 	MAIN LOOP: God Save The Queen  		*/
	/*             			          	*/

 	for(*iter=0; *iter < maxiter; (*iter)++){

	/* UPDATE BETA  */
	 
	for(i = 0; i<n; i++) StatusMinPsiThetaMu[i] = status_double[i]-PsiThetaMu[i];	
	F77_CALL(dgemm)("N","N",&one_int,&p,&n,&one_double,StatusMinPsiThetaMu,&one_int,X,&n,&zero_double,GradBeta,&one_int);
	for(i=0; i<n; i++){
		sum = sqrt(PsiThetaMu[i]);
		for(j=0; j<p; j++) Z[i+j*n]=X[i+j*n]*sum; 
	}
	F77_CALL(dsyrk)("U","T",&p,&n,&one_double,Z,&n,&zero_double,HessianBeta,&p);	
	F77_CALL(dpotrf)("U",&p,HessianBeta,&p,chol_inf);
	F77_CALL(dpotri)("U",&p,HessianBeta,&p,inv_inf);
 	for (j=1; j<p; j++){
		for (jj=0; jj<j; jj++) HessianBeta[j+jj*p] = HessianBeta[jj+j*p]; 
        }
   
	F77_CALL(dgemm)("N","N",&p,&one_int,&p,&one_double,HessianBeta,&p,GradBeta,&p,&zero_double,StepBeta,&p);         
	for(j = 0; j<p; j++){
		Beta_old[j] = Beta[j];
		Beta[j] = Beta_old[j] + StepBeta[j];
		}
		
	ploglik[1] = ploglik[0];  	
	F77_CALL(dgemm)("N","N",&n,&one_int,&p,&one_double,X,&n,Beta,&p,&zero_double,Mu,&n); 
	loglik = 0;
	for(i = 0; i<n; i++){
		Mu[i] = exp(Mu[i]);
		PsiThetaMu[i] = PsiTheta[i]*Mu[i];
		if(status[i] == 1)
			loglik += log(psiTheta[i]*Mu[i]);
		loglik -= PsiThetaMu[i];
	}
	ploglik[0] = (1-gamma)*loglik - gamma*pen;
	
	/* ADAPT NEWTON STEP if needed */
	if(ploglik[0]<ploglik[1]){
		r = 0;
		while(ploglik[0]<ploglik[1]){
			r++;
			for(j = 0; j<p; j++){
				StepBeta[j] /= kappa;
				Beta[j] = Beta_old[j] + StepBeta[j];
				}
			F77_CALL(dgemm)("N","N",&n,&one_int,&p,&one_double,X,&n,Beta,&p,&zero_double,Mu,&n); 
			loglik = 0;
			for(i = 0; i<n; i++){
				Mu[i] = exp(Mu[i]);
				PsiThetaMu[i] = PsiTheta[i]*Mu[i];
				if(status[i] == 1) loglik += log(psiTheta[i]*Mu[i]);
				loglik -= PsiThetaMu[i];
			}
			ploglik[0] = (1-gamma)*loglik - gamma*pen;
			if(r > 500) break;
		}
	}
	
	/* UPDATE THETA */
	
  	for(i=0; i<n; i++){
		for(u=0; u<m; u++) W[i+u*n]=psi[i+u*n]/psiTheta[i];
	}
	F77_CALL(dgemm)("N","N",&one_int,&m,&n,&one_double,status_double,&one_int,W  ,&n,&zero_double,GradTheta,&one_int);
	F77_CALL(dgemm)("N","N",&one_int,&m,&n,&one_double,Mu           ,&one_int,PSI,&n,&zero_double,PsiTMu   ,&one_int);
        for(u=0; u<m; u++){
        	GradTheta[u] = (1-gamma)*(GradTheta[u]-PsiTMu[u])- 2*gamma*RTheta[u];
		sTheta[u] = Theta[u];
		if(RTheta[u]>0)
			sTheta[u] /= (1-gamma)*PsiTMu[u] + 2*gamma*RTheta[u] + 0.3;
		else
			sTheta[u] /= (1-gamma)*PsiTMu[u] + 0.3;
		StepTheta[u] = GradTheta[u]*sTheta[u];
        }
        
       	for(u = 0; u<m; u++){
		Theta_old[u] = Theta[u];
		Theta[u] = Theta_old[u]+StepTheta[u];
		if(Theta[u]<minTheta) Theta[u]= minTheta;
	}
        
   	ploglik[1] = ploglik[0];
	F77_CALL(dgemm)("N","N",&n,&one_int,&m,&one_double,PSI,&n,Theta,&m,&zero_double,PsiTheta,&n); 
	F77_CALL(dgemm)("N","N",&n,&one_int,&m,&one_double,psi,&n,Theta,&m,&zero_double,psiTheta,&n); 
	F77_CALL(dgemm)("N","N",&m,&one_int,&m,&one_double,R  ,&m,Theta,&m,&zero_double,RTheta  ,&m); 
	loglik = 0;
	for(i = 0; i<n; i++){
		PsiThetaMu[i] = PsiTheta[i]*Mu[i];
		if(status[i] == 1) loglik += log(psiTheta[i]*Mu[i]);
		loglik -= PsiThetaMu[i];
	}	
	pen = 0; 
	for(u = 0; u<m; u++) pen += Theta[u] * RTheta[u];
	ploglik[0] = (1-gamma)*loglik - gamma*pen;	

	/* ADAPT NEWTON STEP if needed */
	if(ploglik[0]<ploglik[1]){
		r = 0;
		while(ploglik[0]<ploglik[1]){
			r++;
			for(u = 0; u<m; u++){
				StepTheta[u] /= kappa;
				Theta[u] = Theta_old[u]+StepTheta[u];
				if(Theta[u]<minTheta) Theta[u]= minTheta;
				}
			F77_CALL(dgemm)("N","N",&n,&one_int,&m,&one_double,PSI,&n,Theta,&m,&zero_double,PsiTheta,&n); 
			F77_CALL(dgemm)("N","N",&n,&one_int,&m,&one_double,psi,&n,Theta,&m,&zero_double,psiTheta,&n); 
			F77_CALL(dgemm)("N","N",&m,&one_int,&m,&one_double,R  ,&m,Theta,&m,&zero_double,RTheta  ,&m); 
			loglik = 0;
			for(i = 0; i<n; i++){
				PsiThetaMu[i] = PsiTheta[i]*Mu[i];
				if(status[i] == 1) loglik += log(psiTheta[i]*Mu[i]);
				loglik -= PsiThetaMu[i];
			}	
			pen = 0; 
			for(u = 0; u<m; u++) pen += Theta[u] * RTheta[u];
			ploglik[0] = (1-gamma)*loglik - gamma*pen;	
			if(r > 500) break;
		}
	}
	        
	/*                       */	
	/* Check for convergence */
	/*                       */
	varepsilon = 0;
	for(j=0; j<p; j++){
		sum = fabs(Beta[j]-Beta_old[j]);
		if(sum>varepsilon) varepsilon=sum;
		}
	for(u=0; u<m; u++){
		sum = fabs(Theta[u]-Theta_old[u]);
		if(sum>varepsilon) varepsilon=sum;
		}
	if(varepsilon<convVal) break;
	}

	/* correction factor for Theta */
	correction = 0;
	for(j=0; j<p; j++) correction += -meanX[j]*Beta[j];
	correction = exp(correction);
	ploglik[1] = correction;
	for(u=0; u<m; u++) Theta[u] *= correction;	
	if(*chol_inf!=0) flag[0] += 1;
	if(*inv_inf !=0) flag[0] += 1;
	
	/*             		       			*/
	/* 	INFERENCE: God Save The King  		*/
	/*             			          	*/
	
	SEXP Q2, Ma2, Mb2, H2, Sp2, S2, P2, B2; 
	double *Q, *Ma, *Mb, *H, *Sp, *S, *P, *B;
	PROTECT(Q2 = allocMatrix(REALSXP,l,l)); nprotect++;
	Q = REAL(Q2);
	PROTECT(Ma2 = allocMatrix(REALSXP,l,l)); nprotect++;
	Ma = REAL(Ma2);	
	PROTECT(H2 = allocMatrix(REALSXP,l,l)); nprotect++;	
	H = REAL(H2);	
	PROTECT(Sp2 = allocVector(REALSXP,n*l)); nprotect++;
	Sp = REAL(Sp2);			
	PROTECT(S2 = allocVector(REALSXP,n*l)); nprotect++;
	S = REAL(S2);		
	PROTECT(P2 = allocVector(REALSXP,p*m)); nprotect++;	
	P = REAL(P2);
	PROTECT(B2 = allocVector(REALSXP,m*m)); nprotect++;
	B = REAL(B2);
	
	F77_CALL(dsyrk)("U","T",&p,&n,&one_double,Z,&n,&zero_double,HessianBeta,&p);
 	for (j=1; j<p; j++){
		for (jj=0; jj<j; jj++) HessianBeta[j+jj*p] = 0; 
        }	
	for(i=0; i<n; i++){
		for(j=0; j<p; j++) Z[i+j*n]=X[i+j*n]*Mu[i]; 
	}
	F77_CALL(dgemm)("T","N",&p,&m,&n,&one_double,Z,&n,PSI,&n,&zero_double,P,&p);     
	for(i=0; i<n; i++){
		sum = sqrt(status_double[i]/(psiTheta[i]*psiTheta[i]));
		for(u=0; u<m; u++) W[i+u*n]=psi[i+u*n]*sum; 
	}
	F77_CALL(dsyrk)("U","T",&m,&n,&one_double,W,&n,&zero_double,B,&m);	
	for(j=0;j<p;j++){
		for(jj=j; jj<p; jj++) H[j+jj*l] = HessianBeta[j+jj*p];
		for( u=0;  u<m;  u++) H[j+(p+u)*l] = P[j+u*p];
	}
	for(u=0;u<m;u++){
		for( r=u;  r<m;  r++) H[p+u+(p+r)*l] = B[u+r*m];
	}
	PROTECT(Mb2=duplicate(H2));nprotect++;
	Mb = REAL(Mb2);
	for(u=0;u<m;u++){
		for( r=u;  r<m;  r++) Mb[p+u+(p+r)*l] += 2*gamma/(1-gamma)*R[u+r*m];
	} 
	for(i=0; i<n; i++){
		for(j=0; j<p; j++){
			S[i+j*n]  = -StatusMinPsiThetaMu[i]*X[i+j*n];
			Sp[i+j*n] =  S[i+j*n];
		}
		for(u=0; u<m; u++){
			S[i+(u+p)*n]  = PSI[i+u*n]*Mu[i]-psi[i+u*n]*(status_double[i]/psiTheta[i]);
			Sp[i+(u+p)*n] = S[i+(u+p)*n]+RTheta[u]*2*gamma/(n*(1-gamma));
		}			  
        }
     	F77_CALL(dsyrk)("U","T",&l,&n,&one_double,S,&n,&zero_double,Q,&l);
	F77_CALL(dgemm)("T","N",&l,&l,&n,&one_double,Sp,&n,S,&n,&zero_double,Ma,&l);
	for(i=1; i<l; i++){
		for(ii=0; ii<i; ii++){
			H[i+ii*l]  = H[ii+i*l];
			Q[i+ii*l]  = Q[ii+i*l];
			S[i+ii*l]  = S[ii+i*l];
			Sp[i+ii*l] = Sp[ii+i*l];
			Mb[i+ii*l] = Mb[ii+i*l];
		}
	}
	
	/*             		*/
	/* 	OUTPOUT		*/
	/*                     	*/
	
	PROTECT(list_coef = allocVector(VECSXP, 2)); nprotect++;
	SET_VECTOR_ELT(list_coef, 0, Beta2);
	SET_VECTOR_ELT(list_coef, 1, Theta2);	
	PROTECT(list_coef_names = allocVector(STRSXP, 2)); nprotect++;
	SET_STRING_ELT(list_coef_names, 0, mkChar("Beta"));
	SET_STRING_ELT(list_coef_names, 1, mkChar("Theta"));	
	setAttrib(list_coef, R_NamesSymbol, list_coef_names);
	PROTECT(list_matrices = allocVector(VECSXP,5)); nprotect++;
	SET_VECTOR_ELT(list_matrices, 0, Ma2);
	SET_VECTOR_ELT(list_matrices, 1, Mb2);
	SET_VECTOR_ELT(list_matrices, 2, Q2);
	SET_VECTOR_ELT(list_matrices, 3, H2);	
	SET_VECTOR_ELT(list_matrices, 4, HessianBeta2);			
	PROTECT(list_matrices_names = allocVector(STRSXP, 5)); nprotect++;
	SET_STRING_ELT(list_matrices_names, 0, mkChar("M1"));
	SET_STRING_ELT(list_matrices_names, 1, mkChar("M2"));	
	SET_STRING_ELT(list_matrices_names, 2, mkChar("Q"));	
	SET_STRING_ELT(list_matrices_names, 3, mkChar("H"));
	SET_STRING_ELT(list_matrices_names, 4, mkChar("HessianBeta"));	
	setAttrib(list_matrices, R_NamesSymbol, list_matrices_names);		
	PROTECT(list_out= allocVector(VECSXP, 5)); nprotect++;
	SET_VECTOR_ELT(list_out, 0, list_coef);
	SET_VECTOR_ELT(list_out, 1, iter2);	
	SET_VECTOR_ELT(list_out, 2, ploglik2);	
	SET_VECTOR_ELT(list_out, 3, flag2);	
	SET_VECTOR_ELT(list_out, 4, list_matrices);			
	PROTECT(list_out_names = allocVector(STRSXP, 5)); nprotect++;
	SET_STRING_ELT(list_out_names, 0, mkChar("coef"));
	SET_STRING_ELT(list_out_names, 1, mkChar("iter"));
	SET_STRING_ELT(list_out_names, 2, mkChar("ploglik"));
	SET_STRING_ELT(list_out_names, 3, mkChar("flag"));
	SET_STRING_ELT(list_out_names, 4, mkChar("matrices"));
	setAttrib(list_out, R_NamesSymbol, list_out_names);
	unprotect(nprotect);
	return(list_out);
}





































