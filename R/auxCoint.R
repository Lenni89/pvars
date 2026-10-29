

#############################
###  AUXILIARY FUNCTIONS  ###
#############################
#
# The following functions serve as modules nested 
# in the calling functions. Notation within each 
# function mostly corresponds to the cited literature.


# check, define and stack matrices for RRR, from Johansen 1995:90, ch.6
aux_stackRRR <- function(y, dim_p, x=NULL, dim_q=0, type=NULL, t_D1=list(), t_D2=list(), D1=NULL, D2=NULL){
  ### see also Juselius 2007:115, ch.7 / Luetkepohl,Kraetzig 2004:93, ch.3.3.2 / Luetkepohl 2005:294, ch.7.2.3 (ch.10.3.2 for conditional VECM)
  # endogenous variables
  y = aux_asDataMatrix(y, "y")  # named matrix y is KxT regardless of input 
  y_diff = t(diff(t(y)))
  dim_K  = nrow(y)  # number of endogenous variables
  
  dim_pq = max(dim_p, dim_q)   # overall lag-order
  dim_T  = ncol(y) - dim_pq    # number of observations without presample
  idx_t  = dim_pq:(ncol(y)-1)  # index for excluding time periods t=1,..,(p-1) and T from Z0 and Z1
  
  # weakly exogenous variables
  if(!is.null(x)){
    x = aux_asDataMatrix(x, "x") # named matrix x is LxT regardless of input
    x_diff = t(diff(t(x)))
    dim_L  = nrow(x) # number of weakly exogenous variables
  }else{
    x_diff = NULL
    dim_L  = 0
  }
  
  # deterministic term
  if(!is.null(type)){
    # determin. term:  --restricted--   ; --unrestricted-- , from Johansen 1995:96, ch.6.2 / Juselius 2007, ch.6.3
    if(type=="Case1"){ type_D1 = "none" ; type_D2 = "none"  }else
    if(type=="Case2"){ type_D1 = "const"; type_D2 = "none"  }else
    if(type=="Case3"){ type_D1 = "none" ; type_D2 = "const" }else
    if(type=="Case4"){ type_D1 = "trend"; type_D2 = "const" }else
    if(type=="Case5"){ type_D1 = "none" ; type_D2 = "both"  }else{ 
      stop("Incorrect specification of the deterministic case!") }
    
    t2_shift = c(t_D1$t_break, t_D2$t_shift)
    if(!is.null(t2_shift)){
      t2_shift = unique(sort( t2_shift ))  # prevent duplicated dummies, which would lead to perfect multicollinearity
    }
    ### lagged shift dummies for t_break are passed to the impulse dummies instead, see Trenkler et al 2008:335
    t2_impulse = c(t_D1$t_break, t_D1$t_shift)
    if(!is.null(t2_impulse)){
      t2_impulse = sapply(t2_impulse, FUN=function(t) t + 0:(dim_pq-1))  # sacrifice degrees-of-freedom against non-linear LS estimation
      t2_impulse = unique(sort( c(t2_impulse, t_D2$t_impulse) ))  # prevent duplicated dummies, which would lead to perfect multicollinearity
    }else{
      t2_impulse = unique(sort( t_D2$t_impulse ))
    }
    D1 = aux_dummy(dim_T=ncol(y), type=type_D1, dummies=D1, t_break=t_D1$t_break, t_shift=t_D1$t_shift, t_impulse=t_D1$t_impulse)
    D2 = aux_dummy(dim_T=ncol(y), type=type_D2, dummies=D2, t_shift=t2_shift, t_impulse=t2_impulse, t_blip=t_D2$t_blip, n.season=t_D2$n.season, center=TRUE)
    ### see Saikkonen,Luetkepohl 2001:6, Eq.2.11 / Luetkepohl et al 2004:650 / Trenkler et al 2008:335 with prior-detrending
    ### and Johansen et al 2001:219 / Johansen 2016:7, Eq.17 for generalization of the additive and innovative representation
  }
  
  # stack variables for full or partial VECM, from Johansen 1992:393
  Z0 = y_diff[ , idx_t, drop=FALSE]  # regressand matrix, first-differenced and without presample  
  Z1 = rbind(y, x, D1)[ , idx_t, drop=FALSE]  # regressors in levels for the lagged cointegrating term
  Z2 = aux_stack(y=y_diff, dim_p=dim_p-1, x=x_diff, dim_q=dim_q-1, D=D2[ ,-1, drop=FALSE])  # regressors in first differences for short-run effects
  
  # return result
  rownames(Z1)[1:(dim_K+dim_L)] = paste0(rownames(Z1)[1:(dim_K+dim_L)], ".l1")
  result = list(Z0=Z0, Z1=Z1, Z2=Z2, D1=D1, D2=D2, y=y, x=x, 
                dim_p=dim_p, dim_q=dim_q, dim_T=dim_T, dim_K=dim_K, dim_L=dim_L)
  return(result)
}


# reduced rank ML-estimation, from Johansen 1995:90, ch.6
aux_RRR <- function(Z0, Z1, Z2, via_R0_R1=FALSE){
  # define
  dim_T = ncol(Z0) # number of observations without presample
  
  # product moment matrices, from Johansen 1995:90, Eq.6.4
  M00 = tcrossprod(Z0) / dim_T  # = Z0%*%t(Z0)/dim_T
  M11 = tcrossprod(Z1) / dim_T
  M22 = tcrossprod(Z2) / dim_T
  M01 = tcrossprod(Z0, Z1) / dim_T
  M02 = tcrossprod(Z0, Z2) / dim_T
  M12 = tcrossprod(Z1, Z2) / dim_T
  M10 = t(M01)  # use transpose instead of Z1%*%t(Z0)/dim_T, from Johansen 1995:90, Eq.6.4(II)
  M20 = t(M02)
  M21 = t(M12)
  if(nrow(M22)==0){ M22inv = M22 }else{  # if matrix Z2 is 0xT (i.e. VECM has neither first-differenced regressors nor unrestricted deterministic term)
    M22inv = solve(M22)
  }
  
  # concentrate out short-run effects by multivariate OLS-annihilator matrix, from Johansen 1995:90, Eq.6.6 and Eq.6.7
  if(via_R0_R1){
    R0 = Z0 - M02 %*% M22inv %*% Z2  # I(0)-residuals using first differences Z0 as regressand
    R1 = Z1 - M12 %*% M22inv %*% Z2  # I(1)-residuals using lagged levels Z1 as regressand
  }else{ R0 = R1 = NULL }
  
  # product moment matrices for reduced rank regression, from Johansen 1995:90, Eq.6.10
  S00 = M00 - M02 %*% M22inv %*% M20  # = tcrossprod(R0, R0) / dim_T
  S01 = M01 - M02 %*% M22inv %*% M21  # = tcrossprod(R0, R1) / dim_T
  S10 = t(S01)  # = M10 - M12 %*% M22inv %*% M20  # = tcrossprod(R1, R0) / dim_T
  S11 = M11 - M12 %*% M22inv %*% M21  # = tcrossprod(R1, R1) / dim_T
  S00inv = solve(S00)
  
  # decompose S11 for solving the generalized eigenvalue problem, from Pfaff 2008:81 / Johansen 1995:95
  Ctemp <- chol(S11, pivot = TRUE)
  pivot <- attr(Ctemp, "pivot")
  oo <- order(pivot)
  C <- t(Ctemp[, oo])
  Cinv <- solve(C)
  
  # eigenvectors as cointegrating vectors, from Johansen 1995:92,95
  cc_eigen = eigen(Cinv %*% S10 %*% S00inv %*% S01 %*% t(Cinv), symmetric=TRUE) # eigenvalue problem, from Pfaff 2008:81, Eq.4.13 / Johansen 1995:95
  lambda = cc_eigen$values  # squared canonical correlation coefficients, see Juselius 2007:132
  e = cc_eigen$vector  # orthonormal eigenvectors, see Luetkepohl 2005:295 ### unit vectors with length of sqrt(colSums(e^2)) = sqrt(rowSums(e^2)) = (1,...,1) 
  V = t(Cinv) %*% e  # retransformed eigenvectors as cointegrating vectors ### already normalized to V'S_{11}V=I due to mathematical convention of e'Ie=I
  
  # return result
  result = list(lambda=lambda, V=V, 
                S00=S00, S01=S01, S10=S10, S11=S11, S00inv=S00inv,
                M02=M02, M22inv=M22inv, M12=M12, 
                R0=R0, R1=R1, Z0=Z0, Z1=Z1, Z2=Z2)
  return(result)
}


# check and stack deterministic regressors for additive and innovative formulation
aux_stackDtr <- function(type_SL, t_D=NULL, dim_p, dim_T){
  # deterministic regressors for prior ML-RR-Regression (innovative) #
  # deterministic term:    --restricted--   ; --unrestricted-- , see Luetkepohl et al 2004:650
  if(type_SL=="SL_mean" ){ type_D1 = "const"; type_D2 = "none"  }else
  if(type_SL=="SL_ortho"){ stop("Specification 'SL_ortho' not implemented.") }else
  if(type_SL=="SL_trend"){ type_D1 = "trend"; type_D2 = "const" }else{
    stop("Incorrect specification of the deterministic case!") }
  
  # lagged shift dummies for t_break are passed to the impulse dummies instead, see Trenkler et al 2008:335
  t2_impulse = c(t_D$t_break, t_D$t_shift)
  if(!is.null(t2_impulse)){
    t2_impulse = sapply(t2_impulse, FUN=function(t) t + 0:(dim_p-1))  # sacrifice degrees-of-freedom against non-linear LS estimation
    t2_impulse = unique(sort( c(t2_impulse, t_D$t_impulse) ))
  }else{
    t2_impulse = unique(sort( t_D$t_impulse ))
  }
  
  # unrestricted D2 already contains effects from restricted D1
  t1_impulse = t_D$t_impulse
  if(!is.null(t1_impulse)){
    idx_mcol   = t1_impulse %in% t2_impulse  # detect and ...
    t1_impulse = t1_impulse[!idx_mcol]  # prevent perfect multicollinearity
  }
  
  t1_shift = t_D$t_shift
  t2_shift = t_D$t_break
  if(!is.null(t1_shift)){
    idx_mcol = t1_shift %in% t2_shift  # detect and ...
    t1_shift = t1_shift[!idx_mcol]  # prevent perfect multicollinearity
  }
  
  D1 = aux_dummy(dim_T=dim_T, type=type_D1, t_break=t_D$t_break, t_shift=t1_shift, t_impulse=t1_impulse)
  D2 = aux_dummy(dim_T=dim_T, type=type_D2, t_shift=t2_shift, t_impulse=t2_impulse, t_blip=t_D$t_blip, n.season=t_D$n.season, center=TRUE)
  ### see Saikkonen,Luetkepohl 2001:6, Eq.2.11 / Luetkepohl et al 2004:650 / Trenkler et al 2008:335 with prior-detrending
  ### and Johansen et al 2001:219 / Johansen 2016:7, Eq.17 for generalization of the additive and innovative representation.
  
  # deterministic regressors for GLS-detrending (additive)
  type    = switch(type_SL, "SL_mean"="const", "SL_ortho"="both",  "SL_trend"="both")
  t_shift = c(t_D$t_break, t_D$t_shift)
  if(!is.null(t_shift)){
    t_shift = unique(sort( t_shift ))  # prevent duplicated dummies, which would cause perfect multicollinearity
  }
  D = aux_dummy(dim_T=dim_T, t_break=t_D$t_break, t_shift=t_shift, t_impulse=t_D$t_impulse, n.season=t_D$n.season, type=type)
  ### lagged impulse dummies controlling for non-linearities do not enter the VAR in additive det.term formulation, from Trenkler et al 2008:336, Eq.8
  
  # return result
  result = list(D=D, D1=D1, D2=D2)
  return(result)
}


# GLS-detrended series for VECM, from Saikkonen,Luetkepohl (2000) / Trenkler (2008)
aux_GLStrend <- function(y, dim_p, OMEGA, A, D, alpha=NULL, alpha_oc=NULL){
  # define
  y = aux_asDataMatrix(y, "y")  # named matrix y is KxT regardless of input
  D = aux_asDataMatrix(D, "d")  # named matrix D is nxT regardless of input
  dim_T = ncol(y)  # number of total observations incl. presample
  dim_K = nrow(y)  # number of endogenous variables
  dim_n = nrow(D)  # number of deterministic regressors
  
  OMEGAinv = solve(OMEGA)  # inverse covariance matrix of residuals
  idx_slp  = ncol(A) - (dim_K*dim_p - 1):0  # column indices of the slope coefficients in matrix A
  A = A[ , idx_slp, drop=FALSE]  # slope coefficients of the VAR in levels
  
  # remove dynamic effects A(L)
  Dz = kronecker(t(D), diag(dim_K))
  Yz = y
  for(j in 1:dim_p){ 
    Yj = cbind(matrix(0, ncol=j, nrow=dim_K), y[ ,1:(dim_T-j), drop=FALSE]) # initial value assumption x_t=0, from Saikkonen,Luetkepohl 2000:438, Eq.3.1
    Dj = cbind(matrix(0, ncol=j, nrow=dim_n), D[ ,1:(dim_T-j), drop=FALSE])
    Aj = A[ ,(j-1)*dim_K + 1:dim_K]
    Dz = Dz - kronecker(t(Dj), Aj)
    Yz = Yz - Aj %*% Yj
  }
  
  if(is.null(alpha)){
    # OPTION (1): GLS estimation of the original model, from Saikkonen,Luetkepohl 2000:438, Eq.3.1 / Hubrich et al 2001:265
    Q  = NULL  # no transformation matrix needed
    OD = crossprod(kronecker(diag(dim_T), OMEGAinv), Dz)    # Oz = kronecker(diag(dim_T), OMEGAinv)
    MU = solve(crossprod(OD, Dz)) %*% crossprod(OD, c(Yz))  # MU = solve(t(Dz) %*% Oz %*% Dz) %*% (t(Dz) %*% Oz %*% c(Yz))
  
  }else{
    # matrix for GLS-transformation, from Saikkonen,Luetkepohl 2000:438, Eq.3.2
    dim_r = ncol(alpha) # cointegration rank
    if(dim_r==0){ Q1_left = NULL }else{
      Q1_eigen = eigen( t(alpha) %*% OMEGAinv %*% alpha )  # eigendecomposition of a real symmetric matrix, whose eigenvectors are orthogonal
      Q1_eival = diag(sqrt(Q1_eigen$value)^-1, nrow=dim_r) # square root and inverse via eigendecomposition
      Q1_left  = OMEGAinv %*% alpha %*% (Q1_eigen$vector %*% Q1_eival %*% t(Q1_eigen$vector)) 
    }
    if(dim_K-dim_r==0){ Q2_right = NULL }else{
      Q2_eigen = eigen( t(alpha_oc) %*% OMEGA %*% alpha_oc )
      Q2_eival = diag(sqrt(Q2_eigen$value)^-1, nrow=dim_K-dim_r)
      Q2_right = alpha_oc %*% (Q2_eigen$vector %*% Q2_eival %*% t(Q2_eigen$vector))
    }
    Q = cbind(Q1_left, Q2_right) # transformation matrix
    
    # OPTION (2): OLS estimation of the transformed model, from Saikkonen,Luetkepohl 2000:439, Eq.3.4
    QD = kronecker(diag(dim_T), t(Q)) %*% Dz  # (TK x nK) regressor matrix transformed by Q'A(L)
    QY = c(t(Q) %*% Yz)  # vectorized (TK x 1) regressand matrix transformed by Q'A(L)
    MU = solve(crossprod(QD)) %*% crossprod(QD, QY)
  }
  
  # detrended series
  MU = matrix(MU, nrow=dim_K, ncol=dim_n, byrow=FALSE, dimnames=list(rownames(A), rownames(D))) # coefficients for the deterministic term in VMA
  x  = y - MU%*%D
  ### TODO direct result on Z0, Z1, Z2
  ### TODO: if(SL_ort) different adjustments on series Z0, Z1 and Z2, see Trenkler 2008:23, Eq.2.8
  
  # return result
  result = list(MU=MU, x=x, Q=Q)
  return(result)
}


# moments of asymptotic LR-test distribution via response surface approximation
aux_CointMoments <- function(dim_K, dim_L=0, dim_T=NULL, r_H0=0:(dim_K-1), t_D1=NULL, type="", rs_coef=coint_rscoef[[type]]){
  # define
  is_SL = (type %in% c("SL_mean", "SL_orth", "SL_trend","TSL_trend"))
  t_sub = if(is_SL){ t_D1$t_break }else{ c(t_D1$t_break, t_D1$t_shift) }  # TSL (2008) is only necessary for trend breaks.
  t_sub = unique(t_sub)    # time periods of structural breaks
  dim_q = length(t_sub)+1  # number of sub-samples (separated by the breaks)
  d_H0  = dim_K+dim_L - r_H0  # vector for the numbers of stochastic trends
  
  if(dim_q == 1){
    # basic trace and max.eigenvalue test, from Doornik (1998) / Trenkler (2008) #
    rs_term = cbind(d_H0^2, d_H0, sqrt(d_H0), 1, ifelse(d_H0==1, 1, 0), ifelse(d_H0==2, 1, 0))
    moments = rs_term %*% rs_coef  # first and second moments
    
    if(dim_L != 0){
      # trace test correction for conditional VECM: Doornik 1998:586, Eq.11 / Harbo et al. (1998) #
      rs_cond = (dim_K-r_H0) / d_H0
      cov_TR  = dim_L * (dim_K-r_H0) * c(Case1=-1.27, Case2=-1.066, Case3=NA, Case4=-1.35, Case5=NA)[type]
      moments = moments * cbind(rs_cond, rs_cond, NA, NA) - cbind(0, cov_TR, NA, NA)  # No coefficients for ME_EZ and ME_VZ available.
    }
    
  }else if(dim_q <= 3){
    # trace test with up to two breaks resp. three sub-samples #
    zeros = rep(0, 4-dim_q)
    ### TSL (2008:349) use \tau = [T \lambda] with \lambda as the relative break point
    ### ... resp. (T-t_sub)/T as in their example with (184-125)/184 = 0.321. 
    ### Using '+is_SL', the replications with pvars are closer to the original results.
    vs = sort(c(zeros, t_sub-1+is_SL, dim_T))
    ls = sort(diff(vs)) / dim_T
    l1 = ls[1]  # smallest and second smallest relative sub-sample length ...
    l2 = ls[2]  # l1=l2=0 if there is no break resp. a single sub-sample, see Johansen et al. 2000:226 / Trenkler et al. 2008:350
    
    if(type %in% "TSL_trend"){
      # Trenkler et al. 2008:349, Tab.A1 #
      rs_term = cbind(1, d_H0, l1, l2, d_H0^2, d_H0*l1, d_H0*l2, l1^2, l1*l2, l2^2, 
                      d_H0^3, d_H0^2*l1, d_H0^2*l2, d_H0*l1^2, d_H0*l1*l2, d_H0*l2^2, l1^3, l1^2*l2, l1*l2^2, l2^3,
                      1/d_H0,   l1/d_H0,   l2/d_H0,   l1^2/d_H0, (l1*l2)/d_H0,  l2^2/d_H0,  l1^3/d_H0, l1^2*l2/d_H0, l1*l2^2/d_H0, l2^3/d_H0,
                      1/d_H0^2, l1/d_H0^2, l2/d_H0^2, l1^2/d_H0^2, l2^2/d_H0^2, l1^3/d_H0^2, (l1^2*l2)/d_H0^2, (l1*l2^2)/d_H0^2, l2^3/d_H0^2)
      moments = exp(rs_term %*% rs_coef)
      moments = cbind(moments, ME_EZ=NA, ME_VZ=NA)
      
    }else if(type %in% c("JMN_Case2", "JMN_Case4")){
      # Johansen et al. 2000:229, Tab.4 #
      rs_subs = (3-dim_q) * d_H0
      rs_term = cbind(1, d_H0, l1, l2, d_H0^2, d_H0*l1, d_H0*l2, l1^2, l1*l2, l2^2,
                      d_H0^3, d_H0*l1^2, l1^3, l1*l2^2, l1^2*l2, l2^3,
                      1/d_H0, l1/d_H0, l2/d_H0, l1^2/d_H0, l1*l2/d_H0, l2^2/d_H0, l1^3/d_H0, l1*l2^2/d_H0, l2^3/d_H0,
                      1/d_H0^2, l2/d_H0^2, l1^2/d_H0^2, l2^2/d_H0^2, l1^3/d_H0^2, l2^3/d_H0^3)
      moments = exp(rs_term %*% rs_coef) - cbind(rs_subs, 2*rs_subs)  # from Johansen et al. 2000:228,Eq.3.12/3.13
      moments = cbind(moments, ME_EZ=NA, ME_VZ=NA)
    
    }else if(type %in% c("KN_Case2", "KN_Case4")){
      # Kurita,Nielsen 2019:21/22, Tab.A1/A2 #
      rs_cond = (dim_K-r_H0) / d_H0
      rs_subs = (3-dim_q) * (dim_K-r_H0)
      rs_term = cbind(1, d_H0, d_H0^2, d_H0^3, 1/d_H0, 1/d_H0^2, 1/d_H0^3, 
                      l1, l1^2, l1^3, l2, l2^2, l2^3, l1*l2,  l1^2*l2, l1*l2^2, 
                      d_H0*l1, d_H0*l1^2, d_H0*l2, d_H0*l2^2, d_H0*l1*l2, d_H0^2*l2, 
                      l1/d_H0, l2/d_H0, l1^2/d_H0, l1*l2/d_H0, l2^2/d_H0, l1^3/d_H0, (l1*l2^2)/d_H0, l2^3/d_H0,
                      l1/d_H0^2, l2/d_H0^2, l1^2/d_H0^2, l1*l2/d_H0^2, l2^2/d_H0^2,  l1^3/d_H0^2,  (l1^2*l2)/d_H0^2,  (l1*l2^2)/d_H0^2,  l2^3/d_H0^2,
                      ifelse(d_H0==2,  1, 0), ifelse(d_H0==4,    1, 0), ifelse(d_H0==1, d_H0, 0), ifelse(d_H0==3, d_H0, 0), ifelse(d_H0==2, d_H0^3, 0),
                      ifelse(d_H0==1, l1, 0), ifelse(d_H0==1, l1^2, 0), ifelse(d_H0==1, l1^3, 0), ifelse(d_H0==2, l1^3, 0),  
                      ifelse(d_H0==1, l2, 0), ifelse(d_H0==1, l2^2, 0), ifelse(d_H0==2, l2^2, 0), ifelse(d_H0==3, l2^2, 0), 
                      ifelse(d_H0==1, l2^3, 0), ifelse(d_H0==2, l2^3, 0), l1*l2^2,
                      ifelse(d_H0==3,   d_H0*l1, 0), ifelse(d_H0==3,   d_H0*l2, 0), ifelse(d_H0==2, d_H0*l1^2, 0), ifelse(d_H0==2, d_H0*l2^2, 0),
                      ifelse(d_H0==2, d_H0^2*l1, 0), ifelse(d_H0==3, d_H0^2*l1, 0), ifelse(d_H0==2, d_H0^2*l2, 0))
      params  = exp(rs_term %*% rs_coef[ , c("TR_lam", "TR_del")])
      cov_TR  = dim_L*(dim_K-r_H0) * (rs_term %*% rs_coef[ , "TR_cov"])
      moments = cbind(TR_EZ=params[ , "TR_del"], TR_VZ=params[ , "TR_del"]^2) * params[ , "TR_lam"]
      moments = moments*cbind(rs_cond, rs_cond) - cbind(0, cov_TR) - cbind(rs_subs, 2*rs_subs)  # from  Kurita,Nielsen 2019:17,Eq.33/34
      moments = cbind(moments, ME_EZ=NA, ME_VZ=NA)
      
    }else{
      warning("p-values are NA because response surface coefficients are not available for these specifications.")
      moments = matrix(NA, nrow=length(d_H0), ncol=4, dimnames=list(NULL, c("TR_EZ", "TR_VZ", "ME_EZ", "ME_VZ")))
    }
    
  }else{
    warning("p-values are NA because JMN (2000), TSL (2008), resp. KN (2019) 
    provide response surface coefficients  for up to 
    two structural breaks resp. three sub-samples only.")
    moments = matrix(NA, nrow=length(d_H0), ncol=4, dimnames=list(NULL, c("TR_EZ", "TR_VZ", "ME_EZ", "ME_VZ")))
  }
  
  # return result
  rownames(moments) = paste0("d_", d_H0)
  return(moments)
}


# cointegration-rank by nested LR-test procedure
aux_LRrank <- function(lambda, dim_T, dim_K, r_H0=0:(dim_K-1), moments){
  # LR-test statistics
  LR.stats_ME = sapply(r_H0, FUN=function(k) -dim_T*log(1 - lambda[k+1]))  # maximum eigenvalue test, from Johansen 1995:92, Eq.6.18
  LR.stats_TR = sapply(r_H0, FUN=function(k) -dim_T*sum(log(1 - lambda[(k+1):dim_K])))  # trace test, from Johansen 1995:92, Eq.6.14
  LR.stats = cbind(TR=LR.stats_TR, ME=LR.stats_ME)
  
  # approximated p-values from Gamma distribution, from Doornik 1998:576, Eq.4 / Trenkler 2008:26, Eq.3.1
  m = moments[ ,c(1,3), drop=FALSE]  # means for TR and ME
  v = moments[ ,c(2,4), drop=FALSE]  # variances for TR and ME
  LR.pvals = 1-pgamma(LR.stats, shape=m^2/v, rate=m/v)
  
  # return result
  result = data.frame(r_H0, LR.stats, LR.pvals, row.names=NULL)
  colnames(result) = c("r_H0", "stats_TR", "stats_ME", "pvals_TR", "pvals_ME")
  return(result)
}
### For small samples, consider also: ###
# LR.stats_TRs = LRrank$stats_TR * (dim_T - dim_p*dim_K)/dim_T  # small-sample correction by Reinsel,Ahn (1992), from Hubrich et al. 2001:253
# LR.pvals_TRs = tsDyn:::gamma_doornik_all(LRrank$stats_TR, nmp=dim_K:1, test=type_Doornik, type="trace", smallSamp=TRUE, T=adjT) # gamma-approximation for small samples


# cointegration-rank by nested LM-test procedure, from Saikkonen (1999) / Breitung (2005)
aux_LMrank <- function(R0, R1, Z1=R1, L.alpha_oc, L.beta_oc, moments){
  # define
  dim_T = ncol(R0)  # number of observations without presample 
  dim_K = nrow(R0)  # number of endogenous variables
  idx_k = 1:dim_K   # index for the endogenous variables
  r_H0  = sapply(L.alpha_oc, FUN=function(x) dim_K-ncol(x))  # cointegration ranks under H0
  ### For respecting weakly exogenous variables, see LM-test construction in Luetkepohl 2005:696. 
  
  # restricted deterministic term in auxiliary regression, from Breitung 2005:160
  plus = Z1[-idx_k, , drop=FALSE]  # using Z1=R1 (the default), short-run effects Z2 have been concentrated out from deterministic term D1
  
  # LM-test statistic for each cointegration rank
  LM.stats = NULL
  for(r in 1:length(L.alpha_oc)){
    # variables for auxiliary regression, from Breitung 2005:158, Eq.13 / 2005:160, Eq.18
    U = t(L.alpha_oc[[r]]) %*% R0 
    W = rbind((t(L.beta_oc[[r]]) %*% R1[idx_k, , drop=FALSE]), plus)
    
    # LM-statistic, from Breitung 2005:158, Eq.13(II) / Saikkonen 1999:198, Eq.14
    UW     = tcrossprod(U, W)  # = U%*%t(W) = t(W%*%t(U))
    WWinv  = solve(tcrossprod(W))
    UUinv  = solve(tcrossprod(U))
    lambda = dim_T * sum(diag( UW %*% WWinv %*% t(UW) %*% UUinv ))
    
    LM.stats = c(LM.stats, lambda)
    rm(U, W, UW, WWinv, UUinv, lambda)
  }
  
  # approximated p-values, see Doornik 1998:578 for equivalence of simulated LM- and LR-trace test distribution
  m = moments[ , 1, drop=FALSE]  # means for TR resp. LM
  v = moments[ , 2, drop=FALSE]  # variances for TR resp. LM
  LM.pvals = 1-pgamma(LM.stats, shape=m^2/v, rate=m/v)  # p-values from Gamma distribution, from Doornik 1998:576, Eq.4
  
  # return result
  LM.rank = data.frame(r_H0=r_H0, stats_LM=LM.stats, pvals_LM=LM.pvals, row.names=NULL)
  return(LM.rank)
}


# orthogonal complement, upgraded from MASS::Null()
aux_oc <- function(M){
  tmp <- qr(M)
  set <- if(tmp$rank == 0L) seq_len(max(dim(M))) else -seq_len(tmp$rank)
  qr.Q(tmp, complete = TRUE)[, set, drop = FALSE]
}


# test restriction of weak exogeneity, from Johansen 1995:126
aux_LRrest <- function(lambda, lambda_rest, dim_r, dim_T, dim_L){
  ### see also Johansen 1990:199 / Juselius 2007:195
  # define
  #H # restriction matrix #### TODO test for weak exogeneity of each variable k=1,..,K 
  #dim_r = ncol(alpha) # cointegration rank
  #dim_L = ncol(H)     # number of weakly exogenous variables 
  #dim_T  # number of observations without presample
  
  ### TODO: aux_LRrest() test restrictions on 
  ### cointegrating vectors (some known, different restrictions), see Luetkepohl,Kraetzig 2004:105
  ### deterministic case / breaks, see Johansen et al. 2000:230, Ch.4
  ### beta and alpha jointly
  
  # LR-test statistics and p-values
  LR.stats = dim_T*sum(log((1 - lambda_rest[1:dim_r])/(1 - lambda[1:dim_r])))  # from Johansen 1995:126, Eq.8.17
  LR.pvals = 1-pchisq(LR.stats, df=dim_r*dim_L)
  
  # return result
  result = data.frame(stats=LR.stats, pvals=LR.pvals)
  return(result)
}


# restrict cointegration space to cointegrating vectors
aux_beta <- function(V, dim_r, normalize=FALSE){
  ### see also Johansen 1991:1556, Ch.3
  if(dim_r==0){
    beta = V[ ,0, drop=FALSE] # matrix(0, ncol=ncol(V), nrow=nrow(V))
    return(beta)
  }else if(normalize==FALSE){  # preserve normalization such that \beta'S_{11}\beta=I
    beta = V[ ,1:dim_r, drop=FALSE]  # cointegrating matrix with restricted rank, from Johansen 1995:93, Eq.6.16
  }else if(normalize=="first element"){  # normalize each eigenvector to its first element, see Juselius 2007:122
    beta = apply(V[ ,1:dim_r, drop=FALSE], MARGIN=2, FUN=function(v_k) v_k/v_k[1])
  }else if(normalize=="natural"){  # "natural" normalization, from Johansen 1995:179, Eq.13.3(I)
    beta = V[ ,1:dim_r, drop=FALSE] %*% solve(V[1:dim_r, 1:dim_r]) # restrict to \beta = [I_r:\beta_1']' with \beta_1 being (K-r x r)
  }
  ### TODO: argument for normalization where a variable is I(0), i.e. cointegrated on its own

  # return result
  colnames(beta) = paste0("ect.", 1:dim_r)
  return(beta) 
}



#######################################
###  RESPONSE SURFACE COEFFICIENTS  ###
#######################################

# response surface for moments of asymptotic LR test distributions under r_{H0 \perp} = d
coint_TSL = list(
  # coefficients from Trenkler et al 2008:349, Tab.AI for log-moments of TR test statistics' distribution 
  ###### 1(const.)      d       l1       l2      d^2     d*l1     d*l2     l1^2    l1*l2     l2^2
  TR_EZa=c(2.4402, 0.5664,  1.6881, -0.1674, -0.0367, -0.1265,  0.0286, -7.2613, -1.9837, -1.6794),
  TR_VZa=c(2.2377, 0.6725, -1.8646,  1.5842, -0.0440,  0.0000, -0.2485, 12.0954,  5.0822, -1.5583),
  ######      d^3  d^2*l1   d^2*l2   d*l1^2  d*l1*l2  d*l2^2      l1^3  l1^2*l2  l1*l2^2     l2^3
  TR_EZb=c(0.0012, 0.0044, -0.0014,  0.1830,  0.0293, 0.0303,  11.8030, -2.4871,  4.0200,  2.1430),
  TR_VZb=c(0.0013, 0.0105,  0.0135, -0.4765, -0.2405, 0.0898, -22.1045,  7.7659, -8.7651, -0.3356),
  ######       1/d     l1/d     l2/d    l1^2/d (l1*l2)/d   l2^2/d     l1^3/d  l1^2*l2/d  l1*l2^2/d  l2^3/d
  TR_EZc=c(-3.0135,  1.1124,  5.1272,   4.3452,   3.5022, -8.6823,  -16.7672,    5.9728,   -7.0978, 5.7110),
  TR_VZc=c(-1.6753, 11.7097, -1.8672, -60.2299, -10.1422,  4.5029,  129.7558,  -58.2770,   32.3138, 0.0000),
  ######    1/d^2   l1/d^2   l2/d^2  l1^2/d^2  l2^2/d^2  l1^3/d^2 (l1^2*l2)/d^2  (l1*l2^2)/d^2  l2^3/d^2
  TR_EZd=c(1.0331, -0.6479, -2.9655,   0.0000,   7.6083,   5.7696,      -6.5948,        0.0000,  -6.9392),
  TR_VZd=c(0.2956, -4.9776,  4.3265,  30.9656, -14.4186, -82.5994,      48.3167,      -15.3335,  10.8817)
)


# response surface for moments of asymptotic LR test distributions under r_{H0 \perp} = d
coint_JMN = list(
  # coefficients from Johansen et al. 2000:229, Tab.4 for log-moments of TR test statistics' distribution
  ####### 1(const.)      d      l1     l2       d^2     d*l1     d*l2     l1^2    l1*l2     l2^2
  TR_EZc1=c(2.8000, 0.5010, 1.4300, 0.399, -0.03090, -0.0600,  0.0000, -5.7200, -1.1200, -1.7000),
  TR_VZc1=c(3.7800, 0.3460, 0.8590, 0.000, -0.01060, -0.0339,  0.0000, -2.3500,  0.0000,  0.0000),
  TR_EZl1=c(3.0600, 0.4560, 1.4700, 0.993, -0.02690, -0.0363, -0.0195, -4.2100,  0.0000, -2.3500),
  TR_VZl1=c(3.9700, 0.3140, 1.7900, 0.256, -0.00898, -0.0688,  0.0000, -4.0800,  0.0000,  0.0000),
  #######        d^3 d*l1^2  l1^3 l1*l2^2 l1^2*l2    l2^3    
  TR_EZc2=c(0.000974, 0.168, 6.34,   1.89,   0.00,  1.850),
  TR_VZc2=c(0.000000, 0.000, 3.95,   0.00,   0.00, -0.282),
  TR_EZl2=c(0.000840, 0.000, 6.01,   0.00,  -1.33,  2.040),
  TR_VZl2=c(0.000000, 0.000, 4.75,   0.00,   0.00, -0.587),
  #######     1/d    l1/d  l2/d l1^2/d  l1*l2/d l2^2/d  l1^3/d  l1*l2^2/d  l2^3/d
  TR_EZc3=c(-2.19, -0.438, 1.79,  6.03,    3.08, -1.97,  -8.08,     -5.79,   0.00),
  TR_VZc3=c(-2.73,  0.874, 2.36, -2.88,    0.00, -4.44,   0.00,      0.00,   4.31),
  TR_EZl3=c(-2.05, -0.304, 1.06,  9.35,    3.82,  2.12, -22.80,     -7.15,  -4.95),
  TR_VZl3=c(-2.47,  1.620, 3.13, -4.52,   -1.21, -5.87,   0.00,      0.00,   4.89),
  #######   1/d^2  l2/d^2  l1^2/d^2 l2^2/d^2  l1^3/d^2  l2^3/d^3
  TR_EZc4=c(0.717, -1.290,    -1.52,    2.87,      0.0,    -2.03),
  TR_VZc4=c(1.020, -0.807,     0.00,    0.00,      0.0,     0.00),
  TR_EZl4=c(0.681, -0.828,    -5.43,    0.00,     13.1,     1.50),
  TR_VZl4=c(0.874, -0.865,     0.00,    0.00,      0.0,     0.00)
)


# response surface for moments of asymptotic LR test distributions under r_{H0 \perp} = d
coint_KN = list(
  # coefficients from Kurita,Nielsen 2019:21/22, Tab.A1/A2 for moments of TR test statistics' distribution
  #######   1(const.)        d     d0^2     d_H0^3       1/d    1/d^2     1/d^3
  TR_lamc1=c( 4.95486,  0.0000, 0.01738, -0.000840, -9.2630, 9.16200, -3.66200),
  TR_delc1=c( 0.44720,  0.0000, 0.00000,  0.000000,  0.0000, 1.17564, -1.52940),
  TR_covc1=c(-1.53100,  0.0000, 0.01579, -0.001300,  0.9029, 0.00000,  0.00000),
  
  TR_laml1=c( 4.14000,  0.1700, 0.00000, -0.000124, -6.3010, 5.88420, -2.32576),
  TR_dell1=c( 0.59870, -0.0538, 0.00686, -0.000330,  0.0000, 0.00000,  0.00000),
  TR_covl1=c(-1.29800,  0.0000, 0.00000,  0.000000,  0.0000, 0.00000, -2.02200),
  #######         l1    l1^2      l1^3       l2     l2^2    l2^3    l1*l2 l1^2*l2  l1*l2^2
  TR_lamc2=c( 3.0500, -14.61,   21.560,  0.3315,  -2.419,  3.030,  -4.140,  0.00,   5.560),
  TR_delc2=c( 0.0000,  0.000,   -2.084,  0.8286,   0.000, -0.788,   1.750,  0.00,  -3.698),
  TR_covc2=c( 4.1640,  0.000,  -19.650,  0.0000, -14.150, 17.430, -27.160, 14.03,  42.200),
  
  TR_laml2=c( 2.6165, -7.550,   10.400,  2.5245,  -7.412,  5.851,  -5.323,  0.00,   6.096),
  TR_dell2=c(-1.0390,  5.547,  -10.420, -0.3900,   1.841, -2.553,   2.331,  0.00,  -4.325),
  TR_covl2=c(-8.6890, 59.770, -133.500,  2.2250,  -5.156,  0.000,  24.310,  0.00, -59.050),
  #######       d*l1  d*l1^2     d*l2   d*l2^2  d*l1*l2  d^2*l2
  TR_lamc3=c(-0.1280, 0.3264,  0.0000, 0.02660, 0.1302,  0.0000),
  TR_delc3=c( 0.0000, 0.0000, -0.0646, 0.04051, 0.0000,  0.0000),
  TR_covc3=c( 0.0000, 0.0000,  0.3388, 0.00000, 0.0000, -0.0167),
  
  TR_laml3=c(-0.0572, 0.0000, -0.0971, 0.17900, 0.1610,  0.0000),
  TR_dell3=c( 0.0000, 0.0000,  0.0000, 0.00000, 0.0000,  0.0000),
  TR_covl3=c( 0.0000, 0.0000,  0.0000, 0.00000, 0.0000,  0.0000),
  #######       l1/d     l2/d  l1^2/d  l1*l2/d   l2^2/d   l1^3/d  (l1*l2^2)/d   l2^3/d
  TR_lamc4=c( -5.742,   3.339,  44.20,   9.660,  -4.440,  -81.67,      -15.20,    0.00),
  TR_delc4=c( -4.819,  -3.897,  30.49,  -5.108,   2.273,  -40.90,       13.37,    0.00),
  TR_covc4=c(-77.720, -20.520, 278.70, 313.600, 169.100, -461.70,     -562.90, -221.20),
  
  TR_laml4=c( -8.860,  -4.948,  46.15,  31.850,  26.120,  -86.58,      -50.50,  -28.78),
  TR_dell4=c(  9.905,   1.862, -61.09, -17.090, -11.480,  117.68,       35.19,   18.60),
  TR_covl4=c(-29.550, -66.580,   0.00,   0.000, 255.300,  280.50,      155.30, -240.00),
  #######    l1/d^2  l2/d^2 l1^2/d^2  l1*l2/d^2  l2^2/d^2  l1^3/d^2  (l1^2*l2)/d^2  (l1*l2^2)/d^2   l2^3/d^2 
  TR_lamc5=c( 2.410, -3.440,  -24.23,     0.00,      9.60,   47.34,         0.000,       0.000,       -7.22),
  TR_delc5=c(16.000,  3.795, -110.50,     0.00,      0.00,  184.80,         0.000,      -4.478,        0.00),
  TR_covc5=c(81.640,  0.000, -315.00,  -384.80,   -114.60,  804.00,      -290.000,     860.700,      205.20),
  
  TR_laml5=c( 5.296,  2.386,  -29.03,   -19.46,    -13.42,   62.00,        -5.880,     34.590,        15.93),
  TR_dell5=c(-8.836,  1.033,   66.94,    10.84,      0.00, -140.88,         0.000,    -30.160,       -10.05),
  TR_covl5=c(21.320, 71.680,    0.00,     0.00,   -305.70,    0.00,      -321.100,      0.000,       332.10),
  #######      1(2)    1(4)   d*1(1)  d*1(3)  d^3*1(2)  l1*1(1) l1^2*1(1) l1^3*1(1) l1^3*1(2)  l2*1(1)  l2^2*1(1) l2^2*1(2) l2^2*1(3) l2^3*1(1)  l2^3*1(2)  l1*l2^2
  TR_lamc6=c(0.00000,  0.000, 0.0000, 0.000,  0.00000,   0.000,      0.00,     0.00,     0.00,   0.000,     0.000,    0.00,   0.000,    0.000,     0.000,   0.000),
  TR_delc6=c(0.00000,  0.000, 0.5014, 0.000,  0.00000,  -9.833,     73.02,  -130.20,   -14.06,   0.000,    -5.835,    0.00,   0.000,    4.743,     1.944,   0.000),
  TR_covc6=c(0.00000,  0.000, 0.0000, 0.000, -0.00017,   0.000,      0.00,     0.00,     0.00,   0.000,     0.000,    0.18,   0.000,    0.000,     0.000,   0.000),
  
  TR_laml6=c(0.00000,  0.000, 0.0000, 0.000,  0.00000,   0.000,      0.00,     0.00,     0.00,   0.000,     0.000,    0.00,   0.000,    0.000,     0.000,   0.000),
  TR_dell6=c(0.00000,  0.000, 0.0000, 0.000,  0.00000,   2.107,    -20.63,    45.85,     0.00,  -1.029,     3.511,    0.00,   0.000,    0.000,     0.000,   4.267),
  TR_covl6=c(0.03616, -0.027, 0.0000, 0.038,  0.00000,   0.000,      0.00,     0.00,     0.00,   0.000,     0.000,    0.00,  -0.184,    0.000,     0.000,   0.000),
  #######   d*l1*1(3)  d*l2*1(3)  d*l1^2*1(2) d*l2^2*1(2)  d^2*l1*1(2)  d^2*l1*1(3) d^2*l2*1(2)
  TR_lamc7=c(  0.000,    0.0000,       0.000,      0.000,      0.0000,       0.000,    0.00000),
  TR_delc7=c(  0.000,    0.0000,       3.765,     -0.884,     -0.2472,       0.000,    0.06919),
  TR_covc7=c(  1.337,   -0.0215,       0.000,      0.000,      0.0000,      -0.408,    0.00000),
  
  TR_laml7=c(  0.000,    0.0000,       0.000,      0.000,      0.0000,       0.000,    0.00000),
  TR_dell7=c(  0.000,    0.0000,       0.000,      0.062,      0.0000,       0.000,    0.00000),
  TR_covl7=c(  0.000,    0.0000,       0.000,      0.000,      0.0000,       0.000,    0.00000)
)


# response surface for moments of asymptotic LR test distributions under r_{H0 \perp} = d
coint_rscoef = list(
  # coefficients from Doornik 1998:591, Tab.7 for TR and Tab.8 for ME test statistics' distribution
  # (Mistake in Tab.7 for TR_VZ: The row for "const. 1" is actually given at second position. Subsequent rows are aligned accordingly.)
  ######## d^2    d      sqrt(d)  1(const.)  dum(d=1)  dum(d=2)  
  Case1 = cbind(  # H_z
    TR_EZ=c(2, -1.0000,  0.00000,  0.07000,  0.07000,  0.000000),
    TR_VZ=c(3, -0.3300,  0.00000, -0.55000,  0.00000,  0.000000),
    ME_EZ=c(0,  6.0019, -2.77640, -2.75580,  0.67185,  0.114900),
    ME_VZ=c(0,  1.8806, 14.71400, -15.4990,  1.11360,  0.070508)),
  Case2 = cbind(  # H_c
    TR_EZ=c(2,  2.0100,  0.00000,  0.00000,  0.06000,  0.050000),
    TR_VZ=c(3,  3.6000,  0.00000,  0.75000, -0.40000, -0.300000),
    ME_EZ=c(0,  5.9498, -2.36690,  0.43402,  0.04836,  0.018198),
    ME_VZ=c(0,  2.2231, 12.05800, -7.90640,  0.58592, -0.034324)),
  Case3 = cbind(  # H_lc
    TR_EZ=c(2,  1.0500,  0.00000, -1.55000, -0.50000, -0.230000),
    TR_VZ=c(3,  1.8000,  0.00000,  0.00000, -2.80000, -1.100000),
    ME_EZ=c(0,  5.8271, -1.56660, -1.64870, -1.61180, -0.259490),
    ME_VZ=c(0,  2.0785, 13.07400, -9.78460, -3.36800, -0.245280)),
  Case4 = cbind(  # H_l
    TR_EZ=c(2,  4.0500,  0.00000,  0.50000, -0.23000, -0.070000),
    TR_VZ=c(3,  5.7000,  0.00000,  3.20000, -1.30000, -0.500000),
    ME_EZ=c(0,  5.8658, -1.75520,  2.55950, -0.34443, -0.077991),
    ME_VZ=c(0,  1.9955, 12.84100, -5.54280,  1.24250,  0.419490)),
  Case5 = cbind(  # H_ql
    TR_EZ=c(2,  2.8500,  1.35000, -5.10000, -0.10000, -0.060000),
    TR_VZ=c(3,  4.0000,  0.00000,  0.80000, -5.80000, -2.660000),
    ME_EZ=c(0,  5.6364, -0.21447, -0.90531, -3.51660, -0.479660),
    ME_VZ=c(0,  2.0899, 12.39300, -5.33030, -7.15230, -0.252600)),
  
  # coefficients from Johansen et al. 2000:229, Tab.4 for log-moments of TR test statistics' distribution
  JMN_Case2 = cbind(  # H_c
    TR_EZ=c(coint_JMN$TR_EZc1, coint_JMN$TR_EZc2, coint_JMN$TR_EZc3, coint_JMN$TR_EZc4),
    TR_VZ=c(coint_JMN$TR_VZc1, coint_JMN$TR_VZc2, coint_JMN$TR_VZc3, coint_JMN$TR_VZc4)),
  JMN_Case4 = cbind(  # H_l
    TR_EZ=c(coint_JMN$TR_EZl1, coint_JMN$TR_EZl2, coint_JMN$TR_EZl3, coint_JMN$TR_EZl4),
    TR_VZ=c(coint_JMN$TR_VZl1, coint_JMN$TR_VZl2, coint_JMN$TR_VZl3, coint_JMN$TR_VZl4)),
  
  # coefficients from Kurita,Nielsen 2019:21/22, Tab.A1/A2 for moments of TR test statistics' distribution
  KN_Case2 = cbind(  # H_c
    TR_lam=c(coint_KN$TR_lamc1, coint_KN$TR_lamc2, coint_KN$TR_lamc3, coint_KN$TR_lamc4, coint_KN$TR_lamc5, coint_KN$TR_lamc6, coint_KN$TR_lamc7),
    TR_del=c(coint_KN$TR_delc1, coint_KN$TR_delc2, coint_KN$TR_delc3, coint_KN$TR_delc4, coint_KN$TR_delc5, coint_KN$TR_delc6, coint_KN$TR_delc7),
    TR_cov=c(coint_KN$TR_covc1, coint_KN$TR_covc2, coint_KN$TR_covc3, coint_KN$TR_covc4, coint_KN$TR_covc5, coint_KN$TR_covc6, coint_KN$TR_covc7)),
  KN_Case4 = cbind(  # H_l
    TR_lam=c(coint_KN$TR_laml1, coint_KN$TR_laml2, coint_KN$TR_laml3, coint_KN$TR_laml4, coint_KN$TR_laml5, coint_KN$TR_laml6, coint_KN$TR_laml7),
    TR_del=c(coint_KN$TR_dell1, coint_KN$TR_dell2, coint_KN$TR_dell3, coint_KN$TR_dell4, coint_KN$TR_dell5, coint_KN$TR_dell6, coint_KN$TR_dell7),
    TR_cov=c(coint_KN$TR_covl1, coint_KN$TR_covl2, coint_KN$TR_covl3, coint_KN$TR_covl4, coint_KN$TR_covl5, coint_KN$TR_covl6, coint_KN$TR_covl7)),
  
  # coefficients from Trenkler 2008:30, Tab.6 for TR and Tab.7 for ME test statistics' distribution 
  ########    d^2       d    sqrt(d)  1(const.)  dum(d=1) dum(d=2)  
  SL_trend = cbind(
    TR_EZ=c( 1.9996,  0.0000,  0.0000,   1.0365, -0.3469, -0.1112),
    TR_VZ=c( 2.9715,  0.0000,  0.0000,   1.4089,  0.0000,  0.4297),
    ME_EZ=c(-0.0039,  6.1600, -3.3281,  -0.5071,  0.3725,  0.0850),
    ME_VZ=c(-0.0418,  3.4915,  9.2061,  -8.9114,  0.6652,  0.0000)),
  SL_ortho = cbind(
    TR_EZ=c( 2.0008, -2.0990,  0.4463,   0.0000,  0.0000, -0.0503),
    TR_VZ=c( 3.0152, -3.0099,  2.1117,   0.0000,  0.0000, -0.8004),
    ME_EZ=c( 0.0000,  5.8766, -1.9791,  -4.8042,  0.0000,  0.0000),
    ME_VZ=c( 0.0000,  1.3279, 17.6880,   1.3279,  0.0000,  0.0000)),
  SL_mean = cbind(
    TR_EZ=c( 2.0000, -1.0134,  0.0000,   0.1309,  0.0218,  0.0000),
    TR_VZ=c( 2.9778,  0.0000,  0.0000,  -1.7144,  0.9507,  0.4259),
    ME_EZ=c(-0.0035,  6.1365, -3.2161,  -2.3701,  0.5970,  0.1007),
    ME_VZ=c(-0.0258,  2.6655, 12.4462, -13.6992,  0.8563,  0.0000)),
  
  # coefficients from Trenkler et al 2008:349, Tab.AI for log-moments of TR test statistics' distribution 
  TSL_trend = cbind(
    TR_EZ=c(coint_TSL$TR_EZa, coint_TSL$TR_EZb, coint_TSL$TR_EZc, coint_TSL$TR_EZd),
    TR_VZ=c(coint_TSL$TR_VZa, coint_TSL$TR_VZb, coint_TSL$TR_VZc, coint_TSL$TR_VZd))
)


