

#############################
###  AUXILIARY FUNCTIONS  ###
#############################
#
# The following functions serve as modules nested 
# in the calling functions. Notation within each 
# function mostly corresponds to the cited literature.


# stack regressors of seasonal, transitory blip, impulse and shift/step dummy variables as well as trend breaks
aux_dummy <- function(dim_T, dummies=NULL, type=NULL, t_break=NULL, t_shift=NULL, t_impulse=NULL, t_blip=NULL, n.season=NULL, center=FALSE){
  # check and define
  if(any(dim_T < c(t_break, t_shift, t_impulse, t_blip, n.season))){ 
    stop("A period-specific t_* exceeds sample size T.") }
  idx_t = 1:dim_T
  zeros = rep(0, dim_T)
  
  # given dummies
  if(!is.null(dummies)){
    dummies = aux_asDataMatrix(dummies, "d.d")  # matrix D is nxT regardless of input
    idx_td  = seq.int(to=ncol(dummies), length.out=dim_T)
    dummies = dummies[ , idx_td, drop=FALSE]  
    ### remove potential presample, i.e. last observation is reference period across time series!
  }
  
  # deterministic term
  if(!is.null(type)){
    const = trend = NULL
    if(type %in% c("const", "both")){ const = rep(1, dim_T) }
    if(type %in% c("trend", "both")){ trend = seq(1, dim_T) }
    dummies = rbind(const, trend, dummies)
    ### in vars::VAR and urca::ca.jo, the trend begins in the presample!
  }
  names_d = rownames(dummies)
  
  # trend breaks
  for(t in t_break){
    dum = zeros
    dum[t <= idx_t] = 1:(dim_T-t+1)
    dummies = rbind(dummies, dum)
    names_d = c(names_d, paste0("b.bk", t))
  }
  
  # shift/step dummies
  for(t in t_shift){
    dum = zeros
    dum[t <= idx_t] = 1
    dummies = rbind(dummies, dum)
    names_d = c(names_d, paste0("d.sh", t))
  }
  
  # impulse dummies
  for(t in t_impulse){
    dum = zeros
    dum[t] = 1
    dummies = rbind(dummies, dum)
    names_d = c(names_d, paste0("d.im", t))
  }
  
  # blip dummies (transitory in non-stationary data)
  for(t in t_blip){
    dum = zeros
    dum[t] = 1
    dum[t+1] = -1
    dummies = rbind(dummies, dum)
    names_d = c(names_d, paste0("d.bl", t))
  }
  
  # seasonal dummies
  if(!is.null(n.season)){ if(n.season > 1){
    sea = diag(n.season)[-1, , drop=FALSE]
    if(center){
      sea = sea - 1/n.season
    }
    dum = sea
    for(i in 1:ceiling(dim_T/n.season)){
      dum = cbind(dum, sea)
    }
    dummies = rbind(dummies, dum[ ,1:dim_T, drop=FALSE])
    names_d = c(names_d, paste0("d.se", 2:n.season))
  }}
  
  # return result
  rownames(dummies) = names_d
  return(dummies)
}


# stack regressors for VARX(p,q) model
aux_stack <- function(y=NULL, dim_p=0, x=NULL, dim_q=0, D=NULL, Dt.first=TRUE, Z=NULL){
  ### Note that the last period in matrices x, y and D must be congruent!
  # define
  y = aux_asDataMatrix(y, "y")  # named matrix y is KxT regardless of input
  dim_pq = max(dim_p, dim_q) # overall lag-order
  dim_K  = nrow(y)  # number of endogenous variables
  dim_T  = ncol(y) - dim_pq  # number of observations without presample
  
  # endogenous variables
  y_lag = matrix(NA, nrow=dim_K*dim_p, ncol=dim_T)
  if(dim_p!=0){
    idx_ty = seq.int(to=ncol(y), length.out=dim_T)
    rownames(y_lag) = c(sapply(1:dim_p, FUN=function(j) paste0(rownames(y), ".l", j)))
    for (j in 1:dim_p){
      idx_row = 1:dim_K + dim_K*(j-1)
      y_lag[idx_row, ] = y[ , idx_ty-j, drop=FALSE]
    }
  }
  Z = rbind(Z, y_lag)
  
  # exogenous variables, see Luetkepohl 2005:396, Eq.10.3.3
  if(!is.null(x)){
    x = aux_asDataMatrix(x, "x")  # named matrix x is LxT regardless of input
    dim_L  = nrow(x)  # number of exogenous variables
    idx_tx = seq.int(to=ncol(x), length.out=dim_T)
    x_lag  = matrix(NA, nrow=dim_L*(dim_q+1), ncol=dim_T)
    rownames(x_lag) = c(sapply(0:dim_q, FUN=function(j) paste0(rownames(x), ".l", j)))
    for (j in 0:dim_q){
      idx_row = 1:dim_L + dim_L*j 
      x_lag[idx_row, ] = x[ , idx_tx-j, drop=FALSE]
    }
    Z = rbind(Z, x_lag)
  }
  
  # deterministic regressors
  if(!is.null(D)){
    D = aux_asDataMatrix(D, "d.d")  # named matrix D is nxT regardless of input
    idx_td = seq.int(to=ncol(D), length.out=dim_T)
    D = D[ , idx_td, drop=FALSE]  # remove potential presample, i.e. last observation is reference period across time series!
    Z = if(Dt.first){ rbind(D, Z) }else{ rbind(Z, D) }
  }
  
  # return result
  return(Z)
}


# check, define and stack matrices for OLS, from Luetkepohl 2005:396, Eq.10.3.3
aux_stackOLS <- function(y, dim_p, x=NULL, dim_q=0, type="none", t_D=list(), D=NULL){
  # endogenous variables
  y = aux_asDataMatrix(y, "y")  # named matrix y is KxT regardless of input
  dim_K  = nrow(y)  # number of endogenous variables
  
  dim_pq = max(dim_p, dim_q)  # overall lag-order
  dim_T  = ncol(y) - dim_pq   # number of observations without presample
  idx_t  = (dim_pq+1):ncol(y) # index for excluding time periods t=1,..,p from regressands
  
  # exogenous variables
  if(!is.null(x)){
    x = aux_asDataMatrix(x, "x")  # named matrix x is LxT regardless of input
    dim_L = nrow(x)  # number of exogenous variables
  }else{
    dim_L = 0
  }
  
  # deterministic regressors
  D = aux_dummy(dim_T=ncol(y), type=type, dummies=D, t_break=t_D$t_break, t_shift=t_D$t_shift, 
                t_impulse=t_D$t_impulse, t_blip=t_D$t_blip, n.season=t_D$n.season, center=FALSE)
  
  # stack variables
  Y = y[ , idx_t, drop=FALSE]  # regressand matrix, i.e. without presample  
  Z = aux_stack(y=y, dim_p=dim_p, x=x, dim_q=dim_q, D=D)  # regressors
  
  # return result
  result = list(Y=Y, Z=Z, D=D, y=y, x=x, dim_p=dim_p, dim_q=dim_q, 
                dim_T=dim_T, dim_K=dim_K, dim_L=dim_L, dim_Kpn=nrow(Z))
  return(result)
}


# estimate VARX(p,q) model, from Luetkepohl 2005:396, Eq.10.3.3
aux_VARX <- function(y=NULL, dim_p=0, x=NULL, dim_q=0, type="none", t_D=list(), D=NULL, method="solve"){
  # define
  def = aux_stackOLS(y=y, dim_p=dim_p, x=x, dim_q=dim_q, type=type, t_D=t_D, D=D)
  dim_T   = def$dim_T    # number of observations without presample
  dim_Kpn = def$dim_Kpn  # number of regressors
  Y = def$Y  # regressands
  Z = def$Z  # regressors
  
  # estimate by multivariate OLS
  if(method=="solve"){
    # ... using inverse of product moment matrix
    YZ    = tcrossprod(Y, Z)  # = Y%*%t(Z)
    ZZinv = solve(tcrossprod(Z))  # inverse moment matrix
    A     = YZ %*% ZZinv  # multivariate OLS estimator
    resid = Y - A %*% Z   # residuals
    OMEGA = tcrossprod(resid) / dim_T        # MLE covariance matrix of residuals
    SIGMA = OMEGA * (dim_T/(dim_T-dim_Kpn))  # OLS covariance matrix of residuals
    
  }else if(method=="qr"){
    # ... using QR decomposition with column pivoting
    R.lm  = .lm.fit(x=t(Z), y=t(Y))
    A     = t(R.lm$coefficients);  dimnames(A) = list(rownames(Y), rownames(Z))
    resid = t(R.lm$residuals); dimnames(resid) = dimnames(Y)
    OMEGA = tcrossprod(resid) / dim_T        # MLE covariance matrix of residuals
    SIGMA = OMEGA * (dim_T/(dim_T-dim_Kpn))  # OLS covariance matrix of residuals
    
  }else{ stop("Argument 'method' must be either 'solve' or 'qr' to call .lm.fit().") }
  
  # return result
  result = c(list(A=A, OMEGA=OMEGA, SIGMA=SIGMA, resid=resid), def)
  return(result)
}


# (modified) information criteria of VAR model, from Qu,Perron 2007:651 / Kilian,Luetkepohl 2017:53
aux_MIC <- function(Ucov, COEF, dim_T, lambda_H0=NULL){
  # define
  detCov  = det(Ucov)   # determinant of MLE covariance matrix
  dim_K   = nrow(COEF)  # number of equations in the system
  dim_Kpn = ncol(COEF)  # number of coefficients per equation
  n_coef  = dim_K * dim_Kpn  # number of coefficients in the system
  
  # select ordinary or Qu,Perron's (2007) modified information criteria
  if(is.null(lambda_H0)){
    TR = 0
  }else{
    TR = -dim_T*sum(log(1 - lambda_H0))  # trace-test statistic under r_H0
  }
  
  # calculate information criteria
  FPE = detCov * ((dim_T+dim_Kpn)/(dim_T-dim_Kpn))^dim_K  # (ordinary) Final Prediction Error, from Luetkepohl 2005:147
  AIC = 2           # for Akaike (1973) Information Criterion
  SIC = log(dim_T)  # for Schwarz (1978) Bayesian Information Criterion
  HQC = 2*log(SIC)  # for Hannan,Quinn (1979) Criterion
  MIC = log(detCov) + c(AIC, HQC, SIC) * (TR+n_coef)/dim_T  # modification with TR, from Qu,Perron 2007:651, Eq.10
  
  # return result
  result = c(MIC, FPE)
  names(result) = c("AIC", "HQC", "SIC", "FPE")
  return(result)
}


# estimate rank-restricted VECM
aux_VECM <- function(beta, RRR){
  # define
  M02 = RRR$M02  # product moment matrices 
  M12 = RRR$M12
  M22inv = RRR$M22inv
  S00 = RRR$S00
  S01 = RRR$S01
  S10 = RRR$S10
  S11 = RRR$S11
  
  Z0 = RRR$Z0  # regressand matrix, first-differenced and without presample
  Z1 = RRR$Z1  # regressor matrix in levels for the cointegrating term
  Z2 = RRR$Z2  # regressor matrix in first differences for short-run effects
  
  dim_r   = ncol(beta)  # cointegration rank
  dim_T   = ncol(Z0)  # number of observations without presample
  dim_Kpn = nrow(Z1) + nrow(Z2)  # number of deterministic and endogenous regressors
  
  # respect normalization of the cointegrating vectors
  if(ncol(beta) == 0){ C = diag(dim_r) }else{  # also normalization \beta'S_{11}\beta=I cancels (\beta'S_{11}\beta)^{-1} ...
    C = solve(t(beta) %*% S11 %*% beta) }      # ... out from Johansen's formulas, see Pfaff 2008:82
  
  # calculate coefficient matrices
  # ... from Johansen 1996:93 (after Eq.6.17): For given \hat{\beta}, the remaining estimators are just OLS-estimators.
  alpha = S01 %*% beta %*% C  # matrix of loading weights, from Johansen 1996:91, Eq.6.11
  PI    = alpha %*% t(beta)  # coefficient matrix of the lagged variables in levels, from Johansen 1996:96, Eq. 6.23
  GAMMA = M02 %*% M22inv - PI %*% M12 %*% M22inv  # short-run effects coefficient matrix, PSI from Johansen 1996:90, Eq.6.5
  OMEGA = S00 - PI %*% S10  # MLE covariance matrix of residuals, from Johansen 1996:91, Eq.6.12 with PI = S01 %*% beta %*% C %*% t(beta)
  SIGMA = OMEGA * (dim_T/(dim_T-dim_Kpn))  # OLS covariance matrix of residuals
  resid = Z0 - PI%*%Z1 - GAMMA%*%Z2  # residuals
  
  # return result
  result = list(alpha=alpha, PI=PI, GAMMA=GAMMA, OMEGA=OMEGA, SIGMA=SIGMA, resid=resid)
  return(result)
}


# transform conditional VECM into full-system VECM, from Johansen 1996:122
aux_con2vec <- function(beta, VECM_c, dim_p, VECM_x, dim_q){
  ### see also Jacobs,Wallis 2010:109, Eq.6 / Kurita,Nielsen 2019:4, Eq.3/4
  ### restrictions and stacking operations also used by Banerjee et al 2017:1075, Eq.23/24
  # define
  alpha_y = VECM_c$alpha
  GAMMA_c = VECM_c$GAMMA
  resid_c = VECM_c$resid
  OMEGA_c = VECM_c$OMEGA
  SIGMA_c = VECM_c$SIGMA
  
  GAMMA_xx = VECM_x$GAMMA
  resid_x  = VECM_x$resid
  OMEGA_xx = VECM_x$OMEGA
  SIGMA_xx = VECM_x$SIGMA
  
  names_k = rownames(GAMMA_c)
  names_l = rownames(GAMMA_xx)
  
  dim_r  = ncol(beta)     # cointegration rank
  dim_K  = nrow(resid_c)  # number of endogenous variables
  dim_L  = nrow(resid_x)  # number of weakly exogenous variables
  dim_T  = min(ncol(resid_c), ncol(resid_x))  # number of observations without overall presample
  dim_pq = max(dim_p, dim_q)  # overall lag-order
  dim_n  = ncol(GAMMA_xx)- dim_L*(dim_q-1)  # number of deterministic regressors
  dim_nc = ncol(GAMMA_c) - dim_K*(dim_p-1) - dim_L*dim_q  # for lags from 0 to q-1 of the first-differences
  if(dim_n != dim_nc){ warning("Deterministic regressors must be the same in partial and marginal mode") }
  
  idx_k    = 1:dim_K
  idx_l    = 1:dim_L + dim_K
  idx_dtr  = 0:dim_n  # index for deterministic term in both models
  idx_tx   = 1:dim_T + (dim_pq-dim_q)  # index for periods in the marginal model belonging to the overall total sample
  idx_x.l0 = dim_n + dim_K*(dim_p-1) + 1:dim_L  # index for the coefficients of instantaneous effects
  
  # rearrange short-run and deterministic effects
  GAMMA_x = GAMMA_xx[ , idx_dtr, drop=FALSE]
  GAMMA_y = GAMMA_c[  , idx_dtr, drop=FALSE]
  Null_LK = matrix(0, nrow=dim_L, ncol=dim_K) # Under strong exogeneity, short-run effects of y on x are GAMMA_xyj = 0.
  if(dim_pq > 1){ for(j in 2:dim_pq){  # order according to lag-structure of first-difference regressors
    idx_xxj = dim_n + (j-2)*dim_L + 1:dim_L
    idx_yyj = dim_n + (j-2)*dim_K + 1:dim_K
    idx_yxj = idx_xxj + (dim_p-1)*dim_K + dim_L
    if(j > dim_p){
      GAMMA_xxj = GAMMA_xx[ , idx_xxj, drop=FALSE]
      GAMMA_yxj = GAMMA_c[  , idx_yxj, drop=FALSE]
      GAMMA_yyj = matrix(0, nrow=dim_K, ncol=dim_K)
    }else if(j > dim_q){
      names_yxj = paste0(names_l, ".l", j-1)
      GAMMA_xxj = matrix(0, nrow=dim_L, ncol=dim_L)
      GAMMA_yxj = matrix(0, nrow=dim_K, ncol=dim_L, dimnames=list(NULL, names_yxj))
      GAMMA_yyj = GAMMA_c[ , idx_yyj, drop=FALSE]
    }else{
      GAMMA_xxj = GAMMA_xx[ , idx_xxj, drop=FALSE]
      GAMMA_yxj = GAMMA_c[  , idx_yxj, drop=FALSE]
      GAMMA_yyj = GAMMA_c[  , idx_yyj, drop=FALSE]
    }
    GAMMA_x = cbind(GAMMA_x, Null_LK,   GAMMA_xxj)
    GAMMA_y = cbind(GAMMA_y, GAMMA_yyj, GAMMA_yxj)
  }}
  
  # transform using partial and marginal model
  LAMBDA  = GAMMA_c[ , idx_x.l0, drop=FALSE]  # coefficient matrix (KxL) for instantaneous effects from x to y
  GAMMA_y = GAMMA_y + LAMBDA %*% GAMMA_x  # transform short-run and deterministic effects, see Johansen 1996:122, Th.8.3
  resid_y = resid_c + LAMBDA %*% resid_x[ , idx_tx, drop=FALSE]
  alpha_x = matrix(0, nrow=dim_L, ncol=dim_r, dimnames=list(rownames(resid_x), NULL)) # under weak exogeneity, see Johansen 1996:122, Th.8.1
  # alpha_y = alpha_c + LAMBDA %*% alpha_x
  
  OMEGA_yx = LAMBDA %*% OMEGA_xx
  OMEGA_xy = t(OMEGA_yx)
  OMEGA_yy = OMEGA_c + LAMBDA %*% OMEGA_xy
  
  SIGMA_yx = LAMBDA %*% SIGMA_xx  # Note that, if lag-orders p!=q, SIGMA_yx and SIGMA_xx are corrected for different number of regressors.
  SIGMA_xy = t(SIGMA_yx)
  SIGMA_yy = SIGMA_c + LAMBDA %*% SIGMA_xy
  
  # stack to full-system VECM
  resid = rbind(resid_y, resid_x[ , idx_tx, drop=FALSE])  # residuals
  alpha = rbind(alpha_y, alpha_x)  # loading weights
  PI    = alpha %*% t(beta)        # coefficient matrix of the lagged variables in levels
  GAMMA = rbind(GAMMA_y, GAMMA_x)  # short-run effects coefficient matrix (incl. unrestricted deterministic effects)
  OMEGA = rbind(cbind(OMEGA_yy, OMEGA_yx), cbind(OMEGA_xy, OMEGA_xx))  # MLE covariance matrix of residuals
  SIGMA = rbind(cbind(SIGMA_yy, SIGMA_yx), cbind(SIGMA_xy, SIGMA_xx))  # OLS covariance matrix of residuals 
  
  # stack instantaneous structural parameter matrix, see Lanne,Luetkepohl 2008:1133
  B = diag(dim_K+dim_L)
  dimnames(B) = list(c(names_k, names_l), paste0("u[ ", c(names_k, names_l), " ]"))
  B[idx_k, idx_l] = LAMBDA
  
  # return result
  result = list(alpha=alpha, PI=PI, GAMMA=GAMMA, OMEGA=OMEGA, SIGMA=SIGMA, resid=resid, LAMBDA=LAMBDA)
  return(result)
}


# accumulate coefficient matrices, e.g. A(1) or GAMMA(1)
aux_accum <- function(coef, dim_p){
  # define
  names_k = rownames(coef)
  dim_K   = nrow(coef)  # number of endogenous variables
  dim_Kpn = ncol(coef)  # number of coefficients per equation
  
  # collect slope coefficient matrices in array and accumulate
  if(dim_p==0){ 
    result = diag(dim_K)
  }else{ 
    idx_slp = dim_Kpn - (dim_K*dim_p - 1):0  # column indices of the slope coefficients
    A.coef  = array(coef[ , idx_slp], dim=c(dim_K, dim_K, dim_p))  # array from KxK blocks j=1,...,p
    result  = diag(dim_K) - apply(A.coef, MARGIN=1:2, sum)
  }
  
  # return result
  dimnames(result) = list(names_k, names_k)
  return(result)
}


# TODO: transform VECM into VMA of Granger Representation Theorem
aux_vec2vma <- function(GAMMA, alpha_oc, beta_oc, dim_p, B=diag(dim_K), n.ahead=20){
  ### see also Johansen 1990:49, Eq.4.6 / Juselius 2007:276 / MellanderEtAl (1992)
  # define
  names_k = rownames(GAMMA)
  dim_K   = nrow(GAMMA)  # number of endogenous variables
  GAMMA_1 = aux_accum(GAMMA, dim_p=dim_p-1)  # cumulative GAMMA(1)

  # long-run multiplier matrix, from Luetkepohl 2005:252, Eq.6.3.12
  XI = beta_oc %*% solve( t(alpha_oc) %*% GAMMA_1 %*% beta_oc ) %*% t(alpha_oc)
  dimnames(XI) = list(names_k, NULL)
  if(n.ahead=="XI"){ return(XI) }
  
  # short-run multiplier matrices, from Hansen (2005) / Luetkepohl 2005:252
  XI_star = NA
  dimnames(XI_star) = list(names_k, NULL, NULL)
  ### # optional normalization of 'B', see Jentsch, Lunsford 2021:5
  ### if( !is.null(normf) ){ B = normf(B) }
  ### ?use aux_var2vma() for transforming the transitory dynamics of the GRT?
  
  # multiplier matrix of the deterministic term, from Johansen,Nielsen (2018)
  ### TODO for IRF in consequence of an impulse in a dummy
  
  # structural multiplier matrices, from Hansen 2005:26, Rm.1
  UPSILON = XI %*% B  # long-run
  UP_star = apply(XI_star, MARGIN=3, FUN=function(n) n %*% B)  # short-run
  THETA   = apply(UP_star, MARGIN=3, FUN=function(n) n + UPSILON)  # IRF
  
  # return result
  result = list(XI=XI, THETA=THETA)
  return(result)
}


# transform VECM into rank-restricted VAR in level-representation, from Luetkepohl 2005:289
aux_vec2var <- function(PI, GAMMA, dim_p, type_VECM=NULL){
  # define
  dim_K = nrow(GAMMA)  # number of  endogenous variables
  dim_L = 0  # number of exogenous variables
  names_k = rownames(GAMMA)
  names_slp = c(sapply(1:dim_p, function(j) paste0(names_k, ".l", j)))
  
  if(dim_p==1){
    idx_slp = 0  # column indices of the slope coefficients in GAMMA
    idx_dtr = 0:ncol(GAMMA)  # column indices of the coefficients for the deterministic term in GAMMA
  }else{
    idx_slp = ncol(GAMMA) - (dim_K*(dim_p-1) - 1):0  # column indices of the slope coefficients in GAMMA
    idx_dtr = -idx_slp  # column indices of the coefficients for the deterministic term in GAMMA
  }
  A.GAMMA = array(0, dim=c(dim_K, dim_K+dim_L, dim_p))  # array for GAMMA_j with j=1,...,p incl. zero-matrix GAMMA_p
  A.GAMMA[ , , -dim_p] = GAMMA[ , idx_slp, drop=FALSE]  # block matrices of slope coefficients
  
  # slope coefficients A_j with j=1,...,p
  A = PI[ , 1:dim_K, drop=FALSE] + diag(dim_K) + A.GAMMA[ , , 1]  # coefficients A_1
  if(dim_p > 1){ for(j in 2:dim_p){
    A = cbind(A, A.GAMMA[ , , j] - A.GAMMA[ , , j-1]) }}  # coefficients A_j for j=2,...,p
  colnames(A) = names_slp
  
  # coefficients of the deterministic term
  A = cbind(GAMMA[ ,idx_dtr, drop=FALSE],  # unrestricted (first column vectors in GAMMA)
            PI[ ,-(1:dim_K), drop=FALSE],  # restricted (last column vectors in PI)
            A)                             # add slope coefficients
  type_VAR = if(is.null(type_VECM)){ NULL }else{ switch(type_VECM, "Case1"="none", "Case2"="const", "Case3"="const", "Case4"="both", "Case5"="both") }
  
  # return result
  result = list(A=A, type=type_VAR)
  return(result)
}


# transform conditional VAR into full-system VAR, from Luetkepohl 2005:392
aux_con2var <- function(VAR_c, dim_p, VAR_x, dim_q){
  # define
  A_c     = VAR_c$A
  resid_c = VAR_c$resid
  OMEGA_c = VAR_c$OMEGA
  SIGMA_c = VAR_c$SIGMA
  
  A_xx     = VAR_x$A
  resid_x  = VAR_x$resid
  OMEGA_xx = VAR_x$OMEGA
  SIGMA_xx = VAR_x$SIGMA
  
  names_k = rownames(A_c)
  names_l = rownames(A_xx)
  
  dim_K  = nrow(resid_c)  # number of endogenous variables
  dim_L  = nrow(resid_x)  # number of exogenous variables
  dim_T  = min(ncol(resid_c), ncol(resid_x))  # number of observations without overall presample
  dim_pq = max(dim_p, dim_q)  # overall lag-order
  dim_n  = ncol(A_xx)- dim_L*dim_q  # number of deterministic regressors
  dim_nc = ncol(A_c) - dim_K*dim_p - dim_L*dim_q
  if(dim_n != dim_nc){ warning("Deterministic regressors must be the same in partial and marginal mode") }
  
  idx_k    = 1:dim_K
  idx_l    = 1:dim_L + dim_K
  idx_dtr  = 0:dim_n  # index for deterministic term in both models
  idx_tx   = 1:dim_T + (dim_pq-dim_q)  # index for periods in the marginal model belonging to the overall total sample
  idx_x.l0 = dim_n + dim_K*dim_p + 1:dim_L  # index for the coefficients of instantaneous effects
  
  # rearrange lagged and deterministic effects
  A_x = A_xx[ , idx_dtr, drop=FALSE]
  A_y = A_c[  , idx_dtr, drop=FALSE]
  Null_LK = matrix(0, nrow=dim_L, ncol=dim_K) # Under strong exogeneity, effects of y on x are A_xyj = 0.
  for(j in 2:dim_pq){  # order according to lag-structure of first-difference regressors
    idx_xxj = dim_n + (j-2)*dim_L + 1:dim_L
    idx_yyj = dim_n + (j-2)*dim_K + 1:dim_K
    idx_yxj = idx_xxj + (dim_p-1)*dim_K + dim_L
    if(j > dim_p){
      A_xxj = A_xx[ , idx_xxj, drop=FALSE]
      A_yxj = A_c[  , idx_yxj, drop=FALSE]
      A_yyj = matrix(0, nrow=dim_K, ncol=dim_K)
    }else if(j > dim_q){
      names_yxj = paste0(names_l, ".l", j-1)
      A_xxj = matrix(0, nrow=dim_L, ncol=dim_L)
      A_yxj = matrix(0, nrow=dim_K, ncol=dim_L, dimnames=list(NULL, names_yxj))
      A_yyj = A_c[ , idx_yyj, drop=FALSE]
    }else{
      A_xxj = A_xx[ , idx_xxj, drop=FALSE]
      A_yxj = A_c[  , idx_yxj, drop=FALSE]
      A_yyj = A_c[  , idx_yyj, drop=FALSE]
    }
    A_x = cbind(A_x, Null_LK, A_xxj)
    A_y = cbind(A_y, A_yyj,   A_yxj)
  }
  
  # transform using partial and marginal model
  LAMBDA  = A_c[ , idx_x.l0, drop=FALSE]  # coefficient matrix (KxL) for instantaneous effects from x to y
  A_y     = A_y + LAMBDA %*% A_x  # transform coefficient matrix
  resid_y = resid_c + LAMBDA %*% resid_x[ , idx_tx, drop=FALSE]
  
  OMEGA_yx = LAMBDA %*% OMEGA_xx
  OMEGA_xy = t(OMEGA_yx)
  OMEGA_yy = OMEGA_c + LAMBDA %*% OMEGA_xy
  
  SIGMA_yx = LAMBDA %*% SIGMA_xx  # Note that, if lag-orders p!=q, SIGMA_yx and SIGMA_xx are corrected for different number of regressors.
  SIGMA_xy = t(SIGMA_yx)
  SIGMA_yy = SIGMA_c + LAMBDA %*% SIGMA_xy
  
  # stack to full-system VAR
  A     = rbind(A_y, A_x)  # coefficient matrix of the lagged variables and deterministic effects
  resid = rbind(resid_y, resid_x[ , idx_tx, drop=FALSE])  # residuals  
  ### TODO: Do not trash residuals, but include Null_mat instead? Bootstrap? idx_ty?
  OMEGA = rbind(cbind(OMEGA_yy, OMEGA_yx), cbind(OMEGA_xy, OMEGA_xx))  # MLE covariance matrix of residuals
  SIGMA = rbind(cbind(SIGMA_yy, SIGMA_yx), cbind(SIGMA_xy, SIGMA_xx))  # OLS covariance matrix of residuals 
  
  # stack instantaneous structural parameter matrix, see Lanne,Luetkepohl 2008:1133
  B = diag(dim_K+dim_L)
  dimnames(B) = list(c(names_k, names_l), paste0("u[ ", c(names_k, names_l), " ]"))
  B[idx_k, idx_l] = LAMBDA
  
  # return result
  result = list(A=A, B=B, OMEGA=OMEGA, SIGMA=SIGMA, resid=resid, LAMBDA=LAMBDA)
  return(result)
}


# transform VAR(p) into VAR(1) companion form representation, from Luetkepohl 2005:15, Eq.2.1.8
aux_var2companion <- function(A, dim_p){
  # define
  dim_K   = nrow(A)  # number of endogenous variables
  dim_Kpn = ncol(A)  # number of coefficients per equation
  idx_slp = dim_Kpn - (dim_K*dim_p - 1):0  # column indices of the slope coefficients in matrix A
  
  # construct companion matrix
  result = matrix(0, ncol=dim_K*dim_p, nrow=dim_K*dim_p)
  result[1:dim_K, 1:(dim_K*dim_p)] = A[ , idx_slp, drop=FALSE]
  if(dim_p>1){
    result[(dim_K+1):(dim_K*dim_p), 1:(dim_K*(dim_p-1))] = diag(dim_K*(dim_p-1))}
  
  # return result
  return(result)
}


# transform VAR into (structural) VMA, from Luetkepohl 2005:52 / Kilian,Luetkepohl 2017:109
aux_var2vma <- function(A, B=diag(dim_K), dim_p, n.ahead=20, normf=NULL){
  # optional normalization of 'B', see Jentsch, Lunsford 2021:5
  if( !is.null(normf) ){ B = normf(B) }
  
  # define
  names_k = rownames(A)  # names of endogenous variables
  names_s = colnames(B)  # names of structural shocks
  dim_K = nrow(A)  # number of endogenous variables
  dim_S = ncol(B)  # number of structural shocks
  idx_k = 1:dim_K  # equivalent to selection matrix 'J'
  
  # transform
  A.c   = aux_var2companion(A=A, dim_p=dim_p)  # coefficient matrix for the VAR(1)-companion representation
  PHI   = array(NA, dim=c(dim_K, dim_K, n.ahead+1))
  THETA = array(NA, dim=c(dim_K, dim_S, n.ahead+1))
  A.n   = diag(dim_K*dim_p)  # initial A.c^0
  for(n in (0:n.ahead)+1){
    PHI[ , , n]   = A.n[idx_k, idx_k, drop=FALSE]  # responses to forecast errors
    THETA[ , , n] = PHI[ , , n] %*% B  # responses to structural shocks
    A.n = A.n %*% A.c  # A.c^(n+1) for the subsequent run in this loop
  }
  
  # return result
  dimnames(PHI)   = list(names_k, names_k, NULL)
  dimnames(THETA) = list(names_k, names_s, NULL)
  result = list(PHI=PHI, THETA=THETA)
  return(result)
}


# bias-correction in bootstrap-after-bootstrap, from Kilian (1998)
aux_BaB <- function(A_hat, dim_p, PSI, deltas=cumprod((100:0)/100)){
  ### These geometric deltas correspond to the recursive loop in Kilian 1998:220, Step 1b.
  m_hat = Mod(eigen(aux_var2companion(A=A_hat, dim_p=dim_p))$values[1])
  if( m_hat >= 1 ){
    # no bias-correction under non-stationarity
    A_bc = A_hat
  }else{
    # keep bias-corrected A within the range of stationarity
    for(delta_i in deltas){
      A_bc = A_hat - delta_i * PSI  # bias-corrected coefficients, from Kilian 1998:220
      m_bc = Mod(eigen(aux_var2companion(A=A_bc, dim_p=dim_p))$values[1])
      if( m_bc < 1 ){ break }
  }}
  
  # return result
  return(A_bc)
}


