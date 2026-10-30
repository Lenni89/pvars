

###############################
###  COMBINATION FUNCTIONS  ###
###############################
#
# Combination approaches for
# the KxN individual tests.


# panel tests based on the standardized mean of individual test statistics, see Larsson et al. (2001) and Im et al. 2003:59
ptest.STATSbar <- function(STATS, distribution="theoretical", moments){
  # define
  STATS = as.matrix(STATS)  # matrix of individual test statistics
  if(ncol(STATS)==1){ STATS = t(STATS) }  # ... in case there is just one hypothesis test
  dim_K = nrow(STATS)  # number of endogenous variables resp. of panel tests
  dim_N = ncol(STATS)  # number of individuals
  
  # test statistics
  EZ = moments[ , 1]  # expectation values for statistics' distribution Z_d with d = K-r_H0 = K,...,1 for r_H0 = 0,...,K-1
  VZ = moments[ , 2]  # variances for Z_d
  STATSbar = rowMeans(STATS)  # mean of the N individual test statistics for each r_H0, from Larsson et al. 2001:112, Eq.11
  UPSILON  = sqrt(dim_N) * (STATSbar - EZ) / sqrt(VZ) # standardized statistics, from Larsson et al. 2001:112, Eq.12
  
  # p-values
  if(is.matrix(distribution)){  # use empirical distribution provided as argument
    pt.pvals = rowMeans(UPSILON < distribution, na.rm=TRUE)
  }else if(distribution=="none"){  # return test statistics directly
    return(UPSILON)
  }else if(distribution=="theoretical"){  # use theoretical distribution
    pt.pvals = 1-pnorm(q=UPSILON, mean=0, sd=1)  # right-tailed test using standard normal distribution, from Larsson et al. 2001:114
  }else{ stop("Incorrect specification of the test distribution!") }
  
  # return result
  result = list(stats=UPSILON, pvals=pt.pvals)
  return(result)
}


# panel tests based on the meta-analytical combination of individual p-values, see Choi (2001) and Maddala,Wu (1999)
ptest.METApval <- function(PVALS, distribution="theoretical"){
  # define
  PVALS = as.matrix(PVALS)  # matrix of individual p-values
  if(ncol(PVALS)==1){ PVALS = t(PVALS) }  # ... in case there is just one hypothesis test
  dim_K = nrow(PVALS)  # number of endogenous variables resp. of panel tests 
  dim_N = ncol(PVALS)  # number of individuals
  
  # test statistics
  P  = -2*rowSums(log(PVALS))  # inverse chi-square test by Fisher (1932) for finite N, from Choi 2001:253, Eq.8 / Maddala,Wu 1999:636
  Pm = (P - 2*dim_N) / sqrt(4*dim_N)  # modified for infinite N, from Choi 2001:255, Eq.18
  Z  = rowSums(qnorm(PVALS)) / sqrt(dim_N)  # inverse normal test is already invariant to infinite N, from Choi 2001:253, Eq.9
  pt.stats = cbind(P, Pm, Z)
  
  # p-values
  if(is.matrix(distribution)){  # use empirical distribution provided as argument
    pt.pvals = rowMeans(c(P, Pm, Z) < distribution, na.rm=TRUE)
    pt.pvals = matrix(pt.pvals, nrow=dim_K, ncol=3, dimnames=dimnames(pt.stats))
  }else if(distribution=="none"){  # return test statistics directly
    return(pt.stats)
  }else if(distribution=="theoretical"){  # use theoretical distribution
    pt.pvals = cbind(
      P  = 1-pchisq(P, df=2*dim_N),   # from Choi 2001:255, Eq.11 / Maddala,Wu 1999:636
      Pm = 1-pnorm(Pm, mean=0, sd=1), # right-tailed test using standard normal distribution, from Choi 2001:255, Eq.19 / Eq.23
      Z  = pnorm(Z, mean=0, sd=1))    # left-tailed test using standard normal distribution, from Choi 2001:255, Eq.20 / Eq.24
  }else{ stop("Incorrect specification of the test distribution!") } 
  
  # return result
  result = list(stats=pt.stats, pvals=pt.pvals)
  return(result)
}


# panel tests based on the correlation-augmented combination of individual probits, from Hartung (1999) and Oersal,Arsova (2020)
ptest.CAIN <- function(PVALS, distribution="theoretical", dim_K=nrow(PVALS), r_H0=0:(dim_K-1), rho_eps){
  # define
  PVALS = as.matrix(PVALS)  # matrix of individual p-values
  if(ncol(PVALS)==1){ PVALS = t(PVALS) }  # ... in case there is just one hypothesis test
  dim_N = ncol(PVALS)  # number of individuals
  PROBITS = qnorm(PVALS)  # KxN matrix of individual probits
  
  # test statistics according to Hartung (1999), from Oersal,Arsova 2016:6, Eq.6-8
  LAMBDA   = array(1, dim=dim(PROBITS))  ### TODO function argument: KxN matrix of weights
  rho_prbt = c(1 - 1/(dim_N-1) * rowSums((PROBITS - rowSums(PROBITS)/dim_N)^2))
  rho_star = sapply(rho_prbt, FUN=function(x) max(c(-1/(dim_N), x)))
  kappa    = cbind(kappa_1 = 0.2, kappa_2 = 0.1*(1 + 1/(dim_N-1) - rho_star))
  rho_HA   = rho_star + kappa * sqrt(2/(dim_N+1)) * (1-rho_star)  # correction factor
  pt.stdev = sqrt(rowSums(LAMBDA^2) + (rowSums(LAMBDA)^2 - rowSums(LAMBDA^2)) * rho_HA)
  pt.kappa = rowSums(LAMBDA * qnorm(PVALS)) / pt.stdev  # test statistic
  
  # response surface approximation, from Oersal,Arsova 2016:10, Tab.2
  cain_rscoef = c(0.6319575, -0.5193669, 0.2721753, 
                  0.1821374, -0.0856903, 0.0041125,
                  0.0766267, -0.1008678, 0.1874919, 
                  0.1410229, -0.2029126, 0.0052557, -0.0000327)
  d_H0 = dim_K - r_H0  # number of stochastic trends under H0
  cain_term = cbind(rho_eps^2, 
                    sqrt(dim_K)* rho_eps^2, 
                    sqrt(dim_K)* rho_eps^4, 
                    r_H0/dim_K * rho_eps^2, 
                    r_H0/dim_K * rho_eps^4,
                    (r_H0 * rho_eps)^2, 
                    r_H0 * rho_eps^2, 
                    r_H0 * rho_eps^4, 
                    sqrt(d_H0) * rho_eps^2, 
                    1/d_H0 * rho_eps^2, 
                    1/d_H0 * rho_eps^4,
                    d_H0^2 * rho_eps^2, 
                    d_H0^4 * rho_eps^4)
  
  # test statistics, from Oersal,Arsova 2016:6, Eq.11
  rho_CAIN = c(cain_term %*% cain_rscoef)  # correction factor
  CAIN     = rowSums(PROBITS) / sqrt(dim_N + (dim_N^2 - dim_N) * rho_CAIN)
  pt.stats = cbind(pt.kappa, CAIN)
  
  # p-values
  if(is.matrix(distribution)){  # use empirical distribution provided as argument
    pt.pvals = rowMeans(c(pt.stats) < distribution, na.rm=TRUE)
    pt.pvals = matrix(pt.pvals, nrow=dim_K, ncol=3, dimnames=dimnames(pt.stats))
  }else if(distribution=="none"){  # return test statistics directly
    return(pt.stats)
  }else if(distribution=="theoretical"){  # use theoretical distribution
    pt.pvals = pnorm(pt.stats, mean=0, sd=1)  # left-tailed test using standard normal distribution 
  }else{ stop("Incorrect specification of the test distribution!") }
  
  # return result
  result = list(stats=pt.stats, pvals=pt.pvals, rho_tilde=rbind(t(rho_HA), CAIN=rho_CAIN))
  return(result)
}



#############################
###  AUXILIARY FUNCTIONS  ###
#############################
#
# The following functions serve as modules nested 
# in the calling functions. Notation within each 
# function mostly corresponds to the cited literature.


# factor-restricted SUR for VAR via least squares, from Empting,Herwartz (2024, Sec.3.1)
aux_pvar <- function(L.def, n.factors=FALSE, n.iterations=FALSE){
  # define
  L.dim_T = sapply(L.def, FUN=function(i) i$dim_T)
  dim_T = min(L.dim_T)      # number of time periods without presample in total panel
  dim_K = L.def[[1]]$dim_K  # number of endogenous variables
  dim_N = length(L.def)     # number of individuals
  
  names_i = names(L.def)
  names_k = rownames(L.def[[1]]$Y)
  names_u = paste0("u[ ", names_k, " ]")
  names_F = paste0("F[ ", 0:n.factors, " ]")[-1]
  
  # iterative estimation
  L.ZZinv = lapply(L.def, FUN=function(i) solve(tcrossprod(i$Z)))
  L.idx_j = lapply(1:dim_N, FUN=function(i) dim_K*(i-1) + 1:dim_K)
  L.idx_t = lapply(1:dim_N, FUN=function(i) seq.int(to=L.dim_T[i], length.out=dim_T))
  L.zeros = lapply(1:dim_N, FUN=function(i) matrix(0, nrow=n.factors, ncol=L.dim_T[i]-dim_T))
  L.resid = L.A = list()
  
  sqrt_T = sqrt(dim_T)
  LAMBDA = matrix(0, nrow=dim_N*dim_K, ncol=n.factors)  # initialize as OLS
  Ft     = matrix(0, nrow=n.factors, ncol=dim_T)
  U      = matrix(0, nrow=dim_N*dim_K, ncol=dim_T)
  
  for(n in 0:n.iterations){
    # estimate the VARX
    for(i in 1:dim_N){
      Y.idi = L.def[[i]]$Y - LAMBDA[L.idx_j[[i]], , drop=FALSE] %*% cbind(L.zeros[[i]], Ft)
      L.A[[i]] = tcrossprod(Y.idi, L.def[[i]]$Z) %*% L.ZZinv[[i]]  # multivariate OLS estimator
      L.resid[[i]] = L.def[[i]]$Y - L.A[[i]] %*% L.def[[i]]$Z  # residuals (complete MFE)
      U[L.idx_j[[i]], ] = L.resid[[i]][ , L.idx_t[[i]], drop=FALSE]  # transpose 'usd' to (T x K*N) OR use PCA$v if (K*N x T) 
    }
    
    # estimate the residual factors via PCA
    if(n.factors){
      usd = scale(t(U), center=FALSE, scale=TRUE)  # scaling against variables' differing units or magnitudes
      PCA = svd(usd, nu=n.factors)  # PCA via singular value decomposition 
      Ftt = sqrt_T*PCA$u  # normalized eigenvectors ev conforming to (Ft Ft')/T=I by construction of SVD
      LAMBDA = U %*% Ftt / dim_T
      Ft  = t(Ftt)
    }else{
      PCA = NULL
    }
    
    ### TODO: break after convergence of estimates or likelihood
  }
  
  # create individual "varx" objects
  L.varx = list()
  for(i in 1:dim_N){
    dim_T   = L.def[[i]]$dim_T
    dim_Kpn = L.def[[i]]$dim_Kpn
    OMEGA = tcrossprod(L.resid[[i]]) / dim_T  # MLE covariance matrix of residuals
    SIGMA = OMEGA * (dim_T/(dim_T-dim_Kpn))   # OLS covariance matrix of residuals
    B = cbind(diag(dim_K), LAMBDA[L.idx_j[[i]], , drop=FALSE])
    
    # return result
    dimnames(B) = list(names_k, c(names_u, names_F))
    L.varx[[i]] = c(list(A=L.A[[i]], B=B, OMEGA=OMEGA, SIGMA=SIGMA, resid=L.resid[[i]]), L.def[[i]])
    class(L.varx[[i]]) = "varx"
  }
  
  # return result
  names(L.varx) = names_i
  result = list(L.varx=L.varx, PCA=PCA, n.iterations=n)
  return(result)
}


# factor-restricted SUR for VECM, from Empting,Herwartz (2024, Sec.3.2)
aux_pvec <- function(L.def, L.beta=NULL, dim_r, idx_pool=NULL, n.factors=FALSE, n.iterations=FALSE){
  # define
  L.dim_T = sapply(L.def, FUN=function(i) i$dim_T)
  dim_T = min(L.dim_T)      # number of time periods without presample in total panel
  dim_K = L.def[[1]]$dim_K  # number of endogenous variables
  dim_N = length(L.def)     # number of individuals
  idx_k = 1:dim_K
  idx_i = 1:dim_N
  
  names_i = names(L.def)
  names_k = rownames(L.def[[1]]$y)
  names_u = paste0("u[ ", names_k, " ]")
  names_F = paste0("F[ ", 0:n.factors, " ]")[-1]
  
  # iterative estimation
  L.idx_j = lapply(1:dim_N, FUN=function(i) dim_K*(i-1) + 1:dim_K)
  L.idx_t = lapply(1:dim_N, FUN=function(i) seq.int(to=L.dim_T[i], length.out=dim_T))
  L.zeros = lapply(1:dim_N, FUN=function(i) matrix(0, nrow=n.factors, ncol=L.dim_T[i]-dim_T))
  
  sqrt_T = sqrt(dim_T)
  LAMBDA = matrix(0, nrow=dim_N*dim_K, ncol=n.factors)  # initialize
  Ft     = matrix(0, nrow=n.factors, ncol=dim_T)
  X      = matrix(0, nrow=dim_N*dim_K, ncol=dim_T)
  
  # for(n in 0:n.iterations){
    # concentrate out short-run effects
    L.RRR = lapply(L.def, FUN=function(i) aux_RRR(i$Z0, i$Z1, i$Z2, via_R0_R1=TRUE))
    
    # estimate cointegrating vectors
    if( is.null(L.beta) ){
      if( is.null(idx_pool) ){  # ... by RRR
        L.beta = lapply(L.RRR, FUN=function(i) aux_beta(V=i$V, dim_r=dim_r, normalize="natural"))
      }else{  # ... by two-step
        L.beta = aux_2StepBR(L.RRR, r_H0=dim_r, idx_pool=idx_pool, n.iterations=n.iterations)$L.beta  # iterations only via main-loop
    }}
    
  #   # gather errors of the long-run relation, compare Bai et al 2009:86, Eq.11 under r=1
  #   L.error = sapply(idx_i, FUN=function(i) t(L.beta[[i]][idx_k, , drop=FALSE]) %*% L.RRR[[i]]$R1[idx_k, , drop=FALSE], simplify=FALSE)
  #   ### TODO: allow for individual-specific dim_T and dim_p! Maybe fill Ft with zeros?
  #   
  #   # estimate the common factors via PCA
  #   if(n.factors){
  #     X = do.call("rbind", L.error)  # matrix X must be of dimension (T x r*N)
  #     xsd = scale(t(X), center=FALSE, scale=TRUE)  # scaling against variables' differing units or magnitudes
  #     PCA = svd(xsd, nu=n.factors)  # ... PCA via singular value decomposition
  #     Ftt = sqrt_T*PCA$u  # normalized eigenvectors ev conforming to (Ft Ft')/T=I by construction of SVD
  #     LAMBDA = X %*% Ftt / dim_T
  #     Ft  = t(Ftt)
  #   }else{
     PCA = NULL
  #   }
  #   
  #   ### TODO: break after convergence of estimates or likelihood
  # }
  
  # create individual "varx" objects
  L.varx = list()
  for(i in 1:dim_N){
    vecm = aux_VECM(beta=L.beta[[i]], RRR=L.RRR[[i]])
    rvar = aux_vec2var(PI=vecm$PI, GAMMA=vecm$GAMMA, dim_p=L.def[[i]]$dim_p)
    B = cbind(diag(dim_K), LAMBDA[L.idx_j[[i]], , drop=FALSE])
    
    # return result
    dimnames(B) = list(names_k, c(names_u, names_F))
    L.varx[[i]] = list(A=rvar$A, B=B, y=L.def[[i]]$y, x=L.def[[i]]$x, D1=L.def[[i]]$D1, D2=L.def[[i]]$D2, 
           RRR=L.RRR[[i]], beta=L.beta[[i]], VECM=vecm, 
           resid=vecm$resid, OMEGA=vecm$OMEGA, SIGMA=vecm$SIGMA, dim_r=dim_r, dim_K=dim_K, 
           dim_L=L.def[[i]]$dim_L, dim_T=L.def[[i]]$dim_T, dim_p=L.def[[i]]$dim_p, dim_q=NULL)
    class(L.varx[[i]]) = "varx"
  }
  
  # return result
  names(L.varx) = names_i
  result = list(L.varx=L.varx, PCA=PCA, n.iterations=n.iterations)
  return(result)
}


# parametric two-step estimator for homogeneous cointegrating vectors in panels, from Breitung (2005)
aux_2StepBR <- function(L.RRR, r_H0=0:(dim_K-1), idx_pool=1:dim_K, n.iterations=FALSE){
  # define
  dim_N = length(L.RRR)       # number of individuals
  dim_K = nrow(L.RRR[[1]]$R0) # number of endogenous variables
  names_i = names(L.RRR)            # names for individuals
  names_r = paste0("r", r_H0)       # names for cointegration ranks
  names_k = rownames(L.RRR[[1]]$R0) # names for endogenous variables
  idx_pool = unique(sort(idx_pool)) # preserve the row ordering within beta
  
  # pooled regressor matrices for reduced rank regression, \tilde{y}_it from Breitung 2005:160, Eq.17 (II)
  R1 = do.call("cbind", lapply(L.RRR[], FUN=function(x) x$R1[idx_pool, , drop=FALSE]))  # matrix (K+n x NT) of concentrated I(1)-regressors
  RR = tcrossprod(R1, R1)  # moment matrix of concentrated regressors for second step OLS-regression
  
  # estimation
  L.beta  = matrix(list(), nrow=dim_N, ncol=length(r_H0), dimnames=list(names_i, names_r))  # cointegrating vectors for each i and r
  L.alpha = matrix(list(), nrow=dim_N, ncol=length(r_H0), dimnames=list(names_i, names_r))  # orthogonal complements for each i and r
  L.iter  = rep(0, each=length(r_H0)); names(L.iter) = names_r  # number of iterations till convergence
  for(r in r_H0){
    idx_r = paste0("r", r)
    idx_1 = 0:r  # row index for sub-vectors of the r first elements in y_it, from Breitung 2005:156, Eq.6
    idx_2 = setdiff(idx_pool, idx_1)  # row index for second sub-vectors in y_it
    I_r   = diag(r); rownames(I_r) = rownames(L.RRR[[1]]$R1)[idx_1]
    
    if(r==0){
      # first step, from Breitung 2005:155
      alpha = matrix(0, nrow=dim_K, ncol=0, dimnames=list(names_k, NULL))
      L.alpha[ , idx_r] = lapply(1:dim_N, function(i) alpha)
      
      # second step, from Breitung 2005:156/160, Eq.4/17
      L.beta[ , idx_r] = lapply(L.RRR, function(i) matrix(0, nrow=nrow(i$R1), ncol=0, dimnames=list(rownames(i$R1), NULL)))
      rm(alpha)
      
    }else{
      # iterative estimation, from Wagner,Hlouskova 2009:195
      L.beta[ , idx_r] = lapply(L.RRR, function(i) aux_beta(V=i$V, dim_r=r, normalize="natural"))  # first step
      L.vecm = list()  # individual VECM
      L.a    = list()  # individual \alpha_i^{(+)}
      for(n in 0:n.iterations){
        
        # first step, from Breitung 2005:155
        for(i in 1:dim_N){
          L.vecm[[i]] = aux_VECM(beta=L.beta[[i, idx_r]], RRR=L.RRR[[i]])
          
          # calculate \alpha_i^{(+)}, from Breitung 2005:155, Eq.3
          gamma_tr = t(L.vecm[[i]]$alpha) %*% solve(L.vecm[[i]]$SIGMA)  # transposed \gamma_i
          L.a[[i]] = solve(gamma_tr %*% L.vecm[[i]]$alpha) %*% gamma_tr  # \alpha_i^{(+)}
          rm(gamma_tr)
        }
        
        # second step using pooled conditional OLS-estimator, from Breitung 2005:156/160, Eq.6/17
        if(length(idx_2) != 0){
          # ... with homogeneous cointegrating coefficients B_2S
          R0_plus = NULL
          for(i in 1:dim_N){
            beta_13  = L.beta[[i, idx_r]][-idx_2, , drop=FALSE]  # I_r and heterogeneous beta_3 to be partialled out
            R0i_plus = L.a[[i]] %*% L.RRR[[i]]$R0 - t(beta_13) %*% L.RRR[[i]]$R1[-idx_2, , drop=FALSE]
            R0_plus  = cbind(R0_plus, R0i_plus)  # matrix (r x NT) of pooled concentrated regressands
            rm(beta_13, R0i_plus) }
          R1R1inv = solve(RR[idx_2, idx_2, drop=FALSE])
          B_2S    = tcrossprod(R0_plus, R1[idx_2, , drop=FALSE]) %*% R1R1inv
          beta_2S = rbind(I_r, t(B_2S))
        }else{
          # ... without homogeneous coefficients B_2S
          beta_2S = I_r
        }
        
        # second step using individual conditional OLS-estimator, see Ahn,Reinsel (1990)
        for(i in 1:dim_N){
          dim_Kn1 = nrow(L.RRR[[i]]$R1)  # individual-specific K+n_{1i}
          idx_3   = setdiff(1:dim_Kn1, c(idx_1, idx_2))
          if(length(idx_3) != 0){
            # ... with heterogeneous coefficients Bi_2S
            R1i      = L.RRR[[i]]$R1[idx_3, , drop=FALSE]
            R0i_plus = L.a[[i]] %*% L.RRR[[i]]$R0 - t(beta_2S) %*% L.RRR[[i]]$R1[-idx_3, , drop=FALSE]
            R1R1inv  = solve(tcrossprod(R1i, R1i))
            Bi_2S    = tcrossprod(R0i_plus, R1i) %*% R1R1inv
            idx_123  = c(idx_1, idx_2, idx_3)  # index to restore the original ordering of beta_i
            L.beta[[i, idx_r]] = rbind(beta_2S, t(Bi_2S))[idx_123, , drop=FALSE]
          }else{
            # ... without heterogeneous coefficients Bi_2S
            L.beta[[i, idx_r]] = beta_2S
          }
        }
        
        # break loop after convergence of the likelihood
        ### c0 = -(dim_N*dim_K*dim_T)/2 * log(2*pi)   
        ### logLik[n+1] = c0 - sum(sapply(L.vecm, function(i) log(det(i$OMEGA)*(dim_T/2)))  # log-likelihood, from Breitung 2005:154, Eq.2
        ### if(isTRUE(logLik[n+1] - logLik[n] < conv.crit)){ break() }
      }
      
      # collect results for rank 'r'
      L.alpha[ , idx_r] = lapply(L.vecm, function(i) i$alpha)  # previous-step \alpha for LM-test, from Breitung 2005:158
      L.iter[[idx_r]]   = n
    }
  }
  
  # return result
  result = list(L.beta=L.beta, L.alpha=L.alpha, L.iter=L.iter)
  return(result)
}


# mean-group estimator, from Rebucci (2010) / Pesaran 2006:982, Eq.53 / Eq.58
aux_MG <- function(x, idx_par, w=NULL){ 
  # collect individual coefficient matrices in array
  if(is.list(x)){
    A.coef = aux_getPAR(L.varx=x, idx_par=idx_par, A.fill=0)  # fill A_ij=0 for j>p_i
  }else if(is.array(x)){
    A.coef = x
  }
  
  # define
  n.dims = length(dim(A.coef))  # number of dimensions
  idx_mg = 1:(n.dims-1)  # apply on the last margin
  
  if(is.null(w)){
    # panel estimation by MG, from Pesaran 2006:982, Eq.53 / Eq.58
    A.mean = apply(A.coef, MARGIN=idx_mg, FUN=function(x) mean(x))  # group mean
    A.var  = apply(A.coef, MARGIN=idx_mg, FUN=function(x) var(x))   # group (sample) variance
  }else if(is.logical(w) | is.character(w)){
    # subset MG
    L.idx  = c(replicate(n.dims-1, list(TRUE)), list(w))
    A.coef = do.call(`[`, c(list(A.coef), L.idx, list(drop=FALSE))) # subset
    A.mean = apply(A.coef, MARGIN=idx_mg, FUN=function(x) mean(x))  # group mean
    A.var  = apply(A.coef, MARGIN=idx_mg, FUN=function(x) var(x))   # group (sample) variance
  }else if(is.numeric(w)){
    # weighted MG
    w = w/sum(w)  # normalize the weights
    L.idx  = c(replicate(n.dims-1, list(TRUE)), list(w != 0))
    A.coef = apply(A.coef, MARGIN=idx_mg, FUN=function(x) w*x)      # weighting (transposes the array!)
    A.coef = aperm(A.coef, perm=c(idx_mg+1, 1))  # re-transpose
    A.coef = do.call(`[`, c(list(A.coef), L.idx, list(drop=FALSE))) # subset (remove zeros)
    A.mean = apply(A.coef, MARGIN=idx_mg, FUN=function(x) sum(x))   # group mean
    A.var  = array(NA, dim=dim(A.mean))  ### Consider the bootstrap functions instead.
  }else{ 
    stop("Please provide NULL, a numeric, a logical, or a character vector for argument 'w'.") 
  }
  
  # return result
  result = list(coef=A.coef, mean=A.mean, var=A.var)
  return(result)
}


# estimate common factors for PANIC, from Bai,Ng (2004)
aux_ComFact <- function(X, trend=TRUE, scaled=TRUE, SVD=TRUE, n.factors, D=NULL){
  ### see also Silvestre,Surdeanu 2011:5 for VAR process
  # define
  if(!is.null(D)){ warning("The factor estimation by PCA ignores period-specific deterministic regressors.", call.=FALSE) }  # Bai,Ng 2013:23, Ch.5, adopt D in factor models.
  X = as.matrix(X)  # matrix X must be of dimension (T x K*N)
  dim_T = nrow(X)   # number of time periods
  xit = diff(X)  # first-differences against non-stationarity,  from Bai,Ng 2004:1138, Eq.9
  xit = scale(xit, center=trend, scale=FALSE)  # demeaning against a (first-differenced) linear trend
  xsd = scale(xit, center=FALSE, scale=scaled) # scaling against variables' differing units or magnitudes, from Oersal,Arsova 2017:67
  
  # principal component analysis
  if(SVD){
    PCA_svd   = svd(xsd, nu=n.factors)  # ... via singular value decomposition (faster and more stable)
    eigenvecs = PCA_svd$u   ### [ , 0:n.factors , drop=FALSE]  # normalized eigenvectors ev conforming to ev'ev=I by construction of SVD
    eigenvals = PCA_svd$d^2 ### / (dim_T-1) for scaled eigenvalues from PCA on the covariance or correlation matrix
  }else{
    PCA_eigen = eigen(tcrossprod(xsd))  # ... via spectral decomposition, from Bai,Ng 2004:1132,1138
    eigenvecs = PCA_eigen$vectors[ ,0:n.factors , drop=FALSE]  # normalized eigenvectors ev conforming to ev'ev=I by convention
    eigenvals = PCA_eigen$values ### / (dim_T-1) for scaled eigenvalues from PCA on the covariance or correlation matrix
  }
  
  # estimate factors and loadings
  ft = sqrt(dim_T-1) * eigenvecs  # common factors as (T-1 x n.factors) matrix
  Ft = apply(rbind(0, ft),  MARGIN=2, FUN=function(x) cumsum(x))  # cumulate first-differences into levels (without trend if x_it has zero mean)
  LAMBDA = t(xit)%*%ft / (dim_T-1)  # individual loadings as (K*N x n.factors) matrix; by multivariate OLS and normalization ft'ft/(T-1)=ev'ev=I
  
  # estimate idiosyncratic components
  zit = xit - ft%*%t(LAMBDA)  # from Bai,Ng 2004:1133, Eq.6(II) / Silvestre,Surdeanu 2011:5
  eit = apply(rbind(0, zit), MARGIN=2, FUN=function(x) cumsum(x))  # idiosyncratic components in levels (without trend if x_it has zero mean)
  eit2 = X - Ft%*%t(LAMBDA)  # idiosyncratic components in levels as (T x K*N) matrix, from Arsova,Oersal 2018:1036
  
  # return result
  result = list(LAMBDA=LAMBDA, Ft=Ft, eit=eit, eit2=eit2, eigenvals=eigenvals)
  return(result)
}


# determine n.factors by edge distribution, from Onatski (2010)
aux_ONC <- function(eigenvals, r_max=20, n.iterations=4){
  # step 1: difference eigenvalues from the sample covariance matrix XX'/T
  ### Scaling the matrix XX' for PCA and thus all its eigenvalues does not affect ONC.
  eigendiffs = eigenvals[1:r_max] - eigenvals[2:(r_max+1)]
  
  # edge distribution (ED), from Onatski 2010:1008
  j = r_max + 1
  delta = rep(NA, n.iterations)
  r_hat = rep(NA, n.iterations)
  for(i in 1:n.iterations){
    # step 2: auxiliary OLS regression
    y = eigenvals[j:(j+4)]
    x = cbind(1, ((j-1):(j+3))^(2/3))
    b = c(solve(crossprod(x)) %*% crossprod(x, y))
    delta[i] = 2*abs(b[2])
    
    # step 3: estimate number of factors
    idx_pass = c(0, which(eigendiffs >= delta[i]))
    r_hat[i] = max(idx_pass)  # highest rank whose eigendiff still passes the threshold
    
    # step 4: iterate to reach convergence
    j = r_hat[i] + 1
  }
  
  # return result
  idx_cv = c(which(r_hat[-n.iterations] == r_hat[-1]), NA)  # find convergence
  result = list(selection=r_hat[idx_cv][1], converge=rbind(delta, r_hat))
  return(result)
}


# determine n.factors by ratio of eigenvalues, from Ahn,Horenstein (2013)
aux_AHC <- function(eigenvals, r_max=20){
  # define
  dim_m = length(eigenvals)  # = min(dim(X))
  idx_k = 1:(r_max+1)  # for k=0,...,r_max
  
  # eigenvalues from the sample covariance matrix XX'/(TN)
  ### Scaling the matrix XX' for PCA and thus eigenvals does not affect the ratios ER(k) and GR(k)
  ### possibly first-differenced I(1)-data, see Ahn,Horenstein 2013:1210
  eigen_sum  = sum(eigenvals)
  eigen_zero = eigen_sum / (dim_m * log(dim_m))  # for k=0
  eigen_val0 = c(eigen_zero, eigenvals[idx_k])
  eigen_star = eigen_val0 / (eigen_sum - cumsum(c(0, eigenvals[idx_k])))
  
  # eigenvalue ratios, from Ahn,Horenstein 2013:1207 / Corona et al. 2017:357, Eq.9/10
  ER = eigen_val0[idx_k] / eigen_val0[idx_k+1]
  GR = log(1 + eigen_star[idx_k]) / log(1 + eigen_star[idx_k+1])
  
  # return result
  idx_mx = c(ER=which.max(ER), GR=which.max(GR))  # index on 1:(r_max+1) for k=0,...,r_max
  result = list(selection=idx_mx-1, criteria=rbind(k=0:r_max, ER, GR))
  return(result)
}


# panel criteria for a given number of common factors 'k', from Bai,Ng (2002) / Bai (2004)
aux_PIC <- function(X, k, Vk0, Vkmax0, Vk1=Vk0, Vkmax1=Vkmax0){
  ### based on variance of idiosyncratic components: Vk0 (Bai,Ng 2002:201) and Vk1 (Bai 2004:145)
  # define
  X = as.matrix(X)  # data panel as (T x K*N) matrix
  dim_KN = ncol(X)  # number of time series in the (high-dimensional) data panel
  dim_T  = nrow(X)  # number of time periods
  n_vals = dim_KN * dim_T  # total number of values in the panel
  C2_NT  = min(dim_KN, dim_T)  # minimum dimension of the data panel
  
  # calculate stationary panel criteria, from Bai,Ng 2002:201, Eq.9
  ### possibly first-differenced I(1)-data, 
  ### ... but not for factor models with distributed lags of Ft, see Bai 2004:143, Ch.3.1 / 2004:151
  p1 = (dim_KN + dim_T) / n_vals * log(n_vals / (dim_KN + dim_T))
  p2 = (dim_KN + dim_T) / n_vals * log(C2_NT)
  p3 = log(C2_NT) / C2_NT
  PC = Vk0 + k * Vkmax0 * c(p1, p2, p3)  # "Panel C_p Criteria"
  IC = log(Vk0) + k * c(p1, p2, p3)  # "Panel Information Criteria"
  
  # calculate integrated panel criteria, from  Bai 2004:145, Eq.12
  ### not for cointegrated resp. mixed I(1) and I(0) factors, 
  ### ... but allows DFM with distributed lags of Ft, see Bai 2004:151
  alphaT = dim_T / (4*log(log(dim_T)))  # law of iterated logarithm
  ip3 = (dim_KN + dim_T - k) / n_vals * log(n_vals)
  IPC = Vk1 + k * Vkmax1 * alphaT * c(p1, p2, ip3)  # "Integrated Panel Criteria"
  
  # return result
  result = cbind(PC, IC, IPC)
  dimnames(result) = list(c("p1", "p2", "p3"), c("PC", "IC", "IPC"))
  return(result)
}


# simulated moments for the asymptotic distribution Z_d of the trace statistic and MSB-test statistic distributions
coint_moments = list(
  # moments from Larsson et al. 2001:114, Tab.1 
  # ... (they do not state the number of replications and time periods; Johansen 1996:212 uses 100,000 replications with T=400)
  # ... as used in the panel tests by Larsson et al. (2001) and Breitung (2005) 
  Case1 = cbind(
    TR_EZ=c(1.137,  6.086, 14.955, 27.729, 44.392,  64.960,  89.360, 117.519, 149.441, 185.082, 224.450, 267.708),
    TR_VZ=c(2.212, 10.535, 24.733, 45.264, 71.284, 103.452, 139.680, 183.997, 233.053, 286.483, 343.179, 411.679)),
  # moments from Breitung 2005:171, Tab.B1 (20,000 replications with T=500)
  Case2 = cbind(
    TR_EZ=c(3.051,  9.99, 20.88, 35.67,  54.33,  76.94),
    TR_VZ=c(7.003, 18.46, 35.86, 58.07,  85.13, 119.70)),
  Case3 = cbind(
    TR_EZ=c( 0.98,  8.27, 19.35, 34.18,  53.05,  75.61),
    TR_VZ=c( 1.91, 14.28, 31.84, 54.28,  83.50, 116.70)),
  Case4 = cbind(
    TR_EZ=c( 6.27, 16.28, 30.21, 48.01,  69.65,  94.93),
    TR_VZ=c(10.45, 25.50, 45.13, 72.95, 104.07, 139.70)),
  # moments from Oersal,Droge 2014:6 Tab.1 / Oersal 2008:106, Tab.5.1 (20,000 replications with T=1000)
  # ... as used in the panel SL tests by Droge,Oersal (2014) and Arsova,Oersal (2018)
  SL_trd14 = cbind(
    TR_EZ=c(2.69,  8.86, 18.85, 32.78, 50.58,  72.44,  97.91, 127.55, 161.20, 198.43, 239.70, 284.87),
    TR_VZ=c(4.38, 13.37, 28.23, 47.94, 73.74, 105.33, 143.68, 187.28, 238.00, 300.91, 357.05, 424.86)),
  # moments from Arsova,Oersal 2018:Appendix p.19, Tab.A1 (via response surface approximation)
  SL_trend = cbind(
    TR_EZ=c(2.689,  8.924, 19.011, 33.036, 51.023,  73.042,  99.036, 129.025, 163.003, 200.971, 242.960, 289.002),
    TR_VZ=c(4.396, 13.725, 28.501, 48.837, 75.430, 107.953, 147.468, 193.158, 241.215, 297.598, 360.760, 428.035)),
  # moments from Silvestre,Surdeanu 2011:10, Tab.1(II) (10,000 replications with T=1000, T=500 and T=200)
  # ... as used in their panel MSB test
  MSB_mean1000 = cbind(
    EZ=c(0.50445, 0.08580, 0.04048, 0.02578, 0.01884, 0.01475),
    VZ=c(0.34863, 0.00364, 0.00039, 0.00009, 0.00003, 0.00002)),
  MSB_trend1000 = cbind(
    EZ=c(0.16622, 0.05539, 0.03132, 0.02172, 0.01651, 0.01323),
    VZ=c(0.02267, 0.00095, 0.00017, 0.00005, 0.00002, 0.00001)),
  MSB_mean500 = cbind(
    EZ=c(0.49960, 0.08806, 0.04191, 0.02708, 0.01995, 0.01571),
    VZ=c(0.31146, 0.00378, 0.00040, 0.00010, 0.00004, 0.00002)),
  MSB_trend500 = cbind(
    EZ=c(0.16825, 0.05705, 0.03281, 0.02270, 0.01744, 0.01427),
    VZ=c(0.02122, 0.00106, 0.00018, 0.00005, 0.00002, 0.00001)),
  MSB_mean200 = cbind(
    EZ=c(0.50454, 0.09049, 0.04356, 0.02864, 0.02115, 0.01710),
    VZ=c(0.32442, 0.00386, 0.00042, 0.00011, 0.00003, 0.00002)),
  MSB_trend200 = cbind(
    EZ=c(0.17637, 0.05943, 0.03431, 0.02454, 0.01900, 0.01570),
    VZ=c(0.02195, 0.00109, 0.00018, 0.00006, 0.00002, 0.00001)),
  # moments from Silvestre,Surdeanu 2011:supplement, GAUSS code (T=100 and T=50)
  MSB_mean100 = cbind(
    EZ=c(0.50363233, 0.09293090, 0.04653255, 0.03106941, 0.02378946, 0.01969598),
    VZ=c(0.29973853, 0.00387378, 0.00042435, 0.00010149, 0.00003807, 0.00001723)),
  MSB_trend100 = cbind(
    EZ=c(0.17829789, 0.06172200, 0.03708937, 0.02696152, 0.02160916, 0.01831139),
    VZ=c(0.02079461, 0.00108132, 0.00019160, 0.00005735, 0.00002413, 0.00001183)),
  MSB_mean50 = cbind(
    EZ=c(0.49534851, 0.09760915, 0.05269258, 0.03900235, 0.03267331, 0.02924112),
    VZ=c(0.30860402, 0.00380154, 0.00037111, 0.00008941, 0.00002977, 0.00001272)),
  MSB_trend50 = cbind(
    EZ=c(0.17309672, 0.06820676, 0.04504875, 0.03587122, 0.03122986, 0.02864541),
    VZ=c(0.01673810, 0.00094535, 0.00016972, 0.00004706, 0.00001728, 0.00000913))
)

