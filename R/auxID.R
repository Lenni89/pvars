

#############################
###  AUXILIARY FUNCTIONS  ###
#############################
#
# The following functions serve as modules nested 
# in the calling functions. Notation within each 
# function mostly corresponds to the cited literature.


# unit-wise decomposition of covariance matrices for pooling the residuals, from Herwartz, Wang (2024), Herwartz 2017:12
aux_UPool <- function(L.resid, L.Ucov){  
  ### arguments: list of individual residuals and of their covariance matrices ###
  # define
  dim_N = length(L.resid)    # number of individuals
  dim_K = nrow(L.Ucov[[1]])  # number of endogenous variables
  names_i = if( is.null(names(L.resid)) ){ 1:dim_N }else{ names(L.resid) }
  names_s = paste0("epsilon[ ", 1:dim_K, " ]")
  
  # initialize
  R.eps  = NULL     # pooled sample (K x T*N)
  L.Ustd = list()   # matrices with standard deviations on the main diagonal
  L.Uchl = list()   # lower triangular, Cholesky-decomposed correlation matrices
  
  # unit-wise decomposition for pooling the residuals, from Herwartz 2017:12, Eq.20
  for(i in names_i){
    u = aux_asDataMatrix(L.resid[[i]], "u")  # individual residual matrix is KxT
    
    U.cov = L.Ucov[[i]]  # covariance matrix; Herwartz 2017:6, Eq.7(II) uses ML estimates
    U.std = diag(sqrt(diag(U.cov))) # matrix with standard deviations on the main diagonal
    U.std_inv = solve(U.std)  # inverse standard deviation matrix
    
    U.cor = U.std_inv %*% U.cov %*% U.std_inv  # correlation matrix
    U.chl = t(chol(U.cor))  # lower triangular, Cholesky-decomposed correlation matrix
    epsil = solve(U.chl) %*% U.std_inv %*% u  # = solve(t(chol(U.cov))) %*% u  # orthogonalized residuals
    
    R.eps = cbind(R.eps, epsil)  # pooled sample (K x T*N), from Herwartz 2017:12, Eq.21
    L.Ustd[[i]] = U.std
    L.Uchl[[i]] = U.chl
    rm(u, U.cov, U.std, U.std_inv, U.cor, U.chl, epsil)
  }
  
  # return result
  rownames(R.eps) = names_s
  result = list(eps=t(R.eps), L.Ustd=L.Ustd, L.Uchl=L.Uchl)
  return(result)
}


# pre-step PCAs for group ICA, from Risk et al. 2014:227 / Calhoun et al. 2001:151
aux_UGroup <- function(L.resid, n.factors){
  # define
  L.dim_T = sapply(L.resid, FUN=function(x) ncol(x))  # residual matrices are (K x T)
  dim_T = min(L.dim_T)  # number of periods without presample in total panel

  # pre-step PCAs, from Risk et al. 2014:227
  L.pca1 = lapply(L.resid, FUN=function(x) svd(tail(t(x), n=dim_T)))  # first-step PCA on (T x K) # first dim in tail() is T
  Yhat   = sqrt(dim_T) * do.call("cbind", lapply(L.pca1, FUN=function(x) x$u))  # stacked whitened data matrices (T x NK)
  R.pca2 = svd(Yhat, nu=n.factors)  # second-step PCA
  Zhat   = sqrt(dim_T) * R.pca2$u   # whitened data matrix (T x n.factors)
  
  # return result
  colnames(Zhat) = paste0("epsilon[ ", 1:n.factors, " ]")
  result = list(L.pca1=L.pca1, R.pca2=R.pca2, eps=Zhat)
  return(result)
}


# return mixing matrix with unique sign and column ordering
aux_sico <- function(B, B.orig=NULL, names_s=NULL){
  if(is.null(B.orig)){
    ### I) maximize the trace, from Herwartz 2017:25
    # define
    dim_K = nrow(B)
    dim_S = ncol(B)
    if(dim_K > dim_S){ stop("Number of columns in B must be at least number of rows.") }
    perms = aux_permute(dim_K)  # (K! x K) indices for the permutation of the K first columns in s=1,...,S
    
    # set column order: maximize the trace
    traces = apply(perms, MARGIN=1, FUN=function(idx_s) sum(abs(diag(B[ , idx_s]))))
    idx_co = perms[which.max(traces), , drop=TRUE]
    idx_co = c(idx_co, c(dim_K:dim_S)[-1])  # does not touch potential excess shocks
    
    # set column signs: positive elements on the main diagonal
    idx_si = (0 > diag(B[ , idx_co, drop=FALSE]))
    idx_si = c(idx_si, rep(FALSE, dim_S-dim_K))  # does not touch potential excess shocks
    
    # return result
    mat_perm = diag(nrow=dim_S)[ , idx_co] %*% diag(-1*idx_si + 1*!idx_si)  # permutation matrix
    colnames(mat_perm) = names_s
    B.sico = B %*% mat_perm
    result = list(B=B.sico, idx_co=idx_co, idx_si=idx_si, mat_perm=mat_perm)
    return(result)
    
  }else{
    ### II) minimize Forbenius distance, from Herwartz 2018 online supp.:2
    ### SOURCE:
    ### https://github.com/cran/steadyICA/blob/master/R/matchfun.R
    # define
    Sigma_u_hat_old = B.orig %*% t(B.orig)
    Sigma_u_star    = B %*% t(B)
    Pstar1 <- suppressMessages(expm::sqrtm(Sigma_u_hat_old)) %*% solve(suppressMessages(expm::sqrtm(Sigma_u_star))) %*% B
    diag_sigma_root <- diag(diag(suppressMessages(expm::sqrtm(Sigma_u_hat_old))))
    M1 = t(solve(diag_sigma_root)%*%Pstar1)
    M2 = t(solve(diag_sigma_root)%*%B.orig)
    
    standardize = TRUE
    if(standardize) {
      temp = diag((diag(M1%*%t(M1)))^(-1/2))
      M1 = temp%*%M1
      temp = diag((diag(M2%*%t(M2)))^(-1/2))
      M2 = temp%*%M2
    }
    
    d = n.comp = nrow(M1)
    p = ncol(M1)
    if(d!=nrow(M2)) stop("M1 and M2 must have the same number of rows")
    if(p!=ncol(M2)) stop("M1 and M2 must have the same number of columns")
    
    l2.mat1 = l2.mat2 = matrix(NA,nrow=n.comp,ncol=n.comp)  
    for (j in 1:n.comp) {
      for (i in 1:n.comp) {
        #since signs are arbitrary, take min of plus and minus:
        l2.mat1[i,j] = sum((M2[i,]-M1[j,])^2)
        l2.mat2[i,j] = sum((M2[i,]+M1[j,])^2)
      }
    }
    l2.mat1 = sqrt(l2.mat1)
    l2.mat2 = sqrt(l2.mat2)
    
    #take the min of plus/min l2 distances. This is okay because solve_LSAP is one to one
    l2.mat = l2.mat1*(l2.mat1<=l2.mat2) + l2.mat2*(l2.mat2<l2.mat1)
    map = as.vector(clue::solve_LSAP(l2.mat))
    
    #retain relevant l2 distances:
    l2.1 = diag(l2.mat1[,map])
    l2.2 = diag(l2.mat2[,map])
    
    #sign.change is for re-ordered matrix 2
    sign.change = -1*(l2.2<l2.1) + 1*(l2.1<=l2.2) 
    
    # return result
    mat_perm = diag(nrow=n.comp)[,map] %*% diag(sign.change)
    B.sico   = B %*% mat_perm
    dimnames(B.sico) = dimnames(B.orig)
    return(B.sico)
  
  ### TODO: III) Coherent: minimize Forbenius distance, from Bernoth,Herwartz 2018:9
  ### Should be the same for all individuals if "group"?  
  ### perm_i = svars:::frobICA_mod(star_pid$L.varx[[i]]$B, L.B[[i]], standardize=TRUE)$perm
  ### star_pid$L.varx[[i]]$B = aux_sico(star_pid$L.varx[[i]]$B, perm_i)
  }
}


# all 'n!' permutations of a vector with 'n' elements
aux_permute <- function(n){
  if(n==1){
    return(matrix(1))
  }else{
    sp  <- aux_permute(n-1)
    p   <- nrow(sp)
    mat <- matrix(nrow=n*p, ncol=n)
    idx <- 1:p
    for(i in 1:n){
      idx_rows <- idx + (i-1)*p
      mat[idx_rows, ] <- cbind(i, sp+(sp>=i))
    }
    return(mat)
  }
}


# orthogonal factor of a signed LQ decomposition via QR decomposition
aux_signLQ <- function(M, S=NULL){
  # define
  dim_M = ncol(M)  # M is a square matrix
  R.qr = qr(t(M))  # Use the default QR decomposition (LAPACK = FALSE) only!
  Q = qr.Q(R.qr)   # Note that qr.Q() returns 1 if dim_M = 1 and <0 x 0> if dim_M = 0.
  
  # sign-correction such that R_tr = t(M) %*% Q is lower-triangular with positive diagonal elements 
  diag_R = if(is.null(S)){ diag(qr.R(R.qr)) }else{ diag(S %*% Q) }  #  ... or of (S %*% Q) if an anchor is supplied.
  idx_si = (0 > diag_R)
  Q = Q %*% diag(c(-1*idx_si + 1*!idx_si), nrow=dim_M)
  
  # return result
  return(Q)
}


# identify structural shocks with proxy variables / external instruments
aux_idIV <- function(u, m, S2=c("MR", "JL", "LT", "NQ"), SIGMA_uu, dim_K1=NULL){
  # define
  u = aux_asDataMatrix(u, "u")  # named residual matrix is KxT
  m = aux_asDataMatrix(m, "m")  # named proxy matrix is LxT
  dim_T = ncol(u)  # number of observation without presample
  dim_K = nrow(u)  # number of endogenous variables in the VAR
  dim_L = nrow(m)  # number of proxy variables
  dim_S = if( S2 %in% c("MR", "JL") & is.null(dim_K1) ){ dim_L }else{ dim_K }  # number of identified shocks
  if(ncol(m) > dim_T){ m = tail(m, n=c(dim_L, dim_T)) }  # remove presample periods
  names_k = if( !is.null(rownames(u)) ){ rownames(u) }else{ paste0("y[ ", 1:dim_K, " ]") }
  names_s = paste0("epsilon[ ", 1:dim_S, " ]")
  
  # void results for some S2
  R.udv = list()
  F_stats = NULL 
  
  # identify using multiple proxy variables
  if( S2=="NQ" | S2=="LT" | !is.null(dim_K1) ){
    # ... from Empting et al. (2026)
    if( is.null(dim_K1) ){ dim_K1 = dim_L }
    B0    = t(chol(SIGMA_uu))
    eps_u = solve(B0)%*%u
    R.udv = aux_idSVD(eps_u=eps_u, eps_m=m, S2=S2, dim_K1=dim_K1)
    B1    = B0 %*% R.udv$Q  # fully identified
    
    if(S2=="MR"){
      R.udv$w = R.udv$w %*% aux_NQ2MR(B1, dim_K1=dim_K1)
      R.udv$Q = R.udv$u %*% R.udv$w
      B1      = B0 %*% R.udv$Q
    }
    
  }else if( S2=="JL" ){
    # ... from Jentsch, Lunsford WP2019:34, Appendix A
    SIGMA_mm = tcrossprod(m)/dim_T
    SIGMA_um = tcrossprod(u, m)/dim_T
    SIGMA_mu = t(SIGMA_um)
    
    Pi = solve(SIGMA_uu) %*% SIGMA_um
    PSIinv = solve(chol(SIGMA_mu %*% Pi))
    B1 = SIGMA_um %*% PSIinv
    ### With r=1 in Lunsford 2016:8, PSIinv is the inverted scalar square root of "phi^2".
    
    # F-test for a weak instrument, from Lunsford (2016)
    SIGMA_ee = tcrossprod((m - t(Pi) %*% u))
    F_stats  = ((dim_T - dim_K)/dim_K) * (SIGMA_mm - SIGMA_ee) / SIGMA_ee  # statistic, from Lunsford 2016:15
    
  }else if( S2=="MR" ){
    # ... from Mertens, Raven 2013:1244, see Jentsch,Lunsford 2019:2659
    SIGMA_mm = tcrossprod(m)/dim_T
    SIGMA_um = tcrossprod(u, m)/dim_T
    SIGMA_mu = t(SIGMA_um)
    
    idx_s = 1:dim_S
    SIGMA_mu11 = SIGMA_mu[ idx_s,  idx_s, drop=FALSE]
    SIGMA_mu12 = SIGMA_mu[ idx_s, -idx_s, drop=FALSE]
    SIGMA_uu11 = SIGMA_uu[ idx_s,  idx_s, drop=FALSE]
    SIGMA_uu21 = SIGMA_uu[-idx_s,  idx_s, drop=FALSE]
    SIGMA_uu22 = SIGMA_uu[-idx_s, -idx_s, drop=FALSE]
    
    B21B11inv = t(solve(SIGMA_mu11) %*% SIGMA_mu12) # Eq.10
    M1_tmp = B21B11inv %*% t(SIGMA_uu21)
    bigZ   = SIGMA_uu22 - M1_tmp - t(M1_tmp) + B21B11inv %*% SIGMA_uu11 %*% t(B21B11inv)
    M2_tmp = SIGMA_uu21 - B21B11inv %*% SIGMA_uu11
    B12B12 = t(M2_tmp) %*% solve(bigZ) %*% M2_tmp
    B11B11 = SIGMA_uu11 - B12B12
    B22B22 = SIGMA_uu22 - B21B11inv %*% B11B11 %*% t(B21B11inv)
    B12B22inv = (t(SIGMA_uu21) - B11B11 %*% t(B21B11inv)) %*% solve(B22B22)
    
    M3_tmp = diag(nrow=dim_S) - B12B22inv %*% B21B11inv
    S1S1 = M3_tmp %*% B11B11 %*% t(M3_tmp)  # Eq.13
    B11 = solve(M3_tmp) %*% t(chol(S1S1))  # Eq.14
    B1 = rbind(B11, B21B11inv %*% B11)  # Eq.9+10
    
  }else{
    stop("Please select the S2-identification 'MR', 'JL', 'LT', or 'NQ'.")
  }
  
  # return result
  dimnames(B1) = list(names_k, names_s)
  result = list(B1=B1, m=m, F_stats=F_stats, udv=R.udv, Q=R.udv$Q)
  return(result)
}


# identify structural rotation with proxy variables / external instruments, from Empting et al. (2026)
aux_idSVD <- function(eps_u, eps_m, S2=c("MR", "JL", "LT", "NQ"), dim_K1){
  # define
  eps_u = aux_asDataMatrix(eps_u, "y")  # named matrix of whitened residuals is KxT
  eps_m = aux_asDataMatrix(eps_m, "m")  # named matrix of whitened proxies is LxT
  dim_T = ncol(eps_u)
  dim_K = nrow(eps_u)
  dim_L = nrow(eps_m)
  dim_K2 = dim_K - dim_K1
  idx_k1 = 0:dim_K1
  names_s = paste0("epsilon[ ", 1:dim_K, " ]")
  
  # moments and SVD
  SIGMA_em = tcrossprod(eps_u, eps_m)/dim_T  # proxy moment matrix
  ZEROS_2  = matrix(0, nrow=dim_K,  ncol=max(dim_K-dim_L, 0))  # for the orthogonal complement
  ZEROS_12 = matrix(0, nrow=dim_K1, ncol=dim_K2)  # for the off-diagonal padding
  R.udv    = svd(cbind(SIGMA_em, ZEROS_2), nv=dim_K1)
  
  # impose n_K1 zero-restrictions
  if( S2=="NQ" & dim_L==dim_K1 ){
    w11 = t(R.udv$v[idx_k1, idx_k1, drop=FALSE])
  }else if( S2=="NQ" ){
    v11 = R.udv$v[idx_k1, idx_k1, drop=FALSE]
    D1  = diag(R.udv$d[idx_k1], nrow=dim_K1)
    udv = svd(v11 %*% D1)
    w11 = udv$v %*% t(udv$u)
  }else if( S2=="LT" ){
    u11 = R.udv$u[idx_k1, idx_k1, drop=FALSE]
    w11 = aux_signLQ(u11)
  }else if( S2=="JL" ){
    v11 = R.udv$v[idx_k1, idx_k1, drop=FALSE]
    D1  = diag(R.udv$d[idx_k1], nrow=dim_K1)
    w11 = aux_signLQ(v11 %*% D1)
  }else{  # S2=="MR" for example
    w11 = diag(nrow=dim_K1)
  }
  
  # impose n_K2 zero-restrictions by "LT", signed on the main diagonal of Q
  u12 = R.udv$u[0:dim_K2, -idx_k1, drop=FALSE]
  u22 = R.udv$u[-idx_k1,  -idx_k1, drop=FALSE]
  w22 = aux_signLQ(u12, S=u22)
  
  # return result
  R.udv$w = cbind(rbind(w11, t(ZEROS_12)), rbind(ZEROS_12, w22))
  colnames(R.udv$w) = names_s
  R.udv$Q = R.udv$u %*% R.udv$w
  return(R.udv)
}


# impose restrictions of Mertens, Raven (2013) on Q=UV, from Empting et al. (2026), Appendix A
aux_NQ2MR <- function(B, dim_K1){
  # define the blocks of B
  idx_k1 = 0:dim_K1
  B11 = B[ idx_k1, idx_k1, drop=FALSE]
  B12 = B[ idx_k1,-idx_k1, drop=FALSE]
  B21 = B[-idx_k1, idx_k1, drop=FALSE]
  B22 = B[-idx_k1,-idx_k1, drop=FALSE]
  
  # restrict the Schur complement of B to be lower-triangular with positive diagonal elements
  S_1 = if(nrow(B22) == 0){ B11 }else{ B11 - B12 %*% solve(B22, B21) }
  w11 = aux_signLQ(S_1)
  
  # return result
  W = diag(nrow=nrow(B))
  W[idx_k1, idx_k1] = w11
  colnames(W) = colnames(B)
  return(W)
}


# ICA by Cramer-von Mises distance, see Herwartz, Ploedt (2016)
aux_cvmICA <- function(eps, dd, itermax = 500, steptol = 100, iter2 = 75){
  ### Note that the functions have been taken from svars and are further ...
  ### ... adjusted here to perform ICA on the pre-whitened shocks 'eps' directly! 
  ### SOURCES:
  ### https://github.com/cran/svars/blob/master/R/id.cvm.R
  ### https://github.com/cran/svars/blob/master/R/rotation.R
  
  # function to compute the rotation matrix
  rotmat <- function(theta) {
    ts_dim <- (1 + sqrt(1 + 8 * length(theta))) / 2
    combns <- combn(ts_dim, 2, simplify = FALSE)
    rotmat <- diag(nrow=ts_dim)

    for (i in seq_along(theta)) {
      tmp <- diag(nrow=ts_dim)
      tmp[combns[[i]][1], combns[[i]][1]] <-  cos(theta[i])
      tmp[combns[[i]][2], combns[[i]][2]] <-  cos(theta[i])
      tmp[combns[[i]][1], combns[[i]][2]] <- -sin(theta[i])
      tmp[combns[[i]][2], combns[[i]][1]] <-  sin(theta[i])
      rotmat <- rotmat %*% tmp
    }
    
    return(rotmat)
  }
  
  # function which rotates the structural errors and calculates their independence
  # ... like Maxand 2017:117, svars:::testlik uses the package "copula" by Hofert et al. 2015 / Kojadinovic,Yan (2010)
  testlik <- function(theta, eps, dd) {
    ser_low <- eps %*% t(rotmat(theta))  # t(rotmat) == solve(rotmat)
    ddtest  <- copula::indepTest(ser_low, dd)
    ddtest$global.statistic # * 10000000
  }
  
  ########### starting the computations ------------------------------------------------------------------------

  k <- ncol(eps)  # number of endogenous variables resp. residual time series
  lower <- rep(0, k * (k - 1) / 2)
  upper <- rep(pi, k * (k - 1) / 2)
  
  ## First step of optimization with DEoptim
  de_control <- list(itermax = itermax, steptol = steptol,  trace = FALSE)
  
  de_res <- DEoptim::DEoptim(testlik, lower = lower, upper = upper,  
                             control = de_control, eps = eps, dd = dd)
  
  ## Second step of optimization. Creating randomized starting angles around the optimized angles
  theta_rot <- matrix(rnorm(n = (k*(k-1)/2)*iter2, mean = de_res$optim$bestmem, sd = 0.3), (k*(k-1)/2), iter2)
  theta_rot <- cbind(theta_rot, de_res$optim$bestmem)
  # start vectors for iterative optimization approach
  startvec_list <- as.list(as.data.frame(theta_rot))
  
  erg_list <- lapply(X = startvec_list,
                     FUN = optim,
                     fn = testlik,
                     gr = NULL,
                     eps = eps,
                     dd = dd,
                     method = ifelse(k == 2, "Brent", "Nelder-Mead"),
                     lower = ifelse(k == 2, -.Machine$double.xmax/2, -Inf),
                     upper = ifelse(k == 2, .Machine$double.xmax/2, Inf),
                     control = list(maxit = 1000),
                     hessian = FALSE)
  
  
  # print log-likelihood values from local maxima
  logliks <- sapply(erg_list, "[[", "value")
  
  if(min(logliks) < de_res$optim$bestval){
    params <- sapply(erg_list, "[[", "par", simplify = FALSE)
    par_o  <- params[[which.min(logliks)]]
    logs   <- min(logliks)
    inc    <- 1
  }else{
    par_o  <- de_res$optim$bestmem
    logs   <- de_res$optim$bestval
    inc    <- 0
  }
  
  result <- list(theta = par_o,        # optimal rotation angles
                 W = t(rotmat(par_o)), # optimal rotation matrix (the estimated unmixing matrix)
                 inc = inc,         # whether the second optimization improved the estimation
                 test.stats = logs, # minimum test statistic obtained
                 iter1 = de_res$optim$iter,    # number of iterations of first optimization
                 test1 = de_res$optim$bestval, # minimum test statistic from first optimization
                 test2 = min(logliks) # minimum test statistic from second optimization
  )
  return(result)
}


# identify structural shocks in GRT by FIML estimation
aux_idGRT <- function(OMEGA, XI, dim_T, LR = NULL, SR = NULL, start = NULL, max.iter = 100, conv.crit = 1e-07, maxls = 1){
  ### scoring algorithm, from Amisano,Giannini 1997:57, ch.4.2 / Luetkepohl,Kraetzig 2004:173, ch.4.4.1
  ### SOURCE:
  ### https://github.com/cran/vars/blob/master/R/SVEC.R
  # define
  Sigma = OMEGA  # MLE covariance matrix of residuals
  Xi = XI  # long-run multiplier matrix
  obs = dim_T  # number of observations without presample
  K = nrow(OMEGA)  # number of endogenous variables
  names_k = rownames(OMEGA)
  names_s = colnames(SR)
  
  K2 <- K^2
  IK <- diag(K)
  IK2 <- diag(K2)
  Kkk <- diag(K2)[, c(sapply(1:K, function(i) seq(i, K2, K)))]
  
  ##
  ## S-Matrix for explicit form
  ##
  Lrres <- sum(is.na(LR))
  SRres <- sum(is.na(SR))
  R0 <- diag(K2)
  select <- c(apply(SR, c(1, 2), function(x) ifelse(identical(x, 0.0), TRUE, FALSE)))
  R.B <- R0[select, ]
  select <- c(apply(LR, c(1, 2), function(x) ifelse(identical(x, 0.0), TRUE, FALSE)))
  R.C1 <- R0[select, ]
  if(identical(nrow(R.C1), as.integer(0))){
    R.C1 <- matrix(0, nrow = K2, ncol = K2)
    nResC1 <- 0
  }
  R.C1 <- R.C1 %*% kronecker(IK, Xi)
  nResC1 <- qr(R.C1)$rank
  ##
  ## Setting up the R matrix (implicit form)
  ##
  if(identical(nrow(R.B), as.integer(0))){
    R <- R.C1
    nResB <- 0
  } else {
    R <- rbind(R.C1, R.B)
    nResB <- qr(R.B)$rank
  }
  ##
  ## Obtaining the S matrix and s vector (explicit form)
  ##
  Sb <- aux_oc(t(R))
  S <- rbind(matrix(0, nrow = K2, ncol = ncol(Sb)), Sb)
  l <- ncol(S)
  s <- c(c(diag(K)), rep(0, K2))
  ##
  ## Test of unidentification
  ##
  if((nResB + nResC1) < (K * (K - 1) / 2)){
    stop("The model is not identified. Use less free parameters.")
  }
  ##
  ## Test identification numerically
  ##
  ifelse(is.null(start), gamma <- start <- rnorm(l), gamma <- start)
  vecab <- S %*% gamma + s
  A <- matrix(vecab[1:K2], nrow = K, ncol = K)
  B <- matrix(vecab[(K2 + 1):(2 * K2)], nrow = K, ncol = K)
  v1 <- (IK2 + Kkk) %*% kronecker(t(solve(A) %*% B), solve(B))
  v2 <- -1 * (IK2 + Kkk) %*% kronecker(IK, solve(B))
  v <- cbind(v1, v2)
  idmat <- v %*% S
  ms <- t(v) %*% v
  auto <- eigen(ms)$values
  rni <- 0
  for (i in 1:l) {
    if (auto[i] < 1e-11)
      rni <- rni + 1
  }
  if (identical(rni, 0)) {
    if (identical(l, as.integer(K * (K + 1)/2))) {
      ident <- paste("The model is just identified.")
    } else {
      ident <- paste("The model is over identified.")
    }
  } else {
    ident <- paste("The model is unidentified. The non-identification rank is", rni, ".")
    stop(ident)
  }      
  ##
  ## Scoring Algorithm
  ##
  iters <- 0
  cvcrit <- conv.crit + 1
  while (cvcrit > conv.crit) {
    z <- gamma
    vecab <- S %*% gamma + s
    A <- matrix(vecab[1:K2], nrow = K, ncol = K)
    B <- matrix(vecab[(K2 + 1):(2 * K2)], nrow = K, ncol = K)
    Binv <- solve(B)
    Btinv <- solve(t(B))
    BinvA <- Binv %*% A
    
    infvecab.mat1 <- rbind(kronecker(solve(BinvA), Btinv), -1 * kronecker(IK, Btinv))
    infvecab.mat2 <- IK2 + Kkk
    infvecab.mat3 <- cbind(kronecker(t(solve(BinvA)), Binv), -1 * kronecker(IK, Binv))
    infvecab <- obs * (infvecab.mat1 %*% infvecab.mat2 %*% infvecab.mat3)
    infgamma <- t(S) %*% infvecab %*% S
    infgammainv <- solve(infgamma)
    
    scorevecBinvA <- obs * c(solve(t(BinvA))) - obs * (kronecker(Sigma, IK) %*% c(BinvA))
    scorevecAB.mat <- rbind(kronecker(IK, Btinv), -1 * kronecker(BinvA, Btinv))
    scorevecAB <- scorevecAB.mat %*% scorevecBinvA
    scoregamma <- t(S) %*% scorevecAB
    
    direction <- infgammainv %*% scoregamma
    length <- max(abs(direction))
    ifelse(length > maxls, lambda <- maxls/length, lambda <- 1)
    gamma <- gamma + lambda * direction
    iters <- iters + 1
    z <- z - gamma
    cvcrit <- max(abs(z))
    if (iters >= max.iter) {
      warning(paste("Convergence not achieved after", iters, "iterations. Convergence value:", cvcrit, "."))
      break
    }
  }
  iter <- iters - 1
  vecab <- S %*% gamma + s
  SR <- B
  ##
  ## Normalizing the sign of SR
  ##
  select <- which(diag(solve(A) %*% B) < 0)
  SR[, select] <- -1 * SR[, select]
  ##
  ## Computing LR and Sigma.U
  ##
  LR <- Xi %*% SR
  Sigma.U <- solve(A) %*% B %*% t(B) %*% t(solve(A))  # equals the reduced-form estimate if the SVAR is exactly identified
  dimnames(SR) = list(names_k, names_s)
  dimnames(LR) = list(names_k, names_s)
  dimnames(Sigma.U) = list(names_k, names_s)
  
  ##
  ## Setting near zero elements to zero
  ##
  idxnull <- which(abs(LR) < 0.1e-8, arr.ind = TRUE)
  LR[idxnull] <- 0.0
  idxnull <- which(abs(SR) < 0.1e-8, arr.ind = TRUE)
  SR[idxnull] <- 0.0
  
  # return result
  result = list(SR = SR, LR = LR, Sigma.U = Sigma.U * 100, Restrictions = c(nResC1, nResB), start = start, iter = iter)
  return(result)
}


