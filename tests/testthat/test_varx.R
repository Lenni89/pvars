

####################
###  UNIT TESTS  ###
####################
#
# For exported varx functions 
# and their auxiliary modules.


test_that("the residuals from the rank-restricted VAR and VECM representation 
          in VECM() are identical.", {
  # prepare data #
  data("PCAP")
  names_k = c("g", "k", "l", "y") # variable names
  data_i  = ts(PCAP[PCAP$id_i=="DEU", names_k], start=1960, end=2019, frequency=1)
  dim_p = 3
  
  for(type in c("Case1", "Case2", "Case3", "Case4", "Case5")){
    # estimate VECM and VAR residuals #
    R.vecm = VECM(y=data_i, dim_r=2, dim_p=dim_p, type=type)
    type_A = switch(type, "Case1"="none", "Case2"="const", "Case3"="const", "Case4"="both", "Case5"="both")
    Y = t(data_i)[ ,-(1:dim_p)]
    D = aux_dummy(type=type_A, dim_T=nrow(data_i))
    Z = aux_stack(data_i, dim_p=dim_p, D=D)
    if(type == "Case4"){ Z["trend", ] = Z["trend", ] - 1 }
    e = Y - R.vecm$A %*% Z
    
    # check identity of VECM and level-VAR representation #
    our_VAR  = e  # ... via residuals
    our_VECM = R.vecm$resid
    expect_equal(our_VAR, our_VECM)
  }
})


test_that("VECM() can reproduce basic examples from urca package 
          at different lag-orders and determinstic regressors D2.", {
  # prepare data #
  library("urca")
  data(denmark)
  sjd = denmark[, c("LRM", "LRY", "IBO", "IDE")]
  tolerant_u = 4e-07  # on residuals
  tolerant_a = 2e-06  # on LR-test stats and pvals for adjustment coefficients
  dim_K = ncol(sjd)
  dim_r = 2
  
  for(dim_p in 2:9){ # lag-order of the endogenous variables in the level-VAR
    # estimate VECM #
    R.cajo = urca::ca.jo(sjd, ecdet="trend", type="trace", K=dim_p, spec="transitory")
    ### R.vecm_urca = urca::cajorls(R.cajo, r=dim_r)
    R.vars = vars::vec2var(R.cajo, r=dim_r)
    R.vecm = VECM(y=sjd, dim_r=dim_r, dim_p=dim_p, type="Case4")
    
    # check via residuals #
    our_VECM  = R.vecm$resid
    urca_VECM = t(R.vars$resid)
    expect_equivalent(our_VECM, urca_VECM, tolerance = tolerant_u)
    
    # perform LR-test on weak exogeneity in conditional VECM #
    DA    = diag(dim_K)[ ,-4, drop=FALSE]  ### matrix(c(1,0,0,0), c(4,1))
    R.alr = urca::alrtest(R.cajo, r=dim_r, A=DA)
    #### R.cvec_urca = urca::cajorls(R.alr, r=dim_r)
    R.cvec = VECM(y=sjd[ ,1:3], dim_p=dim_p,
                 x=sjd[ ,   4], dim_q=dim_p, 
                 dim_r=dim_r, type="Case4")
    our_LRalpha = aux_LRrest(lambda      = R.vecm$RRR$lambda, 
                             lambda_rest = R.cvec$RRR$lambda, 
                             dim_r=dim_r, dim_T=R.cvec$dim_T, dim_L=R.cvec$dim_L)
    
    # check #
    urca_LRalpha = data.frame(stats=R.alr@teststat, pvals=R.alr@pval[1])
    expect_equal(our_LRalpha, urca_LRalpha, tolerance = tolerant_a)
  }
  
  # estimate VECM with customized regressors in D2 #
  R.D2  = t(aux_dummy(dim_T=nrow(sjd), t_blip=40, t_impulse=c(10+0:3, 30), t_shift=10, n.season=4))
  dim_p = 4
  
  R.cajo_D2 = urca::ca.jo(sjd, ecdet="trend", type="trace", K=dim_p, spec="transitory", dumvar=R.D2) # uses lm()
  ### R.vecm_urca_D2 = urca::cajorls(R.cajo_D2, r=dim_r)
  R.vars_D2 = vars::vec2var(R.cajo_D2, r=dim_r)
  R.vecm_D2 = VECM(y=sjd, dim_r=dim_r, dim_p=dim_p, type="Case4", D2=R.D2)
  
  # check via residuals #
  our_VECM_D2  = R.vecm_D2$resid
  urca_VECM_D2 = t(R.vars_D2$resid)
  expect_equivalent(our_VECM_D2, urca_VECM_D2, tolerance = tolerant_u)
})


test_that("auxVAR modules can reproduce OLS examples from vars package.", {
  # prepare data #
  library("vars")
  data(Canada)
  dim_T = nrow(Canada)  # number of observations including presample
  tolerant = 3.8e-07
  
  for(dim_p in 1:10){ # lag-order of the endogenous variables
    # stack operators (for "const") #
    Y = t(Canada)[c(1,3,4),-(1:dim_p)]
    D = aux_dummy(dim_T=dim_T, center=TRUE, n.season=4, type="const")
    Z = aux_stack(y=Canada[ ,c(1,3,4)], dim_p=dim_p, x=Canada[ , 2, drop=FALSE], dim_q=0, D=D)
    B = tcrossprod(Y, Z) %*% solve(tcrossprod(Z)) # Y%*%t(Z) %*% solve(Z%*%t(Z))
    
    R.def  = aux_stackOLS(Canada[, c(1,3,4)], dim_p=dim_p, x=Canada[,  2, drop=FALSE], dim_q=0, D=D)
    R.vars = vars::VAR(y=Canada[,-2], p=dim_p, type="const", exogen=Canada[ , 2, drop=FALSE], season=4)
    
    names_slp = paste0(rownames(Y), ".l1")
    theirs    = t(sapply(R.vars$varresult, FUN = function(k) k$coefficients))[ , c(names_slp, "prod")]
    expect_equal(B[ , names_slp], theirs[ , names_slp], tolerance = tolerant)  # Block of seasonal dummies is shifted.
    expect_identical(Y, R.def$Y)
    expect_identical(Z, R.def$Z)
    
    # multivariate OLS estimation using QR decomposition with column pivoting (for "both") #
    R.sens = aux_dummy(dim_T=dim_T, center=TRUE, n.season=4)
    R.varx = aux_VARX(y=Canada[, -2], dim_p=dim_p, x=Canada[ , 2, drop=FALSE], dim_q=0, type="both", D=R.sens, method="qr")
    R.vars = vars::VAR(y=Canada[,-2], p=dim_p, type="both", exogen=Canada[ , 2, drop=FALSE], season=4)
    
    our_resid  = R.varx$resid
    varx_resid = as.varx(R.vars)$resid  # coerce the VAR to varx object first
    vars_resid = t(resid(R.vars))
    expect_equivalent(our_resid, vars_resid)
    expect_equivalent(our_resid, varx_resid)
    
    our_slp  = R.varx$A[ , c(names_slp, "prod.l0"), drop=FALSE]
    vars_slp = t(sapply(R.vars$varresult, FUN = function(k) k$coefficients))[ , c(names_slp, "prod")]
    expect_equivalent(our_slp, vars_slp)  # Differing base season has no effect on slope coefficients!
  }
  
  # roots of the companion matrix #
  var.2c    = vars::VAR(Canada, p = dim_p, type = "const")
  vars_root = vars::roots(var.2c)
  R.compani = aux_var2companion(A=do.call("cbind", Acoef(var.2c)), dim_p=var.2c$p)
  our_root  = Mod(eigen( R.compani )$values)
  expect_equal(our_root, vars_root, tolerance = tolerant)
  
  # ordinary information criteria #
  R.criteria = NULL
  lag_max = 4
  Y = t(Canada)[,-(1:lag_max)]
  dim_T = ncol(Y)
  for(p in 1:lag_max){
    idx_t = (lag_max+1-p):nrow(Canada)
    R.var = aux_VARX(t(Canada)[ ,idx_t], dim_p=p, type="const", method="qr")
    R.ICp = aux_MIC(R.var$OMEGA, COEF=R.var$A, dim_T, lambda_H0=NULL)
    
    R.criteria = cbind(R.criteria, R.ICp)
    rm(idx_t, R.var, R.ICp)
  }
  R.selection = apply(R.criteria, MARGIN=1, function(x) which.min(x))
  
  ours   = list(selection=R.selection, criteria=R.criteria)
  theirs = vars::VARselect(y=Canada, lag.max=4, type="const")
  expect_equivalent(ours$selection, theirs$selection)
  expect_equivalent(ours$criteria,  theirs$criteria)
})


test_that("aux_vec2vma() provides the correct long-run multiplier matrix XI 
          at different cointegration ranks", {
  # prepare data #
  data("PCAP")
  names_k = c("g", "k", "l", "y") # variable names
  data_i  = ts(PCAP[PCAP$id_i=="DEU", names_k], start=1960, end=2019, frequency=1)
  x2 = urca::ca.jo(data_i, ecdet="trend", type="trace", K=2, spec="transitory")
  tolerant = 5e-05
  
  for(dim_r in 1:3){
    # define #
    x = VECM(y=data_i, dim_r=dim_r, dim_p=2, type="Case4")
    dim_K = x$dim_K  # number of endogenous variables
    dim_p = x$dim_p  # number of lags
    dim_T = x$dim_T  # number of observations without presample
    
    # calculate long-run multiplier matrix #
    alpha_oc = MASS::Null(x$VECM$alpha)
    beta_oc  = MASS::Null(x$beta[1:dim_K, ])
    our_XI   = aux_vec2vma(GAMMA=x$VECM$GAMMA, alpha_oc=alpha_oc, beta_oc=beta_oc, dim_p=dim_p, n.ahead="XI")
    
    # define #
    r = dim_r
    K <- x2@P
    P <- x2@lag - 1
    IK <- diag(K)
    vecr <- urca::cajorls(z = x2, r = r)  # uses lm()
    
    # calculate long-run multiplier matrix, code chunk from vars::SVEC() #
    Coef.vecr <- coef(vecr$rlm)
    alpha <- t(Coef.vecr[1:r, ])
    ifelse(r == 1, alpha.orth <- Null(t(alpha)), alpha.orth <- Null(alpha))
    beta <- vecr$beta[1:K, ]
    beta.orth <- Null(beta)
    Gamma.array <- array(c(t(tail(coef(vecr$rlm), K*P))), c(K, K, P))
    Gamma <- apply(Gamma.array, c(1, 2), sum)
    vars_Xi <- beta.orth %*% solve(t(alpha.orth) %*% (IK - Gamma) %*% beta.orth) %*% t(alpha.orth)  
    
    # check #
    expect_equivalent(our_XI, vars_Xi, tolerance = tolerant)
  }
})


