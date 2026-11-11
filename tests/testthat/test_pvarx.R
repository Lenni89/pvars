

####################
###  UNIT TESTS  ###
####################
#
# For exported "pvarx" functions 
# and their subordinated modules.


test_that("pvarx.VAR() and pvarx.VEC() return identical MG estimates under full-rank r=K", {
  # prepare data #
  data("ERPT")
  names_k = c("lpm5", "lfp5", "llcusd")  # variable names for "Chemicals and related products"
  names_i = levels(ERPT$id_i)[c(1,6,2,5,4,3,7)]  # ordered country names
  L.data  = sapply(names_i, FUN=function(i) ts(ERPT[ERPT$id_i==i, names_k], start=c(1995, 1), frequency=12), simplify=FALSE)
  
  # estimate #
  dim_T = nrow(L.data[[1]])
  dim_K = length(names_k)
  R.D   = aux_dummy(dim_T, t_blip=40, t_impulse=c(10+0:3, 30), t_shift=10, n.season=4)
  R.lags = c(3, 3, 3, 4, 3, 3, 3); names(R.lags)=names_i
  R.pvar = pvarx.VAR(L.data, lags=R.lags, type="none", D=R.D)
  R.pvec = pvarx.VEC(L.data, lags=R.lags, type="Case1", dim_r=dim_K, D2=R.D)
  R.mgA  = aux_MG(R.pvec$L.varx, idx_par="A")  # This MG is only valid under r=K, i.e the full-rank case!
  
  # get eigenvalues of the VAR systems #
  #R.pvar$eigenvals = eigen(pvars:::aux_var2companion(R.pvec$MG_A$mean, dim=2))$values
  #R.pvec$eigenvals = eigen(pvars:::aux_var2companion(R.pvec$A, dim=2))$values
  #R.mgA$eigenvals  = eigen(pvars:::aux_var2companion(R.mgA$mean, dim=2))$values
  
  # check #
  expect_equal(R.pvar$A, R.pvec$A)
  expect_equal(R.pvar$A, R.mgA$mean)
})


test_that("aux_pvec() and VECM() return identical RRR estimates when aux_2StepBR() and SUR-GLS are deactivated", {
  # prepare data #
  data("ERPT")
  names_k = c("lpm5", "lfp5", "llcusd")  # variable names for "Chemicals and related products"
  names_i = levels(ERPT$id_i)[c(1,6,2,5,4,3,7)]  # ordered country names
  L.data  = sapply(names_i, FUN=function(i) ts(ERPT[ERPT$id_i==i, names_k], start=c(1995, 1), frequency=12), simplify=FALSE)
  
  # arguments #
  R.lags = c(3, 3, 3, 4, 3, 3, 3); names(R.lags)=names_i
  dim_r = 2
  type  = "Case4"
  t_0 = as.t_D(list())  # aux_stackRRR() is internal and uses t_D that is already checked.
  
  # estimate #
  L.var = sapply(names_i, FUN=function(i) VECM(L.data[[i]], dim_p=R.lags[i], dim_r=dim_r, type=type), simplify=FALSE)
  L.def = sapply(names_i, FUN=function(i) aux_stackRRR(L.data[[i]], dim_p=R.lags[i], type=type, t_D1=t_0, t_D2=t_0), simplify=FALSE)
  R.est = aux_pvec(L.def, dim_r=dim_r, n.factors=FALSE, n.iterations=FALSE)
  
  # check #
  for(i in names_i){  # not considered:
    R.est$L.varx[[i]]$RRR = L.var[[i]]$RRR = NULL  # via_R0_R1=TRUE/FALSE
    R.est$L.varx[[i]]$RRR = L.var[[i]]$RRR = NULL
    L.var[[i]]$args_varx = NULL
    L.var[[i]]$PARTIAL   = NULL
    L.var[[i]]$MARGINAL  = NULL
  }
  expect_equal(ignore_attr=TRUE, R.est$L.varx, L.var)
})


test_that("aux_pvar() and aux_VARX() return identical OLS estimates when SUR-GLS is deactivated", {
  # prepare data #
  data("ERPT")
  names_k = c("lpm5", "lfp5", "llcusd")  # variable names for "Chemicals and related products"
  names_i = levels(ERPT$id_i)[c(1,6,2,5,4,3,7)]  # ordered country names
  L.data  = sapply(names_i, FUN=function(i) ts(ERPT[ERPT$id_i==i, names_k], start=c(1995, 1), frequency=12), simplify=FALSE)
  
  # estimate #
  R.lags = c(3, 3, 3, 4, 3, 3, 3); names(R.lags)=names_i
  type  = "both"
  L.var = sapply(names_i, FUN=function(i) aux_VARX(L.data[[i]], dim_p=R.lags[i], type=type), simplify=FALSE)
  L.def = sapply(names_i, FUN=function(i) aux_stackOLS(L.data[[i]], dim_p=R.lags[i], type=type), simplify=FALSE)
  R.est = aux_pvar(L.def, n.factors=FALSE, n.iterations=FALSE)
  
  # check #
  for(i in names_i){ R.est$L.varx[[i]]$B = NULL }  # not considered
  expect_equivalent(R.est$L.varx, L.var)
})


test_that("as.pvarx() provides the same 'L.varx' from 'pvarx' and list of 'vars' objects", {
  tolerant_u = 25e-07  # vars::VAR() is based on equation-wise lm().
  tolerant_A = 38e-09
  
  # prepare data #
  data("PCAP")
  names_k = c("g", "k", "l", "y")  # variable names
  names_i = levels(PCAP$id_i)      # country names 
  L.data  = sapply(names_i, FUN=function(i) ts(PCAP[PCAP$id_i==i, names_k], start=1960, end=2019, frequency=1), simplify=FALSE)
  D   = data.frame(trend=1:60)  # customized linear trend, which is kept the same in both functions
  L.D = sapply(names_i, FUN=function(i) D, simplify=FALSE)
  
  # calculate #
  L.vars  = lapply(L.data, FUN=function(x) vars::VAR(x, p=2, type="const", exogen=D))
  R.pvarx = pvarx.VAR(L.data, lags=2, type="const", D=D)
  
  # compare coefficients and residuals #
  L.varx1 = as.pvarx(R.pvarx)$L.varx
  L.varx2 = as.pvarx(L.vars)$L.varx
  i_pvarx = unname(L.varx1[[2]]$resid)
  i_vars  = unname(L.varx2[[2]]$resid)
  expect_equal(i_pvarx, i_vars, tolerance=tolerant_u)
  expect_equal(L.varx1[[3]]$A, L.varx2[[3]]$A, tolerance=tolerant_A)
  expect_equal(R.pvarx$A, as.pvarx(L.vars)$A, tolerance=tolerant_A)
  expect_equal(R.pvarx$A, pvarx.VAR(L.data, lags=2, type="const", D=L.D)$A)
  ### Note that type="trend" would have different starting points 
  ### such that the coefficients of the deterministic term are not equal.
})


test_that("as.pvarx() stops in case of heterogeneous arguments", {
  # prepare data #
  data("PCAP")
  names_k = c("g", "k", "l", "y")  # variable names
  names_i = levels(PCAP$id_i)      # country names 
  L.data  = sapply(names_i, FUN=function(i) ts(PCAP[PCAP$id_i==i, names_k], start=1960, end=2019, frequency=1), simplify=FALSE)
  
  # calculate #
  L.vars_K = L.vars_r = L.svars = lapply(L.data, FUN=function(x) vars::VAR(x, p=2, type="const"))
  L.vars_K[[2]] = vars::VAR(L.data[[2]][, 1:3], p=2, type="const")
  L.vars_r[[2]] = VECM(L.data[[2]], dim_p=2, dim_r=1, type="Case3")
  L.svars[[2]]  = svars::id.chol(L.svars[[2]])
  
  # test the stops #
  expect_error(as.pvarx(L.vars_K),
               "The number of variables 'K' must be the same for all individuals.")
  expect_error(as.pvarx(L.vars_r),
               "The cointegration rank restriction 'r' must be the same for all individuals.")
  expect_error(as.pvarx(L.svars),
               "The identification procedure must be the same for all individuals.")
})


