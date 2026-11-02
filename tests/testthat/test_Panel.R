

####################
###  UNIT TESTS  ###
####################
#
# For panel functions 
# and their auxiliary modules.


test_that("ptest.STATSbar() can reproduce Silvestre,Surdeanu 2011:27, Ch.6,Tab.11(A)", {
  # individual test statistics #
  MSB_stats = rbind(
    m3=c(0.008, 0.009, 0.005, 0.009, 0.004, 0.011, 0.004, 0.004, 0.014,  # for r_H0=0
         0.004, 0.005, 0.007, 0.013, 0.004, 0.010, 0.008, 0.013, 0.017, 0.006),
    m1=c(0.351, 0.121, 0.465, 0.212, 0.220, 0.097, 0.267, 0.108, 0.154,  # for r_H0=2
         0.747, 0.408, 0.515, 0.282, 0.305, 0.105, 0.687, 0.266, 0.275, 0.105))
  
  # combination in panel test #
  moments = coint_moments[["MSB_trend50"]][c(3,1), ]
  ours    = ptest.STATSbar(STATS=MSB_stats, moments=moments)
  ours$pvals = 1 - ours$pvals  # MSB-test is left-tailed based on normal distribution!
  theirs  = cbind(m3=-12.38, m2=3.40, m1=4.26)
  expect_equal(t(ours$stats), theirs[ , -2, drop=FALSE], tolerance = 0.005)
})


test_that("ptest.METApval() can reproduce Oersal,Arsova 2017:68, Ch.5,Tab.5", {
  # individual p-values #
  SL_pvals = rbind(
    r1=c(0.829, 0.404, 0.878, 0.565, 0.478, 0.621, 0.776, 0.897, 0.284,  # for r_H0=1
         0.736, 0.837, 0.351, 0.417, 0.049, 0.975, 0.998, 0.845, 0.075, 0.153),
    r3=c(0.256, 0.940, 0.557, 0.900, 0.442, 0.828, 0.563, 0.612, 0.906,  # for r_H0=3
         0.543, 0.374, 0.875, 0.694, 0.907, 0.147, 0.851, 0.609, 0.958, 0.338))
  
  # combination in panel test #
  ours   = ptest.METApval(PVALS=SL_pvals)
  theirs = cbind(r1=c(P=NA, Pm=-0.98, Z=1.49),
                 r3=c(P=NA, Pm=-2.02, Z=2.12))
  expect_equal(t(ours$stats)[-1, ], theirs[-1, ], tolerance = 0.005)
})


test_that("ptest.CAIN() can reproduce Oersal,Arsova 2016:22, Tab.7/8", {
  ### Oersal,Arsova 2016:22, Tab.7.0/8 ###
  TSL0_pvals = rbind(
    r0=c(0.23, 0.36, 0.22, 0.48, 0.02, 0.52, 0.76),
    r1=c(0.69, 0.93, 0.96, 0.63, 0.62, 0.82, 0.86))
  rho_eps = c(r0=0.69, r1=0.69)
  
  ours   = ptest.CAIN(TSL0_pvals, dim_K=3, rho_eps=rho_eps, r_H0=0:1)
  theirs = c(r0=-0.85, r1=NA, r2=NA)
  expect_equivalent(t(ours$stats)["CAIN", 1], theirs[1], tolerance = 0.015)
  
  ### Oersal,Arsova 2016:22, Tab.7.5/8 ###
  TSL5_pvals = rbind(
    r0=c(0.02, 0.04, 0.02, 0.01, 0.04, 0.03, 0.15),
    r1=c(0.25, 0.40, 0.03, 0.06, 0.04, 0.10, 0.60),
    r2=c(0.43, 0.22, 0.67, 0.89, 0.79, 0.41, 0.73))
  
  ours   = ptest.CAIN(TSL5_pvals, dim_K=3, rho_eps=0.63, r_H0=0:2)
  theirs = cbind(
    kappa_1 = c(r0=0.02, r1=0.07, r2=0.63), 
    kappa_2 = c(r0=0.02, r1=0.06, r2=0.63), 
    CAIN    = c(r0=0.00, r1=0.02, r2=0.70))
  expect_equal(ours$pvals, theirs, tolerance = 0.008)
})


test_that("aux_2StepBR() and aux_RRR() return identical cointegration vectors for a duplicated individual", {
  # prepare data #
  data("ERPT")
  names_k = c("lpm5", "lfp5", "llcusd")  # variable names for "Chemicals and related products"
  names_i = levels(ERPT$id_i)[c(1,6,2,5,4,3,7)]  # ordered country names
  L.data  = sapply(names_i, FUN=function(i) ts(ERPT[ERPT$id_i==i, names_k], start=c(1995, 1), frequency=12), simplify=FALSE)
  
  # estimate #
  L.def = lapply(c(1,1), FUN=function(i) aux_stackRRR(L.data[[i]], dim_p=2, type="Case4")) ### , t_D1=list(t_break=89)
  L.RRR = lapply(L.def,  FUN=function(x) aux_RRR(x$Z0, x$Z1, x$Z2, via_R0_R1=TRUE))
  R.hom = aux_2StepBR(L.RRR, r_H0=0:(4-1), idx_pool=1:4, n.iterations=0)  # all coefficients are homogeneous
  R.trd = aux_2StepBR(L.RRR, r_H0=0:(4-1), idx_pool=1:3, n.iterations=0)  # heterogeneous trend coefficients
  R.het = aux_2StepBR(L.RRR, r_H0=0:(4-1), idx_pool=0,   n.iterations=0)  # all coefficients are heterogeneous
  
  # test # 
  R.pcbr = pcoint.BR(L.data[c(1,1)], lags=2, type="Case4")
  
  # check #
  L.beta = lapply(0:3, FUN=function(r) aux_beta(V=L.RRR[[1]]$V, dim_r=r, normalize="natural"))
  names(L.beta) = paste0("r", 0:3)
  expect_equal(R.hom$L.beta[1, ], L.beta)
  expect_equal(R.trd$L.beta[1, ], L.beta)
  expect_equal(R.het$L.beta[1, ], L.beta)
  expect_equal(R.trd$L.beta[1, 1:3], R.pcbr$beta_H0[1, ])  # exclude r=K
})


test_that("speci.factors() and aux_ComFact() based on EVD and SVD are identical and can reproduce Oersal,Arsova 2017:67, Ch.5", {
  # prepare data #
  data("MERM")
  names_k = c("s", "m", "y", "p") # variable names
  names_i = levels(MERM$id_i)     # country names
  L.data = sapply(names_i, FUN=function(i) ts(MERM[MERM$id_i==i, names_k], start=c(1995, 1), frequency=12), simplify=FALSE)
  
  # perform PCA #
  Y = do.call("cbind", L.data)  # data matrix of dimension (T x K*N)
  dim_T   = nrow(Y)
  PCA_evd = aux_ComFact(X=Y, SVD=FALSE, n.factors=8)
  PCA_svd = aux_ComFact(X=Y, SVD=TRUE,  n.factors=8)
  
  # compare PCA based on eigenvalue decomposition and singular value decomposition #
  dim_NK = min(dim(Y))  # maximum number of singular values
  expect_equal(PCA_evd$eigenvals[1:dim_NK], PCA_svd$eigenvals[1:dim_NK])
  expect_equal(PCA_evd$eit2, PCA_svd$eit2)
  expect_equal(abs(PCA_evd$Ft), abs(PCA_svd$Ft))
  
  # Onatski (2010) ED criterion for Oersal,Arsova 2017:67, Ch.5 #
  speci_lvl = speci.factors(L.data, k_max=20, n.iterations=6)
  expect_equivalent(speci_lvl$selection[[2]]["ED"], 8)
  
  # identity between speci.factors() with selected data transformation and aux_ComFact() #
  speci_diff = speci.factors(L.data, k_max=20, n.iterations=6, differenced=TRUE, centered=TRUE, scaled=TRUE, n.factors=8)
  expect_equal(speci_diff$eigenvals$value, PCA_svd$eigenvals/(dim_T-1))
  expect_equal(speci_diff$Ft, PCA_svd$Ft)
})


