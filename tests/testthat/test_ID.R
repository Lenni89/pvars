

####################
###  UNIT TESTS  ###
####################
#
# For exported "id" functions, 
# their auxiliary modules, 
# and residual diagnostics.


test_that("the residual diagnostics reproduce example from vars package", {
  tolerant = 0.0001  # for comparing p-values
  
  # VAR model #
  library("vars")
  data(Canada)
  var.2c <- vars::VAR(Canada, p = 2, type = "const") ### p=3 generates computationally singular matrices in test.serial(lag.h>2)
  
  # vars package #
  R.norm = vars::normality.test(var.2c, multivariate.only = FALSE)
  R.arch = vars::arch.test(var.2c, lags.multi = 5, lags.single = 5, multivariate.only = FALSE)
  R.seBG = vars::serial.test(var.2c, lags.bg = 2, type = "BG")
  R.seES = vars::serial.test(var.2c, lags.bg = 2, type = "ES")
  
  vars_norm = sapply(R.norm$jb.mul, function(x) x$p.value)
  vars_arch = R.arch$arch.mul$p.value
  vars_seri = c(R.seBG$serial$p.value, R.seES$serial$p.value)
  
  # pvars package: p-values from multivariate tests #
  resids = resid(var.2c)
  our_norm = test.normality(u=resids)
  our_arch = test.arch(u=resids, lag.h=5)
  our_seri = test.serial(u=resids, lag.h=2, x=var.2c$datamat[ ,-(1:ncol(var.2c$y))])
  
  # check #
  expect_equivalent(our_norm, vars_norm, tolerance=tolerant)
  expect_equivalent(our_arch, vars_arch, tolerance=tolerant)
  expect_equivalent(our_seri, vars_seri, tolerance=tolerant)
  
  # compare test statistics of exported functions of vars and pvars #
  our_norm  = rboot.normality(var.2c, n.boot=5)$stats[ , "MULTI"]
  vars_norm = sapply(R.norm$jb.mul, function(x) x$statistic)
  expect_equivalent(our_norm, vars_norm)
  
  # pvars package: p-values from univariate tests #
  # our_unorm = apply(resids, MARGIN=2, FUN=function(x) test.normality(u=x))
  # our_uarch = apply(resids, MARGIN=2, FUN=function(x) test.arch(u=x, lag.h=5))
  # our_useri = apply(resids, MARGIN=2, FUN=function(x) test.serial(u=x, lag.h=2, x=var.2c$datamat[ ,-(1:ncol(var.2c$y))]))
})


test_that("id.iv() returns the same B for all S2=c('MR', 'JL', 'NQ') if there is L=1 proxy", {
  # prepare data #
  data("PCIT")
  names_k = c("APITR", "ACITR", "PITB", "CITB", "GOV", "RGDP", "DEBT")
  names_l = c("m_PI")  # name of single proxy
  
  # estimate and identify #
  R.vars = vars::VAR(PCIT[ , names_k], p=4, type="const")
  R.ivMR = id.iv(R.vars, iv=PCIT[-(1:4), names_l], cov_u="OMEGA", S2="MR")
  R.ivJL = id.iv(R.vars, iv=PCIT[-(1:4), names_l], cov_u="OMEGA", S2="JL")
  R.ivNQ = id.iv(R.vars, iv=PCIT[-(1:4), names_l], cov_u="OMEGA", S2="NQ")
  
  # check #
  expect_equal(R.ivJL$B, R.ivMR$B)
  expect_equal(R.ivJL$B, R.ivNQ$B[ , 1, drop=FALSE])
})


test_that("id.iv(S2='MR') can reproduce the point estimates of 'PCIT' in Jentsch,Lunsford 2019:2668, Ch.III", {
  ### Bootstrap results in the example of id.iv() are congruent with Fig.1 and 2. 
  tolerant = 0.1
  
  # prepare data #
  data("PCIT")
  R.norm = function(B) B / matrix(-diag(B), nrow(B), ncol(B), byrow=TRUE)
  
  # estimate and identify under ordering "blue" in Fig.1 and 2 #
  names_k1 = c("APITR", "ACITR", "PITB", "CITB", "GOV", "RGDP", "DEBT")
  names_l1 = c("m_PI", "m_CI")  # proxy names
  R.var1 = vars::VAR(PCIT[ , names_k1], p=4, type="const")
  R.iv1  = id.iv(R.var1, iv=PCIT[-(1:4), names_l1], S2="MR", cov_u="OMEGA")
  R.irf1 = irf.varx(R.iv1, n.ahead=20, normf=R.norm)
  
  # check "blue" #
  our_B   = unlist(R.irf1$irf[1, -1])
  their_B = c(#1 | #2 Figure
    APITR = c(-1, 0.05),
    ACITR = c(0.6,  -1),
    PITB  = c(0.6,  NA), 
    CITB  = c(NA,  3.2), 
    GOV   = c(0.1, 0.6), 
    RGDP  = c(1.3, 0.4),
    DEBT  = c(NA,   NA))
  idx_na = is.na(their_B)
  expect_equivalent(our_B[!idx_na], their_B[!idx_na], tolerance=tolerant)
  ### cbind(unname(our_B[!idx_na]), unname(their_B[!idx_na]))
  
  # estimate and identify under ordering "red" in Fig.1 and 2 #
  names_k2 = c("ACITR", "APITR", "PITB", "CITB", "GOV", "RGDP", "DEBT")
  names_l2 = c("m_CI", "m_PI")
  R.var2 = vars::VAR(PCIT[ , names_k2], p=4, type="const")
  R.iv2  = id.iv(R.var2, iv=PCIT[-(1:4), names_l2], S2="MR", cov_u="OMEGA")
  R.irf2 = irf.varx(R.iv2, n.ahead=20, normf=R.norm)
  
  # check "red" #
  our_B   = unlist(R.irf2$irf[1, -1])
  their_B = c(#2 | #1 Figure
    ACITR = c(-1,  0.03),
    APITR = c(0.07, -1),
    PITB  = c(NA,  0.7), 
    CITB  = c(3.2,  NA), 
    GOV   = c(0.6, 0.2), 
    RGDP  = c(0.3, 1.4),
    DEBT  = c(NA,   NA))
  idx_na = is.na(their_B)
  expect_equivalent(our_B[!idx_na], their_B[!idx_na], tolerance=tolerant)
  ### cbind(unname(our_B[!idx_na]), unname(their_B[!idx_na]))
})


test_that("pid.iv(c('pool','indiv')) can reproduce the id.iv(S2='NQ') example for 'PCIT' under duplicated N=1", {
  # prepare data #
  data("PCIT")
  names_k = c("APITR", "ACITR", "PITB", "CITB", "GOV", "RGDP", "DEBT")
  names_l = c("m_PI", "m_CI")  # proxy names
  L.iv = list(USA=PCIT[-(1:4), names_l], USB=PCIT[-(1:4), names_l])
  
  # estimate and identify #
  R.vars = vars::VAR(PCIT[ , names_k], p=4, type="const")
  L.vars = list(USA=R.vars, USB=R.vars )  # list of N=2 duplication of a single individual
  R.iv1  = id.iv(L.vars[[1]], iv=L.iv[[1]], S2="NQ")
  R.iv2  = pid.iv(L.vars, S2="NQ", cov_u="OMEGA", combine="indiv", iv=L.iv[[1]])  # duplicates the common proxy
  R.iv3  = pid.iv(L.vars, S2="NQ", cov_u="OMEGA", combine="pool", iv=L.iv)
  R.iv4  = as.pvarx(list(USA=R.iv1, USB=R.iv1))
  
  # normalized IRF #
  R.norm = function(B) B / matrix(-diag(B), nrow(B), ncol(B), byrow=TRUE)
  R.irf1 = irf(R.iv1, normf=R.norm)
  R.irf2 = irf(R.iv2, normf=R.norm)
  R.irf3 = irf(R.iv3, normf=R.norm)
  R.irf4 = irf(R.iv4, normf=R.norm)
  
  # check #
  expect_equal(R.iv2$args_pid$L.iv, R.iv3$args_pid$L.iv)
  expect_equal(R.iv2$args_pid$L.iv, R.iv4$args_pid$L.iv)
  
  expect_equal(R.iv1$B, R.iv2$B)
  expect_equal(R.iv1$B, R.iv3$B)
  expect_equal(R.iv1$B, R.iv4$B)
  
  expect_equal(R.irf1, R.irf2)
  expect_equal(R.irf1, R.irf3)
  expect_equal(R.irf1, R.irf4)
})


test_that("aux_idSVD() can recover 'MR' and 'JL' from 'NQ' respectively", {
  # prepare data and estimate #
  data("PCIT")
  names_k = c("APITR", "ACITR", "PITB", "CITB", "GOV", "RGDP", "DEBT")
  names_l = c("m_PI", "m_CI")  # proxy names
  dim_K1  = length(names_l)  # number of identified shocks 
  dim_p   = 4  # lag-order
  R.vars  = vars::VAR(PCIT[ , names_k], p=dim_p, type="const")
  
  # identify under "MR" and check #
  R.idMR = id.iv(R.vars, iv=PCIT[-(1:dim_p), names_l], S2="MR", cov_u="OMEGA")
  R.idNQ = id.iv(R.vars, iv=PCIT[-(1:dim_p), names_l], S2="MR", cov_u="OMEGA", dim_K1=dim_K1)
  expect_equal(R.idMR$B, R.idNQ$B[ , 1:dim_K1])
  
  # identify under "JL" and check #
  R.idJL = id.iv(R.vars, iv=PCIT[-(1:dim_p), names_l], S2="JL", cov_u="OMEGA")
  R.idNQ = id.iv(R.vars, iv=PCIT[-(1:dim_p), names_l], S2="JL", cov_u="OMEGA", dim_K1=dim_K1)
  expect_equal(R.idJL$B, R.idNQ$B[ , 1:dim_K1])
})


test_that("the impact matrix from id.grt() fulfills BB'=OMEGA in a just-identified SVECM with period-pecific deterministic term", {
  # prepare data #
  data("PCAP")
  names_k = c("g", "k", "l", "y")  # variable names
  data_i  = ts(PCAP[PCAP$id_i=="DEU", names_k], start=1960, end=2019, frequency=1)
  
  # specifications #
  LR = matrix(NA, nrow=4, ncol=4); LR[  , 1:2] = 0
  SR = matrix(NA, nrow=4, ncol=4); SR[ 1, 2] = 0; SR[ 3, 4] = 0
  startv = c(-0.363, -1.448, -0.365, 0.696, -0.090, 0.957, 1.044, 0.083, -0.169, 0.044)
  
  # estimate and identify #
  R.vecm = VECM(y=data_i, dim_p=2, dim_r=2, type="Case4", t_D1=list(t_break=c(23, 49)))
  R.svec = id.grt(R.vecm, LR=LR, SR=SR, start=startv)
  
  # check #
  B = unname(R.svec$B)
  vecm_OMEGA = unname(R.vecm$VECM$OMEGA)
  svec_OMEGA = unname(R.svec$OMEGA)
  expect_equivalent(B%*%t(B), vecm_OMEGA)
  expect_equivalent(B%*%t(B), svec_OMEGA)
})


test_that("VECM(), id.grt(), and pid.grt() can reproduce the basic example from vars package under r=1 and N=1", {
  # prepare data #
  library(vars)
  data(Canada)
  names_k = c("prod", "e", "U", "rw")
  CAN = Canada[ , names_k]
  
  # joint specifications #
  dim_p = 3  # lag-order
  dim_r = 1  # cointegration rank
  SR = matrix(NA, nrow=4, ncol=4, dimnames=list(names_k, names_k))
  SR[4, 2] = 0
  LR = matrix(NA, nrow=4, ncol=4, dimnames=list(names_k, names_k))
  LR[1, 2:4] = 0
  LR[2:4, 4] = 0
  startv = c(-0.150, 1.011, -1.457, 0.044, 0.085, 0.467, 0.053, -1.073, 2.073, -1.888)
  
  # pvars: estimate and identify #
  set.seed(4679)
  R.vecm = VECM(y=CAN, dim_p=dim_p, dim_r=dim_r, type="Case4")
  R.grt  = id.grt(R.vecm, LR=LR, SR=SR, start=startv)
  
  set.seed(4679)
  L.vecm = list(A1=R.vecm)  # list of N=1 individual
  R.pgrt = pid.grt(L.vecm, LR=LR, SR=SR, start=startv)
  
  # urca and vars: estimate and identify #
  set.seed(4679)
  R.cajo = urca::ca.jo(CAN, type="trace", ecdet="trend", K=dim_p, spec="transitory")
  R.svec = vars::SVEC(R.cajo, LR=LR, SR=SR, r=dim_r, start=startv, lrtest=FALSE, boot=FALSE)
  
  # check #
  expect_equal(R.pgrt$B, as.varx(R.svec)$B)
  expect_equal(R.grt$B,  as.varx(R.svec)$B)
  expect_equal(R.grt$SR, R.svec$SR)
  expect_equal(R.grt$LR, R.svec$LR)
})


test_that("aux_idGRT() is correctly grafted under r=2", {
  # prepare data #
  data("PCAP")
  names_k = c("g", "k", "l", "y")  # variable names
  data_i  = ts(PCAP[PCAP$id_i=="AUT", names_k], start=1960, end=2019, frequency=1)
  
  # joint specifications #
  dim_p = 2  # lag-order
  dim_r = 2  # cointegration rank
  LR = matrix(NA, nrow=4, ncol=4, dimnames=list(names_k, names_k)); LR[  , 1:2] = 0
  SR = matrix(NA, nrow=4, ncol=4, dimnames=list(names_k, names_k)); SR[ 1, 2] = 0; SR[ 3, 4] = 0
  startv  = c(-0.363, -1.448, -0.365, 0.696, -0.090, 0.957, 1.044, 0.083, -0.169, 0.044)
  n.ahead = 10
  n.boot  = 10
  
  # pvars: estimate, identify, and bootstrap #
  set.seed(8349)
  R.vecm = VECM(y=data_i, dim_p=dim_p, dim_r=dim_r, type="Case4")
  R.grt  = id.grt(R.vecm, LR=LR, SR=SR, start=startv)
  R.boot = sboot.mb(R.grt, b.length=1, n.ahead=n.ahead, n.boot=n.boot, n.cores=1)
  
  # urca and vars: estimate, identify, and bootstrap #
  set.seed(8349)
  R.cajo = urca::ca.jo(data_i, type="eigen", ecdet="trend", K=dim_p, spec="transitory")
  R.svec = vars::SVEC(R.cajo, LR=LR, SR=SR, r=dim_r, start=startv, lrtest=FALSE, boot=FALSE)
  R.irf  = vars::irf(R.svec, n.ahead=n.ahead, boot=TRUE, ci=0.9, runs=n.boot)
  
  # align confidence bands, from svars::plot.sboot #
  kk <- ncol(R.boot$true$irf)
  intervals <- array(0, dim = c(n.ahead+1, kk, n.boot))
  for(i in 1:n.boot){
    intervals[,,i] <- as.matrix(R.boot$bootstrap[[i]]$irf)
  }
  R.bootLW <- array(0, dim = c(n.ahead, kk, 1))
  R.bootUP <- array(0, dim = c(n.ahead, kk, 1))
  for(i in 1:n.ahead){
    for(j in 1:kk){
      R.bootLW[i,j, ] <- quantile(intervals[i,j, ], probs = 0.05)
      R.bootUP[i,j, ] <- quantile(intervals[i,j, ], probs = 0.95)
    }
  }
  
  # check #
  expect_equal(R.grt$SR, R.svec$SR)
  expect_equal(R.grt$LR, R.svec$LR)
  #expect_equal(R.bootLW, R.irf$Lower)
  #expect_equal(R.bootUP, R.irf$Upper)
})


test_that("aux_sico() returns the correct permutation with names", {
  ### This is rather a reminder not to change the output. ###
  # define #
  names_s = paste0("eps_", 1:10)
  L.mat = list(A=cbind(2:3, -c(2:3)),
               B=cbind(-c(3:5), c(5, 0, 1), c(1:3), -c(4, 2, 1)),
               C=diag(5))
  
  for(mat_B in L.mat){
    dim_K = nrow(mat_B)
    dim_S = ncol(mat_B)
    
    # apply permutation matrix and its alternative #
    R.sico = aux_sico(mat_B, names_s=names_s[1:dim_S])
    mat_si = matrix(1, nrow=dim_K, ncol=dim_S)
    mat_si[ , R.sico$idx_si] = mat_si[ , R.sico$idx_si, drop=FALSE] * -1
    B.sico = mat_B[ , R.sico$idx_co, drop=FALSE] * mat_si
    
    # check #
    expect_equivalent(R.sico$B, B.sico)
    expect_equal(colnames(R.sico$B), names_s[1:dim_S])
    expect_equal(colnames(R.sico$mat_perm), names_s[1:dim_S])
  }
})


