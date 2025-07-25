

#' @title Bootstrap procedure for JB-test
#' @description Bootstraps the distribution of the Jarque-Bera test 
#'   for individual VAR and VECM as described by Kilian, Demiroglu (2000).
#' @param x VAR object of class '\code{varx}' or any other 
#'   that will be \link[=as.varx]{coerced} to '\code{varx}'.
#' @param n.boot Integer. Number of bootstrap iterations.
#' @param n.cores Integer. Number of allocated processor cores.
#' @param fix_beta Logical. If \code{TRUE}, the cointegrating vectors \eqn{\beta} 
#'   are fixed over all bootstrap iterations. Kilian and Demiroglu (2000:43) 
#'   suggest this for VECM with known \eqn{\beta}. 
#'   Ignored in case of rank-unrestricted VAR. 
#' 
#' @return An array of dimension (\eqn{3 \times (1+K) \times} \code{n.boot}) 
#'   containing the \code{n.boot} iteration results for 
#'   \emph{(i)} the JB, skewness and kurtosis test and for
#'   \emph{(ii)} the multivariate and each univariate test 
#'   for the \eqn{K} residual time series.
#' 
#' @examples 
#' \donttest{
#' # select minimal or full example #
#' is_min = TRUE
#' n.boot = ifelse(is_min, 50, 5000)
#' 
#' # prepare data and estimate VAR model #
#' library("vars")
#' data("Canada")
#' exogen  = cbind(qtrend=(1:nrow(Canada))^2)  # quadratic trend
#' R.vars  = VAR(Canada, p=2, type="both", exogen=exogen)
#' R.resid = resid(R.vars)
#' 
#' # bootstrapped JB-test: p-values #
#' S.jb = boot_JB(x=R.vars, n.boot=n.boot, n.cores=1)
#' pvars:::test.normality(u=R.resid, distribution=S.jb[ ,"JB.stats_mlt", ])
#' pvars:::test.normality(u=R.resid, distribution="theoretical")
#' apply(R.resid, MARGIN=2, FUN=function(e_k) 
#'   pvars:::test.normality(u=e_k, distribution="theoretical"))
#' sapply(colnames(R.resid), FUN=function(k) 
#'   pvars:::test.normality(u=R.resid[ , k], distribution=S.jb[ , k, ]))
#' 
#' # density plot #
#' library("ggplot2")
#' R.data = data.frame(t(S.jb[ , "JB.stats_mlt", ]))
#' R.args = list(df=2*R.vars$K)
#' F.density = ggplot() +
#'   stat_density(data=R.data, aes(x=JB, color="bootstrap"), geom="line") +
#'   stat_function(fun=dchisq, args=R.args, n=500, aes(color="theoretical")) +
#'   labs(x="JB statistic", y="Density", color="Distribution", title=NULL) +
#'   theme_bw()
#' plot(F.density)
#' }
#' 
#' @references Kilian, L. and Demiroglu, U. (2000): 
#'   "Residual-Based Tests for Normality in Autoregressions: 
#'   Asymptotic Theory and Simulation Evidence", 
#'   \emph{Journal of Business and Economic Statistics}, 32, pp. 40-50.
#' 
#' @export
#' 
boot_JB <- function(x, n.boot=1000, n.cores=1, fix_beta=FALSE){
  # define
  y = D = D1 = D2 = A = beta = SIGMA = dim_K = dim_T = dim_p = dim_r = NULL
  R.varx  = aux_assign_varx(x)
  u.cov   = SIGMA  # Kilian,Demiroglu 2000:42, Eq.3(2), use OLS covariance matrix.
  names_k = rownames(u.cov)  # names of residuals
  
  # fixed objects in the bootstrap function
  Ystar = matrix(NA, nrow=dim_K, ncol=dim_T+dim_p)
  if(!is.null(dim_r)){
    D = rbind(D2, D1)
    D = aux_rm_Dnl(D, dim_p, x$t_D1, x$t_D2, MARGIN=1)
    A = aux_rm_Dnl(A, dim_p, x$t_D1, x$t_D2, MARGIN=2)
  }
  
  # resampling in the bootstrap function
  L.star = lapply(1:n.boot, FUN=function(n) list(
    U_H0  = t(MASS::mvrnorm(n=dim_T, mu=rep(0, dim_K), Sigma=u.cov, empirical=FALSE)),  # (K x T) matrix
    idx_p = 1:dim_p + sample(0:dim_T, 1)))
  
  # bootstrap function, from Kilian,Demiroglu 2000:43
  bootf = function(star){
    # generate bootstrap time series
    Ustar = star$U_H0  # residuals drawn from H0-distribution of multivariate normality
    Ystar[ ,1:dim_p] = y[ , star$idx_p]  # presample by a randomly drawn block of length dim_p from {y_t}
    for(t in (1:dim_T)+dim_p){   # VAR process
      Ystar[ ,t] = A %*% c(D[ ,t], Ystar[ ,t-(1:dim_p)]) + Ustar[ ,t-dim_p] }
    
    # re-estimate the model
    if(is.null(dim_r)){  # ... by least squares
      Y = Ystar[ ,-(1:dim_p)]   # matrix of regressands, where presample observations are excluded
      Z = aux_stack(Ystar, dim_p=dim_p, D=D)   # matrix of regressors
      B = Y%*%t(Z) %*% solve(Z%*%t(Z))   # OLS estimation, from Luetkepohl 2005:70 / Luetkepohl,Kraetzig 2004:93, Eq.3.9
      e = Y - B%*%Z  # VAR residuals
      
    }else{  # ... by reduced-rank regression
      def  = aux_stackRRR(Ystar, dim_p=dim_p, D1=D1, D2=D2)
      RRR  = aux_RRR(def$Z0, def$Z1, def$Z2)
      beta = if(fix_beta){ beta }else{ RRR$V[ , 0:dim_r, drop=FALSE] }
      e    = aux_VECM(beta=beta, RRR=RRR)$resid  # VECM residuals
    }
    
    # calculate JB statistics
    e = t(e)
    colnames(e)  = names_k
    JB.stats_mlt = test.normality(u=e, distribution="none")
    JB.stats_uni = apply(e, MARGIN=2, FUN=function(e_k) test.normality(u=e_k, distribution="none"))
    return(cbind(JB.stats_mlt, JB.stats_uni))
  }
  
  # run bootstrap-loop and return result
  result = pbapply::pbsapply(L.star, FUN=function(n) bootf(n), cl=n.cores, simplify="array")
  return(result)
}



##############################
###  RESIDUAL DIAGNOSTICS  ###
##############################
#
# Uni- and multivariate tests on 
# normality, ARCH, and serial correlation. 
# Fast and flexible matrix implementation as 
# core modules, e.g. test.normality() is also 
# used for bootstrapping the test distribution.


# Jarque-Bera normality test, from Luetkepohl 2005:175 / Luetkepohl,Kraetzig 2004:129
test.normality <- function(u, distribution="theoretical"){
  # define
  u = aux_asDataMatrixTr(u, "u")  # named matrix u is TxK for scale()
  u = t(scale(u, center=TRUE, scale=FALSE))  # demean matrix  
  dim_K = nrow(u)  # number of endogenous variables to consider 
  dim_T = ncol(u)  # number of vectors in time series
  
  # test statistics
  u.cov = u%*%t(u)/dim_T # MLE covariance matrix 
  ### ...as proposed for this test in Luetkepohl,Kraetzig 2004:129; Luetkepohl 2005:175 uses SSR/dim_T-1 ###
  u.P = t(chol(u.cov))   # lower triangular, Choleski-decomposed covariance matrix
  u.s = solve(u.P)%*%u   # standardized residuals  
  b1 = rowMeans(u.s^3)   # third non-central moment vector (skewness)
  b2 = rowMeans(u.s^4)   # fourth non-central moment vector (kurtosis)
  
  s23 = c(dim_T*t(b1) %*% b1/6)        # test statistic for multivariate skewness
  s24 = c(dim_T*t(b2-3) %*% (b2-3)/24) # test statistic for multivariate kurtosis
  JB  = s23 + s24                      # Jarque-Bera test statistic
  JB.stats = c(JB=JB, Skew=s23, Kurt=s24)
  
  # p-values
  if(is.matrix(distribution)){  # use empirical distribution provided as argument
    result = rowMeans(JB.stats < distribution, na.rm=TRUE)
  }else if(distribution=="none"){  # return test statistics directly
    return(JB.stats)
  }else if(distribution=="theoretical"){  # use theoretical distribution
    result = 1-pchisq(JB.stats, df=c(2, 1, 1)*dim_K)
  }else{ stop("Incorrect specification of the test distribution!") }
  
  # return result
  return(result)
}


# ARCH-LM test, from Luetkepohl 2005:576 / Luetkepohl,Kraetzig 2004:130
test.arch <- function(u, lag.h){
  # define
  u = aux_asDataMatrix(u, "u") # named matrix of residuals from the original VAR is KxT
  dim_K = nrow(u)  # number of endogenous variables to consider
  vech = lower.tri(matrix(NA, dim_K, dim_K), diag = TRUE)      # column-stacking operator
  u.sq = apply(u, MARGIN=2, FUN=function(x) (x%*%t(x))[vech])  # squared residual vectors
  if(dim_K==1){u.sq = t(u.sq)}
  
  # auxiliary VAR models for u.sq_t, from Luetkepohl 2005:70
  Y = u.sq[ , -(1:lag.h), drop=FALSE]
  dim_T = ncol(Y)  # number of time periods without presample
  
  Zh = aux_stack(u.sq, dim_p=lag.h, D=rep(1, dim_T))
  Bh = Y%*%t(Zh) %*% solve(Zh%*%t(Zh)) # OLS estimation of VAR(h) model
  eh = Y - Bh%*%Zh    # residuals of VAR(h)
  eh.cov = cov(t(eh)) # covariance matrix
  
  Z0 = t(rep(1, dim_T))
  B0 = Y%*%t(Z0) %*% solve(Z0%*%t(Z0))  # OLS estimation of VAR(0) model
  e0 = Y - B0%*%Z0    # residuals of VAR(0)
  e0.cov = cov(t(e0)) # covariance matrix
  
  # LM-test
  TR = sum(diag(eh.cov %*% solve(e0.cov))) # trace of covariance product
  LM = 0.5*dim_T*dim_K*(dim_K+1) - dim_T*TR
  p.value = 1-pchisq(LM, df=lag.h*dim_K^2*(dim_K+1)^2/4)
  
  # return result
  result = p.value
  names(result) = paste0("ARCH(", lag.h ,")")
  return(result)
}


# LM test for serial correlation, from Luetkepohl 2005:173,347 / Luetkepohl,Kraetzig 2004:129
test.serial <- function(u, lag.h, x){
  # define
  u = aux_asDataMatrix(u, "u")  # named matrix of residuals from the original VAR is KxT
  x = aux_asDataMatrix(x, "x")  # named matrix of all regressors from the original VAR is KxT
  dim_K = nrow(u)  # number of endogenous variables to consider
  dim_T = ncol(u)  # number of vectors in time series
  
  # auxiliary VAR models for u_t, from Luetkepohl 2005:70
  Y = u
  Z = x
  for(i in 1:lag.h){
    u.i = u[ , 1:(dim_T-i), drop=FALSE]
    u.i = cbind(matrix(0, ncol=i, nrow=dim_K), u.i)
    Z   = rbind(Z, u.i)
  }
  B = Y%*%t(Z) %*% solve(Z%*%t(Z))  # OLS estimation of unrestricted model
  e = Y - B%*%Z             # residuals of unrestricted model
  e.cov = e%*%t(e)/dim_T    # MLE covariance matrix
  
  BR = Y%*%t(x) %*% solve(x%*%t(x))  # OLS estimation of restricted model
  eR = Y - BR%*%x           # residuals of restricted model
  R.cov = eR%*%t(eR)/dim_T  # MLE covariance matrix
  
  # Breusch-Godfrey LM-test
  TR = sum(diag(solve(R.cov) %*% e.cov))  # trace of covariances' product
  LM = dim_T*(dim_K-TR)
  p.value_LM = 1-pchisq(LM, df=lag.h*dim_K^2)
  
  # Edgerton,Shukur (1999) F-test, from Pfaff (2008) / Luetkepohl,Kraetzig (2004)
  m = dim_K*lag.h       # number of additional regressors in the auxiliary system
  r = ((dim_K^2 * m^2 - 4) / (dim_K^2 + m^2 - 5))^0.5
  if(dim_K==1){r = 1}  # according to Doornik(1996), avoids 0/0=NaN
  q = 0.5*dim_K*m-1
  n = nrow(x)       # number of regressors in the original system
  N = dim_T - n - m - 0.5*(dim_K-m+1)
  #R2r = 1 - det(e.cov)/det(R.cov)
  #LMF = (1-(1-R2r)^(1/r))/(1-R2r)^(1/r) * (N*r - q )/(dim_K*m)
  #p.value_LMF = 1-pf(LMF, df1=lag.h*dim_K^2, df2=floor(N*r-q))
  R = det(R.cov)/det(e.cov)
  LMF = (R^(1/r)-1) * (N*r - q)/(dim_K*m)
  p.value_LMF = 1-pf(LMF, df1=lag.h*dim_K^2, df2=N*r-q)
  
  # return result
  result = c(p.value_LM, p.value_LMF)
  names(result) = paste0("AC(", lag.h, ")", c("_LM", "_LMF"))
  return(result)
}


