

#' @title Panel cointegration rank tests
#' @description Performs test procedures for the rank of cointegration in a panel of VAR models. 
#'   First, the chosen individual procedure is applied over 
#'   all \eqn{N} individual entities for \eqn{r_{H0}=0,\ldots,K-1}. 
#'   Then, the \eqn{K \times N} individual statistics and \eqn{p}-values 
#'   are combined to panel test results on each \eqn{r_{H0}} 
#'   using all combination approaches available for the chosen procedure.
#' 
#' @param L.data List of '\code{data.frame}' objects for each individual. 
#'   The variables must have the same succession \eqn{k = 1,\ldots,K} 
#'   in each individual '\code{data.frame}'.
#' @param lags Integer or vector of integers. 
#'   Lag-order of the VAR models in levels, which is
#'   either a common \eqn{p} for all individuals or 
#'   individual-specific \eqn{p_i} for each individual.  
#'   In the vector, \eqn{p_i} must have the same succession 
#'   \eqn{i = 1,\ldots,N} as argument \code{L.data}.
#' @param type Character. The conventional case of the \link[=as.t_D]{deterministic term}.
#' @param t_D1 List of vectors. The activating break periods \eqn{\tau} 
#'   for the period-specific \link[=as.t_D]{deterministic regressors} in \eqn{d_{1,it}}, 
#'   which are restricted to the cointegration relations. 
#'   The accompanying lagged regressors are automatically included in \eqn{d_{2,it}}. 
#'   The \eqn{p}-values are calculated for up to two breaks resp. three sub-samples.
#' @param t_D2 List of vectors. The activating break periods \eqn{\tau} 
#'   for the period-specific \link[=as.t_D]{deterministic regressors} in \eqn{d_{2,it}}, 
#'   which are unrestricted.
#' @param t_D  List of vectors. The activating break periods \eqn{\tau} 
#'   for the period-specific \link[=as.t_D]{deterministic regressors} in \eqn{d_{it}} 
#'   of the SL-procedure. 
#'   The accompanying lagged regressors are automatically included in \eqn{d_{it}}. 
#'   The \eqn{p}-values are calculated for up to two breaks resp. three sub-samples.
#'   See \code{\link{speci.VAR}} for the individual procedure of Luetkepohl et al. (2004)
#'   that determines the break periods \eqn{\tau_i} prior to the cointegration analysis.
#' 
#' @return A list of class '\code{pcoint}' with elements: 
#' \item{panel}{List for the panel test results, 
#'   which contains one matrix for the test statistics and one for the \eqn{p}-values.}
#' \item{individual}{List for the individual test results, 
#'   which contains one matrix for the test statistics and one for the \eqn{p}-values.}
#' \item{CSD}{List of measures for cross-sectional dependency. 
#'   \code{NULL} if a first-generation test has been used.}
#' \item{args_pcoint}{List of characters and integers 
#'   indicating the panel cointegration test and specifications that have been used.}
#' \item{beta_H0}{List of matrices, 
#'   which comprise the pooled cointegrating vectors for each rank \eqn{r_{H0}}. 
#' \code{NULL} if a test other than \code{BR} has been used.}
#' 
#' @family panel cointegration rank tests
#' @name pcoint
NULL


#' @describeIn pcoint based on the Johansen procedure.
#' @param n.factors Integer. Number of common factors to be subtracted 
#'   for the PANIC by Arsova and Oersal (2017, 2018). 
#'   Deactivated if \code{FALSE} (the default).
#' @references Larsson, R., Lyhagen, J., and Lothgren, M. (2001): 
#'   "Likelihood-based Cointegration Tests in Heterogeneous Panels",
#'   \emph{Econometrics Journal}, 4, pp. 109-142.
#' @references Choi, I. (2001): 
#'   "Unit Root Tests for Panel Data", 
#'   \emph{Journal of International Money and Finance}, 20, pp. 249-272.
#' @references Arsova, A., and Oersal, D. D. K. (2018): 
#'   "Likelihood-based Panel Cointegration Test in the Presence of 
#'   a Linear Time Trend and Cross-Sectional Dependence", 
#'   \emph{Econometric Reviews}, 37, pp. 1033-1050.
#' 
#' @export
#' 
pcoint.JO <- function(L.data, lags, type=c("Case1", "Case2", "Case3", "Case4"), 
                      t_D1=NULL, t_D2=NULL, n.factors=FALSE){
  # define and check
  L.y   = aux_check(L.data, "L.data", tr=TRUE)
  dim_N = length(L.data)
  dim_K = ncol(L.y[[1]])
  
  names_H0 = paste0("r_H0 = ", 0:(dim_K-1))
  names_pt = c("LRbar", "Choi_P", "Choi_Pm", "Choi_Z")
  names_i  = names(L.data)
  
  L.dim_T = sapply(L.y, FUN=function(x) nrow(x))
  L.dim_p = aux_check(lags, "lags", dim_N, names_i)
  L.t_D1  = aux_check(t_D1, "t_D1", dim_N, names_i)
  L.t_D2  = aux_check(t_D2, "t_D2", dim_N, names_i)
  
  # defactor the time series for PANIC in Johansen procedure, according to Arsova,Oersal 2018:1036
  if(n.factors>0){
    trend = switch(type, "Case1"=FALSE, "Case2"=FALSE, "Case3"=TRUE, "Case4"=TRUE, "Case5"=TRUE)
    Ycd = do.call("cbind", L.y)  # data matrix of dimension (T x K*N), which includes common and deterministic components
    PCA = aux_ComFact(X=Ycd, trend=trend, n.factors=n.factors, D=c(t_D1, t_D2))  # PANIC (2004) decomposition
    Xit = PCA$eit   # panel of idiosyncratic components without deterministic term
    L.y = lapply(1:dim_N, FUN=function(i) Xit[ ,dim_K*(i-1) + 1:dim_K])  # overwrite for subsequent VECM estimation
    
    type_JO  = "Case1"  # idiosyncratic components need no rotation by beta_oc or deterministic term in the Johansen procedure, from Arsova,Oersal 2018:1041 
    type_mom = "SL_trend"  # for defactored series, use moments of panel SL test, see Arsova,Oersal 2018:1041
  }else{
    PCA = NULL
    type_JO  = type
    type_mom = type
  }
  
  # apply Johansen (1996) procedure to each individual
  moments  = aux_CointMoments(dim_K=dim_K, r_H0=0:(dim_K-1), type=type_mom)
  L.define = lapply(1:dim_N,  FUN=function(i) aux_stackRRR(L.y[[i]], dim_p=L.dim_p[i], type=type_JO, t_D1=L.t_D1[[i]], t_D2=L.t_D2[[i]]))
  L.lambda = lapply(L.define, FUN=function(x) aux_RRR(x$Z0, x$Z1, x$Z2)$lambda)
  L.LRrank = lapply(1:dim_N,  FUN=function(i) aux_LRrank(L.lambda[[i]], dim_T=L.dim_T[i]-L.dim_p[i], dim_K, moments=moments))
  LR.stats = sapply(L.LRrank, FUN=function(x) x$stats_TR)  # matrix of trace statistics (K x N) with ordering r_H0 = 0,...,K-1
  LR.pvals = sapply(L.LRrank, FUN=function(x) x$pvals_TR)  # matrix of p-values (K x N) with ordering r_H0 = 0,...,K-1
  
  # apply panel tests
  idx_D = any(unlist(lapply(L.t_D1, FUN=function(x) names(x) %in% c("t_shift", "t_break"))))
  if(idx_D){ warning("'pcoint.JO' does not adjust the test distribution to 't_shift' or 't_break'.") }
  ### An LRbar version of JMN (2000) using common relative break periods is not available in the literature. ###
  moments = coint_moments[[type_mom]][dim_K:1, , drop=FALSE]
  LRbar = ptest.STATSbar(STATS=LR.stats, distribution="theoretical", moments=moments)  # from Larsson et al. 2001:112-114
  Choi  = ptest.METApval(PVALS=LR.pvals, distribution="theoretical")  # from Choi 2001:253-256
  
  # return result
  indiv = list(stats=t(LR.stats), pvals=t(LR.pvals), lags=L.dim_p, t_D1=L.t_D1, t_D2=L.t_D2)
  panel = list(stats=t(cbind(LRbar$stats, Choi$stats)), pvals=t(cbind(LRbar$pvals, Choi$pvals)))
  
  dimnames(indiv$stats) = dimnames(indiv$pvals) = list(names_i,  names_H0)
  dimnames(panel$stats) = dimnames(panel$pvals) = list(names_pt, names_H0)
  
  argues = list(method="Johansen Procedure", type=type_mom, n.factors=n.factors)
  result = list(panel=panel, individual=indiv, CSD=PCA, args_pcoint=argues)
  class(result) = "pcoint"
  return(result)
}


#' @describeIn pcoint based on the pooled two-step estimation of the cointegrating vectors.
#' @param n.iterations Integer. The (maximum) number of iterations for the 
#'   pooled estimation of the cointegrating vectors.
#' @references Breitung, J. (2005): 
#'   "A Parametric Approach to the Estimation of Cointegration Vectors in Panel Data",
#'   \emph{Econometric Reviews}, 24, pp. 151-173.
#' @export
#' 
pcoint.BR <- function(L.data, lags, type=c("Case1", "Case2", "Case3", "Case4"), 
                      t_D1=NULL, t_D2=NULL, n.iterations=FALSE){
  # define and check
  L.y   = aux_check(L.data, "L.data")
  dim_N = length(L.data)
  dim_K = nrow(L.y[[1]])
  
  names_H0 = paste0("r_H0 = ", 0:(dim_K-1))
  names_pt = c("LRbar", "Choi_P", "Choi_Pm", "Choi_Z")
  names_i  = names(L.data)
  
  L.dim_p = aux_check(lags, "lags", dim_N, names_i)
  L.t_D1  = aux_check(t_D1, "t_D1", dim_N, names_i)
  L.t_D2  = aux_check(t_D2, "t_D2", dim_N, names_i)
  
  # estimate individual VECM with homogeneous cointegrating vectors
  L.def = lapply(1:dim_N, FUN=function(i) aux_stackRRR(L.y[[i]], dim_p=L.dim_p[i], type=type, t_D1=L.t_D1[[i]], t_D2=L.t_D2[[i]]))
  L.RRR = lapply(L.def,   FUN=function(i) aux_RRR(i$Z0, i$Z1, i$Z2, via_R0_R1=TRUE))
  R.est = aux_2StepBR(L.RRR, r_H0=0:(dim_K-1), idx_pool=1:dim_K, n.iterations=n.iterations)
  L.alpha_oc = apply(R.est$L.alpha, MARGIN=1:2, FUN=function(alpha_ir) aux_oc(alpha_ir[[1]]), simplify=FALSE)
  L.beta_oc  = apply(R.est$L.beta,  MARGIN=1:2, FUN=function(beta_ir)  aux_oc(beta_ir[[1]][1:dim_K, , drop=FALSE]), simplify=FALSE)
  
  # individual LM-tests with homogeneous cointegrating vectors, see Breitung 2005:158, Eq.13(II)
  moments  = aux_CointMoments(dim_K=dim_K, r_H0=0:(dim_K-1), type=type)
  L.LMrank = lapply(1:dim_N,  FUN=function(i) aux_LMrank(R0=L.RRR[[i]]$R0, 
                                                         R1=L.RRR[[i]]$R1, 
                                                         L.alpha_oc=L.alpha_oc[i, ], 
                                                         L.beta_oc=L.beta_oc[i, ], 
                                                         moments=moments))
  LM.stats = sapply(L.LMrank, FUN=function(x) x$stats_LM)  # matrix of LM statistics (K x N) with ordering r_H0 = 0,...,K-1
  LM.pvals = sapply(L.LMrank, FUN=function(x) x$pvals_LM)  # matrix of p-values (K x N) with ordering r_H0 = 0,...,K-1
  
  # apply panel tests
  idx_D = any(unlist(lapply(L.t_D1, FUN=function(x) names(x) %in% c("t_shift", "t_break"))))
  if(idx_D){ warning("'pcoint.BR' does not adjust the test distribution to 't_shift' or 't_break'.") }
  ### An LRbar version of JMN (2000) using common relative break periods is not available in the literature. ###
  LRbar = ptest.STATSbar(STATS=LM.stats, distribution="theoretical", moments=moments)  # from Breitung 2005:158, Th.2
  Choi  = ptest.METApval(PVALS=LM.pvals, distribution="theoretical")  # panel tests based on combined p-values are not given in Breitung (2005)!
  
  # return result
  indiv = list(stats=t(LM.stats), pvals=t(LM.pvals), lags=L.dim_p, t_D1=L.t_D1, t_D2=L.t_D2)
  panel = list(stats=t(cbind(LRbar$stats, Choi$stats)), pvals=t(cbind(LRbar$pvals, Choi$pvals)))
  
  dimnames(indiv$stats) = dimnames(indiv$pvals) = list(names_i,  names_H0)
  dimnames(panel$stats) = dimnames(panel$pvals) = list(names_pt, names_H0)
  
  argues = list(method="Pooled Two-Step Estimation", type=type, n.iterations=R.est$n)
  result = list(panel=panel, individual=indiv, CSD=NULL, args_pcoint=argues, beta_H0=R.est$L.beta)
  class(result) = "pcoint"
  return(result)
}


#' @describeIn pcoint based on the Saikkonen-Luetkepohl procedure.
#' @param n.factors Integer. Number of common factors to be subtracted 
#'   for the PANIC by Arsova and Oersal (2017, 2018). 
#'   Deactivated if \code{FALSE} (the default).
#' @references Oersal, D. D. K., and Droge, B. (2014): 
#'   "Panel Cointegration Testing in the Presence of a Time Trend", 
#'   \emph{Computational Statistics & Data Analysis}, 76, pp. 377-390.
#' @references Oersal, D. D. K., and Arsova, A. (2017): 
#'   "Meta-Analytic Cointegrating Rank Tests for Dependent Panels", 
#'   \emph{Econometrics and Statistics}, 2, pp. 61-72.
#' @references Arsova, A., and Oersal, D. D. K. (2018): 
#'   "Likelihood-based Panel Cointegration Test in the Presence of 
#'   a Linear Time Trend and Cross-Sectional Dependence", 
#'   \emph{Econometric Reviews}, 37, pp. 1033-1050.
#' 
#' @examples
#' ### reproduce Oersal,Arsova 2017:67, Ch.5 ###
#' data("MERM")
#' names_k = colnames(MERM)[-(1:2)] # variable names
#' names_i = levels(MERM$id_i)      # country names
#' L.data  = sapply(names_i, FUN=function(i) 
#'    ts(MERM[MERM$id_i==i, names_k], start=c(1995, 1), frequency=12), 
#'    simplify=FALSE)
#' 
#' # Oersal,Arsova 2017:67, Tab.5 #
#' R.lags = c(2, 2, 2, 2, 1, 2, 2, 4, 2, 3, 2, 2, 2, 2, 2, 1, 1, 2, 2)
#' names(R.lags) = names_i  # individual lags by AIC (lag_max=4)
#' n.factors = 8  # number of common factors by Onatski's (2010) criterion
#' R.pcsl = pcoint.SL(L.data, lags=R.lags, n.factors=n.factors, type="SL_trend")
#' R.pcjo = pcoint.JO(L.data, lags=R.lags, n.factors=n.factors, type="Case4")
#' 
#' # Oersal,Arsova 2017:67, Tab.6 #
#' R.Ftsl = coint.SL(y=R.pcsl$CSD$Ft, dim_p=2, type_SL="SL_trend")  # lag-order by AIC
#' R.Ftjo = coint.JO(y=R.pcsl$CSD$Ft, dim_p=2, type="Case4")
#' 
#' @export
#' 
pcoint.SL <- function(L.data, lags, type="SL_trend", t_D=NULL, n.factors=FALSE){
  # define and check
  L.y   = aux_check(L.data, "L.data", tr=TRUE)
  dim_N = length(L.data)
  dim_K = ncol(L.y[[1]])
  
  names_H0 = paste0("r_H0 = ", 0:(dim_K-1))
  names_pt = c("LRbar", "Choi_P", "Choi_Pm", "Choi_Z")
  names_i  = names(L.data)
  
  L.dim_p = aux_check(lags, "lags", dim_N, names_i)
  L.t_D   = aux_check(t_D,  "t_D",  dim_N, names_i)
  
  # defactor the time series for PANIC in SL procedure, according to Arsova,Oersal 2018:1038
  if(n.factors>0){
    trend = switch(type, "SL_mean"=FALSE, "SL_trend"=TRUE)
    Ycd = do.call("cbind", L.y)  # data matrix of dimension (T x K*N), which includes common and deterministic components
    PCA = aux_ComFact(X=Ycd, trend=trend, n.factors=n.factors, D=t_D)  # PANIC (2004) decomposition
    Y.star = PCA$eit2  # panel of idiosyncratic components with deterministic term
    L.y = lapply(1:dim_N, FUN=function(i) Y.star[ ,dim_K*(i-1) + 1:dim_K])  # overwrite for subsequent VECM estimation
    type = "SL_trend"  # defactoring introduces deterministic and stochastic trend in VECM, see Arsova,Oersal 2018:1037 / 2017:64
  }else{
    PCA = NULL
  }
  
  # apply SL procedure to each individual, see Arsova,Oersal 2018:1037,1040 / Saikkonen,Luetkepohl (2000) / Trenkler (2008)
  cointf <- function(i){
    # define
    y = aux_asDataMatrix(L.y[[i]], "y")  # named matrix y is KxT regardless of input
    dim_p = L.dim_p[[i]]  # individual lag-order of VAR in levels
    R.dtr = aux_stackDtr(type_SL=type, t_D=L.t_D[[i]], dim_p=dim_p, dim_T=ncol(y))
    R.def = aux_stackRRR(y=y, dim_p=dim_p, D1=R.dtr$D1, D2=R.dtr$D2)
    dim_K = R.def$dim_K  # number of endogenous variables
    dim_T = R.def$dim_T  # number of observations without presample
    
    # rotate, detrend and test under each r_H0
    RRR    = aux_RRR(Z0=R.def$Z0, Z1=R.def$Z1, Z2=R.def$Z2)
    result = NULL
    for(r_H0 in 0:(dim_K-1)){
      
      if(n.factors>0 & r_H0>0){
        # rotate defactored series for testing at higher cointegrating rank, from Arsova,Oersal 2018:1040
        y_H0    = t(MASS::Null(RRR$V[1:dim_K, 0:r_H0, drop=FALSE])) %*% y  # trending series under H0
        def_H0  = aux_stackRRR(y=y_H0, dim_p=dim_p, D1=R.def$D1, D2=R.def$D2)  # use the same deterministic regressors
        RRR_H0  = aux_RRR(Z0=def_H0$Z0, Z1=def_H0$Z1, Z2=def_H0$Z2)
        rank_H0 = 0  # test for no cointegration in the rotated series
        dim_d   = dim_K - r_H0  # number of stochastic trends under H0
        
      }else{ 
        # no rotation for non-defactored series, from Oersal,Droge 2014:5, Eq.11,12
        y_H0    = y
        RRR_H0  = RRR
        rank_H0 = r_H0
        dim_d   = dim_K
      }
      
      # cointegration test with prior adjustment for deterministic term, from Trenkler 2008:22
      vecm = aux_VECM(beta=RRR_H0$V[ , 0:rank_H0, drop=FALSE], RRR=RRR_H0)
      A_H0 = aux_vec2var(PI=vecm$PI, GAMMA=vecm$GAMMA, dim_p=dim_p)$A
      x_H0 = aux_GLStrend(y=y_H0, OMEGA=vecm$OMEGA, A=A_H0, D=R.dtr$D, dim_p=dim_p)$x
      Z_H0 = aux_stackRRR(y=x_H0, dim_p=dim_p)  # no deterministic term in detrended data
      l_H0 = aux_RRR(Z0=Z_H0$Z0, Z1=Z_H0$Z1, Z2=Z_H0$Z2)$lambda 
      rank = aux_LRrank(lambda=l_H0, dim_T=dim_T, dim_K=dim_d, r_H0=rank_H0, moments=moments[r_H0+1, , drop=FALSE])  # moments order r_H0=0,...,K-1
      
      result = rbind(result, rank)
      rm(y_H0, RRR_H0, rank_H0, dim_d, vecm, A_H0, x_H0, Z_H0, l_H0, rank)
    }
    
    # return result
    return(result)
  }
  moments  = aux_CointMoments(dim_K=dim_K, r_H0=0:(dim_K-1), type=type)
  L.SLrank = lapply(1:dim_N,  FUN=function(i) cointf(i))
  LR.stats = sapply(L.SLrank, FUN=function(x) x$stats_TR)  # matrix of trace statistics (K x N) with ordering r_H0 = 0,...,K-1
  LR.pvals = sapply(L.SLrank, FUN=function(x) x$pvals_TR)  # matrix of p-values (K x N) with ordering r_H0 = 0,...,K-1
  
  # apply panel tests
  idx_D = any(unlist(lapply(L.t_D, FUN=function(x) names(x) == "t_break")))
  if(idx_D){ warning("'pcoint.SL' does not adjust the test distribution to 't_break'.") }
  ### An LRbar version of TSL (2008) using common relative break periods is not available in the literature. ###
  if(type=="SL_trend"){ moments = coint_moments[[type]][dim_K:1, , drop=FALSE] }  # use exactly those moments which Arsova,Oersal (2018) have simulated
  LRbar = ptest.STATSbar(STATS=LR.stats, distribution="theoretical", moments=moments)  # from Arsova,Oersal 2018:1039, Eq.12
  Choi  = ptest.METApval(PVALS=LR.pvals, distribution="theoretical")  # from Oersal,Arsova 2017:63, Eq.7-10
  
  # return result
  indiv = list(stats=t(LR.stats), pvals=t(LR.pvals), lags=L.dim_p, t_D=L.t_D)
  panel = list(stats=t(cbind(LRbar$stats, Choi$stats)), pvals=t(cbind(LRbar$pvals, Choi$pvals)))
  
  dimnames(indiv$stats) = dimnames(indiv$pvals) = list(names_i,  names_H0)
  dimnames(panel$stats) = dimnames(panel$pvals) = list(names_pt, names_H0)
  
  argues = list(method="SL-Procedure", type=type, n.factors=n.factors)
  result = list(panel=panel, individual=indiv, CSD=PCA, args_pcoint=argues)
  class(result) = "pcoint"
  return(result)
}


#' @describeIn pcoint accounting for correlated probits between the individual SL-procedures.
#' @references Hartung, J. (1999):
#'   "A Note on Combining Dependent Tests of Significance",
#'   \emph{Biometrical Journal}, 41, pp. 849-855.
#' @references Luetkepohl, H., Saikkonen, P., and Trenkler, C. (2004): 
#'   "Testing for the Cointegrating Rank of a VAR Process with Level Shift at Unknown Time", 
#'   \emph{Econometrica}, 72, pp. 647-662.
#' @references Arsova, A., and Oersal, D. D. K. (2021): 
#'   "A Panel Cointegrating Rank Test with Structural Breaks and Cross-Sectional Dependence", 
#'   \emph{Econometrics and Statistics}, 17, pp. 107-129.
#' 
#' @examples
#' ### reproduce Oersal,Arsova 2016:13, Ch.6 ###
#' data("ERPT")
#' names_k = c("lpm5", "lfp5", "llcusd")  # variable names for "Chemicals and related products"
#' names_i = levels(ERPT$id_i)[c(1,6,2,5,4,3,7)]  # ordered country names
#' L.data  = sapply(names_i, FUN=function(i) 
#'    ts(ERPT[ERPT$id_i==i, names_k], start=c(1995, 1), frequency=12), 
#'    simplify=FALSE)
#' 
#' # Oersal,Arsova 2016:21, Tab.6 (only for individual results) #
#' R.lags = c(3, 3, 3, 4, 3, 3, 3); names(R.lags)=names_i  # lags of VAR model by MAIC
#' R.cain = pcoint.CAIN(L.data, lags=R.lags, type="SL_trend")
#' R.pcsl = pcoint.SL(L.data,   lags=R.lags, type="SL_trend")
#' 
#' # Oersal,Arsova 2016:22, Tab.7/8 #
#' R.lags = c(3, 3, 3, 4, 4, 3, 4); names(R.lags)=names_i  # lags of VAR model by MAIC
#' R.t_D  = list(t_break=89)  # a level shift and trend break in 2002_May for all countries
#' R.cain = pcoint.CAIN(L.data, lags=R.lags, t_D=R.t_D, type="SL_trend")
#' 
#' @export
#' 
pcoint.CAIN <- function(L.data, lags, type="SL_trend", t_D=NULL){
  # define and check
  L.y   = aux_check(L.data, "L.data")
  dim_N = length(L.data)
  dim_K = nrow(L.y[[1]])
  
  names_i  = names(L.data)
  names_H0 = paste0("r_H0 = ", 0:(dim_K-1))
  names_pt = c("Hartung_K1", "Hartung_K2", "CAIN")
  
  L.dim_T = sapply(L.y, FUN=function(x) ncol(x))
  L.dim_p = aux_check(lags, "lags", dim_N, names_i)
  L.t_D   = aux_check(t_D,  "t_D",  dim_N, names_i)
  
  dim_pN = max(L.dim_p)  # maximum lag-order of all individuals
  dim_TN = max(L.dim_T)  # total number of periods incl. presamples
  idx_tN = -(0:dim_pN)   # periods without presample of any individual
  
  # apply (T)SL procedure to each individual, from Oersal,Arsova 2016:3 / Trenkler et al. (2008)
  cointf <- function(i){
    # define
    y = L.y[[i]]  # named matrix y is KxT regardless of input
    dim_p = L.dim_p[[i]]  # individual lag-order of VAR in levels
    R.dtr = aux_stackDtr(type_SL=type, t_D=L.t_D[[i]], dim_p=dim_p, dim_T=ncol(y))
    R.def = aux_stackRRR(y=y, dim_p=dim_p, D1=R.dtr$D1, D2=R.dtr$D2)
    dim_K = R.def$dim_K  # number of endogenous variables
    dim_T = R.def$dim_T  # number of observations without presample
    
    # residuals for the overall sample
    resids = array(NA, dim=c(dim_K, dim_TN, dim_K))
    idx_t  = seq.int(to=dim_TN, length.out=dim_T)
    
    # detrend and test under each r_H0
    type_mom = if(is.null(L.t_D[[i]]$t_break)){ type }else{ "TSL_trend" }
    moments  = aux_CointMoments(dim_K=dim_K, dim_T=dim_T+dim_p, t_D1=L.t_D[[i]], r_H0=0:(dim_K-1), type=type_mom)
    RRR      = aux_RRR(Z0=R.def$Z0, Z1=R.def$Z1, Z2=R.def$Z2)
    SLrank   = NULL
    for(r_H0 in 0:(dim_K-1)){
      # cointegration test with prior adjustment for deterministic term, from Trenkler 2008:22
      vecm = aux_VECM(beta=RRR$V[ , 0:r_H0, drop=FALSE], RRR=RRR)
      A_H0 = aux_vec2var(PI=vecm$PI, GAMMA=vecm$GAMMA, dim_p=dim_p)$A
      x_H0 = aux_GLStrend(y=y, OMEGA=vecm$OMEGA, A=A_H0, D=R.dtr$D, dim_p=dim_p)$x
      Z_H0 = aux_stackRRR(y=x_H0, dim_p=dim_p)  # no deterministic term in detrended data
      l_H0 = aux_RRR(Z0=Z_H0$Z0, Z1=Z_H0$Z1, Z2=Z_H0$Z2)$lambda
      rank = aux_LRrank(lambda=l_H0, dim_T=dim_T, dim_K=dim_K, r_H0=r_H0, moments=moments[r_H0+1, , drop=FALSE])  # moments order r_H0=0,...,K-1
      
      SLrank = rbind(SLrank, rank)
      resids[ , idx_t, r_H0+1] = vecm$resid
      rm(vecm, A_H0, x_H0, Z_H0, l_H0, rank)
    }
    
    # return result
    result = list(SLrank=SLrank, resids=resids)
    return(result)
  }
  L.SLrank = lapply(1:dim_N,  FUN=function(i) cointf(i))
  LR.stats = sapply(L.SLrank, FUN=function(x) x$SLrank$stats_TR)  # matrix of trace statistics (K x N) with ordering r_H0 = 0,...,K-1
  LR.pvals = sapply(L.SLrank, FUN=function(x) x$SLrank$pvals_TR)  # matrix of p-values (K x N) with ordering r_H0 = 0,...,K-1
  resids   = sapply(L.SLrank, FUN=function(x) x$resids[ , idx_tN, 1], simplify="array")  # residuals (K x T-p x N) under r_H0 = 0, from Oersal,Arsova 2016:7 
  
  # apply panel tests
  idx_lower = lower.tri(diag(dim_N))  # gets the lower triangular part of the NxN correlation matrices
  rho_resid = apply(resids, MARGIN=1, FUN=function(k) cor(k, method="pearson")[idx_lower])
  rho_eps   = mean(abs(rho_resid))  # average absolute cross-sectional correlation, from Oersal,Arsova 2016:7, Eq.12
  CAIN      = ptest.CAIN(LR.pvals, distribution="theoretical", dim_K=dim_K, r_H0=0:(dim_K-1), rho_eps=rho_eps)
  
  # return result
  CSD   = list(rho_tilde=CAIN$rho_tilde, rho_eps=rho_eps, rho_resid=rho_resid)
  indiv = list(stats=t(LR.stats), pvals=t(LR.pvals), lags=L.dim_p, t_D=L.t_D)
  panel = list(stats=t(CAIN$stats), pvals=t(CAIN$pvals))
  
  dimnames(indiv$stats) = dimnames(indiv$pvals) = list(names_i,  names_H0)
  dimnames(panel$stats) = dimnames(panel$pvals) = list(names_pt, names_H0)
  dimnames(CSD$rho_tilde) = list(names_pt, names_H0)
  
  argues = list(method="CAIN", type=type)
  result = list(panel=panel, individual=indiv, CSD=CSD, args_pcoint=argues)
  class(result) = "pcoint"
  return(result)
}


#### S3 methods for objects of class 'pcoint' ####
#' @export
print.pcoint <- function(x, ...){
  # define
  dim_N = nrow(x$individual$stats)
  dim_K = ncol(x$individual$stats)
  
  names_pt = rownames(x$panel$stats)
  names_H0 = colnames(x$individual$stats)
  names_i  = rownames(x$individual$stats)
  
  n.char_H0 = sum(nchar(names_H0)) + dim_K-1
  n.char_Pt = max(nchar(names_pt)) + 1
  n.char_i  = if(is.null(names_i)){ nchar(dim_N)+3 }else{ max(nchar(names_i)) }
  
  # create table
  if(n.char_i < 6){ 
    new_line = "\n"
    n.char_x = 0
  }else{
    new_line = NULL
    n.char_x = nchar("Individual:")
  }
  
  if(is.matrix(x$individual$lags)){
    header_lags = c("lags",  rep(".", times=n.char_H0-4), " ")
    n.char_lags = 1
  }else{
    header_lags = NULL
    n.char_lags = nchar("lags")+2
  }
  
  header_args  = c("### Cointegration Rank by ", x$args_pcoint$method, " with ", x$args_pcoint$type, " ###")
  header_base  = c("statistics", rep(".", times=n.char_H0-10), " ", "p-values", rep(".", times=n.char_H0-8))
  header_indiv = c("Individual:", new_line, rep(" ", times=n.char_i-n.char_x+n.char_lags), header_lags, header_base)
  matrix_indiv = cbind(lags=x$individual$lags, 
                       round(x$individual$stats, 3), 
                       round(x$individual$pvals, 3))
  rownames(matrix_indiv) = names_i
  header_panel = c("Panel:", rep(" ", times=n.char_Pt-6), header_base)
  matrix_panel = cbind(round(x$panel$stats, 3), 
                       round(x$panel$pvals, 3))
  
  # print
  cat(header_args,  "\n", sep="")
  cat(header_indiv, "\n", sep="")
  print(matrix_indiv, quote=FALSE)
  cat("\n", header_panel, "\n", sep="")
  print(matrix_panel, quote=FALSE)
}


#' @export
summary.pcoint <- function(object, ...){
  # define
  dim_N = nrow(object$individual$stats)
  dim_K = ncol(object$individual$stats)
  
  # print
  print.pcoint(object, ...)
  cat("\n", paste0(c("Specifications", rep(".", times=20))), sep="", "\n")
  cat("Method: ", object$args_pcoint$method, "\n")
  cat("Deterministic term: ", object$args_pcoint$type, "\n")
  cat("Number of individuals: ", dim_N, "\n")
  
  if(!is.null(object$args_pcoint$n.factors)){
    cat("Number of common factors: ", as.integer(object$args_pcoint$n.factors), "\n") }
  if(object$args_pcoint$method == "Pooled Two-Step Estimation"){
    cat("Number of iterations: ", as.integer(object$args_pcoint$n.iterations), "\n")
    cat("Pooled cointegrating vectors: \n")
    print(object$beta_H0)}
  if(object$args_pcoint$method == "CAIN"){
    cat("Correction factors: ", "\n")
    print(object$CSD$rho_tilde)}
}


#' @export
toLatex.pcoint <- function(object, ..., digits=3){
  # define
  dim_N = nrow(object$individual$stats)
  dim_K = ncol(object$individual$stats)
  
  names_pt = rownames(object$panel$stats)
  names_H0 = colnames(object$individual$stats)
  names_i  = rownames(object$individual$stats)
  
  # sanitize subscripts into math mode
  idx_math = grepl(x=names_pt, pattern="_")
  names_pt = gsub(x=names_pt, pattern="_", replacement=" $")
  names_pt = paste0(names_pt, ifelse(idx_math, "$", ""))
  
  # create tabular
  header_H0    = paste0("$ r_{H0} = ", 0:(dim_K-1), " $")
  header_lags  = "\\multicolumn{2}{r}{ lags }"
  # header_lags2 = if(is.matrix(object$individual$lags)){ names_H0 } 
  ### TODO header for pcoint.MSB resp header_break in header_specs for pcoint.CAIN
  idx_col      = 2 + c(1, dim_K, 1+dim_K, 2*dim_K) + c(0, 0, 1, 1) # second vector adds separating column
  header_base  = paste0("\\multicolumn{", dim_K, "}{c}{ statistics }", " &  & ", 
                        "\\multicolumn{", dim_K, "}{c}{ $ p $-values }")
  header_cline = paste0("\\cline{", idx_col[1], "-", idx_col[2] ,"} ", 
                        "\\cline{", idx_col[3], "-", idx_col[4] ,"}")
  
  header_indiv = c("\\multicolumn{2}{l}{ \\textbf{Individual} }", header_base)
  matrix_indiv = cbind(if(is.null(names_i)){ rep(" ", dim_N) }else{ names_i }, 
                       lags=object$individual$lags, 
                       format(round(object$individual$stats, digits=digits), nsmall=digits), " ",
                       format(round(object$individual$pvals, digits=digits), nsmall=digits))
  
  header_panel = c("\\multicolumn{2}{l}{ \\textbf{Panel} }", header_base)
  matrix_panel = cbind(sapply(names_pt, FUN=function(x) paste0("\\multicolumn{2}{l}{", x, "}")), 
                       format(round(object$panel$stats, digits=digits), nsmall=digits), " ",
                       format(round(object$panel$pvals, digits=digits), nsmall=digits))
  
  # return result
  tabular_cols = paste0(rep("r", each=2*dim_K+1), collapse = "")
  result = c(
    paste0("\\begin{tabular}{lr", tabular_cols, "}"),
    "\t\\hline \\hline",
    
    paste0("\t", paste(header_indiv, collapse=" & "), " \\\\"),
    paste0("\t", header_cline),
    paste0("\t", paste(c(header_lags, header_H0, " ", header_H0), collapse=" & "), " \\\\"),
    "\t\\hline",
    apply(matrix_indiv, MARGIN=1, FUN=function(x) paste0("\t", paste(x, collapse=" & "), " \\\\")),
    
    "\t\\hline",
    paste0("\t", paste(header_panel, collapse=" & "), " \\\\"),
    paste0("\t", header_cline),
    paste0("\t", paste(c(" ", " ", header_H0, " ", header_H0), collapse=" & "), " \\\\"),
    "\t\\hline",
    apply(matrix_panel, MARGIN=1, FUN=function(x) paste0("\t", paste(x, collapse=" & "), " \\\\")),
    
    "\t\\hline \\hline",
    "\\end{tabular}"
  )
  class(result) = "Latex"
  return(result)
}


