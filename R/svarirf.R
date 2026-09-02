

#' @title Impulse Response Functions for panel SVAR models
#' @description Calculates impulse response functions for panel VAR objects.
#' @param x Panel VAR object of class '\code{pid}' or '\code{pvarx}' 
#'   or a list of VAR objects that will be \link[=as.varx]{coerced} to '\code{varx}'.
#' @param ... Currently not used.
#' @param n.ahead Integer. Number of periods to consider after the initial impulse, i.e. the horizon of the IRF.
#' @param normf Function. A given function that normalizes the \eqn{K \times S} input-matrix 
#'   into an output matrix of same dimension. See the example in \code{\link{id.iv}} 
#'   for the normalization of Jentsch and Lunsford (2022) 
#'   that fixes the size of the impact response.
#' @param w Numeric, logical, or character vector. 
#'   \eqn{N} numeric elements weighting the individual coefficients, or 
#'   names or \eqn{N} logical elements selecting a subset from the 
#'   individuals \eqn{i = 1, \ldots, N} for the MG estimation. If \code{NULL} 
#'   (the default), all \eqn{N} individuals are included without weights.
#' @param MG_IRF Logical. If \code{TRUE}, the mean-group of individual 
#'   IRF is calculated in accordance with Gambacorta et al. (2014). 
#'   If \code{FALSE} (the default recommended in Empting et al. 2026), 
#'   the IRF is calculated for the mean-group of individual VAR estimates.
#' 
#' @return A list of class '\code{svarirf}' holding the impulse response functions as a '\code{data.frame}'.
#' 
#' @references Luetkepohl, H. (2005): 
#'   \emph{New Introduction to Multiple Time Series Analysis}, 
#'   Springer, 2nd ed.
#' @references Gambacorta L., Hofmann B., and Peersman G. (2014):
#'   "The Effectiveness of Unconventional Monetary Policy at the Zero Lower Bound: A Cross-Country Analysis",
#'   \emph{Journal of Money, Credit and Banking}, 46, pp. 615-642.
#' @references Jentsch, C., and Lunsford, K. G. (2022):
#'   "Asymptotically Valid Bootstrap Inference for Proxy SVARs",
#'   \emph{Journal of Business and Economic Statistics}, 40, pp. 1876-1891.
#' @references Empting, L. F. T., Maxand, S., Oeztuerk, S., and Wagner, K. (2026): 
#'   "Inference in Panel SVARs with Two-Way Dependence".
#' 
#' @examples
#' data("PCAP")
#' names_k = c("g", "k", "l", "y")  # variable names
#' names_i = levels(PCAP$id_i)      # country names
#' L.data  = sapply(names_i, FUN=function(i) 
#'   ts(PCAP[PCAP$id_i==i, names_k], start=1960, end=2019, frequency=1), 
#'   simplify=FALSE)
#' 
#' # estimate and identify panel SVAR #
#' L.vars = lapply(L.data, FUN=function(x) vars::VAR(x, p=2, type="both"))
#' R.pid  = pid.chol(L.vars, order_k=names_k)
#' 
#' # calculate and plot MG-IRF #
#' library("ggplot2")
#' R.irf = irf(R.pid, n.ahead=60)
#' F.irf = plot(R.irf, selection=list(2:4, 1:2))
#' as.pplot(F.irf=F.irf, color_g="black", n.rows=3)$F.plot + guides(color="none")
#'
#' @import vars
#' @method irf pvarx
#' @export
#' 
irf.pvarx <- function(x, ..., n.ahead=20, normf=NULL, w=NULL, MG_IRF=FALSE){
  # define and check
  x = as.pvarx(x, w=w)
  
  # names for variables, for shocks, and for the header of each IRF
  names_k   = if( !is.null(rownames(x$A)) ){ rownames(x$A) }else{ paste0("y[ ", 1:nrow(x$A), " ]") }
  names_s   = if( !is.null(colnames(x$B)) ){ colnames(x$B) }else{ paste0("epsilon[ ", 1:ncol(x$B), " ]") }
  names_IRF = c(sapply(names_k, FUN=function(k) paste0(names_s, " %->% ", k)))
  
  # calculate structural IRF
  if(MG_IRF){  # ... by MG of individual IRF, from Gambacorta et al. 2014:627
    A.irf = sapply(x$L.varx, simplify="array", FUN=function(x) 
      aux_var2vma(A=x$A, B=x$B, dim_p=x$dim_p, n.ahead=n.ahead, normf=normf)$THETA)
    THETA = aux_MG(A.irf, w=w)$mean  # (optionally weighted) group mean
  
  }else{  # ... by IRF for MG of individual VAR coefficients
    dim_p = max(sapply(x$L.varx, FUN=function(i) i$dim_p))
    THETA = aux_var2vma(A=x$A, B=x$B, dim_p=dim_p, n.ahead=n.ahead, normf=normf)$THETA
  }
  
  # return result
  IRF = aperm(THETA, perm=c(2,1,3))
  IRF = matrix(IRF, nrow=n.ahead+1, byrow=TRUE, dimnames=list(NULL, names_IRF))
  IRF = as.data.frame(cbind(V1=0:n.ahead, IRF))
  result = list(irf=IRF)
  class(result) = "svarirf"
  return(result)
}


#' @title Impulse Response Functions
#' @description Calculates impulse response functions.
#' @param x VAR object of class '\code{varx}', '\code{id}', or any other 
#'   that will be \link[=as.varx]{coerced} to '\code{varx}'.
#' @param ... Currently not used.
#' @param n.ahead Integer. Number of periods to consider after the initial impulse, i.e. the horizon of the IRF.
#' @param normf Function. A given function that normalizes the \eqn{K \times S} input-matrix 
#'   into an output matrix of same dimension. See the example in \code{\link{id.iv}} 
#'   for the normalization of Jentsch and Lunsford (2022) 
#'   that fixes the size of the impact response.
#' 
#' @return A list of class '\code{svarirf}' holding the impulse response functions as a '\code{data.frame}'.
#' 
#' @references Luetkepohl, H. (2005): 
#'   \emph{New Introduction to Multiple Time Series Analysis}, 
#'   Springer, 2nd ed.
#' @references Jentsch, C., and Lunsford, K. G. (2022):
#'   "Asymptotically Valid Bootstrap Inference for Proxy SVARs",
#'   \emph{Journal of Business and Economic Statistics}, 40, pp. 1876-1891.
#' @examples
#' data("PCIT")
#' names_k = c("APITR", "ACITR", "PITB", "CITB", "GOV", "RGDP", "DEBT")
#' names_l = c("m_PI", "m_CI")  # proxy names
#' names_s = paste0("epsilon[ ", c("PI", "CI"), " ]")  # shock names
#' dim_p   = 4  # lag-order
#' 
#' # estimate and identify proxy SVAR #
#' R.vars = vars::VAR(PCIT[ , names_k], p=dim_p, type="const")
#' R.idBL = id.iv(R.vars, iv=PCIT[-(1:dim_p), names_l], S2="MR", cov_u="OMEGA")
#' colnames(R.idBL$B) = names_s  # labeling
#' 
#' # calculate and plot normalized IRF #
#' R.norm = function(B) B / matrix(-diag(B), nrow(B), ncol(B), byrow=TRUE)
#' plot(irf(R.idBL, normf=R.norm))
#' 
#' @import vars
#' @method irf varx
#' @export
#' 
irf.varx <- function(x, ..., n.ahead=20, normf=NULL){
  # define
  x = as.varx(x)
  
  # names for variables, for shocks, and for the header of each IRF
  names_k   = if( !is.null(rownames(x$A)) ){ rownames(x$A) }else{ paste0("y[ ", 1:nrow(x$A), " ]") }
  names_s   = if( !is.null(colnames(x$B)) ){ colnames(x$B) }else{ paste0("epsilon[ ", 1:ncol(x$B), " ]") }
  names_IRF = c(sapply(names_k, FUN=function(k) paste0(names_s, " %->% ", k)))
  
  # calculate structural IRF
  THETA = aux_var2vma(A=x$A, B=x$B, dim_p=x$dim_p, n.ahead=n.ahead, normf=normf)$THETA
  
  # return result
  IRF = aperm(THETA, perm=c(2,1,3))
  IRF = matrix(IRF, nrow=n.ahead+1, byrow=TRUE, dimnames=list(NULL, names_IRF))
  IRF = as.data.frame(cbind(V1=0:n.ahead, IRF))
  result = list(irf=IRF)
  class(result) = "svarirf"
  return(result)
}


#' @title Persistence Profiles
#' @description Calculates persistence profiles for each of the \eqn{r} long-run relationships.
#' @param x Rank-restricted VAR object of class '\code{varx}' or any other 
#'   that can be \link[=as.varx]{coerced} to '\code{varx}', e.g. '\code{\link[vars]{vec2var}}'. 
#'   If the object is also of child class '\code{id}', \code{\link{PP.variable}} calculates 
#'   the persistence profiles which are initiated by the provided structural shocks.
#' @param n.ahead Integer. Number of periods to consider after the initial impulse, i.e. the horizon of the PP.
#' 
#' @return A list of class '\code{svarirf}' holding the persistence profiles as a '\code{data.frame}'.
#' 
#' @references Lee, K., C., Pesaran, M. H. (1993): 
#'   "Persistence Profiles and Business Cycle Fluctuations in a Disaggregated Model of UK Output Growth", 
#'   \emph{Richerche Economiche}, 47, pp. 293-322.
#' @references Pesaran, M. H., and Shin, Y. (1996): 
#'   "Cointegration and Speed of Convergence to Equilibrium", 
#'   \emph{Journal of Econometrics}, 71, pp. 117-143.
#' 
#' @examples
#' data("PCAP")
#' names_k = c("g", "k", "l", "y") # variable names
#' names_i = levels(PCAP$id_i)     # country names
#' L.data  = sapply(names_i, FUN=function(i) 
#'   ts(PCAP[PCAP$id_i==i, names_k], start=1960, end=2019, frequency=1), 
#'   simplify=FALSE)
#'   
#' # estimate VAR for DNK under rank-restriction r=2 #
#' dim_r  = 2  # cointegrataion rank
#' R.t_D1 = list(t_break=c(23, 49))  # trend breaks
#' R.vecm = VECM(y=L.data$DNK, dim_p=2, dim_r=dim_r, type="Case4", t_D1=R.t_D1)
#' 
#' # define shocks #
#' shock1 = diag(4)   # 4 separate shocks
#' shock2 = cbind(c(1, 0,  0, 0),  # positive shock on "g"
#'                c(0, 0, -1, 0),  # negative shock on "l"
#'                c(0, 0,  1, 1))  # simultaneous shocks
#' 
#' # calculate persistence profiles #
#' R.ppv1 = PP.variable(R.vecm, n.ahead=50, shock=shock1)
#' R.ppv2 = PP.variable(R.vecm, n.ahead=50, shock=shock2)
#' R.ppsy = PP.system(R.vecm, n.ahead=50)
#' 
#' # edit plots #
#' library("ggplot2")
#' as.pplot(ppv1=plot(R.ppv1), n.rows=4)$F.plot + guides(color="none")
#' as.pplot(ppv2=plot(R.ppv2), n.rows=3, color_g="black") # reshape facet array
#' plot(R.ppsy, selection=list(1, c(1,4)))  # dismiss cross-term PP
#' 
#' @name PP
NULL


#' @describeIn PP PP due to a system-wide shock
#' @export
#' 
PP.system <- function(x, n.ahead=20){
  # define
  x = as.varx(x)
  A = x$A  # rank-restricted VAR model in levels
  dim_r = ncol(x$beta)  # cointegration rank
  dim_p = x$dim_p  # lag-order of the VAR model in levels
  dim_K = nrow(A)  # number of endogenous variables
  idx_k = 1:dim_K  # equivalent to selection matrix 'J' from Pesaran,Shin 1996:130
  U.cov = x$SIGMA  # OLS covariance matrix of residuals
  beta  = x$beta[idx_k, , drop=FALSE]  # cointegrating matrix without deterministic terms
  
  # names of the long-run relations, with cross-term responses, and for the PP
  names_r  = if( !is.null(colnames(beta)) ){ colnames(beta) }else{ paste0("ect.", 1:dim_r) }
  names_rc = sapply(names_r, FUN=function(r) paste0(r, " %*% ", names_r))
  names_PP = paste0("epsilon %->% ", names_rc)
  idx_diag = 1:dim_r + 0:(dim_r-1)*dim_r  # index for responses on the main diagonal
  names_PP[idx_diag] = paste0("epsilon %->% ", names_r)
  
  # calculate Persistence Profiles, from Pesaran,Shin 1996:129
  PHI = aux_var2companion(A=A, dim_p=dim_p)  # coefficient matrix for the VAR(1)-companion representation
  Hz = hz = array(NA, dim=c(dim_r, dim_r, n.ahead+1))
  PHI.n = diag(dim_K*dim_p)  # initial PHI^0 = I_{K*p}
  for(n in 0:n.ahead+1){
    B.n = PHI.n[idx_k, idx_k, drop=FALSE]  # from Pesaran,Shin 1996:125, Eq.19'
    Hz[ , , n] = t(beta) %*% B.n %*% U.cov %*% t(B.n) %*% beta  # unscaled PP, from Pesaran,Shin 1996:125, Eq.11
    g = diag(as.matrix(Hz[ , , 1]))^-0.5   # scaling matrix' vector of diagonal elements
    G = diag(g, ncol=dim_r, nrow=dim_r)    # scaling matrix from PP's impact coefficients
    hz[ , , n] = G %*% Hz[ , , n] %*% G    # scaled PP from Pesaran,Shin 1996:125, Eq.13
    PHI.n = PHI.n %*% PHI  # PHI^(n+1) for the subsequent run in this loop
  }
  
  # return result
  PP = matrix(hz, nrow=n.ahead+1, byrow=TRUE, dimnames=list(NULL, names_PP))
  PP = as.data.frame(cbind(V1=0:n.ahead, PP))
  result = list(irf=PP, unscaled=Hz, scaled=hz)
  class(result) = "svarirf"
  return(result)
}


#' @describeIn PP PP due to a structural or variable-specific shock
#' @param shock Matrix. Each column vector specifies a set of simultaneous shocks, 
#'   which initiate \eqn{r} persistence profiles. If \code{NULL} (the default),
#'   a separate unit impulse is set for each shock. 
#' @export
#' 
PP.variable <- function(x, n.ahead=20, shock=NULL){
  # define
  x = as.varx(x)
  A = x$A  # rank-restricted VAR model in levels
  dim_r = ncol(x$beta)  # cointegration rank
  dim_p = x$dim_p  # lag-order of the VAR model in levels
  dim_K = nrow(A)  # number of variables
  idx_k = 1:dim_K  # equivalent to selection matrix 'J' from Pesaran,Shin 1996:130
  beta  = x$beta[idx_k, , drop=FALSE]  # cointegrating matrix without deterministic terms
  
  # find mixing matrix
  if( identical(unname(x$B), diag(dim_K)) ){
    U.cov = x$SIGMA  # OLS covariance matrix of residuals
    U.P   = t(chol(U.cov))  # lower triangular, Choleski-decomposed covariance matrix
    names_k = if( !is.null(rownames(A)) ){ rownames(A) }else{ paste0("y", idx_k) }  # variable names
    names_s = paste0("epsilon[ ", names_k, " ]")
  }else{
    U.P = x$B  # structural impact matrix
    names_s = if( !is.null(colnames(x$B)) ){ colnames(x$B) }else{ paste0("epsilon[ ", 1:ncol(x$B), " ]") }
  }
  
  # define shocks
  dim_S = ncol(U.P)
  shock = if( !is.null(shock) ){ as.matrix(shock) }else{ diag(dim_S) }
  
  # names of the long-run relations, of (combined) shocks, and for the PP 
  names_r  = if( !is.null(colnames(beta)) ){ colnames(beta) }else{ paste0("ect.", 1:dim_r) }
  names_sc = apply(shock, MARGIN=2, FUN=function(s) paste0(names_s[s!=0], collapse="+"))
  names_PP = c(sapply(names_sc, FUN=function(s) paste0(s, " %->% ", names_r)))
  
  # calculate Persistence Profiles, from Pesaran,Shin 1996:129
  PHI = aux_var2companion(A=A, dim_p=dim_p)  # coefficient matrix for the VAR(1)-companion representation
  PSI = array(NA, dim=c(dim_r, ncol(shock), n.ahead+1))
  PHI.n = diag(dim_K*dim_p)  # initial PHI^0 = I_{K*p}
  for(n in 0:n.ahead+1){
    B.n = PHI.n[idx_k, idx_k, drop=FALSE]  # from Pesaran,Shin 1996:125, Eq.19'
    PSI[ , , n] = t(beta) %*% B.n %*% U.P %*% shock  # PP, from Pesaran,Shin 1996:122, Eq.8
    PHI.n = PHI.n %*% PHI  # PHI^(n+1) for the subsequent run in this loop
  }
  
  # return result
  PP = matrix(PSI, nrow=n.ahead+1, byrow=TRUE, dimnames=list(NULL, names_PP))
  PP = as.data.frame(cbind(V1=0:n.ahead, PP))
  result = list(irf=PP, PSI=PSI, U.P=U.P)
  class(result) = "svarirf"
  return(result)
}


