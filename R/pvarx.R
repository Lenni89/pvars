

#' @title Estimation of VAR models for heterogeneous panels
#' @description Performs the (pooled) mean-group estimation of a panel VAR model. 
#'   First, VAR models are estimated for all \eqn{N} individual entities. 
#'   Then, their (pooled) mean-group estimate is calculated for each coefficient.
#' 
#' @param L.data List of '\code{data.frame}' objects for each individual. 
#'   The variables must have the same succession \eqn{k = 1,\ldots,K} 
#'   in each individual '\code{data.frame}'.
#' @param lags Integer or vector of integers. 
#'   Lag-order of the VAR models in levels, which is
#'   either a common \eqn{p} for all individuals or 
#'   individual-specific \eqn{p_i} for each individual.  
#'   In the vector, \eqn{p_i} must have the same succession 
#'   \eqn{i = 1,\ldots,N} as argument '\code{L.data}'.
#' @param type Character. The conventional case of the \link[=as.t_D]{deterministic term}.
#' @param n.factors Integer. Number of common factors to be used for SUR. 
#'   Deactivated if \code{FALSE} (the default).
#' @param n.iterations Integer. The (maximum) number of iterations 
#'   for the estimation of SUR resp. the pooled cointegrating vectors.
#' 
#' @return A list of class '\code{\link[=as.pvarx]{pvarx}}' with the elements: 
#' \item{A}{Matrix. The lined-up coefficient matrices \eqn{A_j, j=1,\ldots,p} 
#'         for the lagged variables in the panel VAR estimated by mean-group.}
#' \item{B}{Matrix. Placeholder for the structural impact matrix.}
#' \item{beta}{Matrix. The \eqn{((K+n_{d1}) \times r)} cointegrating matrix
#'   of the VAR model if transformed from a rank-restricted VECM.}
#' \item{L.varx}{List of '\code{varx}' objects for the individual estimation results.}
#' \item{L.data}{List of '\code{data.frame}' objects for each individual.}
#' \item{CSD}{List of measures for cross-sectional dependency. 
#'   \code{NULL} if the individual VAR models have been estimated under independence.}
#' \item{args_pvarx}{List of characters and integers 
#'   indicating the estimator and specifications that have been used.}
#' 
#' @references Canova, F., and Ciccarelli, M. (2013): 
#'   "Panel Vector Autoregressive Models: A Survey",
#'   \emph{Advances in Econometrics}, 32, pp. 205-246.
#' @references Hsiao, C. (2014): 
#'   \emph{Analysis of Panel Data}, 
#'   Econometric Society Monographs, Cambridge University Press, 3rd ed.
#' @references Luetkepohl, H. (2005): 
#'   \emph{New Introduction to Multiple Time Series Analysis}, 
#'   Springer, 2nd ed.
#' 
#' @examples
#' data("PCAP")
#' names_k = c("g", "k", "l", "y")  # variable names
#' names_i = levels(PCAP$id_i)      # country names 
#' L.data  = sapply(names_i, FUN=function(i) 
#'   ts(PCAP[PCAP$id_i==i, names_k], start=1960, end=2019, frequency=1), 
#'   simplify=FALSE)
#' R.lags = c(2, 4, 2, 3, 2, 4, 4, 2, 2, 3, 3, 3, 2, 4, 4, 2, 2, 2, 4, 2, 2, 2, 4)
#' names(R.lags) = names_i
#' 
#' @family panel estimation functions
#' @name pvarx
NULL


#' @describeIn pvarx Mean Group (MG) of VAR models in levels.
#' @param t_D List of vectors. The activating break periods \eqn{\tau} 
#'   for the period-specific \link[=as.t_D]{deterministic regressors} 
#'   in \eqn{d_{it}} of the VAR model in levels.
#' @param D List. A single '\code{data.frame}' of common deterministic 
#'   regressors or a list of \eqn{N} '\code{data.frame}' objects of the 
#'   individual \link[=as.t_D]{deterministic regressors} added to \eqn{d_{it}}. 
#'   These customized regressors correspond to '\code{exogen}' in \strong{vars}' 
#'   \code{\link[vars]{VAR}}, which is fixed over bootstrap iterations. 
#' 
#' @references Pesaran, M. H., and Smith R. J. (1995): 
#'   "Estimating Long-Run Relationships from Dynamic Heterogeneous Panels",
#'   \emph{Journal of Econometrics}, 68, pp. 79-113.
#' @references Rebucci, A. (2010): 
#'   "Estimating VARs with Long Stationary Heterogeneous Panels: 
#'   A Comparison of the Fixed Effect and the Mean Group Estimators",
#'   \emph{Economic Modelling}, 27, pp. 1183-1198.
#' 
#' @examples
#' ### MG of VAR by OLS ###
#' R.t_D  = list(t_shift=10)  # common level shift for all countries 
#' R.pvar = pvarx.VAR(L.data, lags=R.lags, type="both", t_D=R.t_D)
#' R.pirf = irf(R.pvar, n.ahead=50)  # MG of individual forecast-error IRF
#' plot(R.pirf)
#' 
#' @export
#' 
pvarx.VAR <- function(L.data, lags, type=c("const", "trend", "both", "none"), 
                      t_D=NULL, D=NULL, n.factors=FALSE, n.iterations=FALSE){
  # define and check
  dim_N   = length(L.data)
  names_i = names(L.data)
  L.y     = aux_check(L.data, "L.data")
  L.dim_p = aux_check(lags, "lags", dim_N, names_i)
  L.t_D   = aux_check(t_D, "t_D", dim_N, names_i)
  L.D     = aux_check(D, "D", dim_N, names_i)
  if( is.null(names_i) ){ names_i = 1:dim_N }
  
  # estimate individual VAR models and their panel mean-group
  L.def = sapply(names_i, FUN=function(i) aux_stackOLS(L.y[[i]], dim_p=L.dim_p[[i]], type=type, t_D=L.t_D[[i]], D=L.D[[i]]), simplify=FALSE)
  R.est = aux_pvar(L.def, n.factors=n.factors, n.iterations=n.iterations)  # OLS if n.factors==0
  R.mgA = aux_MG(R.est$L.varx, idx_par="A")
  R.mgB = aux_MG(R.est$L.varx, idx_par="B")
  
  # return result
  for(i in names_i){ R.est$L.varx[[i]]$args_varx = list(method="solve", type=type, t_D=L.t_D[[i]]) }
  args_pvarx = list(method="MG of VAR", type=type, dim_K=R.est$L.varx[[1]]$dim_K, 
                    dim_N=dim_N, n.factors=n.factors, n.iterations=R.est$n.iterations, w=NULL)
  result = list(A=R.mgA$mean, B=R.mgB$mean, MG_A=R.mgA, MG_B=R.mgB, 
                L.varx=R.est$L.varx, L.data=L.data, CSD=R.est$PCA, args_pvarx=args_pvarx)
  class(result) = "pvarx"
  return(result)
}


#' @describeIn pvarx (Pooled) Mean Group (PMG) of VECM.
#' @param t_D1 List of vectors. The activating break periods \eqn{\tau} 
#'   for the period-specific \link[=as.t_D]{deterministic regressors} 
#'   in \eqn{d_{1,it}}, which are restricted to the cointegration relations. 
#'   Just as in \code{\link{pcoint}}, the 
#'   accompanying lagged regressors are automatically included in \eqn{d_{2,it}}.
#' @param t_D2 List of vectors. The activating break periods \eqn{\tau} 
#'   for the period-specific \link[=as.t_D]{deterministic regressors} 
#'   in \eqn{d_{2,it}}, which are unrestricted.
#' @param D1 List. A single '\code{data.frame}' of common deterministic regressors  
#'   regressors or a list of \eqn{N} '\code{data.frame}' objects of 
#'   individual \link[=as.t_D]{deterministic regressors} added 
#'   to \eqn{d_{1,it}}, which are restricted to the cointegration relations.
#'   Unlike '\code{t_D1}', these customized regressors require potential
#'   accompanying lagged regressors to be manually included in \eqn{d_{2,it}}
#' @param D2 List. A single '\code{data.frame}' of common deterministic 
#'   regressors or a list of \eqn{N} '\code{data.frame}' objects of 
#'   individual \link[=as.t_D]{deterministic regressors} added 
#'   to \eqn{d_{2,it}}, which are unrestricted. 
#'   These customized regressors correspond to '\code{dumvar}' in \strong{urca}'s 
#'   \code{\link[urca]{ca.jo}}, which is fixed over bootstrap iterations.
#' @param dim_r Integer. Common cointegration-rank \eqn{r} of the VECM.
#' @param idx_pool Vector. Position \eqn{k = 1,...,K+n_1} of each variable 
#'   to be pooled using the two-step estimator by Breitung (2005). 
#'   The integer vector specifies throughout heterogeneous coefficients up to 
#'   the uniform upper block \eqn{I_{r}} estimated with the individual estimator 
#'   by Ahn and Reinsel (1990) if exclusively in the interval \eqn{[0,...,r]}. 
#'   Deactivated if \code{NULL} (the default). 
#' 
#' @references Ahn, S. K., and Reinsel (1990): 
#'   "Estimation for Partially Nonstationary Multivariate Autoregressive Models",
#'   \emph{Journal of the American Statistical Association}, 85, pp. 813-823.
#' @references Breitung, J. (2005): 
#'   "A Parametric Approach to the Estimation of Cointegration Vectors in Panel Data",
#'   \emph{Econometric Reviews}, 24, pp. 151-173.
#' @references Johansen, S. (1996): 
#'   \emph{Likelihood-based Inference in Cointegrated Vector Autoregressive Models}, 
#'   Advanced Texts in Econometrics, Oxford University Press, USA.
#' @references Pesaran, M. H., Shin, Y, and Smith R. J. (1999): 
#'   "Pooled Mean Group Estimation of Dynamic Heterogeneous Panels",
#'   \emph{Journal of the American Statistical Association}, 94, pp. 621-634.
#' 
#' @examples
#' ### Pooled MG of rank-restricted VAR ###
#' R.pvec = pvarx.VEC(L.data, lags=R.lags, dim_r=2, idx_pool=1:4, type="Case4")
#' R.pirf = irf(R.pvec, n.ahead=50)  # MG of individual forecast-error IRF
#' plot(R.pirf)
#' 
#' @export
#' 
pvarx.VEC <- function(L.data, lags, dim_r, 
                      type=c("Case1", "Case2", "Case3", "Case4", "Case5"), 
                      t_D1=NULL, t_D2=NULL, D1=NULL, D2=NULL, 
                      idx_pool=NULL, n.factors=FALSE, n.iterations=FALSE){
  # define and check
  dim_N   = length(L.data)
  names_i = names(L.data)
  L.y     = aux_check(L.data, "L.data")
  L.dim_p = aux_check(lags, "lags", dim_N, names_i)
  L.t_D1  = aux_check(t_D1, "t_D1", dim_N, names_i)
  L.t_D2  = aux_check(t_D2, "t_D2", dim_N, names_i)
  L.D1    = aux_check(D1, "D1", dim_N, names_i)
  L.D2    = aux_check(D2, "D1", dim_N, names_i)
  if( is.null(names_i) ){ names_i = 1:dim_N }
  
  if(is.null(idx_pool)){
    name_MG = "MG of rank-restricted VAR"
    name_ec = "RRR"
  }else{
    name_MG = "PMG of rank-restricted VAR" 
    name_ec = "GLS"
  }
  
  # estimate individual VECM (with option to pool cointegrating vectors) and their panel mean-group
  L.def = sapply(names_i, FUN=function(i) aux_stackRRR(L.y[[i]], dim_p=L.dim_p[i], type=type, t_D1=L.t_D1[[i]], t_D2=L.t_D2[[i]], D1=L.D1[[i]], D2=L.D2[[i]]), simplify=FALSE)
  R.est = aux_pvec(L.def, L.beta=NULL, dim_r=dim_r, idx_pool=idx_pool, n.factors=n.factors, n.iterations=n.iterations)  # RRR if n.factors==0 and is.null(idx_pool)
  R.mgB = aux_MG(R.est$L.varx, idx_par="B")
  
  # (pooled) mean-group estimation, from Pesaran et al. (1999: ch.4.3, p.626)
  vecm = list(
    alpha = aux_MG(R.est$L.varx, idx_par="alpha"),
    beta  = aux_MG(R.est$L.varx, idx_par="beta"),  # flexible towards idx_pool 
    GAMMA = aux_MG(R.est$L.varx, idx_par="GAMMA"))
  vecm$PI = vecm$alpha$mean %*% t(vecm$beta$mean)
  rvar = aux_vec2var(PI=vecm$PI, GAMMA=vecm$GAMMA$mean, dim_p=max(L.dim_p))  # rank-restricted VAR model in levels
  
  # return result
  for(i in names_i){ R.est$L.varx[[i]]$args_varx = list(method=name_ec, type=type, t_D1=L.t_D1[[i]], t_D2=L.t_D2[[i]]) }
  args_pvarx = list(method=name_MG, type=type, dim_r=dim_r, dim_K=R.est$L.varx[[1]]$dim_K, dim_N=dim_N, 
                    idx_pool=idx_pool, n.factors=n.factors, n.iterations=n.iterations, w=NULL)
  result = list(A=rvar$A, B=R.mgB$mean, beta=vecm$beta$mean, MG_VECM=vecm, MG_B=R.mgB, 
                L.varx=R.est$L.varx, L.data=L.data, CSD=R.est$PCA, args_pvarx=args_pvarx)
  class(result) = "pvarx"
  return(result)
}


#### S3 methods for objects of class 'pvarx' ####
#' @export
print.pvarx <- function(x, ...){
  # define
  x = as.pvarx(x)
  
  # create headers
  header_0 = "### Panel VAR Model ### \n"
  header_A = "\nEstimated mean-group coefficients: \n"
  header_B = "\nEstimated mean-group B Matrix (unique decomposition of the covariance matrix): \n"
  
  # print
  cat(header_0)
  cat(header_A)
  print(x$A)
  
  if(identical(unname(x$B), diag(x$args_pvarx$dim_K))){
    cat("\nModel: Reduced-form VAR \n")
  }else{
    cat(header_B)
    print(x$B)
  }
}


#' @export
summary.pvarx <- function(object, ..., modulus=TRUE, digits=3){
  # define
  object  = as.pvarx(object)
  names_i = names(object$L.varx)
  dim_N   = length(object$L.varx)
  dim_pN  = max(sapply(object$L.varx, FUN=function(i) i$dim_p))  # maximum lag-order of all individuals
  n.char_max = max(c(nchar(names_i), 4)) + 1
  
  L.stab = sapply(object$L.varx, simplify=FALSE, FUN=function(i) 
    eigen(aux_var2companion(i$A, dim_p=i$dim_p))$values)
  R.stab = if( any(is.na(object$A)) ){ NULL }else{ 
    eigen(aux_var2companion(object$A, dim_p=dim_pN))$values }
  
  # print
  print.pvarx(object, ...)
  cat("\n", paste0(c("Specifications", rep(".", times=20))), sep="", "\n")
  cat("Method: ", object$args_pvarx$method, "\n")
  cat("Number of iterations: ", as.integer(object$args_pvarx$n.iterations), "\n")
  cat("Number of individuals: ", object$args_pvarx$dim_N, "\n")
  
  if(!is.null(object$args_pvarx$n.factors)){
    cat("Number of common factors: ", as.integer(object$args_pvarx$n.factors), "\n") }
  if(!is.null(object$args_pvarx$type)){
    cat("Deterministic term: ", object$args_pvarx$type, "\n") }
  if(!is.null(object$beta)){
    cat("Cointegrating vectors: \n")
    print(object$beta)}
  
  cat("\n", paste0(c("Eigenvalues of the companion matrices", rep(".", times=5))), sep="", "\n")
  for(i in 1:dim_N){
    R.cat = if(modulus) Mod(L.stab[[i]]) else L.stab[[i]]
    R.cat = format(round(R.cat, digits=digits), nsmall=digits)
    R.cat = c(names_i[i], rep("", times=n.char_max-nchar(names_i[i])), R.cat)
    cat(R.cat, "\n", sep=" ")
    rm(R.cat) }
  if(!is.null(R.stab)){
    R.cat = if(modulus) Mod(R.stab) else R.stab
    R.cat = format(round(R.cat, digits=digits), nsmall=digits)
    R.cat = c("-MG-", rep("", times=n.char_max-4), R.cat)
    cat(R.cat, "\n", sep=" ") }
}



#### Generic functions and their S3 methods ####
# CRAN manual: https://cran.r-project.org/doc/manuals/R-exts.html#Generic-functions-and-methods
# Roxygen: https://r-pkgs.org/man.html#man-s3

#' @title Coerce into a '\code{pvarx}' object
#' @description Coerce into a '\code{pvarx}' object. On top of the parent class 
#'   '\code{pvarx}', the child class '\code{pid}' is imposed if the input object 
#'   to be transformed contains a panel SVAR model. 
#' @details \code{\link{as.pvarx}} is used as an intermediary in the \strong{pvars} 
#'   functions to achieve compatibility with different classes of panel VAR objects.
#'   If the user wishes to extend this compatibility with further classes, she 
#'   may simply specify accordant \code{\link{as.pvarx}}-methods instead of 
#'   altering the original \strong{pvars} function.
#'   
#' @param x A panel VAR object to be transformed.
#' @param ... Additional arguments to be passed to or from methods.
#' @param w Numeric, logical, or character vector. 
#'   \eqn{N} numeric elements weighting the individual coefficients, or 
#'   names or \eqn{N} logical elements selecting a subset from the 
#'   individuals \eqn{i = 1, \ldots, N} for the MG estimation. If \code{NULL} 
#'   (the default), all \eqn{N} individuals are included without weights.
#' 
#' @return  A list of class '\code{pvarx}'. Objects of this class contain the elements:
#' \item{A}{Matrix. The lined-up coefficient matrices \eqn{A_j, j=1,\ldots,p} 
#'   for the lagged variables in the panel VAR.}
#' \item{B}{Matrix. The \eqn{(K \times S)} structural impact matrix of the panel SVAR model 
#'   or an identity matrix \eqn{I_K} as a placeholder for the unidentified VAR model.}
#' \item{beta}{Matrix. The \eqn{((K+n_{d1}) \times r)} cointegrating matrix of the VAR model 
#'   if transformed from a rank-restricted VECM.}
#' \item{L.varx}{List of \code{varx} objects for the individual estimation results.}
#' \item{args_pvarx}{List of characters and integers 
#'   indicating the estimator and specifications that have been used.}
#' \item{args_pid}{List of characters and integers 
#'   indicating the identification methods and specifications that have been used. 
#'   This element is specific to the child-class 'pid' for panel SVAR models, 
#'   that inherit from parent-class 'pvarx' for any panel VAR model.}
#' 
#' @examples
#' data("PCAP")
#' names_k = c("g", "k", "l", "y")  # variable names
#' names_i = levels(PCAP$id_i)      # country names 
#' L.data  = sapply(names_i, FUN=function(i) 
#'   ts(PCAP[PCAP$id_i==i, names_k], start=1960, end=2019, frequency=1), 
#'   simplify=FALSE)
#'   
#' L.vars = lapply(L.data, FUN=function(x) vars::VAR(x, p=2, type="both"))
#' as.pvarx(L.vars)
#' 
#' @export
#' 
as.pvarx <- function(x, w=NULL, ...) UseMethod("as.pvarx")


#' @method as.pvarx default
#' @export
as.pvarx.default <- function(x, w=NULL, ...){
  if(is.list(x)){
    # define
    L.varx  = lapply(x, FUN=function(x_i) as.varx(x_i))
    L.data  = lapply(L.varx, FUN=function(x_i) t(x_i$y))
    L.argID = lapply(L.varx, FUN=function(x_i) x_i$args_id)
    L.dim_r = lapply(L.varx, FUN=function(x_i) x_i$dim_r)
    L.dim_p = sapply(L.varx, FUN=function(x_i) x_i$dim_p)
    L.dim_K = sapply(L.varx, FUN=function(x_i) x_i$dim_K)
    
    dim_N = length(L.varx)
    dim_K = L.dim_K[[1]]
    dim_r = L.dim_r[[1]]
    is_r  = any(sapply(L.dim_r, FUN=function(x_i) !is.null(x_i)))
    is_id = any(sapply(L.varx,  FUN=function(x_i) inherits(x_i, "id")))
    is_iv = FALSE
    
    # gather arguments of identification procedure
    if(is_id){
      args_pid = L.argID[[1]]
      args_pid$combine = "indiv"
      is_iv = all(sapply(L.argID, FUN=function(x_i) isTRUE(x_i$method == "Proxy")))  # isTRUE: FALSE if NULL
      if(is_iv){
        args_pid$iv = NULL
        args_pid$L.iv = lapply(L.argID, FUN=function(x_i) x_i$iv)
        for(i in 1:dim_N){ L.argID[[i]]$iv = NULL }
        L.dim_L = sapply(args_pid$L.iv, FUN=function(x_i) nrow(x_i))
        dim_L   = L.dim_L[[1]]
      }}
    
    # check homogeneity
    if( !all(dim_K == L.dim_K) ){ 
      stop("The number of variables 'K' must be the same for all individuals.") }
    if(is_r){  if( !all(sapply(L.dim_r, FUN=function(x_i) identical(x_i, dim_r))) ){ 
      stop("The cointegration rank-restriction 'r' must be the same for all individuals.") }}
    if(is_id){ if( !all(sapply(L.argID, FUN=function(x_i) identical(x_i, L.argID[[1]]))) ){ 
      stop("The identification procedure must be the same for all individuals.") }}
    if(is_iv){ if( !all(dim_L == L.dim_L) ){ 
      stop("The number of proxy variables 'L' must be the same for all individuals.") }}
    
    # add panel estimates
    if(is.null(dim_r)){
      # estimate panel mean-group for VAR
      R.mgA = aux_MG(L.varx, w=w, idx_par="A")
      R.mgB = aux_MG(L.varx, w=w, idx_par="B")
      
      # construct object corresponding to pvarx.VAR()
      args_pvarx = list(method="MG of VAR", dim_K=dim_K, dim_N=dim_N, 
                        n.factors=FALSE, n.iterations=FALSE, w=w)
      result = list(A=R.mgA$mean, B=R.mgB$mean, MG_A=R.mgA, MG_B=R.mgB, 
                    L.varx=L.varx, L.data=L.data, args_pvarx=args_pvarx)
      
    }else{
      # estimate panel mean-group for rank-restricted VAR
      vecm = list(
        alpha = aux_MG(L.varx, w=w, idx_par="alpha"),
        beta  = aux_MG(L.varx, w=w, idx_par="beta"),
        GAMMA = aux_MG(L.varx, w=w, idx_par="GAMMA"))
      vecm$PI = vecm$alpha$mean %*% t(vecm$beta$mean)
      rvar  = aux_vec2var(PI=vecm$PI, GAMMA=vecm$GAMMA$mean, dim_p=max(L.dim_p))
      R.mgB = aux_MG(L.varx, w=w, idx_par="B")
      
      # construct objects corresponding to pvarx.VEC()
      args_pvarx = list(method="MG of rank-restricted VAR", 
                        dim_r=dim_r, dim_K=dim_K, dim_N=dim_N, 
                        n.factors=FALSE, n.iterations=FALSE, w=w)
      result = list(A=rvar$A, B=R.mgB$mean, beta=vecm$beta$mean, MG_VECM=vecm, MG_B=R.mgB,
                    L.varx=L.varx, L.data=L.data, args_pvarx=args_pvarx)
    }
    
    # return result
    if(is_id){
      result$args_pid = args_pid
      class(result) = c("pid", "pvarx")
    }else{
      class(result) = "pvarx" }
    return(result)
    
  }else{
    stop("Argument is not a panel VAR object of suitable class!")
  }
}


#' @method as.pvarx pvarx
#' @export
as.pvarx.pvarx <- function(x, w=NULL, ...){
  if(is.null(w)){
    return(x)
  
  }else if(identical(w, x$args_pvarx$w)){
    return(x)
   
  # add panel mean-group for VAR under new 'w' 
  }else if(is.null(x$args_pvarx$dim_r)){
    x$MG_A = aux_MG(x$L.varx, w=w, idx_par="A")
    x$MG_B = aux_MG(x$L.varx, w=w, idx_par="B")
    x$A = x$MG_A$mean
    x$B = x$MG_B$mean
    x$args_pvarx$w = w
    return(x)
  
  # add panel mean-group for rank-restricted VAR under new 'w'
  }else{
    dim_p = max(sapply(x$L.varx, FUN=function(x_i) x_i$dim_p))
    vecm  = list(
      alpha = aux_MG(x$L.varx, w=w, idx_par="alpha"),
      beta  = aux_MG(x$L.varx, w=w, idx_par="beta"),
      GAMMA = aux_MG(x$L.varx, w=w, idx_par="GAMMA"))
    vecm$PI = vecm$alpha$mean %*% t(vecm$beta$mean)
    x$MG_B  = aux_MG(x$L.varx, w=w, idx_par="B")
    x$MG_VECM = vecm
    x$A = aux_vec2var(PI=vecm$PI, GAMMA=vecm$GAMMA$mean, dim_p=dim_p)$A
    x$B = x$MG_B$mean
    x$beta = vecm$beta$mean
    x$args_pvarx$w = w
    return(x)
  }
}


#' @method as.pvarx sboot2
#' @export
as.pvarx.sboot2 <- function(x, w=NULL, ...){
  result = as.pvarx(x$pvarx, w=w, ...)
  return(result)
}


### For example TODO: 
### #' @method as.pvarx pvarfeols, pvargmm, pvarhk 
### #' @importFrom panelvar ###


