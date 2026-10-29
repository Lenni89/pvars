

#### Generic functions and their S3 methods ####
# CRAN manual: https://cran.r-project.org/doc/manuals/R-exts.html#Generic-functions-and-methods
# Roxygen: https://r-pkgs.org/man.html#man-s3

#' @title Deterministic regressors in \strong{pvars}
#' @description Deterministic regressors can be specified via the arguments of
#'   the \strong{\emph{conventional}} '\code{type}', \strong{\emph{customized}} 
#'   '\code{D}', and \strong{\emph{period-specific}} '\code{t_D}'. 
#'   While '\code{type}' is a single character and 
#'   '\code{D}' a data matrix of dimension \eqn{(n_{\bullet} \times (p+T))}, 
#'   the specifications for \eqn{\tau} in the list '\code{t_D}' are more complex  
#'   and therefore preventively checked by \code{\link{as.t_D}}.
#'   
#' @param x A list of vectors for \eqn{\tau} to be checked. Since \code{'x'} is
#'   only checked, Section "Value" describes both input and output.
#' @param ... Additional arguments to be passed to or from methods.
#' 
#' @return A list of class '\code{t_D}' specifying \eqn{\tau}. 
#'   Objects of this class can exclusively contain the elements:
#' \item{t_break}{Vector of integers. The activating periods for 
#'   trend breaks \eqn{d = [\ldots, 0, 0, 1, 2, 3, \ldots}].}
#' \item{t_shift}{Vector of integers. The activating periods for 
#'   shifts in the constant \eqn{d = [\ldots, 0, 0, 1, 1, 1, \ldots}].}
#' \item{t_impulse}{Vector of integers. The activating periods for 
#'   single impulses \eqn{d = [\ldots, 0, 0, 1, 0, 0, \ldots}].}
#' \item{t_blip}{Vector of integers. The activating periods for 
#'   blips \eqn{d = [\ldots, 0, 0, 1, -1, 0, \ldots}].}
#' \item{n.season}{Integer. The number of seasons.}
#' 
#' @section Reference Time Interval:
#'   The complete time series (i.e. including the presample) constitutes 
#'   the reference time interval. Accordingly, '\code{D}' contains \eqn{p+T} 
#'   observations, and '\code{t_D}' contains the positions of activating 
#'   periods \eqn{\tau} in \eqn{1,\ldots,(p+T)}. In a balanced panel 
#'   \eqn{p_i+T_i = T^*}, the same date implies the same \eqn{\tau} in 
#'   \eqn{1,\ldots,T^*}, as shown in the example for \code{\link{pcoint.CAIN}}. 
#'   However, in an unbalanced panel, the same date can imply different 
#'   \eqn{\tau} across \eqn{i} in accordance with the individual time interval 
#'   \eqn{1,\ldots,(p_i+T_i)}. Note that across the time series in '\code{L.data}', 
#'   it is the last observation in each data matrix that refers to the same date.
#' 
#' @section Conventional Type: 
#'   An overview is given here and a detailed explanation in the package vignette. 
#' \itemize{
#' \item{\strong{type (VAR)} is specified in VAR models just as in \strong{vars}' \code{\link[vars]{VAR}}, 
#'       namely by a '\code{const}', a linear '\code{trend}', '\code{both}', or '\code{none}' of those.}
#' \item{\strong{type_SL} is used in the 'additive' SL procedure for testing the cointegration rank only,
#'   which removes the mean ('\code{SL_mean}') or mean and linear trend 
#'   ('\code{SL_trend}') by GLS-detrending.}
#' \item{\strong{type (VECM)} is used in the 'innovative' Johansen procedure 
#'   for testing the cointegration rank and estimating the VECM. In accordance 
#'   with Juselius (2007, Ch.6.3), the available model specifications are:
#'   '\code{Case1}' for none, 
#'   '\code{Case2}' for a constant in the cointegration relation,
#'   '\code{Case3}' for an unrestricted constant, 
#'   '\code{Case4}' for a  linear trend in the cointegration relation and an unrestricted constant, or 
#'   '\code{Case5}' for an unrestricted constant and linear trend.}
#' }
#' @references Juselius, K. (2007): 
#'   \emph{The Cointegrated VAR Model: Methodology and Applications}, 
#'   Advanced Texts in Econometrics, Oxford University Press, USA, 2nd ed.
#' 
#' @examples
#' t_D = list(t_impulse=c(10, 20, 35), t_shift=10)
#' as.t_D(t_D)
#' 
#' @export
#' 
as.t_D <- function(x, ...) UseMethod("as.t_D")


#' @method as.t_D default
#' @export
as.t_D.default <- function(x, ...){
  if(inherits(x, "t_D")){
    return(x)
  }else if(is.null(x)){
    x = list()
  }else if(is.list(x)){
    # check
    names_x = names(x)
    names_D = c("t_break", "t_shift", "t_impulse", "t_blip", "n.season")
    idx_rm  = !(names_x %in% names_D)
    if(any(idx_rm)){
      warn_rm = paste0(names_x[idx_rm], collapse=", ")
      warning("Unrecognized elements in 't_D' have been removed: ", warn_rm)
      x = x[idx_rm]
    }
  }
    
  # return
  class(x) = "t_D"
  return(x)
}


### TODO: use dates for t_D 
### #' @method as.t_D Date 
### #' @importFrom svars getStructuralBreak



#############################
###  AUXILIARY FUNCTIONS  ###
#############################
#
# The following functions serve as supporting modules 
# nested in the calling functions. Notation within 
# each function originates from the considered packages.


# check arguments of data matrices
aux_asDataMatrix <- function(x, names_x=NULL){
  # check dimensions
  result = as.matrix(x)
  if(nrow(result) > ncol(result)){
    result = t(result)  # supposed to be (nrow x T) regardless of input
  } ### Note that a transposed time-series object loses its ts-specifications.
  
  # check names
  if(is.null(rownames(result))){
    rownames(result) = paste0(names_x, 0:nrow(result))[-1] 
  }  ### Note that [-1] lets (0 x T) matrices pass through.
  
  # return result
  return(result)
}


# check arguments of transposed data matrices (for scale() in PANIC and test.normality)
aux_asDataMatrixTr <- function(x, names_x=NULL){
  # check dimensions
  result = as.matrix(x)
  if(nrow(result) < ncol(result)){
    result = t(result)  # supposed to be (T x ncol) regardless of input
  } ### Note that a transposed time-series object loses its ts-specifications.
  
  # check names
  if(is.null(colnames(result))){
    colnames(result) = paste0(names_x, 0:ncol(result))[-1] 
  }  ### Note that [-1] lets (T x 0) matrices pass through.
  
  # return result
  return(result)
}


# check arguments of the panel functions
aux_check <- function(x, type_arg, dim_N=NULL, names_i=NULL, tr=FALSE){
  ### This function checks initial data- and specification-arguments. 
  ### Higher-level input objects are checked by as.varx() and as.pvarx().
  
  # check data argument "L.data"
  if(type_arg == "L.data"){
    if(is.data.frame(x)){
      stop("Argument 'L.data' must be a list of N individual data.frame objects.")
      ### TODO: accept objects of class pdata.frame and automatically construct L.data
    }else if(is.list(x)){
      if(tr){
        result = lapply(x, FUN=function(x_i) aux_asDataMatrixTr(x_i, "y"))
        L.ncol = sapply(result, FUN=function(x_i) ncol(x_i))
        if( !all(L.ncol[1] == L.ncol) ){ 
          stop("The number of variables in 'L.data' must be the same for all individuals.") 
        }
      }else{
        result = lapply(x, FUN=function(x_i) aux_asDataMatrix(x_i, "y"))
        L.nrow = sapply(result, FUN=function(x_i) nrow(x_i))
        if( !all(L.nrow[1] == L.nrow) ){ 
          stop("The number of variables in 'L.data' must be the same for all individuals.") 
        }  
      }
    }else{
      stop("Argument 'L.data' must be a list of N individual data.frame objects.")
    }
  
  # check data arguments for variables
  }else if(type_arg %in% c("y", "x", "iv")){
    if(is.matrix(x) | is.data.frame(x)){
      result = aux_asDataMatrix(x, type_arg)
      result = replicate(dim_N, result, simplify=FALSE)
    }else if(length(x) == dim_N){
      result = lapply(x, FUN=function(x_i) aux_asDataMatrix(x_i, type_arg))
      L.nrow = sapply(result, FUN=function(x_i) nrow(x_i))
      if( !all(L.nrow[1] == L.nrow) ){ 
        tmp = switch(type_arg, y="endogenous variables 'K'", x="exogenous variables 'L'", iv="proxies 'L'")
        stop("The number of ",  tmp, " must be the same for all individuals.") 
      }
    }else{
      stop("Argument '", type_arg, "' must be either a single data.frame object or a list of N individual data.frame objects.")
    }
    
  # check data arguments customized for deterministic regressors
  }else if(type_arg %in% c("D", "D1", "D2")){
    if(is.null(x)){
      return(x)
    }else if(is.matrix(x) | is.data.frame(x)){
      result = aux_asDataMatrix(x, "d.d")
      result = replicate(dim_N, result, simplify=FALSE)
    }else if(length(x) == dim_N){
      idx_i  = 1:dim_N; names(idx_i) = names(x)
      result = lapply(idx_i, FUN=function(i) aux_asDataMatrix(x[[i]], paste0("d_", i, ".d")))
    }else{   
      stop("Argument '", type_arg, "' must be either a single data.frame object or a list of N individual data.frame objects.")
    }
    
  # check specification argument "lags"
  }else if(type_arg == "lags"){
    if(length(x) == dim_N){
      result = x
    }else if(length(x) == 1){
      result = rep(x, dim_N)
    }else{
      stop("Argument 'lags' must be either a single integer or a vector of N integers specific to each individual.")
    }
    
  # check specification argument for period-specific deterministic regressors
  }else if(type_arg %in% c("t_D", "t_D1", "t_D2")){
    if(is.null(x)){
      result = replicate(dim_N, list(), simplify=FALSE)
    }else if(length(x) == dim_N){
      result = lapply(x, FUN=function(x_i) as.t_D(x_i))
    }else{
      result = replicate(dim_N, as.t_D(x), simplify=FALSE)
    }
  }
    
  # check individual names
  names_x = names(result)
  if(is.null(names_x)){
    names(result) = names_i
    
  }else if(!is.null(names_i)){
    if(!identical(names_i, names_x)){
      warning("Arguments 'L.data' and ", type_arg, " have mismatching names for individuals.")
    }
  }
  
  # return result
  return(result)
}


# assign objects from "pvarx" to function environment
aux_assign_pvarx <- function(object, w=NULL){
  ### Note: For extending this function by other classes, just specify "as.pvarx" methods for your panel VAR objects. ###
  # define
  object = as.pvarx(object, w=w)
  L.varx = object$L.varx
  
  # gather individual results from "varx"
  L.dim_p = sapply(L.varx, FUN=function(i) i$dim_p)
  L.dim_T = sapply(L.varx, FUN=function(i) i$dim_T)  # number of observations without presample
  L.resid = lapply(L.varx, FUN=function(i) i$resid)
  L.beta  = lapply(L.varx, FUN=function(i) i$beta)
  
  L.D1 = lapply(L.varx, FUN=function(i) i$D1)
  L.D2 = lapply(L.varx, FUN=function(i) i$D2)
  L.D  = lapply(L.varx, FUN=function(i) i$D)
  L.A  = lapply(L.varx, FUN=function(i) i$A)
  L.B  = lapply(L.varx, FUN=function(i) i$B)
  
  # assign arguments
  assign("args_pvarx", object$args_pvarx, envir = parent.frame())
  assign("args_pid",   object$args_pid,   envir = parent.frame())
  
  # assign MG results
  assign("A_MG", object$A,    envir = parent.frame())  # matrix of VAR coefficients (K x (K*p))
  assign("B_MG", object$B,    envir = parent.frame())  # structural impact matrix (K x S)
  assign("beta", object$beta, envir = parent.frame())  # cointegrating vectors (r x K+L+n1)
  
  # assign dimensions
  assign("dim_r", object$args_pvarx$dim_r, envir = parent.frame())  # cointegration rank
  assign("dim_N", object$args_pvarx$dim_N, envir = parent.frame())  # number of individuals
  assign("dim_K", object$args_pvarx$dim_K, envir = parent.frame())  # number of endogenous variables
  assign("dim_S", ncol(object$B), envir = parent.frame())  # number of shocks
  
  # assign individual results
  assign("L.varx",  L.varx,  envir = parent.frame())  # individual "varx" objects
  assign("L.dim_p", L.dim_p, envir = parent.frame())  # lag-order of the VAR model in levels
  assign("L.dim_T", L.dim_T, envir = parent.frame())  # number of observations without presample
  assign("L.resid", L.resid, envir = parent.frame())  # residual matrix (K x T)
  assign("L.beta",  L.beta,  envir = parent.frame())  # cointegrating vectors (r x K+L)
  
  assign("L.A",  L.A,  envir = parent.frame())  # matrix of VAR coefficients (K x (n+K*p)), incl. deterministic 
  assign("L.B",  L.B,  envir = parent.frame())  # structural impact matrix (K x S)
  assign("L.D",  L.D,  envir = parent.frame())  # matrix of deterministic variables (n x (p+T))
  assign("L.D1", L.D1, envir = parent.frame())  # matrix of restricted deterministic variables (n x (p+T))
  assign("L.D2", L.D2, envir = parent.frame())  # matrix of unrestricted deterministic variables (n x (p+T))
  
  # return 'pvarx' object
  return(object)
}


# assign objects from "varx" to function environment
aux_assign_varx <- function(object){
  ### Note: For extending this function by other classes, just specify "as.varx" methods for your VAR objects. ###
  # define
  R.varx = as.varx(object)
  
  # assign results
  assign("A", R.varx$A, envir = parent.frame())  # matrix of VAR coefficients (K x (n+K*p)), incl. deterministic 
  assign("B", R.varx$B, envir = parent.frame())  # structural impact matrix (K x S)
  assign("y", R.varx$y, envir = parent.frame())  # matrix of endogenous variables (K x (p+T))
  assign("x", R.varx$x, envir = parent.frame())  # matrix of exogenous variables (L x (p+T))
  assign("D",  R.varx$D,  envir = parent.frame())  # matrix of deterministic variables (n x (p+T))
  assign("D1", R.varx$D1, envir = parent.frame())  # matrix of restricted deterministic variables (n x (p+T))
  assign("D2", R.varx$D2, envir = parent.frame())  # matrix of unrestricted deterministic variables (n x (p+T))
  
  assign("resid", R.varx$resid, envir = parent.frame())  # residual matrix (K x T)
  assign("OMEGA", R.varx$OMEGA, envir = parent.frame())  # MLE covariance matrix of residuals
  assign("SIGMA", R.varx$SIGMA, envir = parent.frame())  # OLS covariance matrix of residuals
  assign("beta",  R.varx$beta,  envir = parent.frame())  # cointegrating vectors (r x K+L)
  
  # assign dimensions
  assign("dim_r", R.varx$dim_r, envir = parent.frame())  # cointegration rank
  assign("dim_p", R.varx$dim_p, envir = parent.frame())  # lag-order of the VAR model in levels
  assign("dim_T", R.varx$dim_T, envir = parent.frame())  # number of observations without presample
  assign("dim_K", R.varx$dim_K, envir = parent.frame())  # number of endogenous variables
  assign("dim_S", ncol(R.varx$B), envir = parent.frame())  # number of shocks
  
  # assign arguments
  assign("t_D1", R.varx$t_D1, envir = parent.frame())  # period-specific deterministic variables (restricted)
  assign("t_D2", R.varx$t_D2, envir = parent.frame())  # period-specific deterministic variables (unrestricted)
  assign("args_varx", R.varx$args_varx, envir = parent.frame())
  assign("args_id",   R.varx$args_id,   envir = parent.frame())
  
  # return 'varx' object
  return(R.varx)
}


# remove all lagged impulse dummies which control for non-linearities
aux_rm_Dnl <- function(x, dim_p, t_D1, t_D2, MARGIN=1){
  # define
  t_Dnl = c(t_D1$t_break, t_D1$t_shift)
  
  if(is.null(t_Dnl)){
    return(x)  # skip all VAR objects without these dummies
    
  }else{
    t_Dnl = sapply(t_Dnl, FUN=function(x) x + 0:(dim_p-1))
    names_Dnl = paste0("d.im", t_Dnl)
    
    if(!is.null(t_D2$t_imp)){
      names_imp = paste0("d.im", t_D2$t_imp)
      names_Dnl = names_Dnl[names_Dnl != names_imp]  # keep predefined impulse dummies!
    }
    
    if(MARGIN==1){
      idx_Dnl = which(rownames(x) %in% names_Dnl)
      return(x[-idx_Dnl, , drop=FALSE])
    }else if(MARGIN==2){
      idx_Dnl = which(colnames(x) %in% names_Dnl)
      return(x[ ,-idx_Dnl, drop=FALSE])
    }else{
      stop("Argument 'MARGIN' must be an integer of either 1 for rows or 2 for columns.")
    }
  }
}


# collect individual parameters in array
aux_getPAR <- function(L.varx, idx_par, idx_i=TRUE, A.fill=NA){
  if(idx_par %in% c("A", "GAMMA") ){
    # VAR coefficient matrices with different lag-orders and dummy effects
    if(idx_par == "A"){  # ... in levels
      L.coeff = lapply(L.varx[idx_i], FUN=function(i) i[[idx_par]])
      L.dim_p = sapply(L.varx[idx_i], FUN=function(i) i$dim_p)
      
    }else if(idx_par == "GAMMA"){  # ... in first differences
      L.coeff = lapply(L.varx[idx_i], FUN=function(i) i$VECM[[idx_par]])
      L.dim_p = sapply(L.varx[idx_i], FUN=function(i) i$dim_p)-1 
    }
    
    dim_N  = length(L.coeff)  # group size
    dim_K  = nrow(L.coeff[[1]])
    dim_pN = max(L.dim_p)  # maximum lag-order of all individuals
    idx_pN = which.max(L.dim_p)
    
    names_i = names(L.coeff)
    names_k = rownames(L.coeff[[1]])
    names_p = aux_tail(colnames(L.coeff[[idx_pN]]), n=dim_K*dim_pN)
    
    # find common deterministic variables
    L.names = lapply(1:dim_N, FUN=function(i) aux_head(colnames(L.coeff[[i]]), dim_K*L.dim_p[i], rm=TRUE))
    names_d = Reduce(intersect, L.names)
    dim_n   = length(names_d)
    
    # collect in array
    A.dims = c(dim_K, (dim_n + dim_K*dim_pN), dim_N)
    A.name = list(names_k, c(names_d, names_p), names_i)
    result = array(A.fill, dim=A.dims, dimnames=A.name)
    for(i in 1:dim_N){
      coeff_i = L.coeff[[i]]
      idx_slp = aux_tail(1:ncol(coeff_i), dim_K*L.dim_p[i])  # indices of the slope coefficients (without relying on names)
      idx_dtr = match(names_d, colnames(coeff_i))  # indices of the coefficients for common deterministic variables
      idx_mat = c(idx_dtr, idx_slp)  # column indices for the matrix
      idx_arr = 0:(dim_n + dim_K*L.dim_p[i])  # column indices for the array
      result[ , idx_arr, i] = coeff_i[ , idx_mat, drop=FALSE]
    }
    
  }else if(idx_par == "beta"){
    # cointegrating vectors with K variables Y_{i,t-1} and heterogeneous nrow(D1)
    L.coeff = lapply(L.varx[idx_i], FUN=function(i) i[[idx_par]])
    dim_N   = length(L.coeff)  # group size
    dim_K   = nrow(L.varx[[1]]$A)
    idx_k   = 1:dim_K  # indices of the slope coefficients (without relying on names)
    names_i = if( !is.null(names(L.coeff)) ){ names(L.coeff) }else{ 1:dim_N }
    
    # find common variables
    L.names = lapply(L.coeff, FUN=function(i) rownames(i))
    names_d = Reduce(intersect, L.names)
    L.idx_d = lapply(L.names, FUN=function(i) unique(c(idx_k, match(names_d, i))))  # indices of the coefficients
    
    # collect in array
    result = sapply(names_i, FUN=function(i) L.coeff[[i]][L.idx_d[[i]], , drop=FALSE], simplify="array")
    
  }else if(idx_par == "alpha"){
    # loading matrices of equal dimensions in VECM
    result = sapply(L.varx[idx_i], FUN=function(i) i$VECM[[idx_par]], simplify="array")
    
  }else{
    # parameter matrices of equal dimensions, e.g., structural impact coefficients
    result = sapply(L.varx[idx_i], FUN=function(i) i[[idx_par]], simplify="array")
  }
  
  # return result
  return(result)
}


# keep or remove the head or tail of a vector
### Note: In contrast to head() and tail(), these functions are capable ... 
### of removing nothing ("-0") and do not trim exceeding n > length(x) silently.
aux_head <- function(x, n, rm=FALSE){
  # define
  if(rm){
    dim_old = length(x)
    idx = seq_len(dim_old-n)
  }else{
    idx = seq_len(n)
  }
  
  # return result
  x[idx]
}


aux_tail <- function(x, n, rm=FALSE){
  # define
  dim_old = length(x)
  if(rm){
    idx = seq.int(to=dim_old, length.out=dim_old-n)
  }else{
    idx = seq.int(to=dim_old, length.out=n)
  }
  
  # return result
  x[idx]
}


