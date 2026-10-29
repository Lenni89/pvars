

#' @title Estimation of a Vector Error Correction Model
#' @description Estimates a VECM under a given cointegration-rank restriction or cointegrating vectors.
#' 
#' @param y Matrix. A \eqn{(K \times (p+T))} data matrix of the \eqn{K} endogenous time series variables.
#' @param dim_p Integer. Lag-order \eqn{p} for the endogenous variables \code{y}.
#' @param x Matrix. A \eqn{(L \times (p+T))} data matrix of the \eqn{L} weakly exogenous time series variables.
#' @param dim_q Integer. Lag-order \eqn{q} for the distributed lag of the weakly exogenous variables \code{x}.
#' @param dim_r Integer. Cointegration-rank \eqn{r} of the VECM.
#' @param beta Matrix. A \eqn{((K+L+n_{d1}) \times r)} cointegrating matrix to be imposed -- 
#'   or estimated by the reduced-rank regression if \code{NULL} (the default).
#' @param type Character. The conventional case of the 
#'   \link[=as.t_D]{deterministic term} in the Johansen procedure.  
#' @param t_D1 List of vectors. The activating break periods \eqn{\tau} 
#'   for the period-specific \link[=as.t_D]{deterministic regressors} in \eqn{d_{1,t}}, 
#'   which are restricted to the cointegration relations.
#'   Just as in \code{\link{coint}}, the 
#'   accompanying lagged regressors are automatically included in \eqn{d_{2,t}}.
#' @param t_D2 List of vectors. The activating break periods \eqn{\tau} 
#'   for the period-specific \link[=as.t_D]{deterministic regressors} in \eqn{d_{2,t}}, 
#'   which are unrestricted.
#' @param D1 Matrix. A \eqn{(n_{\bullet} \times (p+T))} data matrix of customized 
#'   \link[=as.t_D]{deterministic regressors} added to \eqn{d_{1,t}}, 
#'   which are restricted to the cointegration relations. 
#'   Unlike '\code{t_D1}', these customized regressors require potential
#'   accompanying regressors to be manually included in \eqn{d_{2,t}}.
#' @param D2 Matrix. A \eqn{(n_{\bullet} \times (p+T))} data matrix of customized 
#'   \link[=as.t_D]{deterministic regressors} added to \eqn{d_{2,t}}, 
#'   which are unrestricted. 
#'   These additional regressors correspond to '\code{dumvar}' in \strong{urca}'s 
#'   \code{\link[urca]{ca.jo}}, which is fixed over bootstrap iterations.
#' 
#' @return A list of class '\code{\link[=as.varx]{varx}}'.
#' 
#' @references Johansen, S. (1996): 
#'   \emph{Likelihood-Based Inference in Cointegrated Vector Autoregressive Models}, 
#'   Advanced Texts in Econometrics, Oxford University Press, USA.
#' @references Luetkepohl, H. (2005): 
#'   \emph{New Introduction to Multiple Time Series Analysis}, 
#'   Springer, 2nd ed.
#' 
#' @examples
#' ### extend basic example in "vars" ###
#' library(vars)
#' data(Canada)
#' names_k = c("e", "U", "rw")  # names of endogenous variables
#' names_l = c("prod")  # names of exogenous variables
#' names_s = NULL  # optional shock names
#' x = Canada[ , names_l, drop=FALSE]
#' y = Canada[ , names_k, drop=FALSE]
#' 
#' # colnames of the restriction matrices are passed as shock names #
#' SR = matrix(NA, nrow=4, ncol=4, dimnames=list(NULL, names_s))
#' SR[4, 2] = 0
#' LR = matrix(NA, nrow=4, ncol=4, dimnames=list(NULL, names_s))
#' LR[1, 2:4] = 0
#' LR[2:4, 4] = 0
#' 
#' # estimate, identify, and plot the IRF #
#' R.vecm = VECM(y=y, dim_p=3, x=x, dim_q=3, dim_r=1, type="Case4")
#' R.grt  = id.grt(R.vecm, LR=LR, SR=SR)
#' R.irf  = irf(R.grt, n.ahead=50)
#' plot(R.irf)
#' 
#' @export
#' 
VECM <- function(y, dim_p, x=NULL, dim_q=dim_p, dim_r=NULL, beta=NULL, 
                 type=c("Case1", "Case2", "Case3", "Case4", "Case5"), 
                 t_D1=list(), t_D2=list(), D1=NULL, D2=NULL){
  # define and check
  t_D1  = as.t_D(t_D1)
  t_D2  = as.t_D(t_D2)
  def   = aux_stackRRR(y=y, dim_p=dim_p, x=x, dim_q=dim_q, type=type, t_D1=t_D1, t_D2=t_D2, D1=D1, D2=D2)
  RRR   = aux_RRR(def$Z0, def$Z1, def$Z2)
  dim_T = def$dim_T  # number of observations without presample
  dim_K = def$dim_K  # number of endogenous variables
  dim_L = def$dim_L  # number of weakly exogenous variables
  
  # set cointegrating vectors
  if(is.null(beta)){
    # ... estimated
    beta = aux_beta(V=RRR$V, dim_r=dim_r, normalize="natural")
  }else{
    # ... imposed, see Kilian,Luetkepohl 2017:88
    beta = as.matrix(beta)
    if(nrow(beta) != nrow(RRR$V)){ stop("Matrix 'beta' must have ", nrow(RRR$V), " rows!") }
    dim_r = ncol(beta)  # cointegration rank
    names_r = if( !is.null(colnames(beta)) ){ colnames(beta) }else{ paste0("ect.", 1:dim_r) }  # names of the long-run relations
    dimnames(beta) = list(rownames(RRR$V), names_r)
  }
  
  # estimate short-run dynamics
  vecm = aux_VECM(beta=beta, RRR=RRR)
  if(is.null(x)){
    # ... as full-system VECM
    PARTIAL  = NULL
    MARGINAL = NULL
    dim_q    = NULL
    rvar     = aux_vec2var(PI=vecm$PI, GAMMA=vecm$GAMMA, type_VECM=type, dim_p=dim_p)  # rank-restricted VAR model in levels
    names_k  = rownames(rvar$A)  # names of the endogenous variables
    B = diag(dim_K); dimnames(B) = list(names_k, paste0("u[ ", names_k, " ]"))
  }else{
    # ... as conditional VECM
    PARTIAL  = list(VECM=vecm, dim_K=dim_K)
    MARGINAL = VECM(y=x, dim_p=dim_q, dim_r=0, type=type, t_D1=t_D1, t_D2=t_D2, D1=D1, D2=D2)  # estimate via recursive call
    dim_T = min(def$dim_T, MARGINAL$dim_T)
    dim_K = dim_K + dim_L
    vecm  = aux_con2vec(beta=beta, VECM_c=PARTIAL$VECM, dim_p=dim_p, VECM_x=MARGINAL$VECM, dim_q=dim_q)
    rvar  = aux_vec2var(PI=vecm$PI, GAMMA=vecm$GAMMA, dim_p=max(dim_p, dim_q))  # rank-restricted VAR model in levels
    B     = vecm$B
  }
  
  # return result
  argues = list(method="RRR", type=type, t_D1=t_D1, t_D2=t_D2)
  result = list(A=rvar$A, B=B, y=def$y, x=def$x, D1=def$D1, D2=def$D2, 
                RRR=RRR, beta=beta, VECM=vecm, resid=vecm$resid, OMEGA=vecm$OMEGA, SIGMA=vecm$SIGMA, 
                dim_r=dim_r, dim_K=dim_K, dim_L=dim_L, dim_T=dim_T, dim_p=dim_p, dim_q=dim_q,
                PARTIAL=PARTIAL, MARGINAL=MARGINAL, args_varx=argues)
  class(result) = "varx"
  return(result)
}


#### S3 methods for objects of class 'varx' ####
#' @export
print.varx <- function(x, ...){
  # define
  x = as.varx(x)
  
  # create headers
  header_0 = "### VAR model ### \n"
  header_A = "\nEstimated coefficients: \n"
  header_B = "\nEstimated B matrix (unique decomposition of the covariance matrix): \n"
  
  # print
  cat(header_0)
  cat(header_A)
  print(x$A)
  
  if(identical(unname(x$B), diag(x$dim_K))){
    cat("\nModel: Reduced-form VAR \n")
  }else{
    cat(header_B)
    print(x$B)
  }
}


#' @export
summary.varx <- function(object, ...){
  # define
  object = as.varx(object)
  
  # print
  print.varx(object, ...)
  cat("\n", paste0(c("Specifications", rep(".", times=20))), sep="", "\n")
  cat("Deterministic term: ", object$args_varx$type, "\n")
  if(!is.null(object$args_id)){
    cat("Method: ", object$args_id$method, "\n")
  if(!is.null(object$beta)){
    cat("Cointegrating vectors: \n")
    print(object$beta) }
  }
}



#### Generic functions and their S3 methods ####
# CRAN manual: https://cran.r-project.org/doc/manuals/R-exts.html#Generic-functions-and-methods
# Roxygen: https://r-pkgs.org/man.html#man-s3

#' @title Coerce into a '\code{varx}' object
#' @description Coerce into a '\code{varx}' object. On top of the parent class 
#'   '\code{varx}', the child class '\code{id}' is imposed if the input object 
#'   to be transformed contains an SVAR model. 
#' @details \code{\link{as.varx}} is used as an intermediary in the \strong{pvars} 
#'   functions to achieve compatibility with different classes of VAR objects.
#'   If the user wishes to extend this compatibility with further classes, she 
#'   may simply specify accordant \code{\link{as.varx}}-methods instead of altering 
#'   the original \strong{pvars} function. Classes already covered by \strong{pvars} 
#'   are those of the \strong{vars} ecosystem, in particular the classes
#' \itemize{
#' \item '\code{varest}'  for reduced-form VAR estimates from \code{\link[vars]{VAR}},
#' \item '\code{vec2var}' for reduced-form VECM estimates from \code{\link[vars]{vec2var}},
#' \item '\code{svarest}' for structural VAR estimates from \code{\link[vars]{BQ}},
#' \item '\code{svecest}' for structural VECM estimates from \code{\link[vars]{SVEC}}, and
#' \item '\code{svars}'   for structural VAR estimates from \strong{svars}' 
#'   \code{\link[svars]{id.chol}}, \code{\link[svars]{id.cvm}}, or \code{\link[svars]{id.dc}}.
#' }
#'   By transformation to '\code{varx}', these VAR estimates can thus be subjected 
#'   to \strong{pvars}' bootstrap procedure \code{\link{sboot.mb}} and S3 methods 
#'   such as \code{\link[base]{summary}} and \code{\link[utils]{toLatex}}.
#' 
#' @param x A VAR object to be transformed.
#' @param ... Additional arguments to be passed to or from methods.
#' 
#' @return  A list of class '\code{varx}'. Objects of this class contain the elements:
#' \item{A}{Matrix. The lined-up VAR coefficient matrices \eqn{A_j, j=1,\ldots,p} for the 
#'   lagged variables.}
#' \item{B}{Matrix. The \eqn{(K \times S)} structural impact matrix of the SVAR model 
#'   or an identity matrix \eqn{I_K} as a placeholder for the unidentified VAR model.}
#' \item{beta}{Matrix. The \eqn{((K+n_{d1}) \times r)} cointegrating matrix of the VAR model 
#'   if transformed from a rank-restricted VECM.}
#' \item{SIGMA}{Matrix. The \eqn{(K \times K)} residual covariance matrix estimated by least-squares.}
#'   The following integers indicate the size of dimensions:
#' \item{dim_K}{Integer. The number of endogenous variables \eqn{K} in the full system.}
#' \item{dim_S}{Integer. The number of identified shocks \eqn{S} in the SVAR model.}
#' \item{dim_T}{Integer. The number of time periods \eqn{T} without presample.}
#' \item{dim_p}{Integer. The lag-order \eqn{p} of the VAR model in levels.}
#' \item{dim_r}{Integer. The cointegration rank \eqn{r} of the VAR model
#'   if transformed from a rank-restricted VECM.}
#'   Some further elements required for the bootstrap functions are: 
#' \item{y}{Matrix. The \eqn{(K \times (p+T))} endogenous variables.}
#' \item{D,D1,D2}{Matrices. The \eqn{(n_{\bullet} \times (p+T))}
#'   deterministic variables, fixed over bootstrap iterations, 
#'   (un)restricted to the cointegration relations of the VAR model 
#'   if transformed from a rank-restricted VECM.}
#' \item{resid}{Matrix. The \eqn{(K \times T)} residual matrix.}
#' \item{args_varx}{List of characters and integers indicating the estimator 
#'   and specifications that have been used.} 
#' \item{args_id}{List of characters and integers indicating the identification 
#'   methods and specifications that have been used. This element is specific 
#'   to the child class '\code{id}' for SVAR models, which inherits from
#'   the parent class '\code{varx}' for any VAR model.}
#' 
#' @examples
#' data("PCIT")
#' names_k = c("APITR", "ACITR", "PITB", "CITB", "GOV", "RGDP", "DEBT")
#' 
#' # estimate reduced-form VAR and coerce into 'varx' object #
#' R.vars = vars::VAR(PCIT[ , names_k], p=4, type="const")
#' as.varx(R.vars)
#' 
#' # identify structural VAR and coerce into 'id' object #
#' R.svar = svars::id.chol(R.vars, order_k=names_k)
#' as.varx(R.svar)
#' 
#' @export
#' 
as.varx <- function(x, ...) UseMethod("as.varx")


#' @method as.varx default
#' @export
as.varx.default <- function(x, ...){
  if(inherits(x, "varx")){
    return(x)
  }else{
    stop("Argument is not a VAR object of suitable class!")
  }
}


#' @method as.varx sboot2
#' @export
as.varx.sboot2 <- function(x, ...){
    return(x$varx)
}


#' @method as.varx varest
#' @importFrom vars Bcoef
#' @export
as.varx.varest <- function(x, ...){
  # define
  y = t(x$y)  # Kx(T+p) matrix of endogenous variables
  dim_K = x$K     # number of endogenous variables
  dim_T = x$obs   # number of observations without presample
  dim_p = x$p     # number of lags
  type  = x$type  # type of the VAR model, e.g. 'const'
  
  # matrix of coefficients (K x (n + K*p)) and of deterministic variables (n x T+p)
  A_raw   = vars::Bcoef(x)
  D_raw   = t(x$datamat)[-(1:dim_K), , drop=FALSE]  # remove regressands
  idx_slp = 1:(dim_K*dim_p)
  idx_dtr = ((dim_K*dim_p):ncol(A_raw))[-1]
  dim_n   = length(idx_dtr)  # number of deterministic regressors
  D_pre   = matrix(NA, nrow=dim_n, ncol=dim_p)  # presample of deterministic variables
  
  A = A_raw[ , c(idx_dtr, idx_slp), drop=FALSE]   # estimated VAR parameters
  D = cbind(D_pre, D_raw[idx_dtr, , drop=FALSE])  # deterministic variables
  ### Note that the ordering 'idx_dtr' for 'A_raw' and 'D_raw' can differ from aux_dummy()-outputs! 
  ### See 'season' and 'exogen' in 'D', which is then fixed in pvars' bootstrap functions. 
  ### 'exogen' is treated as deterministic in accordance with vars:::.boot(). 
  ### https://github.com/bpfaff/vars/blob/master/R/internal.R#L139
  
  # residuals and their covariance matrix
  resid = t(resid(x))  # residual matrix (K x T)
  rownames(resid) = rownames(y)
  OMEGA = tcrossprod(resid) / dim_T        # MLE covariance matrix of residuals
  SIGMA = OMEGA * (dim_T/(dim_T-ncol(A)))  # OLS covariance matrix of residuals
  
  # return result
  argues = list(method="qr", type=type)
  result = list(A=A, B=diag(dim_K), y=y, D=D, resid=resid, SIGMA=SIGMA, OMEGA=OMEGA, 
                dim_K=dim_K, dim_T=dim_T, dim_p=dim_p, args_varx=argues)
  class(result) = "varx"
  return(result)
}


#' @method as.varx vec2var
#' @export
as.varx.vec2var <- function(x, ...){
  # define
  y = t(x$y)  # Kx(T+p) matrix of endogenous variables
  dim_r = x$r    # cointegration rank
  dim_K = x$K    # number of endogenous variables
  dim_T = x$obs  # number of observations without presample
  dim_p = x$p    # number of lags
  
  # re-estimate (easier and less error-prone than re-arranging the results)
  type = switch(x$vecm@ecdet, "const"="Case2", "none"="Case3", "trend"="Case4")  # the conventional case
  D2   = if( is.null(x$vecm@dumvar) ){ NULL }else{ t(x$vecm@dumvar) }
  def  = aux_stackRRR(y=y, dim_p=dim_p, type=type, D2=D2, t_D2=list(n.season=x$vecm@season))
  RRR  = aux_RRR(def$Z0, def$Z1, def$Z2)
  beta = aux_beta(V=x$vecm@V, dim_r=dim_r, normalize="natural")
  vecm = aux_VECM(beta=beta, RRR=RRR)
  rvar = aux_vec2var(PI=vecm$PI, GAMMA=vecm$GAMMA, dim_p=dim_p, type_VECM=type)
  
  # return result
  if(x$vecm@spec == "longrun"){ warning("The VECM has been re-specified into 'transitory' form.") }
  argues = list(method="RRR", type=type)
  result = list(A=rvar$A, B=diag(dim_K), y=y, D1=def$D1, D2=def$D2, RRR=RRR, beta=beta, VECM=vecm,
                resid=vecm$resid, SIGMA=vecm$SIGMA, OMEGA=vecm$OMEGA, 
                dim_r=dim_r, dim_K=dim_K, dim_T=dim_T, dim_p=dim_p, args_varx=argues)
  class(result) = "varx"
  return(result)
}


#' @method as.varx svarest
#' @export
as.varx.svarest <- function(x, ...){
  # coerce the reduced-form VAR
  result = as.varx(x$var)
  
  # append "id" results
  result$B = x$B  # structural impact matrix
  
  # append "id" arguments
  if(x$type == "Blanchard-Quah"){
    result$SR = x$B     # structural effects on impact
    result$LR = x$LRIM  # structural long-run effects
    args_id   = list(method=x$type, ...)
    
  ### TODO: }else if(any(x$type == c("A-model", "B-model", "AB-model"))){
    ### args_id = list(method=x$type, ...)
    
  }else{
    warning("svarest type '",  x$type,"' is not fully supported in pvars.")
    args_id = list(method=x$type, ...)
  }
  
  # return result
  result$args_id = args_id
  class(result)  = c("id", "varx")
  return(result)
}


#' @method as.varx svecest
#' @importFrom vars vec2var
#' @export
as.varx.svecest <- function(x, ...){
  # coerce the reduced-form VAR
  R.vars = vars::vec2var(x$var, r=x$r)
  result = as.varx(R.vars)
  
  # append "id" results
  result$B  = x$SR  # estimated B matrix (unique decomposition of the covariance matrix)
  result$SR = x$SR  # structural effects on impact
  result$LR = x$LR  # structural long-run effects
  
  # append "id" arguments
  ### According to the arguments for bootstrapping IRF in 'vars'. 
  ### https://github.com/bpfaff/vars/blob/master/R/internal.R#L380
  args_id = list(method="SVECM", SR=x$SRorig, LR=x$LRorig)
  args_id$start     = if(is.null(x$call$start)){      NULL }else{ x$call$start }
  args_id$max.iter  = if(is.null(x$call$max.iter)){    100 }else{ x$call$max.iter }
  args_id$conv.crit = if(is.null(x$call$conv.crit)){ 1e-07 }else{ x$call$conv.crit }
  args_id$maxls     = if(is.null(x$call$maxls)){       1.0 }else{ x$call$maxls }
  
  # return result
  result$args_id = args_id
  class(result)  = c("id", "varx")
  return(result)
}


#' @method as.varx svars
#' @importFrom copula indepTestSim
#' @export
as.varx.svars <- function(x, dd=NULL, itermax=300, steptol=200, iter2=50, ...){
  # coerce the reduced-form VAR
  result = as.varx(x$VAR)
  ### Note that svars (version 1.3.8) does not consider seasonal dummies or
  ### exogenous variables after 'YLagCr()' or 'get_A_hat()'. 
  ### https://github.com/alexanderlange53/svars/blob/master/R/get_A_hat.R
  
  # append "id" results
  result$B = x$B  # structural impact matrix
  
  # append "id" arguments
  if(x$method == "Cholesky"){
    args_id = list(method=x$method, order_k=x$order_k, ...)
    
  }else if(x$method == "Distance covariances"){
    args_id = list(method=x$method, PIT=x$PIT, n.iterations=100, ...)
    
  }else if(x$method == "Cramer-von Mises distance"){
    if(is.null(dd)){ dd = copula::indepTestSim(n=result$dim_T, p=result$dim_K, verbose=F) }
    args_id = list(method=x$method, dd=dd, itermax=itermax, steptol=steptol, iter2=iter2, ...)
    ### Same defaults as in https://github.com/alexanderlange53/svars/blob/master/R/bootstrap_mb.R#L87
    
  }else{
    warning("svars method '",  x$method,"' is not fully supported in pvars.")
    args_id = list(method=x$method, ...)
  }
  
  # return result
  result$args_id = args_id
  class(result)  = c("id", "varx")
  return(result)
}


