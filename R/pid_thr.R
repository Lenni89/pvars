

#' @title Recursive identification of panel SVAR models via Cholesky decomposition
#' @description Given an estimated panel VAR model, this function uses the Cholesky decomposition to identify  
#'   the structural impact matrix \eqn{B_i} of the corresponding SVAR model
#'   \deqn{y_{it} = c_{it} + A_{i1} y_{i,t-1} + ... + A_{i,p_i} y_{i,t-p_i} + u_{it}}
#'   \deqn{       = c_{it} + A_{i1} y_{i,t-1} + ... + A_{i,p_i} y_{i,t-p_i} + B_i \epsilon_{it}.}
#'   Matrix \eqn{B_i} corresponds to the decomposition of the least squares covariance matrix \eqn{\Sigma_{u,i} = B_i B_i'}.
#'
#' @param x An object of class '\code{pvarx}' or a list of VAR objects 
#'   that will be \link[=as.varx]{coerced} to '\code{varx}'. 
#'   Estimated panel of VAR objects.
#' @param order_k Vector. Vector of characters or integers specifying the assumed structure of the recursive causality. 
#'   Change the causal ordering in the instantaneous effects without permuting variables and re-estimating the VAR model.
#' 
#' @return List of class '\code{pid}' with elements:
#' \item{A}{Matrix. The lined-up coefficient matrices \eqn{A_j, j=1,\ldots,p} 
#'    for the lagged variables in the panel VAR.}
#' \item{B}{Matrix. Mean group of the estimated structural impact matrices \eqn{B_i}, 
#'    i.e. the unique decomposition of the covariance matrices of reduced-form errors.}
#' \item{L.varx}{List of '\code{varx}' objects for the individual estimation results
#'   to which the structural impact matrices \eqn{B_i} have been added.}
#' \item{args_pid}{List of characters and integers indicating the identification methods and specifications that have been used.}
#' \item{args_pvarx}{List of characters and integers indicating the estimator and specifications that have been used.}
#' 
#' @references Luetkepohl, H. (2005): 
#'   \emph{New Introduction to Multiple Time Series Analysis}, 
#'   Springer, 2nd ed.
#' @references Sims, C. A. (2008): 
#'   "Macroeconomics and Reality", 
#'   \emph{Econometrica}, 48, pp. 1-48.
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
#' @family panel identification functions
#' @export
#' 
pid.chol <- function(x, order_k=NULL){ 
  # define
  x = as.pvarx(x)
  dim_K = nrow(x$A)
  
  # names for variables and for shocks
  names_k = if( !is.null(rownames(x$A)) ){ rownames(x$A) }else{ paste0("y", 1:dim_K) }
  names_s = paste0("epsilon[ ", names_k, " ]")
  
  # common permutation matrix for recursive causality
  if(is.null(order_k)){
    order_k = 1:dim_K
  }else if(is.character(order_k)){
    order_k = sapply(order_k, FUN=function(x) which(names_k == x))
    warn_ok = "Check variable names given as characters in 'order_k' or use integers instead!"
    if(is.list(order_k)){ stop(warn_ok) }
  }
  names(order_k) = names_k
  perm  = diag(dim_K)[ , order_k, drop=FALSE]  # permutation matrix
  tperm = t(perm)  # note that solve(perm) == t(perm)
  
  # identify individual SVAR models
  idf <- function(x){
    # define
    SIGMA = x$SIGMA  # OLS covariance matrix of residuals
    
    # Choleski decomposition
    psig = tperm %*% SIGMA %*% perm  # permute in accordance with recursive causality
    P    = t(chol(psig))
    B    = perm %*% P %*% tperm  # repermute via inverse
    
    # add to "varx"-object
    dimnames(B) = list(names_k, names_s)
    x$B = B
    x$args_id = args_id
    class(x) = c("id", "varx")
    return(x)
  }
  
  # calculate structural impact matrices
  args_id  = list(method="Cholesky", order_k=order_k)
  x$L.varx = lapply(x$L.varx, function(x_i) idf(x_i))
  
  # estimate panel mean-group
  x$MG_B = aux_MG(x$L.varx, w=x$args_pvarx$w, idx_par="B")
  x$B    = x$MG_B$mean
  
  # return result
  x$args_pid = args_id
  class(x) = c("pid", "pvarx")
  return(x)
}


#' @title Identification of panel SVEC models by imposing long- and short-run restrictions
#' @description Identifies a panel SVEC model by utilizing a scoring algorithm 
#'   to impose long- and short-run restrictions. 
#'   See the details of \code{\link[vars]{SVEC}} in \strong{vars}.
#' 
#' @param x An object of class '\code{pvarx}' or a list of VECM objects 
#'   that will be \link[=as.varx]{coerced} to '\code{varx}'. 
#'   Panel of VAR objects estimated under rank-restriction.
#' @param LR Matrix. The restricted long-run impact matrix.
#' @param SR Matrix. The restricted contemporaneous impact matrix.
#' @param start Vector. The starting values for \eqn{\gamma}, 
#'   set by \code{\link[stats]{rnorm}} if \code{NULL} (the default).
#' @param max.iter Integer. The maximum number of iterations.
#' @param conv.crit Real number. Convergence value of algorithm.
#' @param maxls Real number. Maximum movement of the parameters between two iterations of the scoring algorithm.
#' 
#' @return List of class '\code{pid}' with elements:
#' \item{A}{Matrix. The lined-up coefficient matrices \eqn{A_j, j=1,\ldots,p} 
#'    or the lagged variables in the panel VAR.}
#' \item{B}{Matrix. Mean group of the estimated structural impact matrices \eqn{B_i}, 
#'    i.e. the unique decomposition of the covariance matrices of reduced-form errors.}
#' \item{L.varx}{List of '\code{varx}' objects for the individual estimation results
#'   to which the structural impact matrices \eqn{B_i} have been added.}
#' \item{args_pid}{List of characters and integers indicating the identification methods and specifications that have been used.}
#' \item{args_pvarx}{List of characters and integers indicating the estimator and specifications that have been used.}
#' 
#' @references Amisano, G. and Giannini, C. (1997): 
#'   \emph{Topics in Structural VAR Econometrics}, 
#'   Springer, 2nd ed.
#' @references Breitung, J., Brueggemann R., and Luetkepohl, H. (2004): 
#'   "Structural Vector Autoregressive Modeling and Impulse Responses", 
#'   in \emph{Applied Time Series Econometrics}, 
#'   ed. by H. Luetkepohl and M. Kraetzig, 
#'   Cambridge University Press, Cambridge.
#' @references Johansen, S. (1996): 
#'   \emph{Likelihood-Based Inference in Cointegrated Vector Autoregressive Models}, 
#'   Advanced Texts in Econometrics, Oxford University Press, USA.
#' @references Luetkepohl, H. (2005): 
#'   \emph{New Introduction to Multiple Time Series Analysis}, 
#'   Springer, 2nd ed.
#' @references Pfaff, B. (2008):
#'   "VAR, SVAR and SVEC Models: Implementation within R Package \strong{vars}",
#'   \emph{Journal of Statistical Software}, 27, pp. 1-32.
#' @seealso \ldots the original \code{\link[vars]{SVEC}} by Pfaff (2008) in \strong{vars}. 
#'   Note that \code{\link{pid.grt}}  relies on this underlying procedure, 
#'   but allows for the additional model specifications in \code{\link{pvarx.VEC}} 
#'   and for the bootstrap procedures in \code{\link{sboot.pmb}}, 
#'   both provided by the \strong{pvars} package.
#' 
#' @examples
#' data("PCAP")
#' names_k = c("g", "k", "l", "y")  # variable names
#' names_i = levels(PCAP$id_i)      # country names
#' names_s = NULL                   # optional shock names
#' L.data  = sapply(names_i, FUN=function(i) 
#'   ts(PCAP[PCAP$id_i==i, names_k], start=1960, end=2019, frequency=1), 
#'   simplify=FALSE)
#' 
#' # colnames of the restriction matrices are passed as shock names #
#' SR = matrix(NA, nrow=4, ncol=4, dimnames=list(names_k, names_s))
#' SR[1, 2] = 0
#' SR[3, 4] = 0
#' LR = matrix(NA, nrow=4, ncol=4, dimnames=list(names_k, names_s))
#' LR[ , 3:4] = 0
#' 
#' # estimate and identify panel SVECM #
#' R.pvec = pvarx.VEC(L.data, lags=2, dim_r=2, type="Case4")
#' R.pid  = pid.grt(R.pvec, LR=LR, SR=SR)
#' 
#' @family panel identification functions
#' @export
#' 
pid.grt <- function(x, LR=NULL, SR=NULL, start=NULL, max.iter=100, conv.crit=1e-07, maxls=1.0){ 
  # define
  x = as.pvarx(x)
  dim_K = nrow(x$A)
  
  # identify individual SVAR models
  idf <- function(x){
    # define
    dim_p = x$dim_p  # number of lags
    dim_T = x$dim_T  # number of observations without presample
    
    # calculate long-run multiplier matrix
    alpha_oc = aux_oc(x$VECM$alpha)
    beta_oc  = aux_oc(x$beta[1:dim_K, ])
    XI = aux_vec2vma(GAMMA=x$VECM$GAMMA, alpha_oc=alpha_oc, beta_oc=beta_oc, dim_p=dim_p, B=NULL, n.ahead="XI")
    
    # identify via FIML scoring algorithm
    GRT = aux_idGRT(OMEGA=x$VECM$OMEGA, XI=XI, dim_T=dim_T, LR=LR, SR=SR, 
                    start=start, max.iter=max.iter, conv.crit=conv.crit, maxls=maxls)
    
    # add to "varx"-object
    x$B  = GRT$SR  # estimated B matrix (unique decomposition of the covariance matrix)
    x$SR = GRT$SR  # structural effects on impact
    x$LR = GRT$LR  # structural long-run effects
    x$args_id = args_id
    class(x) = c("id", "varx")
    return(x)
  }
  
  # calculate structural impact matrices
  args_id  = list(method="SVECM", SR=SR, LR=LR, start=start, max.iter=max.iter, conv.crit=conv.crit, maxls=maxls)
  x$L.varx = lapply(x$L.varx, function(x_i) idf(x_i))
  
  # estimate panel mean-group
  x$MG_B = aux_MG(x$L.varx, w=x$args_pvarx$w, idx_par="B")
  x$B    = x$MG_B$mean
  
  # return result
  x$args_pid = args_id
  class(x) = c("pid", "pvarx")
  return(x)
}


#### S3 methods for objects of class 'pid' ####
#' @export
summary.pid <- function(object, ..., probs=seq(0, 1, 0.25), digits=3){
  # define
  object = as.pvarx(object)
  B = sapply(object$L.varx, FUN=function(x) x$B)
  dim_N = length(object$L.varx)
  dim_K = nrow(object$B)
  dim_S = ncol(object$B)
  
  idx_min = (probs == 0)
  idx_max = (probs == 1)
  names_q = paste0("q", sub(".", "", format(probs)))
  names_q[idx_min] = "min"
  names_q[idx_max] = "max"
  names_tmp = paste0("b_", 1:dim_K)
  names_b   = c(sapply(1:dim_S, FUN=function(s) paste0(names_tmp, s)))
   
  # calculate mean-group statistics, see Herwartz,Bernoth 2021:10, Tab.2
  R.quant = apply(B, MARGIN=1, FUN=function(x) quantile(x, probs=probs))
  R.mean  = apply(B, MARGIN=1, FUN=function(x) mean(x))  # sample mean, i.e. MG estimate
  R.var   = apply(B, MARGIN=1, FUN=function(x) var(x))   # sample variance
  R.ratio = R.mean / sqrt(R.var)  # t-ratios
  R.stats = cbind(t(R.quant), R.mean, R.ratio)
  
  # create table
  header_args = c("### Structural Identification by ", object$args_pid$method, " ###")
  header_tab  = c(names_q, "mean", "t-ratio")
  matrix_tab  = format(round(R.stats, digits=digits), nsmall=digits)
  dimnames(matrix_tab) = list(names_b, header_tab)
  
  # print
  cat(header_args, "\n", sep="")
  print(matrix_tab, quote=FALSE, ...)
  cat("\n", paste0(c("Specifications", rep(".", times=20))), sep="", "\n")
  cat("Number of individuals: ", dim_N, "\n")
  if(!is.null(object$rotation_angles)){
    cat("Rotation angles: ", object$rotation_angles, "\n") }
}


#' @export
toLatex.pid <- function(object, ..., probs=seq(0, 1, 0.25), digits=3){
  # define
  object = as.pvarx(object)
  B = sapply(object$L.varx, FUN=function(x) x$B)
  dim_N = length(object$L.varx)
  dim_K = nrow(object$B)
  dim_S = ncol(object$B)
  
  idx_min = (probs == 0)
  idx_max = (probs == 1)
  names_q = paste0("$ q", sub(".", "", format(probs)), " $")
  names_q[idx_min] = "min"
  names_q[idx_max] = "max"
  names_tmp = paste0("$ \\widehat{\\mathsf{b}}_{", 1:dim_K)
  names_b   = c(sapply(1:dim_S, FUN=function(s) paste0(names_tmp, s, "} $")))
  
  ### TODO: effect names
  #names_k = if( !is.null(rownames(A)) ){ rownames(A) }else{ paste0("y.", 1:nrow(A)) }  # variable names
  #names_s = if( !is.null(colnames(B)) ){ colnames(B) }else{ paste0("epsilon[ ", 1:ncol(B), " ]") }  # shock names
  #names_e = c(sapply(names_k, FUN=function(k) paste0(names_s, " %->% ", k)))  # names for each effect in B
  
  # calculate mean-group statistics, see Herwartz,Bernoth 2021:10, Tab.2
  R.quant = apply(B, MARGIN=1, FUN=function(x) quantile(x, probs=probs))
  R.mean  = apply(B, MARGIN=1, FUN=function(x) mean(x))  # sample mean, i.e. MG estimate
  R.var   = apply(B, MARGIN=1, FUN=function(x) var(x))   # sample variance
  R.ratio = R.mean / sqrt(R.var)  # t-ratios
  R.stats = cbind(t(R.quant), R.mean, R.ratio)
  idx_NaN = cbind(FALSE, is.nan(R.stats))
  
  # create tabular
  header_tab = c("\t ", names_q, "mean", "$t$-ratio")
  matrix_tab = cbind(names_b, format(round(R.stats, digits=digits), nsmall=digits))
  matrix_tab[idx_NaN] = " -- "
  
  # return result
  tabular_core = apply(matrix_tab, MARGIN=1, FUN=function(x) paste0("\t", paste(x, collapse=" & "), " \\\\"))
  tabular_cols = paste0(rep("r", each=ncol(matrix_tab)), collapse = "")
  tabular_rows = c(rep(c(rep("", each=dim_K-1), "[5pt] \n"), times=dim_S-1), rep("", each=dim_K))
  
  result = c(
    paste0("\\begin{tabular}{", tabular_cols, "}", sep=""),
    "\t\\hline \\hline",
    paste0("\t", paste(header_tab, collapse=" & "), " \\\\"),
    "\t\\hline",
    paste0(tabular_core, tabular_rows),
    "\t\\hline \\hline",
    "\\end{tabular}"
  )
  class(result) = "Latex"
  return(result)
}


