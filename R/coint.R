

#' @title Test procedures for the cointegration rank
#' @description Performs test procedures for the rank of cointegration in a single VAR model.
#'   The \eqn{p}-values are approximated by gamma distributions, 
#'   whose moments are automatically adjusted to 
#'   potential period-specific deterministic regressors 
#'   and weakly exogenous regressors in the partial VECM.
#' @param y Matrix. A \eqn{(K \times (p+T))} data matrix of the \eqn{K} endogenous time series variables.
#' @param dim_p Integer. Lag-order \eqn{p} for the endogenous variables \code{y}.
#' @param x Matrix. A \eqn{(L \times (q+T))} data matrix of the \eqn{L} weakly exogenous time series variables.
#' @param dim_q Integer. Lag-order \eqn{q} for the weakly exogenous variables \code{x}. 
#'   The literature uses \code{dim_p} (the default).
#' @param t_D1 List of vectors. The activating break periods \eqn{\tau} 
#'   for the period-specific \link[=as.t_D]{deterministic regressors} in \eqn{d_{1,t}}, 
#'   which are restricted to the cointegration relations.
#'   The accompanying lagged regressors are automatically included in \eqn{d_{2,t}}. 
#'   The \eqn{p}-values are calculated for up to two breaks resp. three sub-samples.
#' @param t_D2 List of vectors. The activating break periods \eqn{\tau} 
#'   for the period-specific \link[=as.t_D]{deterministic regressors} in \eqn{d_{2,t}}, 
#'   which are unrestricted.
#' 
#' @return A list of class '\code{coint}', 
#'   which contains elements of length \eqn{K} for each \eqn{r_{H0}=0,\ldots,K-1}: 
#' \item{r_H0}{Rank under each null hypothesis.}
#' \item{stats_TR}{Trace (TR) test statistics.}
#' \item{stats_ME}{Maximum eigenvalue (ME) test statistics.}
#' \item{pvals_TR}{\eqn{p}-values of the TR test.}
#' \item{pvals_ME}{\eqn{p}-values of the ME test. 
#'   \code{NA} if moments of the gamma distribution are not available 
#'   for the chosen data generating process.}
#' \item{lambda}{Eigenvalues, the squared canonical correlation coefficients 
#'   (saved only for the Johansen procedure).}
#' \item{args_coint}{List of characters and integers 
#'   indicating the cointegration test and specifications that have been used.}
#' 
#' @family cointegration rank tests
#' @name coint
NULL


#' @describeIn coint Johansen procedure.
#' @param type Character. The conventional case of the 
#'   \link[=as.t_D]{deterministic term} in the Johansen procedure.
#' @references Johansen, S. (1988): 
#'   "Statistical Analysis of Cointegration Vectors", 
#'   \emph{Journal of Economic Dynamics and Control}, 12, pp. 231-254.
#' @references Doornik, J. (1998): 
#'   "Approximations to the Asymptotic Distributions of Cointegration Tests", 
#'   \emph{Journal of Economic Surveys}, 12, pp. 573-593.
#' @references Johansen, S., Mosconi, R., and Nielsen, B. (2000): 
#'   "Cointegration Analysis in the Presence of Structural Breaks in the Deterministic Trend",
#'   \emph{Econometrics Journal}, 3, pp. 216-249.
#' @references Kurita, T., Nielsen, B. (2019):
#'   "Partial Cointegrated Vector Autoregressive Models with Structural Breaks in Deterministic Terms",
#'   \emph{Econometrics}, 7, pp. 1-35.
#' @examples
#' ### reproduce basic example in "urca" ###
#' library("urca")
#' data(denmark)
#' sjd = denmark[ , c("LRM", "LRY", "IBO", "IDE")]
#' 
#' # rank test and estimation of the full VECM as in "urca" #
#' R.JOrank = coint.JO(y=sjd,      dim_p=2, type="Case2", t_D2=list(n.season=4))
#' R.JOvecm = VECM(y=sjd, dim_r=1, dim_p=2, type="Case2", t_D2=list(n.season=4))
#' 
#' # ... and of the partial VECM, i.e. after imposing weak exogeneity #
#' R.KNrank = coint.JO(y=sjd[ , c("LRM"), drop=FALSE],   dim_p=2,
#'                     x=sjd[ , c("LRY", "IBO", "IDE")], dim_q=2, 
#'                     type="Case2", t_D1=list(t_shift=36), t_D2=list(n.season=4))
#' R.KNvecm = VECM(y=sjd[ , c("LRM"), drop=FALSE],   dim_p=2,
#'                 x=sjd[ , c("LRY", "IBO", "IDE")], dim_q=2, dim_r=1, 
#'                 type="Case2", t_D1=list(t_shift=36), t_D2=list(n.season=4))
#' 
#' @export
#' 
coint.JO <- function(y, dim_p, x=NULL, dim_q=dim_p, 
                     type=c("Case1", "Case2", "Case3", "Case4", "Case5"), 
                     t_D1=NULL, t_D2=NULL){
  # define and check
  t_D1  = as.t_D(t_D1)
  t_D2  = as.t_D(t_D2)
  R.def = aux_stackRRR(y=y, dim_p=dim_p, x=x, dim_q=dim_q, type=type, t_D1=t_D1, t_D2=t_D2)
  dim_K = R.def$dim_K  # number of endogenous variables
  dim_L = R.def$dim_L  # number of weakly exogenous variables
  dim_T = R.def$dim_T  # number of observations without presample
  
  is_break = !is.null(c(t_D1$t_break, t_D1$t_shift))
  type_mom = if(is_break){ if(dim_L != 0){ paste0("KN_", type) }else{ paste0("JMN_", type) }}else{ type }
  
  # LR-tests, from Johansen 1996:92
  moments = aux_CointMoments(dim_K=dim_K, dim_L=dim_L, dim_T=dim_T+dim_p, t_D1=t_D1, r_H0=0:(dim_K-1), type=type_mom)
  RRR     = aux_RRR(Z0=R.def$Z0, Z1=R.def$Z1, Z2=R.def$Z2)
  result  = aux_LRrank(lambda=RRR$lambda, dim_T=dim_T, dim_K=dim_K, r_H0=0:(dim_K-1), moments=moments)
  
  # return result
  result$lambda = RRR$lambda[1:dim_K]
  class(result) = "coint"
  result$args_coint = list(method="Johansen Procedure", y=y, x=x, dim_p=dim_p, dim_q=dim_q, type=type_mom, t_D1=t_D1, t_D2=t_D2)
  return(result)
}


#' @describeIn coint (Trenkler)-Saikkonen-Luetkepohl procedure.
#' @param type_SL Character. The conventional case of the 
#'   \link[=as.t_D]{deterministic term} in the Saikkonen-Luetkepohl (SL) procedure. 
#' @param t_D List of vectors. The activation periods \eqn{\tau} 
#'   for the period-specific \link[=as.t_D]{deterministic regressors} in \eqn{d_{t}} of the SL-procedure. 
#'   The accompanying lagged regressors are automatically included in \eqn{d_{t}}. 
#'   The \eqn{p}-values are calculated for up to two breaks resp. three sub-samples. 
#' @references Saikkonen, P., and Luetkepohl, H. (2000):
#'   "Trend Adjustment Prior to Testing for the Cointegrating Rank of a Vector Autoregressive Process",
#'   \emph{Journal of Time Series Analysis}, 21, pp. 435-456.
#' @references Trenkler, C. (2008): 
#'   "Determining \eqn{p}-Values for Systems Cointegration Tests with a Prior Adjustment for Deterministic Terms", 
#'   \emph{Computational Statistics}, 23, pp. 19-39.
#' @references Trenkler, C., Saikkonen, P., and Luetkepohl, H. (2008): 
#'   "Testing for the Cointegrating Rank of a VAR Process with Level Shift and Trend Break", 
#'   \emph{Journal of Time Series Analysis}, 29, pp. 331-358.
#' @examples
#' ### reproduce Oersal,Arsova 2016:22, Tab.7.5 "France" ###
#' data("ERPT")
#' names_k = c("lpm5", "lfp5", "llcusd")  # variable names for "Chemicals and related products"
#' names_i = levels(ERPT$id_i)[c(1,6,2,5,4,3,7)]  # ordered country names
#' L.data  = sapply(names_i, FUN=function(i) 
#'    ts(ERPT[ERPT$id_i==i, names_k], start=c(1995, 1), frequency=12), 
#'    simplify=FALSE)
#' R.TSLrank = coint.SL(y=L.data$France, dim_p=3, type_SL="SL_trend", t_D=list(t_break=89))
#' 
#' @export
#' 
coint.SL <- function(y, dim_p, type_SL=c("SL_mean", "SL_trend"), t_D=NULL){
  # define and check
  t_D = as.t_D(t_D)
  y = aux_asDataMatrix(y, "y")  # named matrix y is supposed to be KxT regardless of input
  
  R.dtr = aux_stackDtr(type_SL=type_SL, t_D=t_D, dim_p=dim_p, dim_T=ncol(y))
  R.def = aux_stackRRR(y=y, dim_p=dim_p, D1=R.dtr$D1, D2=R.dtr$D2)
  dim_K = R.def$dim_K  # number of endogenous variables
  dim_T = R.def$dim_T  # number of observations without presample
  
  # detrend and test under each r_H0  
  type_mom = if(is.null(t_D$t_break)){ type_SL }else{ "TSL_trend" }
  moments  = aux_CointMoments(dim_K=dim_K, dim_T=dim_T+dim_p, t_D1=t_D, r_H0=0:(dim_K-1), type=type_mom)
  RRR      = aux_RRR(Z0=R.def$Z0, Z1=R.def$Z1, Z2=R.def$Z2)
  result   = NULL
  for(r_H0 in 0:(dim_K-1)){
    # cointegration test with prior adjustment for deterministic term, from Trenkler 2008:22
    vecm = aux_VECM(beta=RRR$V[ , 0:r_H0, drop=FALSE], RRR=RRR)
    A_H0 = aux_vec2var(PI=vecm$PI, GAMMA=vecm$GAMMA, dim_p=dim_p)$A
    x_H0 = aux_GLStrend(y=y, OMEGA=vecm$OMEGA, A=A_H0, D=R.dtr$D, dim_p=dim_p)$x
    Z_H0 = aux_stackRRR(y=x_H0, dim_p=dim_p)  # no deterministic term in detrended data
    l_H0 = aux_RRR(Z0=Z_H0$Z0, Z1=Z_H0$Z1, Z2=Z_H0$Z2)$lambda
    rank = aux_LRrank(lambda=l_H0, dim_T=dim_T, dim_K=dim_K, r_H0=r_H0, moments=moments[r_H0+1, , drop=FALSE])  # moments order r_H0=0,...,K-1
    
    result = rbind(result, rank)
    rm(vecm, A_H0, x_H0, Z_H0, l_H0, rank)
  }
  
  # return result  
  result$r_H0 = 0:(dim_K-1)
  class(result) = "coint"
  result$args_coint = list(method="SL-Procedure", y=y, dim_p=dim_p, type=type_mom, t_D=t_D)
  return(result)
}


#### S3 methods for objects of class 'coint' ####
#' @export
print.coint <- function(x, ...){
  # create table
  header_args = c("### Cointegration Rank by ", x$args_coint$method, " with ", x$args_coint$type, " ###")
  x$args_coint = NULL
  class(x) = NULL
  df_table = round(as.data.frame(x), 3)
    
  # print
  cat(header_args, "\n", sep="")
  print(df_table, quote=FALSE, row.names=FALSE)
}


#' @export
summary.coint <- function(object, ...){
  cat_q = if(is.null(object$args_coint$x)){ NULL }else{ paste0(" and q = ", object$args_coint$dim_q) }
  
  # print
  print.coint(object, ...)
  cat("\n", paste0(c("Specifications", rep(".", times=20))), sep="", "\n")
  cat("Method: ", object$args_coint$method, "\n")
  cat("Lag-order:  p =", object$args_coint$dim_p, cat_q, "\n")
  cat("Moments of Z_d: ", object$args_coint$type, "\n")
}


#' @export
toLatex.coint <- function(object, ..., digits=3, write_ME=FALSE, add2header=NULL){
  # check and define
  L.coint = list(object, ...)
  L.hypot = sapply(L.coint, function(i) identical(object$r_H0, i$r_H0))
  if(!all(L.hypot)){ stop("The hypotheses must be the same for all test procedures!") }
  
  n.tables = length(L.coint)
  n.header = length(add2header)
  add2header = c(add2header, rep("Trace test", n.tables-n.header))
  
  # initialize
  matrix_table = cbind(object$r_H0)
  header_tests = c("$ r_{H0} $")
  header_args  = rep(" ", times=ncol(matrix_table))
  header_cline = NULL
  
  for(j in 1:n.tables){
    # define
    xj = L.coint[[j]]
    xj$args_coint = NULL
    class(xj) = NULL
    df_table = round(as.matrix(as.data.frame(xj)), digits=digits)
    idx_col = ncol(matrix_table) + 1:4 + c(0, 0, 1, 1)  # second vector adds separating column
    
    # create tabular
    header_args  = c(header_args, paste0("\\multicolumn{2}{c}{ ", add2header[j], " }"), " ", ### add header without "&"!
                     if(write_ME){ c("\\multicolumn{2}{c}{ Max. eigenvalue test }", " ")})
    
    header_cline = paste0(header_cline, "\\cline{", idx_col[1], "-", idx_col[2] ,"} ", 
                          if(write_ME){ paste0("\\cline{", idx_col[3], "-", idx_col[4] ,"}")})
    
    header_tests = c(header_tests, c("statistic", "$ p $-value", " "),
                     if(write_ME){ c("statistic", "$ p $-value", " ")})
    
    matrix_table = cbind(matrix_table, format(df_table[ , c(2, 4), drop=FALSE], nsmall=digits), " ",
                         if(write_ME){ cbind(format(df_table[ , c(3, 5), drop=FALSE], nsmall=digits), " ")})
    rm(xj, df_table, idx_col)
  }
  
  # return result
  tabular_cols = paste0(rep("r", each=ncol(matrix_table)), collapse = "")
  result = c(
    paste0("\\begin{tabular}{", tabular_cols, "}", sep=""),
    "\t\\hline \\hline",
    
    paste0("\t", paste(header_args,  collapse=" & "), " \\\\"),
    paste0("\t", header_cline),
    paste0("\t", paste(header_tests, collapse=" & "), " \\\\"),
    "\t\\hline",
    apply(matrix_table, MARGIN=1, FUN=function(x) paste0("\t", paste(x, collapse=" & "), " \\\\")),
    
    "\t\\hline \\hline",
    "\\end{tabular}"
  )
  class(result) = "Latex"
  return(result)
}


