

#' @title Identification of SVAR models by means of proxy variables
#' @description Given an estimated VAR model, this function uses proxy variables
#'   to identify the structural impact matrix \eqn{B} 
#'   of the corresponding SVAR model
#'   \deqn{y_{t} = c_{t} + A_{1} y_{t-1} + ... + A_{p} y_{t-p} + u_{t}}
#'   \deqn{      = c_{t} + A_{1} y_{t-1} + ... + A_{p} y_{t-p} + B \epsilon_{t}.}
#'   In general, identification procedures determine \eqn{B} up to 
#'   column ordering, scale, and sign. For a unique solution, \code{id.iv} 
#'   follows the literature on proxy SVAR. 
#'   The \eqn{K_1} columns in \eqn{B = [B_1 : B_2]} of the instrumented 
#'   shocks \eqn{\epsilon_{tk}, k=1,\ldots,K_1,} are ordered first, and their 
#'   variances \eqn{\sigma^2_{\epsilon,k} = 1} are normalized to unity 
#'   (see e.g. Mertens, Ravn 2013:1215 and Lunsford 2015:6, Eq. 9). 
#'   A normalization of the impulsed shock that may fix the size of the impact 
#'   response in the IRF can be imposed subsequently via '\code{normf}' 
#'   in \code{\link{irf.varx}} and \code{\link{sboot.mb}}.
#'   
#' @param x VAR object of class '\code{varx}' or any other 
#'   that will be \link[=as.varx]{coerced} to '\code{varx}'. 
#' @param iv Matrix. A \eqn{(L \times T)} data matrix of the \eqn{L} proxy 
#'   time series \eqn{m_t}. 
#' @param S2 Character. Selection of the exact identification within multiple instrumented shock series, see 'Details'. 
#'   In case of \eqn{K_1=1} instrumented shock series, all four choices provide identical results on \eqn{B_1} up to the sign.
#' @param cov_u Character. Selection of the estimated residual covariance matrix \eqn{\hat{\Sigma}_{u}} 
#'   to be used in the identification procedure. 
#'   Either \code{'OMEGA'} (the default) for \eqn{\hat{U} \hat{U}'/T_i} as used in Mertens, Ravn (2013) and Jentsch, Lunsford (2022)
#'   or \code{'SIGMA'} for \eqn{\hat{U}\hat{U}'/(T-n_{z})}, which corrects for the number of regressors \eqn{n_z}. 
#'   Both character options refer to the name of the respective estimate in the '\code{varx}' object.
#' @param dim_K1 Integer. Number of shock series instrumented by the \eqn{L \ge K_1} proxy variables. 
#'   If \code{NULL} (the default), the first \eqn{L=K_1} shock series are instrumented.
#'   
#' @details For the exact identification within the \eqn{K_1 > 1} shock series, 
#'   the literature on proxy SVAR provides the following templates to impose 
#'   \eqn{K_1(K_1-1)/2} additional restrictions. Those are available via \code{S2} as:
#'   \itemize{
#'   \item{'\code{MR}'} for the lower-triangular Schur complement 
#'   \eqn{S_1 = B_{11} - B_{12} B_{22}^{-1} B_{21}} by Mertens, Ravn (2013), 
#'   \item{'\code{JL}'} for the lower-triangular covariance matrix \eqn{\Sigma_{m\epsilon 1}}
#'   of proxies \eqn{m_{t}} and instrumented shocks \eqn{\epsilon_{1,t}} by Jentsch, Lunsford (2022),
#'   \item{'\code{NQ}'} for the nearest orthogonal matrix by Empting et al. (2026), or
#'   \item{'\code{LT}'} for lower-triangular blocks \eqn{Q_{11}} and \eqn{B_{11}} by Empting et al. (2026, Appendix A).
#'   }
#'   If the argument \code{dim_K1} is not \code{NULL}, 
#'   the identification procedure resorts to \code{\link[base]{svd}} and 
#'   a subsequent recovery of the selected \code{S2} restrictions 
#'   (see Empting et al. 2026, Appendix A). 
#'   The solution by \code{\link[base]{svd}} further provides the orthogonal 
#'   complement \eqn{Q_{2}} for a full identification of \eqn{Q = [Q_{1} : Q_{2}]}, 
#'   where \eqn{K_2(K_2-1)/2} additional restrictions are imposed as lower-triangular 
#'   blocks \eqn{Q_{12}} and \eqn{B_{12}} irrespective of the choice of \code{S2}.
#' 
#' @return List of class '\code{\link[=as.varx]{id}}'.
#' 
#' @references Mertens, K., and Ravn, M. O. (2013):
#'   "The Dynamic Effects of Personal and Corporate Income Tax Changes in the 
#'   United States", \emph{American Economic Review}, 103, pp. 1212-1247.
#' @references Lunsford, K. G. (2015):
#'   "Identifying Structural VARs with a Proxy Variable and a Test for a Weak Proxy",
#'   Working Paper, No 15-28, Federal Reserve Bank of Cleveland.
#' @references Jentsch, C., and Lunsford, K. G. (2019):
#'   "The Dynamic Effects of Personal and Corporate Income Tax Changes in the 
#'   United States: Comment", \emph{American Economic Review}, 109, pp. 2655-2678.
#' @references Jentsch, C., and Lunsford, K. G. (2022):
#'   "Asymptotically Valid Bootstrap Inference for Proxy SVARs",
#'   \emph{Journal of Business and Economic Statistics}, 40, pp. 1876-1891.
#' @references Empting, L. F. T., Oeztuerk, S., Maxand, S., and Wagner, K. (2026): 
#'   "Inference in Panel SVARs with Two-Way Dependence", arXiv:2610.11898 [econ.EM].
#' @seealso \ldots the individual identification approaches 
#'   by Lange et al. (2021) in \strong{svars}.
#' @family identification functions
#' @example inst/examples/id_iv.R
#' @export
#' 
id.iv <- function(x, iv, S2=c("MR", "JL", "NQ", "LT"), cov_u="OMEGA", dim_K1=NULL){
  # define
  x = as.varx(x)
  
  # identify via proxy variables
  R.id = aux_idIV(u=x$resid, m=iv, S2=S2, SIGMA_uu=x[[cov_u]], dim_K1=dim_K1)
  
  # return result
  x$B = R.id$B1  # estimated impact matrix B
  x$Q = R.id$Q   # estimated orthogonal matrix
  x$udv = R.id$udv   # SVD results
  x$F_stats = R.id$F_stats  # F-statistic on proxy relevance
  x$args_id = list(method="Proxy", iv=R.id$m, S2=S2, cov_u=cov_u, dim_K1=dim_K1)
  class(x)  = c("id", "varx")
  return(x)
}


#' @title Identification of SVEC models by imposing long- and short-run restrictions
#' @description Identifies an SVEC model by utilizing a scoring algorithm 
#'   to impose long- and short-run restrictions. 
#'   See the details of \code{\link[vars]{SVEC}} in \strong{vars}.
#' @param x VAR object of class '\code{varx}' estimated under rank-restriction 
#'   or any other that will be \link[=as.varx]{coerced} to '\code{varx}'.
#' @param LR Matrix. The restricted long-run impact matrix.
#' @param SR Matrix. The restricted contemporaneous impact matrix.
#' @param start Vector. The starting values for \eqn{\gamma}, 
#'   set by \code{\link[stats]{rnorm}} if \code{NULL} (the default).
#' @param max.iter Integer. The maximum number of iterations.
#' @param conv.crit Real number. Convergence value of algorithm.
#' @param maxls Real number. Maximum movement of the parameters between two iterations of the scoring algorithm.
#' 
#' @return List of class '\code{\link[=as.varx]{id}}'.
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
#'   Note that \code{\link{id.grt}} is just a graftage, but allows for the additional 
#'   model specifications in \code{\link{VECM}} and for the bootstrap procedures 
#'   in \code{\link{sboot.mb}}, both provided by the \strong{pvars} package.
#' 
#' @examples
#' ### reproduce basic example in "vars" ###
#' library(vars)
#' data("Canada")
#' names_k = c("prod", "e", "U", "rw")  # variable names
#' names_s = NULL  # optional shock names
#' 
#' # colnames of the restriction matrices are passed as shock names #
#' SR = matrix(NA, nrow=4, ncol=4, dimnames=list(names_k, names_s))
#' SR[4, 2] = 0
#' LR = matrix(NA, nrow=4, ncol=4, dimnames=list(names_k, names_s))
#' LR[1, 2:4] = 0
#' LR[2:4, 4] = 0
#' 
#' # estimate and identify SVECM #
#' R.vecm = VECM(y=Canada[ , names_k], dim_p=3, dim_r=1, type="Case4")
#' R.grt  = id.grt(R.vecm, LR=LR, SR=SR)
#' 
#' @family identification functions
#' @export
#' 
id.grt <- function(x, LR=NULL, SR=NULL, start=NULL, max.iter=100, conv.crit=1e-07, maxls=1.0){
  # define
  x = as.varx(x)
  dim_K = x$dim_K  # number of endogenous variables
  dim_p = x$dim_p  # number of lags
  dim_T = x$dim_T  # number of observations without presample
  
  # calculate long-run multiplier matrix
  alpha_oc = MASS::Null(x$VECM$alpha)
  beta_oc  = MASS::Null(x$beta[1:dim_K, ])
  XI = aux_vec2vma(GAMMA=x$VECM$GAMMA, alpha_oc=alpha_oc, beta_oc=beta_oc, dim_p=dim_p, B=NULL, n.ahead="XI")
  
  # identify via FIML scoring algorithm
  GRT = aux_idGRT(OMEGA=x$VECM$OMEGA, XI=XI, dim_T=dim_T, LR=LR, SR=SR, 
                  start=start, max.iter=max.iter, conv.crit=conv.crit, maxls=maxls)
  
  # return result
  x$B  = GRT$SR  # estimated B matrix (unique decomposition of the covariance matrix)
  x$SR = GRT$SR  # structural effects on impact
  x$LR = GRT$LR  # structural long-run effects
  x$args_id = list(method="SVECM", SR=SR, LR=LR, start=start, max.iter=max.iter, conv.crit=conv.crit, maxls=maxls)
  class(x)  = c("id", "varx")
  return(x)
}


