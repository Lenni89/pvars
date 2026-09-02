

#' @title Bootstrap with residual moving blocks for individual SVAR models
#' @description Calculates confidence bands for impulse response functions via recursive-design bootstrap.
#' @param x VAR object of class '\code{id}' or '\code{varx}' or any other 
#'   that can be \link[=as.varx]{coerced} to '\code{varx}', e.g. '\code{svars}'.
#'   If a bias term \code{x$PSI_bc} is available for coefficient matrix \eqn{A} (such as in '\code{sboot2}'), 
#'   the bias-corrected second-step of the bootstrap-after-bootstrap procedure by Kilian (1998) is performed.
#' @param b.length Integer. Length \eqn{b_{(t)}} of each residual time series block, 
#'   which is often set to \eqn{T/10}.
#' @param n.ahead Integer. Number of periods to consider after the initial impulse, i.e. the horizon of the IRF.
#' @param n.boot Integer. Number of bootstrap iterations.
#' @param n.cores Integer. Number of allocated processor cores.
#' @param fix_beta Logical. If \code{TRUE} (the default), the cointegrating vectors \eqn{\beta} 
#'   are fixed over all bootstrap iterations. Ignored in case of rank-unrestricted VAR. 
#'   Use this for VECM with known \eqn{\beta}, too. Note that \eqn{\beta} is fixed in \code{vars:::.bootsvec}, 
#'   but not in \code{vars:::.bootirfsvec} or \code{vars:::.bootirfvec2var}.
#' @param deltas Vector. Numeric weights \eqn{\delta_j} that are successively 
#'   multiplied to the bias estimate \eqn{\hat{\Psi}} for a stationarity correction. 
#'   The default weights \code{deltas = cumprod((100:0)/100)} correspond 
#'   to the iterative correction procedure of Step 1b in Kilian (1998). Choosing 
#'   \code{deltas = 1} inhibits the stationarity correction, and 
#'   \code{deltas = 0} inhibits the bias-correction but still returns a bias estimate. 
#'   However, \code{deltas = NULL} deactivates the bootstrap-after-bootstrap procedure.
#' @param normf Function. A given function that normalizes the \eqn{K \times S} input-matrix 
#'   into an output matrix of same dimension. See the example in \code{\link{id.iv}} 
#'   for the normalization of Jentsch and Lunsford (2022) 
#'   that fixes the size of the impact response in the IRF.
#' 
#' @return A list of class '\code{sboot2}' with elements:
#' \item{true}{Point estimate of impulse response functions.}
#' \item{bootstrap}{List of length '\code{n.boot}' holding bootstrap impulse response functions.}
#' \item{A}{List for the VAR coefficients containing 
#'   the matrix of point estimates '\code{par}' and 
#'   the array of bootstrap results '\code{sim}'.}
#' \item{B}{List for the structural impact matrix containing 
#'   the matrix of point estimates '\code{par}' and 
#'   the array of bootstrap results '\code{sim}'.}
#' \item{PSI_bc}{Matrix of the estimated bias term \eqn{\hat{\Psi}} 
#'   for the VAR coefficients \eqn{\hat{A}} according to Kilian (1998).}
#' \item{varx}{Input VAR object of class '\code{varx}' 
#'   that has been subjected to the first-step bias-correction.}
#' \item{nboot}{Number of correct bootstrap iterations.}
#' \item{b_length}{Length of each block.}
#' \item{design}{Character indicating that the recursive design bootstrap has been performed.}
#' \item{method}{Used bootstrap method.}
#' \item{stars}{Matrix of (\eqn{T \times }\code{n.boot}) integers containing 
#'   the \eqn{T} resampling draws from each bootstrap iteration.}
#' 
#' @seealso \code{\link[svars]{mb.boot}}, \code{\link[vars]{irf}}, 
#'   and the panel counterpart \code{\link{sboot.pmb}}.
#' 
#' @examples 
#' \donttest{
#' # select minimal or full example #
#' is_min = TRUE
#' n.boot = ifelse(is_min, 5, 500)
#' 
#' # use 'b.length=1' to conduct basic "vars" bootstraps #
#' set.seed(23211)
#' data("Canada")
#' R.vars = vars::VAR(Canada, p=2, type="const")
#' R.svar = svars::id.chol(R.vars)
#' R.boot = sboot.mb(R.svar, b.length=1, n.boot=n.boot, n.ahead=30, n.cores=1)
#' summary(R.boot, idx_par="A", level=0.9)  # VAR coefficients with 90%-confidence intervals 
#' plot(R.boot, lowerq = c(0.05, 0.1, 0.16), upperq = c(0.95, 0.9, 0.84))
#' 
#' # second step of bootstrap-after-bootstrap #
#' R.bab = sboot.mb(R.boot, b.length=1, n.boot=n.boot, n.ahead=30, n.cores=1)
#' summary(R.bab, idx_par="A", level=0.9)  # VAR coefficients with 90%-confidence intervals 
#' plot(R.bab, lowerq = c(0.05, 0.1, 0.16), upperq = c(0.95, 0.9, 0.84))
#' 
#' # conduct bootstraps for Blanchard-Quah type SVAR from "vars" #
#' set.seed(23211)
#' data("Canada")
#' R.vars = vars::VAR(Canada, p=2, type="const")
#' R.svar = vars::BQ(R.vars)
#' R.boot = sboot.mb(R.svar, b.length=1, n.boot=n.boot, n.ahead=30, n.cores=1)
#' summary(R.boot, idx_par="B", level=0.9)  # impact matrix with 90%-confidence intervals 
#' plot(R.boot, lowerq = c(0.05, 0.1), upperq = c(0.95, 0.9), cumulative=2:3) 
#' # impulse responses of the second and third variable are accumulated
#' 
#' # set 'args_id' to CvM defaults of "svars" bootstraps #
#' set.seed(23211)
#' data("USA")
#' R.vars = vars::VAR(USA, lag.max=10, ic="AIC")
#' R.cob  = copula::indepTestSim(R.vars$obs, R.vars$K, verbose=FALSE)
#' R.svar = svars::id.cvm(R.vars, dd=R.cob)
#' 
#' R.varx = as.varx(R.svar, dd=R.cob, itermax=300, steptol=200, iter2=50)
#' R.boot = sboot.mb(R.varx, b.length=15, n.boot=n.boot, n.ahead=30, n.cores=1)
#' plot(R.boot, lowerq = c(0.05, 0.1, 0.16), upperq = c(0.95, 0.9, 0.84))
#' }
#' 
#' @references Breitung, J., Brueggemann R., and Luetkepohl, H. (2004): 
#'   "Structural Vector Autoregressive Modeling and Impulse Responses", 
#'   in \emph{Applied Time Series Econometrics}, 
#'   ed. by H. Luetkepohl and M. Kraetzig, 
#'   Cambridge University Press, Cambridge.
#' @references Brueggemann R., Jentsch, C., and Trenkler, C. (2016): 
#'   "Inference in VARs with Conditional Heteroskedasticity of Unknown Form", 
#'   \emph{Journal of Econometrics}, 191, pp. 69-85.
#' @references Jentsch, C., and Lunsford, K. G. (2022): 
#'   "Asymptotically Valid Bootstrap Inference for Proxy SVARs",
#'   \emph{Journal of Business and Economic Statistics}, 40, pp. 1876-1891.
#' @references Kilian, L. (1998): 
#'   "Small-Sample Confidence Intervals for Impulse Response Functions",
#'   \emph{Review of Economics and Statistics}, 80, pp. 218-230.
#' @references Luetkepohl, H. (2005): 
#'   \emph{New Introduction to Multiple Time Series Analysis}, 
#'   Springer, 2nd ed.
#' 
#' @importFrom steadyICA steadyICA 
#' @export
#' 
sboot.mb <- function(x, b.length=1, n.ahead=20, n.boot=500, n.cores=1, fix_beta=TRUE, deltas=cumprod((100:0)/100), normf=NULL){
  # define
  PSI_bc = x$PSI_bc  # bias term (Note that 'x' will be overwritten for exogenous variables.)
  if( is.null(PSI_bc) & inherits(x, "sboot") ){ warning("The bootstrap is performed without bias-correction due to unavailable bias estimates.") }
  y = D = D1 = D2 = A = B = beta = resid = NULL
  dim_K = dim_S = dim_T = dim_p = dim_r = t_D1 = t_D2 = NULL
  args_varx = args_id = NULL
  R.varx = aux_assign_varx(x)
  A_hat = A  # estimated VAR coefficients incl. impulse dummies for Dnl
  
  # names for variables, for shocks, and for the header of each IRF
  names_k   = if( !is.null(rownames(A)) ){ rownames(A) }else{ paste0("y[ ", 1:dim_K, " ]") }
  names_s   = if( !is.null(colnames(B)) ){ colnames(B) }else{ paste0("epsilon[ ", 1:dim_S, " ]") }
  names_IRF = c(sapply(names_k, FUN=function(k) paste0(names_s, " %->% ", k)))
  dimnames(B) = list(names_k, names_s)
  
  # calculate structural IRF of the original model (point estimates)
  vma = aux_var2vma(A=A, B=B, dim_p=dim_p, n.ahead=n.ahead, normf=normf)
  IRF = aperm(vma$THETA, perm=c(2,1,3))
  IRF = matrix(IRF, nrow=n.ahead+1, byrow=TRUE, dimnames=list(NULL, names_IRF))
  IRF = list(irf=as.data.frame(cbind(V1=0:n.ahead, IRF)))
  class(IRF) = "svarirf"
  
  # fixed objects in the bootstrap function
  Ystar = cbind(y[ , 1:dim_p, drop=FALSE], matrix(NA, nrow=dim_K, ncol=dim_T))  # presample for initialization
  if(!is.null(dim_r)){
    beta_oc = aux_oc(beta[1:dim_K, ])
    D = rbind(D2, D1)
    D = aux_rm_Dnl(D, dim_p, t_D1, t_D2, MARGIN=1)
    A = aux_rm_Dnl(A, dim_p, t_D1, t_D2, MARGIN=2)
  }
  
  # resampling in the bootstrap function
  ### moving-block bootstrap (MBB), from Brueggemann et al. 2016:73, Ch.4
  ### residual BB for cointegrated processes, see Jentsch et al. 2015:425, Ch.3.2
  n.block = ceiling(dim_T/b.length)  # number of blocks
  idx_t   = 1:dim_T  # restrict to original sample size via [idx_t]
  idx_r   = 0:(dim_T-b.length)  # "in-sample" indexes to be resampled
  idx1122 = rep(1:n.block,  each = b.length)[idx_t]  # broadcasting index for resampling
  idx1212 = rep(1:b.length, times = n.block)[idx_t]  # broadcasting index for centering
  L.star  = lapply(1:n.boot,   FUN = function(n) sample(idx_r, size=n.block, replace=TRUE)[idx1122] + idx1212)
  Umean   = sapply(1:b.length, FUN = function(s) rowMeans(resid[ , s+idx_r, drop=FALSE]))
  Umean   = Umean[ , idx1212, drop=FALSE]
  
  # bootstrap function
  sbootf <- function(star){
    # generate bootstrap time series
    Ustar = resid[ , star, drop=FALSE]  # resampled residuals
    Ustar = Ustar - Umean               # centered residuals
    idx_t = dim_p + (1:dim_T)           # periods of y to be resampled
    for(t in idx_t){  # recursive design: VAR process
      Ystar[ ,t] = A %*% c(D[ ,t], Ystar[ ,t-(1:dim_p)]) + Ustar[ ,t-dim_p] }
    
    # re-estimate the VAR model
    if(!is.null(args_varx$varxf)){  # ... by user-defined estimator
      star_var  = args_varx$varxf(Ystar, args_varx)
      star_beta = star_var$beta
      
    }else if(is.null(dim_r)){  # ... by least squares
      star_var  = aux_VARX(Ystar, dim_p=dim_p, D=D)
      star_beta = NULL
      if(!is.null(PSI_bc) & !is.null(deltas)){  
        # ... under second-step bias-correction, from Kilian 1998:220 (Step 2b)
        star_var$A = aux_BaB(star_var$A, dim_p=dim_p, PSI=PSI_bc, deltas=deltas)
      }
      
    }else{  # ... by reduced-rank regression
      star_def  = aux_stackRRR(Ystar, dim_p=dim_p, D1=D1, D2=D2)
      star_RRR  = aux_RRR(star_def$Z0, star_def$Z1, star_def$Z2)
      star_beta = if(fix_beta){ beta }else{ aux_beta(V=star_RRR$V, dim_r=dim_r, normalize="natural") }
      star_var  = aux_VECM(RRR=star_RRR, beta=star_beta)
      star_var$A= aux_vec2var(PI=star_var$PI, GAMMA=star_var$GAMMA, dim_p=dim_p)$A
    }
    
    # re-identify the structural impact matrix 'B'
    if(!is.null(args_id$idf)){  # ... by user-defined identification procedure
      star_B = args_id$idf(star_var, args_id)
      
    }else if(is.null(args_id)){
      star_B = diag(dim_K)  # for forecast-error IRF
      ### TODO: accommodate LAMBDA - aux_con2var()$B or aux_con2vec()$B
      
    }else if(args_id$method == "Cholesky"){
      perm      = diag(dim_K)[ , args_id$order_k, drop=FALSE]  # permutation matrix
      star_psig = t(perm) %*% star_var$SIGMA %*% perm
      star_P    = t(chol(star_psig))
      star_B    = perm %*% star_P %*% t(perm)
      
    }else if(args_id$method == "Blanchard-Quah"){
      # see Luetkepohl 2005:376
      star_A1 = aux_accum(star_var$A, dim_p)
      star_Xi = solve(star_A1)
      star_LR = t(chol(star_Xi %*% star_var$SIGMA %*% t(star_Xi)))
      star_B  = star_A1 %*% star_LR
      
    }else if(args_id$method == "SVECM"){
      alpha_oc = aux_oc(star_var$alpha)
      beta_oc  = if(fix_beta){ beta_oc }else{ aux_oc(star_beta[1:dim_K, ]) }
      star_XI  = aux_vec2vma(GAMMA=star_var$GAMMA, alpha_oc=alpha_oc, beta_oc=beta_oc, dim_p=dim_p, B=NULL, n.ahead="XI")
      star_B   = aux_idGRT(OMEGA=star_var$OMEGA, XI=star_XI, 
                           dim_T=dim_T, LR=args_id$LR, SR=args_id$SR, 
                           start=args_id$start, max.iter=args_id$max.iter, 
                           conv.crit=args_id$conv.crit, maxls=args_id$maxls)$SR
      
    }else if(args_id$method == "Distance covariances"){
      star_P   = t(chol(star_var$SIGMA))
      star_eps = solve(star_P) %*% star_var$resid  # pre-whitened shocks form baseline decomposition
      star_ICA = suppressMessages(steadyICA::steadyICA(t(star_eps), symmetric=TRUE, PIT=args_id$PIT, maxit=args_id$n.iterations))
      star_B   = aux_sico(star_P %*% star_ICA$W, B.orig=B)
      
    }else if(args_id$method == "Cramer-von Mises distance"){
      star_P   = t(chol(star_var$SIGMA))
      star_eps = solve(star_P) %*% star_var$resid  # pre-whitened shocks form baseline decomposition
      star_ICA = aux_cvmICA(t(star_eps), dd=args_id$dd, itermax=args_id$itermax, steptol=args_id$steptol, iter2=args_id$iter2)
      star_B   = aux_sico(star_P %*% star_ICA$W, B.orig=B)
      
    }else if(args_id$method == "Proxy"){
      # MBB with multiple proxy IV, see Jentsch, Lunsford 2022:9, Ch.4.1
      star_iv = args_id$iv[ , star, drop=FALSE]  # resampled proxies (not centered anymore, see their ft.26)
      star_B  = aux_idIV(star_var$resid, m=star_iv, S2=args_id$S2, SIGMA_uu=star_var[[args_id$cov_u]], dim_K1=args_id$dim_K1)$B1
      ### TODO. Confidence intervals for the correlations of ... 
      ### ... proxies, structural shocks, and regression errors, see Lunsford 2015:22, Tab.3
      
    }else{
      stop("Identification method is not available in this bootstrap function!")
    }
    
    # return bootstrap result
    star_vma = aux_var2vma(A=star_var$A, B=star_B, dim_p=dim_p, n.ahead=n.ahead, normf=normf)
    star_IRF = aperm(star_vma$THETA, perm=c(2,1,3))
    star_IRF = matrix(star_IRF, nrow=n.ahead+1, byrow=TRUE, dimnames=list(NULL, names_IRF))
    star_IRF = list(irf=as.data.frame(cbind(V1=0:n.ahead, star_IRF))); class(star_IRF) = "svarirf"
    result   = list(star_IRF, Pstar=star_B, star_A=star_var$A, star_beta=star_beta)
    return(result)
  }
  
  # run bootstrap procedure
  L.boot    = pbapply::pblapply(L.star, FUN=function(n) try(sbootf(n), silent=TRUE), cl=n.cores)
  idx_error = sapply(L.boot, FUN=function(n) inherits(n, "try-error"))
  if(any(idx_error)){
    L.star[idx_error] = NULL
    L.boot[idx_error] = NULL
    n.boot = length(L.boot)
    warning("Erroneous bootstrap replications have been deleted: ", sum(idx_error), " run(s)")
  }
  
  ipb   = list()
  Bs    = array(NA, dim=c(dim_K, dim_S, n.boot), dimnames=list(names_k, names_s, NULL))
  Aboot = array(NA, dim=c(dim_K, ncol(A_hat), n.boot), dimnames=list(names_k, colnames(A_hat), NULL))
  for(j in 1:n.boot){
    ipb[[j]]    = L.boot[[j]][[1]]  # IRF
    Bs[,, j]    = L.boot[[j]][[2]]  # structural impact matrix
    Aboot[,, j] = L.boot[[j]][[3]]  # VAR coefficients
  }
  
  if(!is.null(dim_r)){
    betas = array(NA, dim=c(dim(beta), n.boot), dimnames=dimnames(beta))
    for(j in 1:n.boot){ betas[,, j] = L.boot[[j]][[4]] }
    betas = aperm(betas, perm=c(2, 1, 3))
    beta  = t(beta)
  }else{ betas = NULL }
  
  # first-step bias-correction of VAR coefficients, from Kilian 1998:220
  if(is.null(PSI_bc) & !is.null(deltas)){
    # Step 1a: bias estimate
    R.varx$PSI_bc = PSI_bc = apply(Aboot, MARGIN=1:2, mean) - A_hat
    # Step 1b: stationary correction
    R.varx$A = aux_BaB(A_hat=A_hat, dim_p=dim_p, PSI=PSI_bc, deltas=deltas)
  }
  
  # return result
  result = list(true = IRF,
                bootstrap = ipb,
                A = list(par=A_hat, sim=Aboot),
                B = list(par=B, sim=Bs),
                beta = list(par=beta, sim=betas),
                PSI_bc = PSI_bc,
                varx = R.varx,
                nboot = n.boot,
                b_length = b.length,
                rest_mat = matrix(NA, nrow=dim_K, ncol=dim_S),  # for 'bonferroni'
                design = "recursive",
                method = "Moving block bootstrap",
                stars = matrix(unlist(L.star), nrow=dim_T))
  class(result) = c("sboot2", "sboot")
  return(result)
}


#' @title Bootstrap with residual panel blocks for panel SVAR models
#' @description Calculates confidence bands for impulse response functions via recursive-design bootstrap. 
#' @details In case of heterogeneous lag-orders \eqn{p_i} or sample sizes \eqn{T_i},
#'   the initial periods are fixed in accordance with the usage of presamples. 
#'   Only the \eqn{(K \times T_{min} \times N)} array of the \eqn{T_{min} = min(T_1,\ldots,T_N)} 
#'   last residuals is resampled.
#'   
#' @param x Panel VAR object of class '\code{pid}' or '\code{pvarx}' 
#'   or a list of VAR objects that will be \link[=as.varx]{coerced} to '\code{varx}'.
#'   If a list \code{x$L.PSI_bc} of \eqn{N} bias terms are available 
#'   for the \eqn{N} coefficient matrices \eqn{A_i} (such as in \code{sboot2}), 
#'   the bias-corrected second-step of the bootstrap-after-bootstrap procedure 
#'   by Empting et al. (2026) is performed.
#' @param b.dim Vector of two integers. The dimensions \eqn{(b_{(t)}, b_{(i)})} 
#'   of each residual panel block for temporal and cross-sectional resampling. The default 
#'   \code{c(1, 1)} specifies an \eqn{i.i.d.} resampling in both dimensions, 
#'   \code{c(1, FALSE)} a temporal resampling, and 
#'   \code{c(FALSE, 1)} a cross-sectional resampling. 
#'   Integers \eqn{> 1} assemble blocks of consecutive residuals. 
#' @param n.ahead Integer. Number of periods to consider after the initial impulse, i.e. the horizon of the IRF.
#' @param n.boot Integer. Number of bootstrap iterations.
#' @param n.cores Integer. Number of allocated processor cores.
#' @param fix_beta Logical. If \code{TRUE} (the default), the cointegrating vectors \eqn{\beta} 
#'   are fixed over all bootstrap iterations. Ignored in case of rank-unrestricted VAR. 
#'   Use this for VECM with known \eqn{\beta}, too. Note that \eqn{\beta} is fixed in \code{vars:::.bootsvec}, 
#'   but not in \code{vars:::.bootirfsvec} or \code{vars:::.bootirfvec2var}.
#' @param deltas Vector. Numeric weights \eqn{\delta_j} that are successively 
#'   multiplied to each bias estimate \eqn{\hat{\Psi}_i} for a stationarity correction. 
#'   The default weights \code{deltas = cumprod((100:0)/100)} correspond 
#'   to the iterative correction procedure of Step 1b in Kilian (1998). Choosing 
#'   \code{deltas = 1} inhibits the stationarity correction, and 
#'   \code{deltas = 0} inhibits the bias-correction but still returns a bias estimate. 
#'   However, \code{deltas = NULL} deactivates the bootstrap-after-bootstrap procedure.
#' @param normf Function. A given function that normalizes the \eqn{K \times S} input-matrix 
#'   into an output matrix of same dimension. See the example in \code{\link{id.iv}} 
#'   for the normalization of Jentsch and Lunsford (2022) 
#'   that fixes the size of the impact response in the IRF.
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
#' @return A list of class '\code{sboot2}' with elements:
#' \item{true}{Mean group estimate of impulse response functions.}
#' \item{bootstrap}{List of length \code{nboot} holding bootstrap impulse response functions.}
#' \item{A}{List for the VAR coefficients containing 
#'   the matrix of point estimates '\code{par}' and 
#'   the array of bootstrap results '\code{sim}'.}
#' \item{B}{List for the structural impact matrix containing 
#'   the matrix of point estimates '\code{par}' and 
#'   the array of bootstrap results '\code{sim}'.}
#' \item{L.PSI_bc}{List of the \eqn{N} estimated bias terms \eqn{\hat{\Psi}_i} 
#'   for the individual VAR coefficients \eqn{\hat{A}_i} according to Kilian (1998).}
#' \item{pvarx}{Input panel VAR object of class '\code{pvarx}' 
#'   that has been subjected to the first-step bias-correction.}
#' \item{b.dim}{Dimensions of each block.}
#' \item{nboot}{Number of correct bootstrap iterations.}
#' \item{design}{Character indicating that the recursive design bootstrap has been performed.}
#' \item{method}{Used bootstrap method.}
#' \item{stars_t}{Matrix of (\eqn{T \times }\code{n.boot}) integers containing 
#'   the \eqn{T} temporal resampling draws from each bootstrap iteration.}
#' \item{stars_i}{Matrix of (\eqn{N \times }\code{n.boot}) integers containing 
#'   the \eqn{N} cross-sectional resampling draws from each bootstrap iteration.}
#' 
#' @seealso For the the individual counterpart see \code{\link{sboot.mb}}.
#' 
#' @references Brueggemann R., Jentsch, C., and Trenkler, C. (2016): 
#'   "Inference in VARs with Conditional Heteroskedasticity of Unknown Form", 
#'   \emph{Journal of Econometrics}, 191, pp. 69-85.
#' @references Empting, L. F. T., Maxand, S., Oeztuerk, S., and Wagner, K. (2026): 
#'   "Inference in Panel SVARs with Two-Way Dependence".
#' @references Kapetanios, G. (2008): 
#'   "A Bootstrap Procedure for Panel Data Sets with many Cross-sectional Units", 
#'   \emph{The Econometrics Journal}, 11, pp.377-395.
#' @references Kilian, L. (1998): 
#'   "Small-Sample Confidence Intervals for Impulse Response Functions",
#'   \emph{Review of Economics and Statistics}, 80, pp. 218-230.
#' @references Gambacorta L., Hofmann B., and Peersman G. (2014):
#'   "The Effectiveness of Unconventional Monetary Policy at the Zero Lower Bound: A Cross-Country Analysis",
#'   \emph{Journal of Money, Credit and Banking}, 46, pp. 615-642.
#' 
#' @examples
#' \donttest{
#' # select minimal or full example #
#' is_min = TRUE
#' n.boot = ifelse(is_min, 5, 500)
#' 
#' # prepare data panel #
#' data("PCAP")
#' names_k = c("g", "k", "l", "y")  # variable names
#' names_i = levels(PCAP$id_i)      # country names 
#' L.data  = sapply(names_i, FUN=function(i) 
#'   ts(PCAP[PCAP$id_i==i, names_k], start=1960, end=2019, frequency=1), 
#'   simplify=FALSE)
#' R.lags = c(2, 4, 2, 3, 2, 4, 4, 2, 2, 3, 3, 3, 2, 4, 4, 2, 2, 2, 4, 2, 2, 2, 4)
#' names(R.lags) = names_i
#' 
#' # estimate, identify, and bootstrap #
#' R.pvar = pvarx.VAR(L.data, lags=R.lags, type="both")
#' R.pid  = pid.chol(R.pvar)
#' R.boot = sboot.pmb(R.pid, n.boot=n.boot)
#' summary(R.boot, idx_par="A", level=0.95)  # VAR coefficients with 95%-confidence intervals
#' plot(R.boot, lowerq = c(0.05, 0.1, 0.16), upperq = c(0.95, 0.9, 0.84))
#' 
#' # second step of bootstrap-after-bootstrap #
#' R.bab = sboot.pmb(R.boot, n.boot=n.boot)
#' summary(R.bab, idx_par="A", level=0.95)  # VAR coefficients with 95%-confidence intervals
#' plot(R.bab, lowerq = c(0.05, 0.1, 0.16), upperq = c(0.95, 0.9, 0.84))
#' }
#' 
#' @export
#' 
sboot.pmb <- function(x, b.dim=c(1, 1), n.ahead=20, n.boot=500, n.cores=1, fix_beta=TRUE, deltas=cumprod((100:0)/100), normf=NULL, w=NULL, MG_IRF=FALSE){
  # define
  L.PSI_bc = x$L.PSI_bc  # list of individual bias terms
  if( is.null(L.PSI_bc) & inherits(x, "sboot") ){ warning("The bootstrap is performed without bias-correction due to unavailable bias estimates.") }
  A_MG = B_MG = beta = dim_N = dim_K = dim_S = dim_r = NULL
  L.varx = L.dim_T = L.dim_p = L.resid = L.beta = L.A = L.B = L.D = L.D1 = L.D2 = NULL
  args_pvarx = args_pid = NULL
  R.pvarx = aux_assign_pvarx(x, w=w)
  
  # names for variables, for shocks, and for the header of each IRF
  names_k   = if( !is.null(rownames(A_MG)) ){ rownames(A_MG) }else{ paste0("y[ ", 1:dim_K, " ]") }
  names_s   = if( !is.null(colnames(B_MG)) ){ colnames(B_MG) }else{ paste0("epsilon[ ", 1:dim_S, " ]") }
  names_IRF = c(sapply(names_k, FUN=function(k) paste0(names_s, " %->% ", k)))
  dimnames(B_MG) = list(names_k, names_s)
  
  # calculate structural IRF of the original model (point estimates) under 'w'
  IRF = irf.pvarx(R.pvarx, n.ahead=n.ahead, normf=normf, w=w, MG_IRF=MG_IRF)
  
  # fixed objects in the bootstrap function
  dim_T   = min(L.dim_T)  # number of time periods without presample in total panel
  L.Ystar = lapply(L.varx, FUN=function(i)  # presample for initialization
    cbind(head(i$y, n=c(dim_K, -dim_T)), matrix(NA, nrow=dim_K, ncol=dim_T)))
  if(!is.null(dim_r)){
    L.D = lapply(L.varx, FUN=function(i) aux_rm_Dnl(rbind(i$D2, i$D1), i$dim_p, i$t_D1, i$t_D2, MARGIN=1))
    L.A = lapply(L.varx, FUN=function(i) aux_rm_Dnl(i$A, i$dim_p, i$t_D1, i$t_D2, MARGIN=2))
    L.beta = if(fix_beta){ L.beta }else{ NULL }
  }
  
  # resampling in the bootstrap function with ...
  L.star  = lapply(1:n.boot, FUN=function(n) list(star_t=NULL, star_i=NULL))
  L.resid = lapply(L.resid, FUN=function(x) tail(x, n=c(dim_K, dim_T)))
  L.Umean = list()
  idx_t = 1:dim_T  # restricts to original sample size via [idx_t]
  idx_i = 1:dim_N  # ... and via [idx_i]
  
  # ... temporal resampling
  if(b.dim[1] > 0){
    n.blockt = ceiling(dim_T/b.dim[1])  # number of temporal blocks to reach T
    idx_rt   = 0:(dim_T-b.dim[1])  # "in-sample" indexes to be resampled
    idx1122t = rep(1:n.blockt,  each = b.dim[1])[idx_t]  # broadcasting indexes for resampling
    idx1212t = rep(1:b.dim[1], times = n.blockt)[idx_t]  # ... and for centering
    
    for(n in 1:n.boot){
      L.star[[n]][[1]] = sample(idx_rt, size=n.blockt, replace=TRUE)[idx1122t] + idx1212t }
    for(i in 1:dim_N){  # account for draw-probability of end periods
      Umean = sapply(1:b.dim[1], FUN = function(s) rowMeans(L.resid[[i]][ , s+idx_rt, drop=FALSE]))
      L.Umean[[i]] = Umean[ , idx1212t, drop=FALSE]}
  }else{ 
    for(n in 1:n.boot){ L.star[[n]][[1]] = idx_t }
    L.Umean = lapply(L.resid, FUN = function(i) rowMeans(i))}
  
  # ... cross-sectional resampling
  if(b.dim[2] > 0){
    n.blocki = ceiling(dim_N/b.dim[2])  # number of cross-sectional blocks to reach N
    idx_ri   = 0:(dim_N-1)  # "in-sample" indexes to be resampled
    idx1122i = rep(1:n.blocki,  each = b.dim[2])[idx_i]  # broadcasting indexes for resampling
    idx1212i = rep(1:b.dim[2], times = n.blocki)[idx_i]
    
    for(n in 1:n.boot){
      star_i = sample(idx_ri, size=n.blocki, replace=TRUE)[idx1122i] + idx1212i
      idx_N1 = star_i > dim_N
      star_i[idx_N1] = star_i[idx_N1] - dim_N  # assume last and first individuals to be neighbors in block
      L.star[[n]][[2]] = star_i }
  }else{ 
    for(n in 1:n.boot){ L.star[[n]][[2]] = idx_i }}
  
  # bootstrap function
  sbootf <- function(star){
    # generate bootstrap panels
    star_t = star[[1]]  # indexes of temporal resampling
    star_i = star[[2]]  # indexes of cross-sectional resampling
    S.data = list()     # bootstrapped panel data
    
    for(i in 1:dim_N){
      A = L.A[[i]]
      D = L.D[[i]]
      dim_p = L.dim_p[i]
      Ystar = L.Ystar[[i]]
      Ustar = L.resid[[star_i[i]]][ , star_t, drop=FALSE]
      Ustar = Ustar - L.Umean[[star_i[i]]]  # resampled and centered residuals
      dim_0 = ncol(Ystar)-dim_T  # presample size (wrt. total sample)
      idx_t = dim_0 + (1:dim_T)  # periods of y to be resampled
      
      # recursive design: VAR process
      for(t in idx_t){
        Ystar[ ,t] = A %*% c(D[ ,t], Ystar[ ,t-(1:dim_p)]) + Ustar[ ,t-dim_0] }
      S.data[[i]] = Ystar
    }
    
    # re-estimate the VAR models
    if(!is.null(args_pvarx$pvarxf)){  # ... by user-defined estimator
      S.pvarx = args_pvarx$pvarxf(S.data, args_pvarx)
    
    }else if(is.null(dim_r)){  # ... by least squares
      L.def = lapply(1:dim_N, FUN=function(i) aux_stackOLS(S.data[[i]], dim_p=L.dim_p[i], D=L.D[[i]]))
      S.est = aux_pvar(L.def, n.factors=args_pvarx$n.factors, n.iterations=args_pvarx$n.iterations)
      if(!is.null(L.PSI_bc) & !is.null(deltas)){
      for(i in 1:dim_N){  # ... under second-step bias-correction, from Kilian 1998:220 (Step 2b)
        S.est$L.varx[[i]]$A = aux_BaB(S.est$L.varx[[i]]$A, dim_p=L.varx[[i]]$dim_p, PSI=L.PSI_bc[[i]], deltas=deltas)
      }}
      S.mgA = aux_MG(S.est$L.varx, w=w, idx_par="A")$mean
      S.pvarx = list(L.varx=S.est$L.varx, A=S.mgA, args_pvarx=args_pvarx)
      class(S.pvarx) = "pvarx"
      
    }else{  # ... by reduced-rank regression
      L.def = lapply(1:dim_N, FUN=function(i) aux_stackRRR(S.data[[i]], dim_p=L.dim_p[i], D1=L.D1[[i]], D2=L.D2[[i]]))
      S.est = aux_pvec(L.def, L.beta=L.beta, dim_r=dim_r, idx_pool=args_pvarx$idx_pool, 
                       n.factors=args_pvarx$n.factors, n.iterations=args_pvarx$n.iterations)
      S.vec = list(
        alpha = aux_MG(S.est$L.varx, w=w, idx_par="alpha")$mean,
        beta  = aux_MG(S.est$L.varx, w=w, idx_par="beta")$mean,  # flexible towards idx_pool 
        GAMMA = aux_MG(S.est$L.varx, w=w, idx_par="GAMMA")$mean)
      S.vec$PI = S.vec$alpha %*% t(S.vec$beta)
      S.mgA = aux_vec2var(PI=S.vec$PI, GAMMA=S.vec$GAMMA, dim_p=max(L.dim_p))$A  # rank-restricted VAR model in levels
      S.pvarx = list(L.varx=S.est$L.varx, A=S.mgA, beta=S.vec$beta, args_pvarx=args_pvarx)
      class(S.pvarx) = "pvarx"
    }
    
    # re-identify the structural impact matrices 'B'
    if(!is.null(args_pid$pidf)){  # ... by user-defined identification procedure
      S.pid = args_pid$pidf(S.pvarx, args_pid)
      
    }else if(is.null(args_pid$method)){
      S.pid   = S.pvarx
      S.pid$B = B_MG  # diagonal matrix with dimnames
      
    }else if(args_pid$method == "Cholesky"){
      S.pid = pid.chol(S.pvarx, order_k=args_pid$order_k)
      
    }else if(args_pid$method == "Blanchard-Quah"){
      # see Luetkepohl 2005:376
      S.pid = S.pvarx
      for(i in 1:dim_N){
        star_A1 = aux_accum(S.pid$L.varx[[i]]$A, L.dim_p[i])
        star_Xi = solve(star_A1)
        star_LR = t(chol(star_Xi %*% S.pid$L.varx[[i]]$SIGMA %*% t(star_Xi)))
        star_B  = star_A1 %*% star_LR 
        S.pid$L.varx[[i]]$B = star_B
      }
      S.pid$MG_B = aux_MG(S.pid$L.varx, w=w, idx_par="B")
      S.pid$B    = S.pid$MG_B$mean
      
    }else if(args_pid$method == "SVECM"){
      S.pid = pid.grt(S.pvarx, LR=args_pid$LR, SR=args_pid$SR, start=args_pid$start, 
                      max.iter=args_pid$max.iter, conv.crit=args_pid$conv.crit, maxls=args_pid$maxls)
      
    }else if(args_pid$method == "Distance covariances"){
      S.pid = pid.dc(S.pvarx, combine=args_pid$combine, n.factors=args_pid$n.factors, 
                     n.iterations=args_pid$n.iterations)
      # sign and column ordering according to the point estimates
      for(i in 1:dim_N){ S.pid$L.varx[[i]]$B = aux_sico(S.pid$L.varx[[i]]$B, B.orig=L.B[[i]]) }
      S.pid$MG_B = aux_MG(S.pid$L.varx, w=w, idx_par="B")
      S.pid$B    = S.pid$MG_B$mean
      
    }else if(args_pid$method == "Cramer-von Mises distance"){
      S.pid = pid.cvm(S.pvarx, combine=args_pid$combine, n.factors=args_pid$n.factors, 
                      dd=args_pid$dd, itermax=args_pid$itermax, steptol=args_pid$steptol, iter2=args_pid$iter2)
      # sign and column ordering according to the point estimates
      for(i in 1:dim_N){ S.pid$L.varx[[i]]$B = aux_sico(S.pid$L.varx[[i]]$B, B.orig=L.B[[i]]) }
      S.pid$MG_B = aux_MG(S.pid$L.varx, w=w, idx_par="B")
      S.pid$B    = S.pid$MG_B$mean
      
    }else if(args_pid$method == "Proxy"){
      L.iv = lapply(args_pid$L.iv, FUN=function(i) tail(i, c(NA, dim_T)))
      S.iv = args_pid$L.iv
      for(i in 1:dim_N){
        # resample the proxy panel exactly like the residual panel
        dim_0 = ncol(S.iv[[i]]) - dim_T  # match with presample of y (wrt. total sample)
        idx_t = dim_0 + (1:dim_T)  # periods of m_it to be resampled
        S.iv[[i]][ , idx_t] = L.iv[[star_i[i]]][ , star_t, drop=FALSE]
      }
      S.pid = pid.iv(S.pvarx, iv=S.iv, S2=args_pid$S2, cov_u=args_pid$cov_u, 
                     dim_K1=args_pid$dim_K1, combine=args_pid$combine)
        
    }else{
      stop("Identification method is not available in this bootstrap function!")
    }
    
    # return bootstrap result
    S.Ais  = lapply(S.pvarx$L.varx, FUN=function(i) i$A)
    S.irf  = irf.pvarx(S.pid, n.ahead=n.ahead, normf=normf, w=w, MG_IRF=MG_IRF)
    result = list(irf=S.irf, Pstar=S.pid$B, star_A=S.pvarx$A, star_beta=S.pvarx$beta, star_Ais=S.Ais)
    return(result)
  }
  
  # run bootstrap procedure
  L.boot    = pbapply::pblapply(L.star, FUN=function(n) try(sbootf(n), silent=TRUE), cl=n.cores)
  idx_error = sapply(L.boot, FUN=function(n) inherits(n, "try-error"))
  if(any(idx_error)){
    L.star[idx_error] = NULL
    L.boot[idx_error] = NULL
    n.boot = length(L.boot)
    warning("Erroneous bootstrap replications have been deleted: ", sum(idx_error), " run(s)")
  }
  
  ipb   = list()
  Bs    = array(NA, dim=c(dim_K, dim_S, n.boot), dimnames=list(names_k, names_s, NULL))
  Aboot = array(NA, dim=c(dim_K, ncol(A_MG), n.boot), dimnames=list(names_k, colnames(A_MG), NULL))
  for(j in 1:n.boot){
    ipb[[j]]    = L.boot[[j]][[1]]  # IRF
    Bs[,, j]    = L.boot[[j]][[2]]  # structural impact matrix
    Aboot[,, j] = L.boot[[j]][[3]]  # MG-VAR coefficients
  }
  
  if(!is.null(dim_r)){
    betas = array(NA, dim=c(dim(beta), n.boot), dimnames=dimnames(beta))
    for(j in 1:n.boot){ betas[,, j] = L.boot[[j]][[4]] }
    betas = aperm(betas, perm=c(2, 1, 3))
    beta  = t(beta)
  }else{ betas = NULL }
  
  # first-step bias-correction of VAR coefficients, from Kilian 1998:220 / Empting et al. (2026)
  if(is.null(L.PSI_bc) & !is.null(deltas)){
    for(i in 1:dim_N){  
      # Step 1a: bias estimate
      A_hat = R.pvarx$L.varx[[i]]$A
      A_sim = array(NA, dim=c(dim_K, ncol(A_hat), n.boot), dimnames=list(names_k, colnames(A_hat), NULL))
      for(j in 1:n.boot){ A_sim[,, j] = L.boot[[j]][[5]][[i]] }
      R.pvarx$L.varx[[i]]$PSI_bc = L.PSI_bc[[i]] = apply(A_sim, MARGIN=1:2, mean) - A_hat
      
      # Step 1b: stationary correction
      R.pvarx$L.varx[[i]]$A = aux_BaB(A_hat=A_hat, dim_p=R.pvarx$L.varx[[i]]$dim_p, PSI=L.PSI_bc[[i]], deltas=deltas)
    }
    names(L.PSI_bc) = names(R.pvarx$L.varx)
    R.pvarx$A = aux_MG(R.pvarx$L.varx, w=w, idx_par="A")$mean
    R.pvarx$L.PSI_bc = L.PSI_bc
  }
  
  # return result
  result = list(true = IRF,
                bootstrap = ipb,
                A = list(par=A_MG, sim=Aboot),
                B = list(par=B_MG, sim=Bs),
                beta = list(par=beta, sim=betas),
                L.PSI_bc = L.PSI_bc,
                pvarx = R.pvarx,
                nboot = n.boot,
                b.dim = b.dim,
                rest_mat = matrix(NA, nrow=dim_K, ncol=dim_S),  # for 'bonferroni'
                design = "recursive",
                method = "Panel block bootstrap",
                stars_t = sapply(L.star, FUN=function(n) n$star_t),
                stars_i = sapply(L.star, FUN=function(n) n$star_i))
  class(result) = c("sboot2", "sboot")
  return(result)
}


#' @title Mean group inference for panel SVAR models
#' @description Calculates confidence bands for impulse response functions via mean group inference. 
#'   The function does not perform bootstraps, but coerces the panel VAR object to class '\code{sboot}'
#'   and, therewith, gives a distributional overview on the estimated parameter heterogeneity.
#' @details MG inference presumes the individual estimates to be the empirical variation 
#'   around a common parameter. In case of heterogeneous lag-orders \eqn{p_i},
#'   specifically the '\code{summary}' of VAR coefficient matrices fills 
#'   \eqn{\hat{A}_{ij} = 0_{K \times K}} for \eqn{p_i < j \le max(p_1,\ldots,p_N)} 
#'   in accordance with the finite order VAR\eqn{(p_i)}.
#'   
#' @param x Panel VAR object of class '\code{pid}' or '\code{pvarx}' 
#'   or a list of VAR objects that will be \link[=as.varx]{coerced} to '\code{varx}'.
#' @param n.ahead Integer. Number of periods to consider after the initial impulse, i.e. the horizon of the IRF.
#' @param normf Function. A given function that normalizes the \eqn{K \times S} input-matrix 
#'   into an output matrix of same dimension. See the example in \code{\link{id.iv}} 
#'   for the normalization of Jentsch and Lunsford (2022) 
#'   that fixes the size of the impact response in the IRF.
#' @param idx_i Logical or character vector. 
#'   Names or \eqn{N} logical elements selecting a subset from the 
#'   individuals \eqn{i = 1, \ldots, N} for the MG estimation. If \code{NULL} 
#'   (the default), all \eqn{N} individuals are included.
#' 
#' @return A list of class '\code{sboot2}' with elements:
#' \item{true}{Mean group estimate of impulse response functions.}
#' \item{bootstrap}{List of length \eqn{N} holding the individual impulse response functions.}
#' \item{A}{List for the VAR coefficients containing 
#'   the matrix of mean group estimates '\code{par}' and 
#'   the array of individual results '\code{sim}'.}
#' \item{B}{List for the structural impact matrix containing 
#'   the matrix of mean group estimates '\code{par}' and 
#'   the array of individual results '\code{sim}'.}
#' \item{pvarx}{Input panel VAR object of class '\code{pvarx}'.}
#' \item{nboot}{Integer '0' indicating that no bootstrap iteration has been performed.}
#' \item{method}{Method used for inference.}
#' 
#' @seealso For an actual panel bootstrap procedure see \code{\link{sboot.pmb}}.
#' 
#' @references Pesaran, M. H., and Smith R. J. (1995):
#'   "Estimating Long-Run Relationships from Dynamic Heterogeneous Panels",
#'   \emph{Journal of Econometrics}, 68, pp. 79-113.
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
#' idx_nord = c("DNK", "FIN", "ISL", "SWE")
#' 
#' R.pvec = pvarx.VEC(L.data, lags=R.lags, dim_r=2, type="Case4")
#' R.pid  = pid.chol(R.pvec)
#' R.boot = sboot.mg(R.pid, idx_i=idx_nord)
#' plot(R.boot, lowerq=c(0, 0.25), upperq=c(1, 0.75))
#' summary(as.pvarx(R.pid$L.varx[idx_nord]))
#' 
#' # suppress imprecise results of restricted cointegrating coefficients #
#' dim_r = R.pvec$args_pvarx$dim_r
#' R.boot$beta$sim[ , 1:dim_r, ] = diag(dim_r)  # for normalized beta
#' summary(R.boot, idx_par="beta", level=0.95)
#' 
#' @export
#' 
sboot.mg <- function(x, n.ahead=20, normf=NULL, idx_i=NULL){
  # define
  if( is.null(idx_i) ){ idx_i = TRUE; w = NULL }else{ w = idx_i }
  A_MG = B_MG = beta = L.varx = dim_K = dim_S = NULL
  R.pvarx = aux_assign_pvarx(x, w=w)
  
  # collect individual IRF estimates and calculate their MG point estimates under 'w'
  IRF = irf.pvarx(R.pvarx, n.ahead=n.ahead, normf=normf, w=w)
  ipb = lapply(L.varx[idx_i], FUN=function(x) irf.varx(x, n.ahead=n.ahead, normf=normf))
  
  # return result
  result = list(true = IRF,
                bootstrap = ipb,
                A = list(par=A_MG, sim=R.pvarx$MG_A$coef),
                B = list(par=B_MG, sim=R.pvarx$MG_B$coef),
                beta = list(par=beta, sim=R.pvarx$MG_VECM$beta$coef),
                pvarx = R.pvarx,
                nboot = 0,
                rest_mat = matrix(NA, nrow=dim_K, ncol=dim_S),  # for 'bonferroni'
                design = NULL,
                method = "Mean group inference")
  class(result) = c("sboot2", "sboot")
  return(result)
}


#### S3 methods for objects of child class 'sboot2' ####
#' @export
print.sboot2 <- function(x, ...){
  # create header
  header_args = c("### 'sboot2' object from ", x$design, " ", x$method,  " ###")
  
  # print
  cat(header_args, "\n", sep="")
  cat("To get a distributional overview on the parameter matrices, \n")
  cat("select 'summary(object, idx_par='A')' for the VAR coefficients \n")
  cat("or 'summary(object, idx_par='B')' for the structural impact matrix. \n")
  cat("\n", paste0(c("Specifications", rep(".", times=20))), sep="", "\n")
  cat("Number of iterations: ", x$nboot, "\n")
  if(!is.null(x$b_length)){
    cat("Block size: ", paste(x$b_length, collapse=" x "), "\n") }
}


#' @export
summary.sboot2 <- function(object, ..., idx_par="A", level=0.95, digits=3){
  # define
  R.est = object[[idx_par]]$par  # point estimates
  R.sim = object[[idx_par]]$sim  # bootstrap results
  dim_K = nrow(R.est)
  if(is.null(R.sim)){ stop("The 'sboot' object does not provide results for parameter '", idx_par, "'.") }
  
  # names for variables and for all regressors
  names_k = if( !is.null(rownames(R.est)) ){ rownames(R.est) }else{ paste0("y.", 1:dim_K) }
  names_j = colnames(R.est)
  
  # calculate metrics of bootstrapped sampling distribution
  R.lower = apply(R.sim, MARGIN=1:2, FUN=function(x) quantile(x, probs=(1-level)/2))  # confidence intervals
  R.upper = apply(R.sim, MARGIN=1:2, FUN=function(x) quantile(x, probs=(1+level)/2))
  R.se    = apply(R.sim, MARGIN=1:2, FUN=function(x) sqrt(var(x)))  # standard errors
  R.ratio = R.est / R.se  # t-ratios
  
  # create and print header
  header_tab  = c("Estimate", "Std.error", "t-value", paste0("[", level*100,"%-conf."), "interval]")
  header_args = c("### Estimates of '", idx_par,"' with sampling measures from ", object$design, " ", object$method, " ###")
  cat(header_args, "\n", sep="")
  
  # create and print a table for each equation
  for(k in 1:dim_K){
    tab_k = cbind(R.est[k, ], R.se[k, ], R.ratio[k, ], R.lower[k, ], R.upper[k, ])
    tab_k = format(round(tab_k, digits=digits), nsmall=digits)
    dimnames(tab_k) = list(names_j, header_tab)
    
    cat("Parameters of equation '", names_k[k], "':", "\n", sep="")
    print(tab_k, quote=FALSE, ...)
    cat("\n")
  }
}


#' @export
toLatex.sboot2 <- function(object, ..., idx_par="A", measure=c("std.error", "t-value", "confint", NA), level=0.95, digits=3){
  # define
  R.est = object[[idx_par]]$par  # point estimates
  R.sim = object[[idx_par]]$sim  # bootstrap results
  dim_K = nrow(R.est)
  idx_m = unlist(lapply(measure, FUN=function(x) switch(x, "std.error"=1, "t-value"=2, "confint"=3:4, 0)))
  dim_M = sum(idx_m!=0)  # number of valid measures selected
  if(is.null(R.sim)){ stop("The 'sboot' object does not provide results for parameter '", idx_par, "'.") }
  
  # names for variables and for all regressors
  names_k = if( !is.null(rownames(R.est)) ){ rownames(R.est) }else{ paste0("y.", 1:dim_K) }
  names_k = paste0("\\texttt{", names_k, "}")
  names_j = paste0("\\texttt{", colnames(R.est), "}")
  
  # calculate measures of the bootstrapped sampling distributions
  R.lower = apply(R.sim, MARGIN=1:2, FUN=function(x) quantile(x, probs=(1-level)/2))  # confidence intervals
  R.upper = apply(R.sim, MARGIN=1:2, FUN=function(x) quantile(x, probs=(1+level)/2))
  R.se    = apply(R.sim, MARGIN=1:2, FUN=function(x) sqrt(var(x)))  # standard errors
  R.ratio = R.est / R.se  # t-ratios
  
  # create tabular
  est_tab = format(round(R.est, digits=digits), nsmall=digits)
  est_tab = apply(est_tab, MARGIN=1:2, FUN=function(x) paste0("$", x, "$"))
  est_tab = cbind(names_k, est_tab)
  if(dim_M){
    matrix_tab = NULL
    for(k in 1:dim_K){
      meas_k = rbind(R.se[k, ], R.ratio[k, ], R.lower[k, ], R.upper[k, ])
      idx_NaN = is.nan(meas_k)
      meas_k = format(round(meas_k, digits=digits), nsmall=digits)
      meas_k[idx_NaN] = " - "
      #### meas_k = rbind(meas_k[1:2, ], paste(meas_k[3, ], meas_k[4, ], sep=", "))
      meas_k[1:2,] = apply(meas_k[1:2,], MARGIN=1:2, FUN=function(x) paste0("${\\scriptstyle (", x, ")}$"))
      meas_k[3, ] = sapply(meas_k[3, ], FUN=function(x) paste0("${\\scriptstyle (", x, ";}$"))
      meas_k[4, ] = sapply(meas_k[4, ], FUN=function(x) paste0("${\\scriptstyle ", x, ")}$"))
      meas_k = cbind("", meas_k[idx_m, , drop=FALSE])
      matrix_tab = rbind(matrix_tab, est_tab[k, ], meas_k)
    }
  }else{ matrix_tab = est_tab }
  
  # return result
  tabular_core = apply(matrix_tab, MARGIN=1, FUN=function(x) paste0("\t", paste(x, collapse=" & "), " \\\\"))
  tabular_cols = paste0(c("l", rep("r", each=ncol(matrix_tab)-1)), collapse = "")
  tabular_rows = rep(c(rep("[-8pt]", each=dim_M), " \n"), times=dim_K)
  
  result = c(
    paste0("\\begin{tabular}{", tabular_cols, "}", sep=""),
    "\t\\hline \\hline",
    paste0("\t", paste(c("", names_j), collapse=" & "), " \\\\"),
    "\t\\hline",
    paste0(tabular_core, tabular_rows),
    "\t\\hline \\hline",
    "\\end{tabular}"
  )
  class(result) = "Latex"
  return(result)
}


