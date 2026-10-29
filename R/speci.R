

#' @title Criteria on the number of common factors
#' @description Determines the number of factors in an approximate factor model
#'   for a data panel, where both dimensions \eqn{(T \times KN)} are large, 
#'   and calculates the factor time series and corresponding list of \eqn{N} idiosyncratic components.
#'   See Corona et al. (2017) for an overview and further details.
#'   
#' @param L.data List of \eqn{N} \code{data.frame} objects each collecting the \eqn{K_i} time series along the rows \eqn{t=1,\ldots,T}.
#'   The \eqn{\sum_{i=1}^{N} K_i = NK} time series are immediately combined into the \eqn{T \times KN} data panel \code{X}.
#' @param k_max Integer. The maximum number of factors to consider.
#' @param n.iterations Integer. Number of iterations for the Onatski criterion.
#' @param differenced Logical. If \code{TRUE}, each time series of the panel 
#'   \code{X} is first-differenced prior to any further transformation.
#'   Thereby, all criteria are calculated as outlined by Corona et al. (2017).
#' @param centered Logical. If \code{TRUE}, each time series of the panel \code{X} is centered.
#' @param scaled Logical. If \code{TRUE}, each time series of the panel \code{X} is scaled.
#'   Thereby, the PCA is applied via the correlation matrix instead of the covariance matrix of \code{X}.
#' @param n.factors Integer. A presumed number of factors under which the idiosyncratic component \code{L.idio} is calculated. 
#'   Deactivated if \code{NULL} (the default).
#' 
#' @return A list of class '\code{speci}', which contains the elements:
#' \item{eigenvals}{Data frame. The eigenvalues of the PCA, which have been used to calculate the criteria, 
#'      and their respective share of the total variance in the data panel.}
#' \item{Ahn}{Matrix. The eigenvalue ratio  \eqn{ER(k)} and growth rate \eqn{GR(k)} 
#'      by Ahn, Horenstein (2013) for \eqn{k=0,\ldots,}\code{k_max} factors.}
#' \item{Onatski}{Matrix. The calibrated threshold \eqn{\delta} and suggested number of factors \eqn{\hat{r}(\delta)} 
#'      by Onatski (2010) for each iteration.}
#' \item{Bai}{Array. The values of the criteria \eqn{PC(k)}, \eqn{IC(k)}, and \eqn{IPC(k)}
#'      with penalty weights \eqn{p1}, \eqn{p2}, and \eqn{p3} for \eqn{k=0,\ldots,}\code{k_max} factors.}
#' \item{selection}{List of the optimal number of common factors:
#'      (1) A matrix of \eqn{k^*} which minimizes each information criterion with each penalty weight. 
#'      (2) A vector of \eqn{k^*} which maximizes \code{ER} and \code{GR} respectively. 
#'      \code{ED} denotes the result by Onatski's (2010) "edge distribution" after convergence.}
#' \item{Ft}{Matrix. The common factors of dimension \eqn{(T \times} \code{n.factors}) estimated by PCA.}
#' \item{LAMBDA}{Matrix. The loadings of dimension \eqn{(KN \times} \code{n.factors}) estimated by OLS.}
#' \item{L.idio}{List of \eqn{N} \code{data.frame} objects each collecting 
#'       the \eqn{K_i} idiosyncratic series \eqn{\hat{e}_{it}} along the rows \eqn{t=1,\ldots,T}. 
#'       The series \eqn{\hat{e}_{it}} are given in levels and may contain a deterministic component with 
#'       (1) the initial \eqn{\hat{e}_{i1}} being non-zero and (2) re-accumulated means of the first-differenced series.}
#' \item{args_speci}{List of characters and integers indicating the specifications that have been used.}
#' 
#' @details If \code{differenced} is \code{TRUE}, the approximate factor model is estimated as proposed by Bai, Ng (2004).
#'   If all data transformations are selected, the estimation results are identical 
#'   to the objects in \code{$CSD} for PANIC analyses in '\code{pcoint}' objects.
#' 
#' @references Ahn, S., and Horenstein, A. (2013): 
#'   "Eigenvalue Ratio Test for the Number of Factors", 
#'   \emph{Econometrica}, 81, pp. 1203-1227.
#' @references Bai, J. (2004): 
#'   "Estimating Cross-Section Common Stochastic Trends in Nonstationary Panel Data", 
#'   \emph{Journal of Econometrics}, 122, pp. 137-183.
#' @references Bai, J., and Ng, S. (2002): 
#'   "Determining the Number of Factors in Approximate Factor Models", 
#'   \emph{Econometrica}, 70, pp. 191-221.
#' @references Bai, J., and Ng, S. (2004): 
#'   "A PANIC Attack on Unit Roots and Cointegration", 
#'   \emph{Econometrica}, 72, pp. 1127-1173.
#' @references Corona, F., Poncela, P., and Ruiz, E. (2017): 
#'   "Determining the Number of Factors after Stationary Univariate Transformations", 
#'   \emph{Empirical Economics}, 53, pp. 351-372.
#' @references Onatski, A. (2010): 
#'   "Determining the Number of Factors from Empirical Distribution of Eigenvalues", 
#'   \emph{Review of Economics and Statistics}, 92, pp. 1004-1016.
#' @examples
#' ### reproduce Oersal,Arsova 2017:67, Ch.5 ###
#' data("MERM")
#' names_k = colnames(MERM)[-(1:2)] # variable names
#' names_i = levels(MERM$id_i)      # country names
#' L.data  = sapply(names_i, FUN=function(i) 
#'    ts(MERM[MERM$id_i==i, names_k], start=c(1995, 1), frequency=12), 
#'    simplify=FALSE)
#' 
#' R.fac1 = speci.factors(L.data, k_max=20, n.iterations=4)
#' R.fac0 = speci.factors(L.data, k_max=20, n.iterations=4, 
#'    differenced=TRUE, centered=TRUE, scaled=TRUE, n.factors=8)
#'    
#' # scree plot #
#' library("ggplot2")
#' pal = c("#999999", RColorBrewer::brewer.pal(n=8, name="Spectral"))
#' lvl = levels(R.fac0$eigenvals$scree)
#' F.scree = ggplot(R.fac0$eigenvals[1:20, ]) +
#'   geom_col(aes(x=n, y=share, fill=scree), color="black", width=0.75) +
#'   scale_fill_manual(values=pal, breaks=lvl, guide="none") +
#'   labs(x="Component number", y="Share of total variance", title=NULL) +
#'   theme_bw()
#' plot(F.scree)
#' 
#' # factor plot (comp. Oersal,Arsova 2017:71, Fig.4) #
#' library("ggfortify")
#' Ft = ts(R.fac0$Ft, start=c(1995, 1), frequency=12)
#' F.factors = autoplot(Ft, facets=FALSE, size=1.5) + 
#'   scale_color_brewer(palette="Spectral") +
#'   labs(x="Year", y=NULL, color="Factor", title=NULL) +
#'   theme_bw()
#' plot(F.factors)
#' 
#' @family specification functions
#' @export
#' 
speci.factors <- function(L.data, k_max=20, n.iterations=4, differenced=FALSE, centered=FALSE, scaled=FALSE, n.factors=NULL){
  # define
  L.dim_K = sapply(L.data, FUN=function(x) ncol(x))
  L.dim_T = sapply(L.data, FUN=function(x) nrow(x))
  dim_NK = sum(L.dim_K)
  dim_T = min(L.dim_T)
  dim_N = length(L.data)
  idx_d = if(differenced){ 1:2 }else{ 1:3 }
  
  # data matrix of dimension (T_min x sum(K_i))
  X = do.call("cbind", lapply(L.data, FUN=function(x) tail(x, n=dim_T))) # first dim in tail() is T
  
  # data transformation and PCA
  if(differenced){ xit = diff(X) }else{ xit = X }
  xit = scale(xit, center=centered, scale=FALSE)
  xsd = scale(xit, center=FALSE, scale=scaled)
  R.svd = svd(xsd)
  evals = R.svd$d^2 / (dim_T-differenced)  # eigenvalues of the covariance matrix (X'X)/T
  
  # number of principal components
  dim_C2 = length(evals)  # min(dim(X))
  if(k_max > dim_C2){
    warning("k_max supasses the minimum dimension of the data matrix. 'k_max = 0.1*min(dim(X)' has been used instead")
    k_max = round(0.1*dim_C2)
  }
  
  # Onatski (2010) and Ahn,Horenstein (2013)
  R.ahc = aux_AHC(eigenvals=evals, r_max=k_max)
  R.onc = aux_ONC(eigenvals=evals, r_max=k_max, n.iterations=n.iterations)
  
  # Bai,Ng (2002) and Bai (2004)
  V = (sum(evals) - cumsum(c(0, evals[1:k_max]))) / dim_NK  # variances for k=0,...,k_max, from Corona et al. 2017:356, Eq.7
  R.bai = sapply(0:k_max, function(k) aux_PIC(X, k=k, Vk0=V[k+1], Vkmax0=V[k_max+1])[ ,idx_d], simplify="array")
  R.min = apply(R.bai, MARGIN=1:2, FUN=function(k) which.min(k)) - 1  # account for k=0
  
  # approximate factor model, see Bai,Ng 2002:198
  if(length(n.factors) > 0){
    evecs = R.svd$u[ , 0:n.factors , drop=FALSE]  # eigenvectors
    
    # estimate
    ft = sqrt(dim_T-differenced) * evecs  # common factors as (T-1 x n.factors) matrix
    LAMBDA = t(xit)%*%ft / (dim_T-differenced)  # loadings as (K*N x n.factors) matrix
    if(differenced){
      Ft = apply(rbind(0, ft),  MARGIN=2, FUN=function(x) cumsum(x)) # cumulate first-differences into levels, see Bai,Ng (2004)
    }else{ Ft = ft }
    et = X - Ft%*%t(LAMBDA)  # idiosyncratic components as (T x K*N) matrix
    
    # panel of idiosyncratic components in list-format
    idx_K  = cumsum(c(0, L.dim_K))
    L.idio = lapply(1:dim_N, FUN=function(i) et[ , idx_K[i] + 1:L.dim_K[i]])
    scree  = as.factor(c(1:n.factors, rep(0, dim_C2-n.factors)))
  }else{
    Ft     = NULL
    LAMBDA = NULL
    L.idio = NULL
    scree  = NA
  }
  
  # return result
  R.eval = data.frame(n = 1:dim_C2, 
                      scree = scree, 
                      value = evals, 
                      share = evals / sum(evals))  # share of explained variation
  select = list(R.min, c(R.ahc$selection, ED=R.onc$selection))
  argues = list(specifies="number of common factors", n.factors=n.factors, k_max=k_max,
                differenced=differenced, centered=centered, scaled=scaled)
  result = list(eigenvals=R.eval, Ahn=R.ahc$criteria, Onatski=R.onc$converge, Bai=R.bai,
                selection=select, Ft=Ft, LAMBDA=LAMBDA, L.idio=L.idio, args_speci=argues)
  class(result) = "speci"
  return(result)
}


#' @title Criteria on the lag-order and break period(s)
#' @description Determines the lag-order \eqn{p} and break period(s) \eqn{\tau} 
#'   jointly via information criteria on the OLS-estimated VAR model for a given 
#'   number of breaks. These \eqn{m} breaks are common to all \eqn{K} equations 
#'   of the system and partial, pertaining to the 
#'   \link[=as.t_D]{deterministic term} only.
#' 
#' @param x VAR object of class '\code{varx}' or any other 
#'   that will be \link[=as.varx]{coerced} to '\code{varx}'. Specifically for 
#'   \strong{vars}' \code{\link[vars]{VAR}}, use \code{p = min(lag_set)} 
#'   or simply \code{p=1} such that the customized \code{$D} from the coerced 
#'   '\code{varx}' object contains no \code{NA} in the effective sample.
#' @param lag_set Vector. Set of candidates for the lag-order \eqn{p}. 
#'   If only a single integer is provided, the criteria just reflect 
#'   the variation of det\eqn{(\hat{U}_{\tau} \hat{U}_{\tau}')} uniformly and 
#'   determine the break period(s) \eqn{\tau} unanimously as \eqn{\hat{\tau} = }
#'   arg min det\eqn{(\hat{U}_{\tau} \hat{U}_{\tau}')} under the given \eqn{p}.
#' @param dim_m Integer. Number of breaks in the deterministic term to consider.
#'   If \code{FALSE} (the default), the criteria determine only 
#'   the lag-order \eqn{p} just like \strong{vars}' \code{\link[vars]{VARselect}}.
#' @param trim Numeric. Either a numeric value \eqn{h \in (p_{max}/T, 1/m)} that 
#'   defines the minimal fraction relative to the total sample size \eqn{T} or 
#'   an integer that defines the minimal number of observations in each sub-sample. 
#'   For example, \eqn{h=0.15} (the default) specifies the window 
#'   \eqn{[0.15 \cdot T, 0.85 \cdot T]} that is often used 
#'   as the set of candidates for \eqn{m=1} single period \eqn{\tau_1}. 
#' @param type_break Character. Whether the \eqn{m} common breaks pertain to the 
#'   '\code{const}' (the default), the linear '\code{trend}', or '\code{both}'. 
#'   Adds the period-specific \link[=as.t_D]{deterministic term} activated 
#'   during \eqn{\tau}.
#' @param add_dummy Logical. If \code{TRUE} (not the default), accompanying 
#'   impulse dummies activated in \eqn{\tau + (0, \ldots, p-1)} are added to each break.
#' @param n.cores Integer. Number of allocated processor cores.
#'   Note that parallel processing is activated only for the combined 
#'   determination of lag-order \eqn{p} and break period(s) \eqn{\tau}.
#' 
#' @return A list of class '\code{speci}', which contains the elements: 
#' \item{df}{A '\code{data.frame}' of \eqn{(1+m) + 4} columns for all admissible 
#'   combinations of candidate \eqn{(p, \tau)} and their values of 
#'   \eqn{AIC(p, \tau)}, \eqn{HQC(p, \tau)}, \eqn{SIC(p, \tau)}, and \eqn{FPE(p, \tau)}.}
#' \item{selection}{A \eqn{(1+m) \times 4} matrix of the specification pairs 
#'   \eqn{(p^*, \tau^*)} suggested by the global minimum of the AIC (Akaike 1969), 
#'   HQC (Hannan, Quinn 1979), SIC (Schwarz 1978), and FPE respectively.}
#' \item{args_speci}{List of characters and integers indicating the specifications that have been used.}
#' 
#' @details The literature on structural breaks in time series deals mostly with 
#'   the determination of the number \eqn{m} and position \eqn{\tau} of breaks 
#'   (e.g. Bai, Perron 1998 and 2003), but leaves the lag-order \eqn{p} aside. 
#'   For example, under a given \eqn{p}, Luetkepohl et al. (2004) use a full-rank 
#'   VAR in levels to determine \eqn{m=1} common break period \eqn{\tau_1} 
#'   and subsequently perform cointegration analysis with \code{\link{coint.SL}} 
#'   (which actually provides \eqn{p}-values for up to \eqn{m=2}). 
#'   Note, however, that the lag-order of a VECM is likewise usually determined via 
#'   information criteria of a full-rank VAR in levels.
#'   
#'   \code{\link{speci.VAR}} combines Bai, Perron (2003) and Approach 3 of Yang (2002)
#'   into a global minimization of information criteria on the pair \eqn{(p,\tau)}. 
#'   Specifically, Yang (2002:378, Ch.2.2) estimates all candidate VAR models by 
#'   OLS and then determines their optimal lag-order \eqn{p^*} and \eqn{m=1} break 
#'   period \eqn{\tau^*} jointly via the global minimum of the information criteria. 
#'   Bai and Perron (2003, Ch.3) determine 
#'   \eqn{\tau^* = (\tau_1^*, \ldots, \tau_m^*)} of multiple breaks via the 
#'   minimum sum of squared residuals from a single-equation model \eqn{(K=1)}. 
#'   They use dynamic programming to reduce the number of least-squares operations. 
#'   Although adapting their streamlined set of admissible combinations for \eqn{\tau}, 
#'   \code{\link{speci.VAR}} still resorts to (parallelized brute-force) OLS estimation 
#'   of all candidate VAR models and therewith circumvents issues of correct 
#'   initialization and iterative updating for the models with partial breaks.
#'   
#' @references Bai, J., and Perron, P. (1998): 
#'   "Estimating and Testing Linear Models with Multiple Structural Changes", 
#'   \emph{Econometrica}, 66, pp. 47-78.
#' @references Bai, J., and Perron, P. (2003): 
#'   "Computation and Analysis of Multiple Structural Change Models", 
#'   \emph{Journal of Applied Econometrics}, 18, pp. 1-22.
#' @references Luetkepohl, H., Saikkonen, P., and Trenkler, C. (2004): 
#'   "Testing for the Cointegrating Rank of a VAR Process with Level Shift at Unknown Time", 
#'   \emph{Econometrica}, 72, pp. 647-662.
#' @references Yang, M. (2002): 
#'   "Lag Length and Mean Break in Stationary VAR Models", 
#'   \emph{Econometrics Journal}, 5, pp. 374-386.
#' @examples
#' ### extend basic example in "urca" ###
#' library("urca")
#' library("vars")
#' data("denmark")
#' sjd = denmark[, c("LRM", "LRY", "IBO", "IDE")]
#' 
#' # use the single lag-order p=2 to determine only the break period #
#' R.vars  = VAR(sjd, type="both", p=1, season=4)
#' R.speci = speci.VAR(R.vars, lag_set=2, dim_m=1, trim=3, add_dummy=FALSE)
#' 
#' library("ggfortify")
#' autoplot(ts(R.speci$df[3:5], start=1+R.speci$args_speci$trim), 
#'  main="For a single 'p', all IC just reflect the variation of det(UU').")
#' print(R.speci)
#' 
#' # perform cointegration test procedure with detrending #
#' R.t_D   = list(t_shift=8, n.season=4)
#' R.coint = coint.SL(sjd, dim_p=2, type_SL="SL_trend", t_D=R.t_D)
#' summary(R.coint)
#' 
#' # m=1: line plot #
#' library("ggplot2")
#' R.speci1 = speci.VAR(R.vars, lag_set=1:5, dim_m=1, trim=6)
#' R.values = c("#BDD7E7", "#6BAED6", "#3182BD", "#08519C", "#08306B")
#' F.line   = ggplot(R.speci1$df) +
#'   geom_line( aes(x=tau_1, y=HQC, color=as.factor(p), group=as.factor(p))) +
#'   geom_point(aes(x=tau_1, y=HQC, color=as.factor(p), group=as.factor(p))) +
#'   geom_point(x=R.speci1$selection["tau_1", "HQC"], 
#'              y=min(R.speci1$df$HQC), color="red") +
#'   scale_x_continuous(limits=c(1, nrow(sjd))) +
#'   scale_color_manual(values=R.values) +
#'   labs(x=expression(tau), y="HQ Criterion", color="Lag order", title=NULL) +
#'   theme_bw()
#' plot(F.line)
#' 
#' # m=2: heat map on integer lattice #
#' R.speci2 = speci.VAR(R.vars, lag_set=2, dim_m=2, trim=3)
#' dim_T    = nrow(sjd)  # total sample size
#' F.heat   = ggplot(R.speci2$df) +
#'   geom_point(aes(x=tau_1, y=tau_2, color=AIC), size=3) +
#'   geom_abline(intercept=0, slope=-1, color="grey") +
#'   scale_x_continuous(limits=c(1, dim_T), expand=c(0, 0)) +
#'   scale_y_reverse(limits=c(dim_T, 1), expand=c(0, 0)) +
#'   scale_color_continuous(type="viridis") +
#'   labs(x=expression(tau[1]), y=expression(tau[2]), color="AIC", title=NULL) +
#'   theme_bw()
#' plot(F.heat)
#' 
#' @family specification functions
#' @export
#' 
speci.VAR <- function(x, lag_set=1:10, dim_m=FALSE, trim=0.15, type_break="const", add_dummy=FALSE, n.cores=1){
  # define
  x = as.varx(x)
  lag_max  = max(lag_set)  # maximum lag-order p_max in the set of candidates
  names_m  = paste0("tau_", 0:dim_m)[-1]
  names_ic = c("AIC", "HQC", "SIC", "FPE")
  dim_T  = ncol(x$y)  # total sample size
  idx_ic = (1 + dim_m) + 1:length(names_ic)
  idx_pm = 1 + 0:dim_m  # index for model specifications (including lag-orders)
  idx_na = any(is.na(x$D[ , min(lag_set)+1]))
  if(idx_na){ stop("The provided VAR model 'x' has missing values in the presample of its deterministic regressors.") }
  
  # define admissible combinations of tau for grid search
  if(dim_m > 0){
    # check break type
    if(type_break %in% c("const", "both") & !(x$args_varx$type %in% c("const", "both"))){
      warning("'type_break' pertains the constant, which is missing in the provided VAR model 'x'.") }
    if(type_break %in% c("trend", "both") & !(x$args_varx$type %in% c("trend", "both"))){
      warning("'type_break' pertains the linear trend, which is missing in the provided VAR model 'x'.") }
    
    # check trim parameter
    if(trim < 1){
      # relative trimming parameter for candidate break periods w.r.t.
      # ... the total sample by truncation (Luetkepohl et al. 2004:649, Eq.2.2)
      # ... the effective sample by the floor (Yang 2002:376, 378)
      # ... the imposition of the minimal sub-sample size (Bai, Perron 2003:12) [CHOSEN]
      trim = max(trunc(trim * dim_T), 1) }
    if(trim <= lag_max+1 & type_break %in% c("trend", "both")){
      stop("Minimum size of sub-samples must be larger than max(lag_set)+1 to accommodate trend breaks.") }
    if(trim <= lag_max){
      stop("Minimum size of sub-samples must be larger than maximum lag-order.") }
    if(trim*(dim_m+1) > dim_T){
      stop("Minimum sizes of the m+1 sub-samples must, in total, not exceed sample size T.") }
    
    # define admissible combinations of tau, from Bai, Perron 2003:4, Ch.3.1
    all_combos = combn(1:dim_T, m=dim_m, FUN=NULL)  # all combinations
    idx_admiss = apply(all_combos, MARGIN=2, FUN=function(combo){
      # ordered first periods of the m+1 sub-samples and the upper bound
      t_starts = c(1, sort(combo), dim_T+1)
      # each sub-sample must contain at least 'trim' integers
      all(diff(t_starts) >= trim) })
    TAUS = t(all_combos[, idx_admiss, drop = FALSE])  # admissible combinations
    colnames(TAUS) = names_m
  }else{
    TAUS = NULL
  }
  
  # function for ICs over break periods for given lag-order 'p'
  specif <- function(p){
    # define under given lag-order
    result = cbind(p=p, TAUS, AIC=NA, HQC=NA, SIC=NA, FPE=NA)  # initialize
    idx_t  = (lag_max+1):dim_T    # same effective sample for all 'p', from Yang 2002:376
    idx_tp = (lag_max+1-p):dim_T  # pre-sample for VAR regressors increases with 'p' (other regressors will match via aux_stack)
    Y      = x$y[ , idx_t, drop=FALSE]
    Z_orig = aux_stack(y=x$y[ , idx_tp, drop=FALSE], dim_p=p, x=x$x, dim_q=x$dim_q, D=x$D)
    
    # calculate and collect information criteria
    for(j in 1:nrow(result)){
      # period-specific regressors, from Luetkepohl et al. 2004:650, Eq.3.1 resp. Eq.3.3
      if(dim_m > 0){
        D = aux_dummy(dim_T=dim_T,
          t_impulse = if(add_dummy){ sapply(TAUS[j, ], FUN=function(tau_j) tau_j + 0:(p-1)) }else{ NULL },
          t_shift   = switch(type_break, "const"=TAUS[j, ], "trend"=NULL, "both"=TAUS[j, ]),
          t_break   = switch(type_break, "const"=NULL, "trend"=TAUS[j, ], "both"=TAUS[j, ]))
        D = D[ , idx_t, drop=FALSE]
      }else{ D = NULL }
      
      # estimate by OLS
      Z = rbind(Z_orig, D)
      A = tcrossprod(Y, Z) %*% solve(tcrossprod(Z))
      E = Y - A %*% Z
      
      # calculate information criteria
      OMEGA = tcrossprod(E) / ncol(E)  # MLE residual covariance matrix
      result[j, idx_ic] = aux_MIC(Ucov=OMEGA, COEF=A, dim_T=ncol(E))
      ### Yang (2002:378, Eq.4) does not include impulse dummies into the estimated model.
      ### Yang 2002:377 treats tau as an estimated parameter and thus uses n.coef+1 in the penalty weights, 
      ### ... but these uniform adjustments do not change the minimizing order anyway (Luetkepohl 2005:147).
      ### Luetkepohl et al. (2004:651, Eq.3.4) use det(UU') only.
    }
    
    # return result
    return(result)
  }
  
  # run procedure of joint IC
  if(length(lag_set) == 1){
    spec = "break period(s) 'tau'"
    df   = specif(lag_set)
  }else if(dim_m == 0){
    spec = "lag-order 'p'"
    L.df = lapply(lag_set, FUN=function(p) specif(p))
    df   = do.call("rbind", L.df)
  }else{
    spec = "lag-order 'p' and break period(s) 'tau'"
    L.df = pbapply::pblapply(lag_set, FUN=function(p) specif(p), cl=n.cores)
    df   = do.call("rbind", L.df)
  }
  
  # return result
  select = sapply(names_ic, FUN=function(ic) df[which.min(df[, ic]), idx_pm])  # global minimum of each IC
  select = matrix(select, ncol=4, dimnames=list(c("p", names_m), names_ic))  # enforce matrix even if m=0.
  argues = list(specifies=spec, trim=trim, lag_set=lag_set)
  result = list(df=as.data.frame(df), selection=select, args_speci=argues)
  class(result) = "speci"
  return(result)
}


#### S3 methods for objects of class 'speci' ####
#' @export
print.speci <- function(x, ...){
  # create table
  header_args = c("### Optimal ", x$args_speci$specifies, " ###")
  
  # print
  cat(header_args, "\n", sep="")
  print(x$selection, quote=FALSE, row.names=FALSE)
}


