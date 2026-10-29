

#' @title Forecast Error Variance Decomposition
#' @description Calculates the forecast error variance decomposition. Respects SVAR 
#'   models of cases \eqn{S \neq K}, i.e. partially identified or excess shocks, too.
#'
#' @param x SVAR object of class '\code{id}' or any other 
#'   that will be \link[=as.varx]{coerced} to ('\code{id}', '\code{varx}').
#' @param n.ahead Integer specifying the steps ahead, i.e. the horizon of the FEVD.
#' @param ... Currently not used.
#'
#' @return A list of class '\code{svarfevd}' holding the forecast error variance decomposition 
#'   of each variable as a '\code{data.frame}'.
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
#' # calculate and plot FEVD under partial identification #
#' plot(fevd(R.idBL, n.ahead=20))
#'
#' @import vars
#' @method fevd id
#' @export
#' 
fevd.id <- function(x, n.ahead=10, ...){
  # define
  x = as.varx(x)
  A = x$A
  B = x$B
  SIGMA = x$SIGMA
  dim_p = x$dim_p
  dim_K = nrow(A)
  dim_S = ncol(B)
  idx_h = 1:n.ahead
  
  # names for variables and for shocks
  names_k = if( !is.null(rownames(A)) ){ rownames(A) }else{ paste0("y[ ", 1:dim_K, " ]") }
  names_s = if( !is.null(colnames(B)) ){ colnames(B) }else{ paste0("epsilon[ ", 1:dim_S, " ]") }
  
  # calculate the FEVD, from Luetkepohl 2005:64, Eq.2.3.37 / Jentsch, Lunsford 2022:6, Ch.2.3
  R.vma  = aux_var2vma(A=A, B=B, dim_p=dim_p, n.ahead=n.ahead)
  df_k   = as.data.frame(matrix(NA, nrow=n.ahead, ncol=dim_S, dimnames=list(NULL, names_s)))
  result = sapply(names_k, FUN=function(k) df_k, simplify=FALSE)
  for(k in 1:dim_K){
    THETA_sq  = matrix(R.vma$THETA[k, , idx_h]^2, byrow=TRUE, nrow=n.ahead, ncol=dim_S)  # normalization cancels out, Eq.7
    fev_total = 0  # mean square error of variable 'k' at step 'j'
    for(j in idx_h){
      fev_shock = colSums(THETA_sq[1:j, , drop=FALSE])
      fev_total = fev_total + c(R.vma$PHI[k, ,j] %*% SIGMA %*% R.vma$PHI[k, ,j])
      ###fev_total = sum(THETA_sq[1:j, ])  # faster but ignores FEV from unidentified shocks
      result[[k]][j, ] = 100 * fev_shock / fev_total
    }
  }
  
  # return result
  class(result) = "svarfevd"
  return(result)
}



#### S3 methods for objects of class 'svarfevd' ####
#' @import ggplot2
#' @importFrom reshape2 melt
#' @method plot svarfevd
#' @export
#' 
plot.svarfevd <- function(x, ...){
  # define
  n.ahead = nrow(x[[1]])  # horizon
  dim_K   = length(x)     # number of variables
  dim_S   = ncol(x[[1]])  # number of shocks
  names_k = names(x)          # names of variables
  names_s = colnames(x[[1]])  # names of shocks
  
  # factors in data.frame preserve the stacking order in geom_bar() and facet_wrap()
  names_d = list(V1=NULL, shock=names_s, variable=paste0("FEVD~of~", names_k))
  ar_fevd = array(unlist(x), dim=c(n.ahead, dim_S, dim_K), dimnames=names_d)
  df_fevd = reshape2::melt(ar_fevd, varnames=names(names_d))  # Column vectors 'shock' and 'variable' are created as factors.
  
  # function to display integer steps for horizon
  integer_breaks <- function(n = 5, ...){
    fxn <- function(x){
      breaks <- floor(pretty(x, n, ...))
      names(breaks) <- attr(breaks, "labels")
      breaks
    }
    return(fxn)
  }
  
  # stfu R CMD check vs. ggplot2 (common practice, aes_ is deprecated)
  V1 = value = shock = NULL
  
  # plot the FEVD
  ggplot2::ggplot() + 
    geom_bar(data=df_fevd, aes(x=V1, y=value, fill=shock), stat="identity", position="stack") +
    facet_wrap(~variable, ncol=1, labeller=label_parsed) +  # parse like in svars:::plot.svarirf
    labs(x="Horizon", y="Contribution to FEV [in %]", fill="Shock") + 
    scale_fill_grey(labels=parse(text=names_s)) +  # Note: parse() sets a time stamp such that
    scale_x_continuous(breaks=integer_breaks()) +  # ... the ggplot objects cannot be all.equal().
    theme_bw() + 
    theme(legend.title=element_blank())
}


