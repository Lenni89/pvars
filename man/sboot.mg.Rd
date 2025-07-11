% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sboot.R
\name{sboot.mg}
\alias{sboot.mg}
\title{Mean group inference for panel SVAR models}
\usage{
sboot.mg(x, n.ahead = 20, normf = NULL, idx_i = NULL)
}
\arguments{
\item{x}{Panel VAR object of class '\code{pid}' or '\code{pvarx}' 
or a list of VAR objects that will be \link[=as.varx]{coerced} to '\code{varx}'.}

\item{n.ahead}{Integer. Number of periods to consider after the initial impulse, i.e. the horizon of the IRF.}

\item{normf}{Function. A given function that normalizes the \eqn{K \times S} input-matrix 
into an output matrix of same dimension. See the example in \code{\link{id.iv}} 
for the normalization of Jentsch and Lunsford (2021) 
that fixes the size of the impact response in the IRF.}

\item{idx_i}{Logical or character vector. 
Names or \eqn{N} logical elements selecting a subset from the 
individuals \eqn{i = 1, \ldots, N} for the MG estimation. If \code{NULL} 
(the default), all \eqn{N} individuals are included.}
}
\value{
A list of class '\code{sboot2}' with elements:
\item{true}{Mean group estimate of impulse response functions.}
\item{bootstrap}{List of length \eqn{N} holding the individual impulse response functions.}
\item{A}{List for the VAR coefficients containing 
  the matrix of mean group estimates '\code{par}' and 
  the array of individual results '\code{sim}'.}
\item{B}{List for the structural impact matrix containing 
  the matrix of mean group estimates '\code{par}' and 
  the array of individual results '\code{sim}'.}
\item{pvarx}{Input panel VAR object of class '\code{pvarx}'.}
\item{nboot}{Integer '0' indicating that no bootstrap iteration has been performed.}
\item{method}{Method used for inference.}
}
\description{
Calculates confidence bands for impulse response functions via mean group inference. 
  The function does not perform bootstraps, but coerces the panel VAR object to class '\code{sboot}'
  and, therewith, gives a distributional overview on the parameter heterogeneity.
}
\details{
MG inference presumes the individual estimates to be the empirical variation 
  around a common parameter. In case of heterogeneous lag-orders \eqn{p_i},
  specifically the '\code{summary}' of VAR coefficient matrices fills 
  \eqn{\hat{A}_{ij} = 0_{K \times K}} for \eqn{p_i < j \le max(p_1,\ldots,p_N)} 
  in accordance with the finite order VAR\eqn{(p_i)}.
}
\examples{
data("PCAP")
names_k = c("g", "k", "l", "y")  # variable names
names_i = levels(PCAP$id_i)      # country names 
L.data  = sapply(names_i, FUN=function(i) 
  ts(PCAP[PCAP$id_i==i, names_k], start=1960, end=2019, frequency=1), 
  simplify=FALSE)
R.lags = c(2, 4, 2, 3, 2, 4, 4, 2, 2, 3, 3, 3, 2, 4, 4, 2, 2, 2, 4, 2, 2, 2, 4)
names(R.lags) = names_i
idx_nord = c("DNK", "FIN", "ISL", "SWE")

R.pvec = pvarx.VEC(L.data, lags=R.lags, dim_r=2, type="Case4")
R.pid  = pid.chol(R.pvec)
R.boot = sboot.mg(R.pid, idx_i=idx_nord)
plot(R.boot, lowerq=c(0, 0.25), upperq=c(1, 0.75))
summary(as.pvarx(R.pid$L.varx[idx_nord]))

# suppress imprecise results of restricted cointegrating coefficients #
dim_r = R.pvec$args_pvarx$dim_r
R.boot$beta$sim[ , 1:dim_r, ] = diag(dim_r)  # for normalized beta
summary(R.boot, idx_par="beta", level=0.95)

}
\references{
Pesaran, M. H., and Smith R. J. (1995):
  "Estimating Long-Run Relationships from Dynamic Heterogeneous Panels",
  \emph{Journal of Econometrics}, 68, pp. 79-113.
}
\seealso{
For an actual panel bootstrap procedure see \code{\link{sboot.pmb}}.
}
