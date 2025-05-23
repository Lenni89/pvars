% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pvars.R
\docType{data}
\name{MDEM}
\alias{MDEM}
\title{Data set for the \emph{Monetary Demand Model}}
\format{
A long-format data panel of class '\code{data.frame}', 
  where the columns \code{id_i} and \code{id_t} 
  indicate the country and year respectively.
}
\source{
\url{https://www.degruyter.com/view/journals/snde/15/4/article-1558-3708.1825.xml.xml}
}
\usage{
data("MDEM")
}
\description{
The data set \code{MDEM} consists of 
  annual observations for the nominal short-term interest rate \eqn{R} and 
  the logarithm of the real money aggregate \eqn{m1} and real GDP \eqn{gdp}.
  It covers the period 1957 to 1996 \eqn{(T=40)} for \eqn{N=19} countries.
}
\references{
Carrion-i-Silvestre, J. L., and Surdeanu L. (2011): 
  "Panel Cointegration Rank Testing with Cross-Section Dependence", 
  \emph{Studies in Nonlinear Dynamics & Econometrics}, 15 (4), pp. 1-43.

Mark, N. C., and Sul, D. (2003): 
  "Cointegration Vector Estimation by Panel DOLS and Long-Run Money Demand," 
  \emph{Oxford Bulletin of Economics and Statistics}, 65, pp. 655-680.
}
\seealso{
Other data sets: 
\code{\link{ERPT}},
\code{\link{EURO}},
\code{\link{EU_w}},
\code{\link{MERM}},
\code{\link{PCAP}},
\code{\link{PCIT}}
}
\concept{data sets}
\keyword{datasets}
