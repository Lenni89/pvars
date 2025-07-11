% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pvars.R
\docType{data}
\name{ERPT}
\alias{ERPT}
\title{Data set on the \emph{Exchange Rate Pass-Through}}
\format{
A long-format data panel of class '\code{data.frame}', 
  where the columns \code{id_i} and \code{id_t} 
  indicate the country and month respectively.
}
\source{
The prepared \emph{Eurostat} data set is directly obtainable from the 
  \emph{ZBW Journal Data Archive}: \doi{10.15456/jae.2022321.0717881037}. 
  This is open data under the \href{https://creativecommons.org/licenses/by/4.0/}{CC BY 4.0 license}.
}
\usage{
data("ERPT")
}
\description{
The data set \code{ERPT} consists of 
  monthly observations for the logarithm of import prices \eqn{lm^*}, 
  foreign prices \eqn{lf^*}, and the exchange rate against the US dollar \eqn{llcusd}.
  It covers the period January 1995 to March 2005 \eqn{(T=123)} for \eqn{N=7} countries. 
  The asterisk denotes the industry of the variables and can take values from 0 to 8:
  \itemize{
  \item{0:} \emph{Food and live animals chiefly for food}
  \item{1:} \emph{Beverages and tobacco}
  \item{2:} \emph{Crude materials (inedible, except fuels)}
  \item{3:} \emph{Mineral fuels, lubricants and related materials}
  \item{4:} \emph{Animal and vegetable oils, fats and waxes}
  \item{5:} \emph{Chemicals and related products}
  \item{6:} \emph{Manufactured goods classified chiefly by materials}
  \item{7:} \emph{Machines, transport equipment}
  \item{8:} \emph{Manufactured goods}
  }
}
\references{
Banerjee, A., and Carrion-i-Silvestre, J. L. (2015): 
  "Cointegration in Panel Data with Structural Breaks and Cross-Section Dependence", 
  \emph{Journal of Applied Econometrics}, 30 (1), pp. 1-23.
}
\seealso{
Other data sets: 
\code{\link{EURO}},
\code{\link{EU_w}},
\code{\link{ICAP}},
\code{\link{MDEM}},
\code{\link{MERM}},
\code{\link{PCAP}},
\code{\link{PCIT}}
}
\concept{data sets}
\keyword{datasets}
