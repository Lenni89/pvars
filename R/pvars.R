

#########################
### ROXYGEN NAMESPACE ###
#########################

#' @import svars
#' @importFrom stats .lm.fit coef cor cov var quantile optim resid residuals
#' @importFrom stats pchisq pf pgamma pnorm qnorm rnorm
#' @importFrom utils combn head tail toLatex
NULL 



########################
### ROXYGEN .Rd TAGS ###
########################

#' pvars: VAR Modeling for Heterogeneous Panels
#' @name pvars
#' @author Lennart Empting 
#' \email{lennart.empting@vwl.uni-due.de} 
#' (ORCID: \href{https://orcid.org/0009-0004-5068-4639}{0009-0004-5068-4639})
#' @description
#' This package implements (1) panel cointegration rank tests, (2) estimators for 
#' panel vector autoregressive (panel VAR) models, and (3) identification methods 
#' for panel structural vector autoregressive (panel SVAR) models as described 
#' in the accompanying vignette. 
#' The implemented functions allow users to account for cross-sectional dependence 
#' and for structural breaks in the deterministic terms of the VAR processes. 
#' @details
#' \bold{(1)} The panel functions to determine the cointegration rank are:
#' \itemize{
#' \item \code{\link{pcoint.JO}}   panel Johansen procedures,
#' \item \code{\link{pcoint.BR}}   panel test with pooled two-step estimation,
#' \item \code{\link{pcoint.SL}}   panel Saikkonen-Luetkepohl procedures,
#  \item \code{\link{pcoint.MSB}}  panel MSB procedures,
#' \item \code{\link{pcoint.CAIN}} correlation-augmented inverse normal test.
#' }
#' 
#' \bold{(2)} The panel functions to estimate the VAR models are:
#' \itemize{
#' \item \code{\link{pvarx.VAR}} mean-group of a panel of VAR models,
#' \item \code{\link{pvarx.VEC}} pooled cointegrating vectors in a panel VECM.
#' }
#' 
#' \bold{(3)} The panel functions to retrieve structural impact matrices are:
#' \itemize{
#' \item \code{\link{pid.chol}} identification of panel SVAR models using Cholesky decomposition to impose recursive causality,
#' \item \code{\link{pid.grt}}  identification of panel SVEC models by imposing long- and short-run restrictions,
#' \item \code{\link{pid.iv}}   identification of panel SVAR models by means of proxy variables,
#' \item \code{\link{pid.dc}}   independence-based identification of panel SVAR models using distance covariance (DC) statistic,
#' \item \code{\link{pid.cvm}}  independence-based identification of panel SVAR models using Cramer-von Mises (CVM) distance.
#' }
#' 
#' Supporting tools, such as the specification functions (\code{\link{speci.VAR}}, 
#' \code{\link{speci.factors}}) and the panel block bootstrap procedure 
#' (\code{\link{sboot.pmb}}), complement the panel VAR functions and complete 
#' this coherent approach to \strong{\emph{VAR modeling for heterogeneous panels}} 
#' within the \strong{vars} ecosystem. The provided data sets further allow for 
#' the exact replication of the implemented literature. 
"_PACKAGE"


#' Weights for the \emph{Euro Monetary Policy Transmission}
#' @docType data
#' @description The data set \code{EU_w} is a vector of 14 elements.
#'   These are weights for \eqn{N=14} member countries of the Euro area, 
#'   constructed as the average share of their respective real GDP 
#'   over the sample period in Herwartz, Wang (2024).
#' @usage data("EU_w")
#' @format A numeric vector containing 14 named elements.
#' @family data sets
#' @references Herwartz, H., and Wang, S. (2024): 
#'   "Statistical Identification in Panel Structural Vector Autoregressive 
#'   Models Based on Independence Criteria", 
#'   \emph{Journal of Applied Econometrics}, 39 (4), pp. 620-639.
#' @source The prepared \emph{Eurostat} data set is directly obtainable from the 
#'   \emph{ZBW Journal Data Archive}: \doi{10.15456/jae.2024044.1425287131}.
#'   This is open data under the \href{https://creativecommons.org/licenses/by/4.0/}{CC BY 4.0 license} 
#'   in accordance with the deposit license of the \emph{ZBW Journal Data Archive}.  
"EU_w"


#' Data set on the \emph{Euro Monetary Policy Transmission}
#' @docType data
#' @description The data set \code{EURO} is a list of 15 '\code{data.frame}' objects,
#'   each consisting of quarterly observations for 
#'   \itemize{
#'   \item the first-difference of log real GDP 
#'         on national \eqn{dl\_GDP} or aggregated EA-level \eqn{EA\_dl\_GDP}, 
#'   \item the annualized inflation of the (log) GDP deflator
#'         on national \eqn{dl\_deflator} or aggregated EA-level \eqn{EA\_pi}, 
#'   \item the EA-wide short-term interest rate \eqn{IR}, 
#'   \item the EA-wide option-adjusted bond spreads \eqn{BBB}, 
#'   \item the first-difference of log real GDP in the remaining countries \eqn{dl\_GDP\_EA}, 
#'   \item the weighted inflation in the remaining countries \eqn{dl\_deflator\_EA}, 
#'   \item the inflation of a world commodity price index \eqn{WCP}, 
#'   \item the US effective federal funds rate \eqn{US\_FFR}, 
#'   \item the trade volume in percentage of GDP \eqn{trade}, and 
#'   \item the government spending in percentage of GDP \eqn{ge}. 
#'   }
#'   The data covers the period Q1 2001 to Q1 2020 \eqn{(T=77)} for 
#'   the aggregate of the Euro area (EA, first element in list) and 
#'   \eqn{N=14} of its member countries (subsequent 14 elements in list).
#' @usage data("EURO")
#' @format A list-format data panel of class '\code{list}' 
#'   containing 15 '\code{data.frame}' objects with named time series.
#' @family data sets
#' @references Herwartz, H., and Wang, S. (2024): 
#'   "Statistical Identification in Panel Structural Vector Autoregressive 
#'   Models Based on Independence Criteria", 
#'   \emph{Journal of Applied Econometrics}, 39 (4), pp. 620-639.
#' @source The prepared \emph{Eurostat} data set is directly obtainable from the 
#'   \emph{ZBW Journal Data Archive}: \doi{10.15456/jae.2024044.1425287131}. 
#'   This is open data under the \href{https://creativecommons.org/licenses/by/4.0/}{CC BY 4.0 license} 
#'   in accordance with the deposit license of the \emph{ZBW Journal Data Archive}. 
"EURO"


#' Data set on the \emph{Exchange Rate Pass-Through}
#' @docType data
#' @description The data set \code{ERPT} consists of 
#'   monthly observations for the logarithm of import prices \eqn{lm^*}, 
#'   foreign prices \eqn{lf^*}, and the exchange rate against the US dollar \eqn{llcusd}.
#'   It covers the period January 1995 to March 2005 \eqn{(T=123)} for \eqn{N=7} countries. 
#'   The asterisk denotes the industry of the variables and can take values from 0 to 8:
#'   \itemize{
#'   \item{0:} \emph{Food and live animals chiefly for food}
#'   \item{1:} \emph{Beverages and tobacco}
#'   \item{2:} \emph{Crude materials (inedible, except fuels)}
#'   \item{3:} \emph{Mineral fuels, lubricants and related materials}
#'   \item{4:} \emph{Animal and vegetable oils, fats and waxes}
#'   \item{5:} \emph{Chemicals and related products}
#'   \item{6:} \emph{Manufactured goods classified chiefly by materials}
#'   \item{7:} \emph{Machines, transport equipment}
#'   \item{8:} \emph{Manufactured goods}
#'   }
#' @usage data("ERPT")
#' @format A long-format data panel of class '\code{data.frame}', 
#'   where the columns \code{id_i} and \code{id_t} 
#'   indicate the country and month respectively.
#' @family data sets
#' @references Banerjee, A., and Carrion-i-Silvestre, J. L. (2015): 
#'   "Cointegration in Panel Data with Structural Breaks and Cross-Section Dependence", 
#'   \emph{Journal of Applied Econometrics}, 30 (1), pp. 1-23.
#' @source The prepared \emph{Eurostat} data set is directly obtainable from the 
#'   \emph{ZBW Journal Data Archive}: \doi{10.15456/jae.2022321.0717881037}. 
#'   This is open data under the \href{https://creativecommons.org/licenses/by/4.0/}{CC BY 4.0 license}. 
"ERPT"


#' Data set on \emph{Infrastructure Capital Stocks}
#' @docType data
#' @description The data set \code{ICAP} consists of annual observations for
#'   \itemize{
#'   \item the real GDP \eqn{Y} in US dollars, 
#'   \item the physical capital stocks \eqn{K} in US dollars,
#'   \item the physical capital stocks \eqn{KBC30} via "backcasting" in US dollars, 
#'   \item the number of workers \eqn{LWDI} providing the total labor force, and 
#'   \item the average years of secondary education \eqn{secondary} of the population. 
#'   }
#'   It further reports physical measures of infrastructure given by
#'   \itemize{
#'   \item the electricity generation capacity \eqn{EGC} in megawatts, 
#'   \item the number of main phone lines \eqn{mlines}, 
#'   \item the kilometers of total roads \eqn{troads}, 
#'   \item the number of cell phone lines \eqn{cells}, 
#'   \item the kilometers of paved roads \eqn{proads}, and 
#'   \item the kilometers of rails \eqn{rails}.
#'   }
#'   It covers the period 1960 to 2000 \eqn{(T=41)} for \eqn{N=97} countries.
#'   The monetary values are given in US-Dollars at 2000 prices, i.e. constant PPP.
#' @usage data("ICAP")
#' @format A long-format data panel of class '\code{data.frame}', 
#'   where the columns \code{id_i} and \code{id_t} indicate the country and year 
#'   respectively. Column \code{COUNTRY} contains the complete country names.   
#' @family data sets
#' @references Calderon, C., Moral-Benito, E., and Serven, L. (2015): 
#'   "Is Infrastructure Capital Productive? A Dynamic Heterogeneous Approach", 
#'   \emph{Journal of Applied Econometrics}, 30 (2), pp. 177-198.
#' @source The prepared data set is directly obtainable from the 
#'   \emph{ZBW Journal Data Archive}: \doi{10.15456/jae.2022321.0717368489}. 
#'   This is open data under the \href{https://creativecommons.org/licenses/by/4.0/}{CC BY 4.0 license}.
"ICAP"


#' Data set for the \emph{Monetary Demand Model}
#' @docType data
#' @description The data set \code{MDEM} consists of 
#'   annual observations for the nominal short-term interest rate \eqn{R} and 
#'   the logarithm of the real money aggregate \eqn{m1} and real GDP \eqn{gdp}.
#'   It covers the period 1957 to 1996 \eqn{(T=40)} for \eqn{N=19} countries.
#' @usage data("MDEM")
#' @format A long-format data panel of class '\code{data.frame}', 
#'   where the columns \code{id_i} and \code{id_t} 
#'   indicate the country and year respectively.
#' @family data sets
#' @references Carrion-i-Silvestre, J. L., and Surdeanu, L. (2011): 
#'   "Panel Cointegration Rank Testing with Cross-Section Dependence", 
#'   \emph{Studies in Nonlinear Dynamics & Econometrics}, 15 (4), pp. 1-43.
#' @references Mark, N. C., and Sul, D. (1999): 
#'   "A Computationally Simple Cointegration Vector Estimator for Panel Data",
#'   Working Paper, Department of Economics, Ohio State University.
#' @references Mark, N. C., and Sul, D. (2003): 
#'   "Cointegration Vector Estimation by Panel DOLS and Long-Run Money Demand", 
#'   \emph{Oxford Bulletin of Economics and Statistics}, 65, pp. 655-680.
#' @source The prepared data is sourced from OECD and IMF's 
#'   \href{https://legacydata.imf.org/?sk=4c514d48-b6ba-49ed-8ab9-52b0c1a0179b}{\emph{International Financial Statistics}}
#'    of the year 1998, see the open 
#'   \href{https://www.imf.org/en/About/copyright-and-terms#data}{terms of use}. 
#'   Employed by Carrion-i-Silvestre and Surdeanu (2011:24, Ch.6.1), it was 
#'   originally compiled and described in the unpublished appendix of Mark and 
#'   Sul (2003). See the related working paper of Mark and Sul (1999, Appendix B).
"MDEM"


#' Data set for the \emph{Monetary Exchange Rate Model}
#' @docType data
#' @description The data set \code{MERM} consists of 
#'   monthly observations for the log-ratios of 
#'   prices \eqn{p}, money supply \eqn{m}, and industrial production \eqn{y} 
#'   as well as the natural logarithm of nominal exchange rates against the dollar \eqn{s}. 
#'   It covers the period January 1995 to December 2007 \eqn{(T=156)} for \eqn{N=19} countries.
#  # The data set in Mark, Sul (2001) covers the period Q1 1973 to Q1 1997 for N=19 countries.
#' @usage data("MERM")
#' @format A long-format data panel of class '\code{data.frame}', 
#'   where the columns \code{id_i} and \code{id_t} 
#'   indicate the country and month respectively.
#' @family data sets
#' @references 
#'   Oersal, D. D. K., and Arsova, A. (2017): 
#'   "Meta-Analytic Cointegrating Rank Tests for Dependent Panels", 
#'   \emph{Econometrics and Statistics}, 2, pp. 61-72.
#  @references Mark, N.C., and Sul, D. (2001): 
#    "Nominal exchange rates and monetary fundamentals: evidence from a small post-Bretton woods panel,"
#    \emph{Journal of International Economics}, 53, pp. 29-52.
#' @source The prepared data set is directly obtainable from the journal website: 
#'   \doi{10.1016/j.ecosta.2016.10.001}. Supplementary Raw Research Data. 
#'   This is open data under the \href{https://creativecommons.org/licenses/by/4.0/}{CC BY 4.0 license}.
"MERM"


#' Data set on \emph{Public Capital Stocks}
#' @docType data
#' @description The data set \code{PCAP} consists of annual observations for
#'   \itemize{
#'   \item the governmental capital stocks \eqn{G} and their logarithm \eqn{g}, 
#'   \item the private capital stocks \eqn{K} and their logarithm \eqn{k}, 
#'   \item the total hours worked \eqn{L} and their logarithm \eqn{l}, and 
#'   \item the real GDP \eqn{Y} and its logarithm \eqn{y}.
#'   }
#'   It is constructed from the annual observations for
#'   \itemize{
#'   \item the governmental investments \eqn{IG}, 
#'   \item the private non-residential investments and capital stocks \eqn{IB} and \eqn{B}, 
#'   \item the private housing investments and capital stocks \eqn{IH} and \eqn{H}, and 
#'   \item the persons employed \eqn{ET} and hours worked per person \eqn{HRS}.
#'   }
#'   It covers the period 1960 to 2019 \eqn{(T=60)} for \eqn{N=23} OECD countries.
#'   All monetary values are given in US-Dollars at 2005 prices, i.e. constant PPP.
#' @usage data("PCAP")
#' @format A long-format data panel of class '\code{data.frame}', 
#'   where the columns \code{id_i} and \code{id_t} 
#'   indicate the country and year respectively.
#' @family data sets
#' @references Empting, L. F. T., and Herwartz, H. (2025): 
#'   "Revisiting the 'Productivity of Public Capital': 
#'   VAR Evidence on the Heterogeneous Dynamics in a New Panel of 23 OECD Countries".
#' @references Feenstra, R. C., Inklaar, R., and Timmer, M. P. (2015):
#'   "The Next Generation of the Penn World Table", 
#'   \emph{American Economic Review}, 105, pp. 3150-3182.
#' @references Kamps, C. (2006):
#'   "New Estimates of Government Net Capital Stocks for 22 OECD Countries, 1960-2001",
#'   \emph{IMF Staff Papers}, 53, pp. 120-150.
#' @source Own compilation based on data from \emph{PWT}, \emph{Eurostat}, and OECD's \emph{Economic Outlook}. 
#'   Capital stocks are derived by the \emph{Perpetual Inventory Method} as described by Kamps (2006). 
#'   This is open data under the \href{https://creativecommons.org/licenses/by/4.0/}{CC BY 4.0 license}.
#' @example inst/examples/PCAP.R
"PCAP"


#' Data set on \emph{Personal and Corporate Income Tax}
#' @docType data
#' @description The data set \code{PCIT} consists of quarterly observations for
#'   \itemize{
#'   \item the average personal income tax rates \eqn{APITR}, 
#'   \item the average corporate income tax rates \eqn{ACITR}, 
#'   \item the logarithm of the personal income tax base \eqn{PITB}, 
#'   \item the logarithm of the corporate income tax base \eqn{CITB}, 
#'   \item the logarithm of government spending \eqn{GOV}, 
#'   \item the logarithm of GDP divided by population \eqn{RGDP}, and 
#'   \item the logarithm of government debt held by the public 
#'   divided by the GDP deflator and population \eqn{DEBT}.
#'   }
#'   Moreover, the proxies for shocks to personal \eqn{m\_PI} and corporate \eqn{m\_CI} 
#'   income taxes are prepended, where non-zero observations from the related 
#'   narratively identified shock series \eqn{T\_PI} resp. \eqn{T\_CI} have been demeaned. 
#'   The data set covers the period Q1 1950 to Q4 2006 \eqn{(T=228)} for the US.
#' @usage data("PCIT")
#' @format A time series data set of class '\code{data.frame}', 
#'   where the column \code{id_t} indicates the quarter of the year.
#' @family data sets
#' @references Mertens, K., and Ravn, M. O. (2013):
#'   "The Dynamic Effects of Personal and Corporate Income Tax Changes in the 
#'   United States", \emph{American Economic Review}, 103, pp. 1212-1247.
#' @references Jentsch, C., and Lunsford, K. G. (2019):
#'   "The Dynamic Effects of Personal and Corporate Income Tax Changes in the 
#'   United States: Comment", \emph{American Economic Review}, 109, pp. 2655-2678.
#' @references Mertens, K., and Ravn, M. O. (2019):
#'   "The Dynamic Effects of Personal and Corporate Income Tax Changes in the 
#'   United States: Reply", \emph{American Economic Review}, 109, pp. 2679-2691.
#' @source The prepared data set is directly obtainable from \emph{openICPSR}: 
#'   \doi{10.3886/E116190V1}. Supplementary Research Data. 
#'   This is open data under the \href{https://creativecommons.org/licenses/by/4.0/}{CC BY 4.0 license}.
"PCIT"

