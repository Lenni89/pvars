### gallery of merged IRF plots ###
\donttest{
library("ggplot2")
data("PCAP")
names_k = c("g", "k", "l", "y")  # variable names
names_i = levels(PCAP$id_i)      # country names 
L.data  = sapply(names_i, FUN=function(i) 
  ts(PCAP[PCAP$id_i==i, names_k], start=1960, end=2019, frequency=1), 
  simplify=FALSE)
L.vars = lapply(L.data, FUN=function(x) vars::VAR(x, p=2, type="both"))
L.chol = lapply(L.vars, FUN=function(x) svars::id.chol(x))

# overlay all IRF to get an overview on the stability #
L.irf = lapply(L.chol, FUN=function(x) plot(irf(x, n.ahead=30)))
summary(as.pvarx(L.vars))
as.pplot(L.irf)

# overlay IRF of selected countries and quantiles of all countries #
F.mg  = plot(sboot.mg(L.chol, n.ahead=30), lowerq=0.05, upperq=0.95)
R.irf = as.pplot(MG=F.mg, L.irf[c("DEU", "FRA", "ITA", "JPN")])
plot(R.irf)  # emphasize MG-IRF in next step
R.irf = as.pplot(R.irf, color_g="black", shape_g=c(20, rep(NA, 4)))
R.irf$F.plot + labs(color="Country")

# compare two mean-groups and their quantiles #
idx_nord = c(5, 6, 10, 17, 20)  # Nordic countries
R.irf = as.pplot(color_g=c("black", "blue"), 
  Other  = plot(sboot.mg(L.chol[-idx_nord])), 
  Nordic = plot(sboot.mg(L.chol[ idx_nord])),
  width_L=0.2, color_L="black")
R.irf$F.plot + guides(fill=guide_legend(title="Group"))
 
# compare different shock ordering for MG-IRF #
R.pid1 = pid.chol(L.vars)
R.pid2 = pid.chol(L.vars, order_k=4:1)
R.pid3 = pid.chol(L.vars, order_k=c(1,4,2,3))

R.pal = RColorBrewer::brewer.pal(n=8, name="Spectral")[c(8, 1, 4)]
R.irf = as.pplot(color_g=R.pal, shape_g=c(2, 3, 20), 
  GKLY = plot(irf(R.pid1, n.ahead=25)), 
  YLKG = plot(irf(R.pid2, n.ahead=25)), 
  GYKL = plot(irf(R.pid3, n.ahead=25)))
plot(R.irf)

R.mg = as.pplot(color_g=R.pal, shape_g=c(2, 3, 20), width_L=0.3, 
  GKLY = plot(sboot.mg(R.pid1, n.ahead=25), lowerq=0.05, upperq=0.95), 
  YLKG = plot(sboot.mg(R.pid2, n.ahead=25), lowerq=0.05, upperq=0.95), 
  GYKL = plot(sboot.mg(R.pid3, n.ahead=25), lowerq=0.05, upperq=0.95))
R.mg$F.plot + guides(fill=guide_legend()) + 
  labs(color="Ordering", fill="Ordering")
         
# colorize and export a single sub-plot to Latex #
library("tikzDevice")
textwidth = 15.5/2.54  # LaTeX textwidth from "cm" into "inch"
file_fig  = file.path(tempdir(), "Fig_irf.tex")

R.irf = as.pplot(
  DEU = plot(irf(L.chol[["DEU"]], n.ahead=50), selection=list(4, 1)), 
  FRA = plot(irf(L.chol[["FRA"]], n.ahead=50), selection=list(4, 1)),
  color_g = c("black", "blue"),
  names_g = c("Germany", "France"),
  names_k = "y", 
  names_s = "\\epsilon_{ g }", 
  Latex   = TRUE)

tikz(file=file_fig, width=1.2*textwidth, height=0.8*textwidth)
  R.irf$F.plot + labs(color="Country") + theme_minimal()
dev.off()
}

