

############
# SETTINGS #
############

# default path for exporting figures and tables #
path_fig = "~/vignettes/images/"
path_tab = "~/vignettes/tables/"

library("tikzDevice")
textwidth = 15.5/2.54  # LaTeX textwidth from "cm" into "inch"
options(tikzMetricsDictionary = paste0(path_fig, "tikzMetricsDictionary"))  # create MetricsDictionary



############
# 4.1 MERM #
############

library("pvars")
library("ggplot2")
library("ggfortify")
data("MERM")
names_k = colnames(MERM)[-(1:2)] # variable names
names_i = levels(MERM$id_i)      # country names
L.data  = sapply(names_i, FUN=function(i) ts(MERM[MERM$id_i==i, names_k], start=c(1995, 1), frequency=12), simplify=FALSE)

head(MERM, n=3)
# L.data = lapply(names_i, FUN=function(i) MERM[MERM$id_i==i, names_k]) 
# names(L.data) = names_i

# determine model specifications #
R.fac1 = speci.factors(L.data, k_max=20, n.iterations=4)
R.fac0 = speci.factors(L.data, k_max=20, n.iterations=4, differenced=TRUE, centered=TRUE, scaled=TRUE, n.factors=8)
R.lags = sapply(R.fac0$L.idio, FUN=function(i) vars::VARselect(i, lag.max=4, type="both")$selection[1])
R.lags = c(2, 2, 2, 2, 1, 2, 2, 4, 2, 3, 2, 2, 2, 2, 2, 1, 1, 2, 2); names(R.lags)=names_i  # individual lags by AIC (lag_max=4)

# scree plot #
pal = c("#999999", RColorBrewer::brewer.pal(n=8, name="Spectral"))
lvl = levels(R.fac0$eigenvals$scree)
F.scree = ggplot(R.fac0$eigenvals[1:20, ]) +
  geom_col(aes(x=n, y=share, fill=scree), color="black", width=0.75) +
  scale_fill_manual(values=pal, breaks=lvl, guide="none") +
  labs(x="Component number", y="Share of total variance", title=NULL) +
  theme_bw()

tikz(file=paste0(path_fig, "Fig_Scree.tex"), width=textwidth, height=0.35*textwidth)
  plot(F.scree)
dev.off()

# Arsova,Oersal 2017:67, Tab.5 #
R.pcsl = pcoint.SL(L.data, lags=R.lags, type="SL_trend", n.factors=8)  # eight factors by Onatski criterion
R.pcjo = pcoint.JO(L.data, lags=R.lags, type="Case4",    n.factors=8)

sink(file=paste0(path_tab, "Tab_MERM.tex"))
  toLatex(R.pcsl)
sink(file=NULL)

# Arsova,Oersal 2017:71, Fig.4 #
Ft   = ts(R.pcsl$CSD$Ft, start=c(1995, 1), frequency=12)
F.Ft = autoplot(Ft, facets=FALSE, size=1.75) + theme_bw() + 
  scale_color_brewer(palette="Spectral") +
  labs(x=NULL, y=NULL, color="Factor", title=NULL)

tikz(file=paste0(path_fig, "Fig_Factors.tex"), width=textwidth, height=0.4*textwidth)
  plot(F.Ft)
dev.off()

# Arsova,Oersal 2017:67, Tab.6 #
vars::VARselect(Ft, lag.max=4, type="both")$selection
speci.VAR(vars::VAR(Ft, p=1, type="both"), lag_set=1:4, type="both")
R.Ftsl = coint.SL(y=R.pcsl$CSD$Ft, dim_p=2, type_SL="SL_trend")  # lag-order by AIC
R.Ftjo = coint.JO(y=R.pcsl$CSD$Ft, dim_p=2, type="Case4")

sink(file=paste0(path_tab, "Tab_MERMft.tex"))
  toLatex(R.Ftsl, R.Ftjo, write_ME=TRUE, add2header=c("\\textbf{SL:} Trace test", "\\textbf{Johansen:} Trace test"))
sink(file=NULL)



############
# 4.2 ERPT #
############

library("pvars")
data("ERPT")
names_k = c("lpm5", "lfp5", "llcusd")  # variable names for "Chemicals and related products"
names_i = levels(ERPT$id_i)[c(1,6,2,5,4,3,7)]  # ordered country names
L.data  = sapply(names_i, FUN=function(i) ts(ERPT[ERPT$id_i==i, names_k], start=c(1995, 1), frequency=12), simplify=FALSE)

# Oersal,Arsova 2016:21, Tab.6 (only for individual results) #
R.lags = c(3, 3, 3, 4, 3, 3, 3); names(R.lags)=names_i  # lags of VAR model by MAIC
R.cain = pcoint.CAIN(L.data, lags=R.lags, type="SL_trend")
R.pcsl = pcoint.SL(L.data,   lags=R.lags, type="SL_trend")

# Oersal,Arsova 2016:22, Tab.7/8 #
R.lags = c(3, 3, 3, 4, 4, 3, 4); names(R.lags)=names_i  # lags of VAR model by MAIC
R.t_D  = list(t_break=89)  # a level shift and trend break in 2002_May for all countries
# L.t_D = sapply(names_i, function(i) list(t_break=89), simplify=FALSE)
# L.t_D$France$t_shift = 89
R.cain = pcoint.CAIN(L.data, lags=R.lags, t_D=R.t_D, type="SL_trend")

R.cain$CSD$rho_tilde
R.cain$CSD$rho_eps
sink(file=paste0(path_tab, "Tab_ERPT.tex"))
  toLatex(R.cain)
sink(file=NULL)


