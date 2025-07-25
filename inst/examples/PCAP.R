### Latex figure: public capital stocks ###
data(PCAP)
names_i = c("DEU", "FRA", "GRC", "ITA", "NLD")  # select 5 exemplary countries
idx_i  = PCAP$id_i %in% names_i
idx_pl = c(1, 3, 5, 9, 10)
pallet = c(RColorBrewer::brewer.pal(n=11, name="Spectral")[idx_pl], "#000000")
names(pallet) = c(names_i, "\\O23")
breaks = factor(c(names_i, "\\O23"), levels=c(names_i, "\\O23"))
events = data.frame(
  label   = paste0("\\quad \\textbf{ ",
            c("Oil Crisis", "Oil Crisis II", "Early 1990s", "Early 2000s", 
              "Great Recession", "European Debt Crisis", "COVID-19"), " }"), 
  t_start = c(1973, 1979, 1990, 2000, 2007, 2010, 2020),
  t_end   = c(1975, 1982, 1993, 2003, 2010, 2013, 2022))

# plot events #
library("ggplot2")
F.base <- ggplot() + 
  geom_hline(yintercept = 0, color="grey") +
  geom_rect(data=events, aes(xmin=t_start, xmax=t_end, ymin=-Inf, ymax=+Inf), 
            fill="black", alpha=0.2) +
  geom_text(data=events, aes(x=(t_start+t_end)/2, y=-Inf, label=label), 
            hjust=0, angle=90, colour='white') +
  scale_x_continuous(breaks = seq(1960, 2020, 10), limits = c(1960, 2022)) +
  theme_bw(base_size=10)

# add levels #
F.G <- F.base +
  geom_line(data=PCAP[idx_i, ], aes(x=id_t, y=G/1e+12, colour=id_i, group=id_i), 
            linewidth=2) +
  stat_summary(data=PCAP[idx_i, ], aes(x=id_t, y=G/1e+12, color="\\O23"), 
               fun=mean, geom="line", linewidth=0) +
  geom_text(data=events, aes(x=(t_start+t_end)/2, y=-Inf, label=label), 
            hjust=0, angle=90, colour='white', alpha=0.6) +
  scale_colour_manual(values=pallet, breaks=breaks) +
  labs(x=NULL, y="Trillion US-\\$ at 2005 prices", colour="Country", title=NULL) +
  guides(colour=guide_legend(override.aes = list(linewidth=2))) +
  theme(legend.position="inside", legend.position.inside=c(0.01, 0.99), 
        legend.justification = c(0, 1), 
        legend.title=element_text(size=8), legend.text=element_text(size=6), 
        legend.key.width = unit(0.35, "cm"), legend.key.height = unit(0.35, "cm"))

# add growth rates #
PCAP$gG = ave(PCAP$G, PCAP$id_i, FUN=function(x) 
  c(diff(x), NA)/x*100)  # beginning-of-the-year observations!
F.gG <- F.base +
  geom_line(data=PCAP[idx_i, ], aes(x=id_t, y=gG, colour=id_i, group=id_i), 
            linewidth=2) +
  stat_summary(data=PCAP, aes(x=id_t, y=gG, color="\\O23"), 
               fun=mean, geom="line", linewidth=2) +
  geom_text(data=events, aes(x=(t_start+t_end)/2, y=-Inf, label=label), 
            hjust=0, angle=90, colour='white', alpha=0.6) +
  scale_colour_manual(values=pallet, breaks=breaks, guide="none") +
  labs(x="Year", y="Growth in \\%", colour="Country", title=NULL)

# export to Latex #
\donttest{
library(tikzDevice)
textwidth = 15.5/2.54  # LaTeX textwidth from "cm" into "inch"
file_fig  = file.path(tempdir(), "Fig_G.tex")

tikz(file=file_fig, width=1.2*textwidth, height=0.8*textwidth)
 # gridExtra::grid.arrange(grobs=list(F.G, F.gG), layout_matrix=cbind(1:2))
 ggpubr::ggarrange(F.G, F.gG, ncol=1, nrow=2, align="v")
dev.off()
}

