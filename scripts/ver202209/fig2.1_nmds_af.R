# Title: Run an NMDS on allele frequencies
# Author: Meixi Lin
# Date: Thu Aug  4 00:20:37 2022

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet")

date()
library(vegan)
library(ggplot2)
library(dplyr)
sessionInfo()

# def functions --------
plot_

# def variables --------

# load data --------
# distance measures


# metadata ========
meta0922 = readRDS(file = 'data/metadata/ver202209/metacensusclim_0922.rds')

# main --------
dd = dist(sampledt)
# calculate distance
# sampledt = t(af285_nona_maf05[sample(1:nrow(af285_nona_maf05), size = 10000),])
sampledt = t(af285_nona_maf05)
dim(sampledt)

# Note that site x species. here we need to transpose the dataset
af285_NMS <-vegan::metaMDS(sampledt,
                           distance = "bray",
                           k = 3,
                           trymax = 20)

str(af285_NMS)

# plot NMDS
forplot = af285_NMS$points %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = 'id') %>%
    dplyr::left_join(., y = meta0922[,c('id','SITE', 'year')], by = 'id')

pp <- ggplot(forplot, aes(x = MDS1, y = MDS2, color = year)) +
    geom_point() +
    ggpubr::theme_pubr()

ggsave(filename = 'plots/nmds_af285_nona_maf05_s10k.pdf', height = 5, width = 5, plot = pp)

# see rep level
forplot1 = forplot[forplot$SITE %in% c(NA, 1:10),]
mycolors = c(pals::alphabet(n = 26), pals::alphabet2(n = 6)); names(mycolors) = unique(as.character(forplot$SITE))
hull_forplot <- forplot %>%
    group_by(SITE) %>%
    slice(chull(MDS1,MDS2))

pp1 <- ggplot(forplot, aes(x = MDS1, y = MDS2, color = as.character(SITE), fill = as.character(SITE))) +
    geom_point() +
    geom_polygon(data = hull_forplot, alpha = 0.2) +
    labs(color = "") +
    scale_color_manual(values = mycolors) +
    scale_fill_manual(values = mycolors) +
    ggpubr::theme_pubr() +
    theme(legend.position = 'none')

pp1.1 <- pp1 +
    lims(x = c(-0.1,0.1), y = c(-0.1,0.1))
    ()
pp2 <- ggplot(forplot, aes(x = MDS1, y = MDS3, color = as.character(SITE))) +
    geom_point() +
    labs(color = "") +
    scale_color_manual(values = mycolors) +
    ggpubr::theme_pubr() +
    theme(legend.position = 'none')

pp3 <- ggplot(forplot, aes(x = MDS3, y = MDS2, color = as.character(SITE))) +
    geom_point() +
    labs(color = "") +
    scale_color_manual(values = mycolors) +
    ggpubr::theme_pubr() +
    theme(legend.position = 'none')

pp <- ggpubr::ggarrange(plotlist = list(pp1,pp3,pp2),ncol = 2, nrow = 2,align = 'h')

ggsave(filename = 'plots/nmds_af285c_nona_maf05_all.png', height = 8, width = 8, plot = pp)
# output files --------

# cleanup --------
date()
closeAllConnections()
