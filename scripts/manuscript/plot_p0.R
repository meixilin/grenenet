# Title: plot_p0
# Author: Meixi Lin
# Date: Wed Mar  8 10:05:50 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses")

library(dplyr)
library(ggplot2)

date()
sessionInfo()

# def functions --------
source('scripts/manuscript/config_manuscript.R')

# def variables --------
plotdir = 'data/manuscript/plot_p0/'
dir.create(plotdir)

# load data --------
load('data/AF/seedmix_p0_231.rda')

# main --------
pp <- ggplot(data = p0, mapping = aes(x = POS, y = p0)) + 
    geom_point(size = 0.1) + 
    facet_wrap(. ~ CHR, ncol = 1, scales = 'free_x')

# output files --------
ggsave(filename = 'seedmix_p0_231.png', width = 8, height = 8, plot = pp, path = plotdir)
# save.image(file = paste0(plotdir, 'plotdata.rda'))

# cleanup --------
date()
closeAllConnections()
