# Title: plot_ecotype_site
# Author: Meixi Lin
# Date: Tue Aug 29 12:04:16 2023

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
plotdir = 'data/manuscript/plot_ecotype_site/'
dir.create(plotdir)

# load data --------

# main --------
pp <- ggplot()

# output files --------
ggsave(filename = , width = , height = , path = plotdir)
save.image(file = paste0(plotdir, 'plotdata.rda'))

# cleanup --------
date()
closeAllConnections()
