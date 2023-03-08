# Title: Plot LMM for bio1
# Author: Meixi Lin
# Date: Wed Mar  8 09:32:30 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses")

library(ggplot2)
library(dplyr)

date()
sessionInfo()

# def functions --------
source('scripts/manuscript/config_manuscript.R')

# def variables --------
pcutoff = 5*10^(-8)
plotdir = 'data/manuscript/plot_lmm_bio1'
dir.create(plotdir)

# load data --------
load('data/lmm_loo/bio1/lmeout_bio1_all_2018.rda')

lmeresallpass = lmeresallp[pass_p==TRUE, ]
lmeresallsamp = lmeresallp[pass_p==FALSE, ][sample(.N,50000)]

forplot = rbind(lmeresallpass, lmeresallsamp)

# main --------
pp1 <- ggplot(data = forplot, aes(x = POS, y = -log10(LR_nofix_p), color = pass_p)) +
    geom_point(size = 0.1) +
    scale_color_manual(values = c('darkgray','green')) +
    facet_wrap(. ~ CHR, ncol = 1, scales = 'free_x') + 
    geom_hline(yintercept = -log10(pcutoff), color = 'darkgray', linetype = 'dashed') +
    labs(x = 'genomic position (bp)', y = '-log10(P)') +
    theme(legend.position = 'none')

ggsave(filename =  'lmeout_bio1_all_2018.png', plot = pp1, path = plotdir, height = 8, width = 8)

# output files --------

# cleanup --------
date()
closeAllConnections()
