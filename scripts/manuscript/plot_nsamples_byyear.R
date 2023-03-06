# Title: Number of samples by year
# Author: Meixi Lin
# Date: Mon Mar  6 13:47:02 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses")
date()
sessionInfo()

library(ggplot2)
library(dplyr)

# def functions --------

# def variables --------

# load data --------
load('data/lmm_loo/test/merged_delta_p_776_30k.rda')

# main --------
metadt = reshape2::colsplit(colnames(deltapt), '_', names = c('site', 'generation', 'plot')) %>% 
    dplyr::group_by(site, generation) %>%
    dplyr::mutate(nplot = n())

pp <- ggplot(metadt, aes(x = generation, fill = nplot, y = as.factor(site))) +
    geom_tile() +
    scale_fill_viridis_b()

# number of sites with at least 6 plots for three years
nsamplesbyyear = metadt %>% 
    dplyr::filter(nplot > 6) %>% 
    dplyr::group_by(site, plot) %>%
    dplyr::summarise(ngeneration = n()) %>% 
    dplyr::group_by(site) %>% 
    dplyr::summarise(ngen = max(ngeneration))

table(nsamplesbyyear$ngen)
# 13 sites have three generations of data that suffice the conditions
# 1  2  3  6 
# 9  6 13  1 


# output files --------

# cleanup --------
date()
closeAllConnections()
