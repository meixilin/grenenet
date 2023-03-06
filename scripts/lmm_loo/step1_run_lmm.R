# Title: Run the linear mixed model (by batch)
# Author: Meixi Lin
# Date: Mon Mar  6 14:23:46 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses")
date()
sessionInfo()

# def functions --------
split_by_year <- function() {
    
}

# def variables --------
args = commandArgs(trailingOnly=TRUE)
mybatch = as.character(args[1])
myyear = as.character(args[2])

# load data --------
load('data/lmm_loo/test/merged_delta_p_776_30k.rda')
load('data/lmm_loo/test/seedmix_p0_231_30k.rda')

# main --------

# output files --------

# cleanup --------
date()
closeAllConnections()
