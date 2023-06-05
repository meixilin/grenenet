# Title: Eco type clustering with deviation from the source environments
# Author: Meixi Lin
# Date: Fri Jun  2 09:49:43 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses")
date()
sessionInfo()

# def functions --------

# def variables --------

# load data --------
load('../metadata/data/worldclim_ecotypesdata.rda')
load('../metadata/data/worldclim_sitesdata.rda')

load('./data-raw/ecotype_freq/merged_ecotype_deltap.rds')

for (ii in ecotype )

# main --------

# output files --------

# cleanup --------
date()
closeAllConnections()
