# Title: Explore the site SNPs patterns
# Author: Meixi Lin
# Date: Wed May 31 14:15:12 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses")
date()
sessionInfo()

library(data.table)

# def functions --------
run_sitecor <- function(ii,jj) {
    
}


# def variables --------

# load data --------
load(file = 'data/AF/merged_deltap10k_745.rda')

# main --------
# get correlation
deltapt = as.matrix(deltapt)
sitecor = cor(deltapt)
# seq along the diagonal matrix
outdt = data.table()
for (ii in seq_len(ncol(deltap))) {
    for (jj in (ii+1):ncol(deltap)) {
        print(ii)
        print(jj)
    }
}

# output files --------

# cleanup --------
date()
closeAllConnections()
