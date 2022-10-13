# Title: Calculate DeltaAF without merging
# Author: Meixi Lin
# Date: Wed Oct 12 11:40:09 2022

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet")
date()
sessionInfo()

library(ggplot2)
library(dplyr)

# def functions --------

# def variables --------

# load data --------
# get the raw AF
afdt0 = readRDS(file = 'data/AF/ver202209/haplotype/AF284_0922.rds')
afdt0[1:5,1:5]
# some positions were either 1. NA or 2. completely the same
clean_scaledAF = readRDS(file = 'data/AF/ver202209/haplotype/DeltaP/scaled_delta_freq_uniq.rds')
goodcoords = rownames(clean_scaledAF)
goodcoordspos = which(rownames(afdt0) %in% goodcoords)
# check for monotonic
all(goodcoordspos == cummax(goodcoordspos))
afdt = afdt0[goodcoords,]
all(rownames(afdt)==rownames(clean_scaledAF))

# get the starting AF
startAF0 = readRDS(file = 'data/AF/ver202209/haplotype/AFSeed_0922.rds')$seed_mix
any(is.na(startAF0))
startAF = startAF0[goodcoordspos]

# main --------
# calculate the deltaAF
# afdt = afdt[,1:20]
deltaAF = apply(afdt, 2, function(xx) {xx - startAF})
afdt[1:10,1:2]
startAF[1:10]
deltaAF[1:10,1:2]

# calculate the scaledAF
scalep = apply(deltaAF,2,function(xx){xx/(startAF*(1-startAF))})
deltaAF[1:10,1:2]
startAF[1:10]
scalep[1:10,1:2]

# output files --------
saveRDS(scalep, file = './data/AF/ver202209/haplotype/scaled_delta_freq_1012.rds')
saveRDS(deltaAF, file = './data/AF/ver202209/haplotype/delta_freq_1012.rds')

# cleanup --------
date()
closeAllConnections()
