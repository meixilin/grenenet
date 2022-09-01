# Title: Check the scaled AF and the cause of the NaN, Inf values
# Author: Meixi Lin
# Date: Wed Aug 31 13:50:29 2022

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet")
date()
sessionInfo()

# def variables --------
# the variables that were filtered out in fig3
# > setdiff(rownames(lmres0), as.character(lmres$POS))
# [1] "8916"     "9094035"  "9527351"  "13914486" "21327459" "21368627"

# load data --------
# find the rows that contained NA in the start
deltadt = readRDS(file = 'data/AF/ver202209/haplotype/DeltaP/delta_freq.rds')
table(apply(deltadt, 1, function(xx){any(is.na(xx))})) # some is NA, should be FALSE 857365
naid = apply(deltadt, 1, function(xx){any(is.na(xx))})

# 13900058 contained NA for site 10_20180130
head(deltadt[naid,])
table(apply(deltadt, 1, function(xx){all(is.na(xx))})) # not all is NA

deltadt_uniq = readRDS(file = 'data/AF/ver202209/haplotype/DeltaP/delta_freq_uniq.rds')
table(apply(deltadt_uniq, 1, function(xx){any(is.na(xx))}))
table(apply(deltadt_uniq, 1, function(xx){all(is.na(xx))}))

# check the scaled AF if it's been affected in the same ways ========
scalep = readRDS(file = 'data/AF/ver202209/haplotype/DeltaP/scaled_delta_freq.rds')
scalep_uniq = readRDS(file = 'data/AF/ver202209/haplotype/DeltaP/scaled_delta_freq_uniq.rds')
table(apply(scalep_uniq, 1, function(xx){any(is.na(xx))})) # contained three more NA rows
table(apply(scalep_uniq, 1, function(xx){all(is.na(xx))}))
names(which(apply(scalep_uniq, 1, function(xx){all(is.na(xx))}))) # chromosome positions that were all NA
# also find the ones that are either NA or Inf
table(apply(scalep_uniq, 1, function(xx){all(!is.finite(xx))}))
missingids = names(which(apply(scalep_uniq, 1, function(xx){all(!is.finite(xx))}))) # the same as the ones that failed LM
missingids

# check start AF ========
startAF = readRDS(file = 'data/AF/ver202209/haplotype/AFSeed_0922.rds')
head(startAF)
# check if the raw deltaP NA locations start AF is also different
head(startAF[naid,]) # it was normal for the locations with missing deltadt

# get the six missing AF positions ========
# check if the scaled NA locations start AF is different
startAF[startAF$pos %in% missingids,] # either 1 or 0
apply(deltadt[missingids,], 1, function(xx){table(xx, useNA = 'always')}) # the raw deltaP
apply(scalep_uniq[missingids,], 1, function(xx){table(xx, useNA = 'always')})

# cleanup --------
date()
closeAllConnections()
