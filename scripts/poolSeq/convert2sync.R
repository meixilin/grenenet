# Title: Convert files to sync
# Author: Meixi Lin
# Date: Thu Mar 16 10:21:10 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses")

library(poolSeq)
library(dplyr)
library(data.table)
date()
sessionInfo()

# def functions --------

# def variables --------
load(file = '../metadata/data/samples_data.rda')
s4 = samples_data %>% dplyr::filter(site == 4)

plot(s4$date)

# load data --------
# load allele frequency at time t
pt = data.table::fread(file = '/NOBACKUP/scratch/xwu/grenet/haf-pipe/hafpipe_231/merged_allele_frequency.csv', nrows = 10000)
# load starting allele frequency
load('data/AF/seedmix_p0_231.rda')
p0 = p0[1:10000, ]

# main --------
estimateNe(p0 = p0$p0, pt = pt[[56]], cov0 = 100, covt = 100, t = 2, Ncensus = 10000)

plot(p0$p0,pt[[56]])
# output files --------

# cleanup --------
date()
closeAllConnections()
