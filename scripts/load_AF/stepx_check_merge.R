# Title: Double checking that the merging of generation is correct
# Author: Meixi Lin
# Date: Wed Mar  1 17:16:39 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses")
date()
sessionInfo()

library(dplyr)
library(diffdf)
# def functions --------

# def variables --------

# load data --------
# table for merged samples
sampdt = read.csv(file = '/NOBACKUP/scratch/xwu/grenet/haf-pipe/hafpipe_231/merged_sample_table.csv')

# raw sample data
load('../metadata/data/samples_data.rda')

# main --------
osampdt = samples_data %>% 
    dplyr::filter(usesample) %>%
    dplyr::mutate(newvar = paste(site, generation, plot, sep = '_')) %>% 
    dplyr::group_by(newvar) %>%
    dplyr::summarise(nruns = n(),
                     nflowers = as.integer(sum(flowerscollected)),.groups = 'drop') %>% 
    dplyr::rename(sample_name = newvar,
                  total_flower_counts = nflowers,
                  sample_times = nruns) %>% 
    dplyr::arrange(sample_name)

sampdt = sampdt %>% 
    dplyr::arrange(sample_name)

# compare
diffdf(osampdt, sampdt[,c(1,2,4)])

# output files --------
# PERFECT the samples are generated correctly. 




# cleanup --------
date()
closeAllConnections()
