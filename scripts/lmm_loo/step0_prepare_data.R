# Title: Partition the data for LMM runs
# Author: Meixi Lin
# Date: Wed Jun 14 21:48:05 PDT 2023
# ARCHIVED versions: 
# March 2023: /NOBACKUP/scratch/meixilin/grenenet/ARCHIVE/deltap_lmm_202303

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses")
date()
sessionInfo()

library(data.table)
library(dplyr)

# def functions --------

# def variables --------
dir.create('data/lmm_loo/info/')

# load data --------
# load LD pruned, scaled (divided by p(1-p)), change in allele freq SNPs
load('data-raw/snp_freq/merged_delta_p_745_ldpruned.rda')
# load starting allele frequencies
load('data-raw/snp_freq/seedmix_p0_231_ldpruned.rda')
# load climate data
load('../metadata/data/weatherstation_bioclim.rda')

# main --------
# output metadata to be used ========
# year here does not match samples_data. only used to match generations. 
metadt = reshape2::colsplit(colnames(deltap_p), '_', names = c('site', 'generation', 'plot')) %>% 
    dplyr::mutate(mergeid = colnames(deltap_p)) %>% 
    dplyr::group_by(site, generation) %>%
    dplyr::mutate(nplot = n(),
                  year = 2017 + generation) %>% 
    dplyr::ungroup()

# merge by year for now (which can substitute for generations)
metadt = dplyr::left_join(metadt, weatherstation_bioclim, by = c('site', 'year')) %>% 
    dplyr::relocate(mergeid) 
table(metadt[,c('generation', 'year')])
# g/y2018 2019 2020
# 1  326    0    0
# 2    0  226    0
# 3    0    0  193
# group samples by generation
metadt_gen = base::split(metadt, metadt$generation)
samp_gen = sapply(metadt_gen, function(xx){xx$mergeid})

# output files --------
# output data split indices ========
save(metadt_gen, file = 'data/lmm_loo/info/weather_bygen.rda')
save(samp_gen, file = 'data/lmm_loo/info/mergeid_bygen.rda')

# cleanup --------
date()
closeAllConnections()



