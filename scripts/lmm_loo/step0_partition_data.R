# Title: Partition the data for LMM runs
# Author: Meixi Lin
# Date: Mon Mar  6 13:39:01 2023
# batch job id: 176633

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
outdir = '/NOBACKUP/scratch/meixilin/grenenet/deltap_lmm/'
dir.create('data/lmm_loo/info/')

# load data --------
load('data/AF/merged_delta_p_776.rda')
# load('data/lmm_loo/test/merged_delta_p_776_30k.rda')
# deltap = deltapt

load('data/AF/seedmix_p0_231.rda')
# load climate data
load('../metadata/data/weatherstation_bioclim.rda')

# main --------
# set up data split indices (1:10000, within chromosome) ========
p0l = base::split(p0, p0$CHR)
p0len = sapply(p0l, nrow)
chrstart = unname(cumsum(p0len[-length(p0len)]) + 1)
p0[chrstart, ]
#    CHR   POS         p0
# 1:   2 99812 0.01548715
# 2:   3   175 0.00915216
# 3:   4  1056 0.00850109
# 4:   5    53 0.06139028

mysplit_start = sort(c(seq(1,nrow(p0), by = 10000), chrstart))
mysplit_end = c(mysplit_start - 1, nrow(p0))[-1]
mysplitdf = data.frame(start = mysplit_start, end = mysplit_end) # 299 arrays

# output metadata to be used ========
metadt = reshape2::colsplit(colnames(deltap), '_', names = c('site', 'generation', 'plot')) %>% 
    dplyr::mutate(mergeid = colnames(deltap)) %>% 
    dplyr::group_by(site, generation) %>%
    dplyr::mutate(nplot = n(),
                  year = case_when(site == 57 ~ 2017 + ceiling(generation/2),
                                   TRUE ~ 2017 + generation)) %>%
    dplyr::ungroup()

# FIXME: merge by year for now
metadt = dplyr::left_join(metadt, weatherstation_bioclim, by = c('site', 'year')) %>% 
    dplyr::relocate(mergeid) 

# group samples by year
metadt_year = base::split(metadt, metadt$year)
samp_year = sapply(metadt_year, function(xx){xx$mergeid})

# output files --------
# split the data by genomic locations ========
for (ii in 1:nrow(mysplitdf)) {
    mystart = mysplitdf[ii,'start']; myend = mysplitdf[ii,'end']
    if (mystart > nrow(deltap)) {
        stop('Subscript out of bounds.')
    }
    tempdf = deltap[mystart:myend]
    save(tempdf, file = paste0(outdir, 'deltap_', stringr::str_pad(ii, width = 3, side = 'left', pad = '0'), 
                               '_', mystart, '-', myend,'.rda'))
}

# output data split indices ========
write.csv(mysplitdf, file = 'data/lmm_loo/info/coord_splits.csv')
save(metadt_year, file = 'data/lmm_loo/info/weather_byyear.rda')
save(samp_year, file = 'data/lmm_loo/info/mergeid_byyear.rda')

# cleanup --------
date()
closeAllConnections()



