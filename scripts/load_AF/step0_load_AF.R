# Title: Load allele frequency data and convert to an RDS
# Author: Meixi Lin
# Date: Wed Feb 22 11:14:01 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses/")

date()
sessionInfo()

library(data.table)

# def functions --------

# def variables --------
# source data directory
deltapf = '/NOBACKUP/scratch/xwu/grenet/haf-pipe/hafpipe_231/merged_delta_p.csv'
p0f = paste0('/NOBACKUP/scratch/xwu/grenet/seed_mix/haf-pipe/231_vcf/chr', 1:5, '_p0.txt')

# load deltap --------
deltap = data.table::fread(input = deltapf, nThread = 6)
any(is.na(deltap))
duplicated(deltap)

# load p0 --------
p0l = lapply(1:5, function(chrn) {
    p0f = paste0('/NOBACKUP/scratch/xwu/grenet/seed_mix/haf-pipe/231_vcf/chr', chrn, '_p0.txt')
    mydt = data.table::fread(input = p0f, nThread = 6) 
    setnames(mydt, 1:2, c('POS', 'p0'))
    mydt[['CHR']] = chrn
    setcolorder(mydt, c('CHR', 'POS', 'p0'))
    return(mydt)
})

p0 = data.table::rbindlist(p0l)

# output files --------
# 2948475 rows x 776 columns. rows are genomic sites. columns are samples. 
save(deltap, file = 'data/AF/merged_delta_p_776.rda') 
# 2948475 rows x 3 columns. columns 1 and 2 are CHR and POS. column 3 is the p0 from 231 founders. 
save(p0, file = 'data/AF/seedmix_p0_231.rda')

# cleanup --------
date()
closeAllConnections()
