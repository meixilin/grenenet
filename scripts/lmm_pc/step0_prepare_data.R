# Title: Load LMM required data 
# Author: Meixi Lin
# Date: Mon Aug 28 14:25:06 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses")

library(dplyr)

date()
sessionInfo()

# def functions --------
# load and check that the ecotypeid for PCs and delta_ecotype_p are the same
load_ecotypeids <- function() {
    PC_ecotypeids = read.delim("/NOBACKUP/scratch/xwu/grenet/hapFIRE_updatedVCF/plink.eigenvec",header = F,sep = " ")[,1]
    ecop0s1 = read.delim('/NOBACKUP/scratch/xwu/grenet/hapFIRE_frequencies/seed_mix/s1_ecotype_frequency.txt', header = F)
    stopifnot(all(PC_ecotypeids == ecop0s1$V1))
    return(PC_ecotypeids)
}

# load mergeids based on merged_hapFIRE_allele_frequency.csv
load_mergeids <- function() {
    mergeids = colnames(read.csv("/NOBACKUP/scratch/xwu/grenet/merged_frequency/merged_hapFIRE_allele_frequency.csv", nrows = 1)) %>% 
        stringr::str_remove(.,'X') 
    # temporarily use generation + 1 as the year. (may not be accurate)
    metadt = mergeids %>% 
        reshape2::colsplit(., '_', names = c('site', 'generation', 'plot')) %>% 
        dplyr::mutate(mergeid = mergeids) %>% 
        dplyr::group_by(site, generation) %>%
        dplyr::mutate(nplot = n(),
                      year = 2017 + generation) %>% 
        dplyr::ungroup()
    print(table(metadt[,c('generation','year')]))
    return(metadt)
}

# def variables --------
outdir = './data/lmm_pc/inputs/'
dir.create(outdir, recursive = T)

generations = 1:3

# load data --------
# load founder VCF PCs
PCs <- read.delim("/NOBACKUP/scratch/xwu/grenet/hapFIRE_updatedVCF/greneNet_final_v1.1_LDpruned_10PCs.txt",header = F,sep = " ")
# load ecotype frequency shifts 
ecodeltap <- read.delim("/NOBACKUP/scratch/xwu/grenet/merged_frequency/delta_ecotype_freq.txt")
ecotypeids <- load_ecotypeids()

# load deltap mergeid
metadt0 = load_mergeids()
# load worldclim climate data
load('../metadata/data/worldclim_sitesdata.rda')
# load loations data
load('../metadata/data/locations_data.rda')

# merge by site
metadt = dplyr::left_join(metadt0, worldclim_sitesdata, by = 'site') %>% 
    dplyr::left_join(., locations_data[,c('site', 'longitude', 'latitude', 'altitude')], by = 'site')
any(is.na(metadt))

# split by generations
metadt_gen = base::split(metadt, metadt$generation)

# main --------
# organize env_sites and pop_strc for each generations
lapply(metadt_gen, function(env_sites) {
    mygen = unique(env_sites$generation)
    # get mergeids in this generation
    mergeids = env_sites$mergeid
    # sanity check that this is unique
    stopifnot(all(!duplicated(mergeids)))
    # subset ecodeltap for this generation
    myecodeltap = t(as.matrix(ecodeltap[,paste0('X', mergeids)]))
    pop_strc = myecodeltap %*% as.matrix(PCs)
    dimnames(pop_strc)[[2]] = paste0('PC', 1:10)
    save(env_sites, pop_strc, file = paste0(outdir, 'envsites_popstrc_gen', mygen, '.rda'))
    return(invisible())
})

# output files --------
save.image(paste0(outdir, 'step0_prepare_data.RData'))

# cleanup --------
date()
closeAllConnections()
