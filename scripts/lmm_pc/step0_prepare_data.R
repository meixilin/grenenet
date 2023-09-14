# Title: Load LMM required data 
# Author: Meixi Lin
# Date: Mon Aug 28 14:25:06 2023
# Modification: Load deltap as well
# Date: Tue Sep 12 13:43:06 2023


# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses")

library(data.table)
library(dplyr)

date()
sessionInfo()

# def functions --------
loadSomeRData <- function(x, file) {
    E <- new.env()
    load(file=file, envir=E)
    return(get(x, envir=E, inherits=F))
}

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

# column headers are defined in env_sites
# deltaps: normalized and merged change in allele frequnecy [dp/(p0(1-p0))]
# total rows: 3235480 no batching performed
# after removing completely LD sites, total rows: 1243018
load_normdeltap <- function() {
    # this file had a header V1,V2,V3, so you need to skip = 1 to read the first SNP (but fread also automatically find the header)
    # TODO: not able to verify the header, but assuming it's the same as merged_hapFIRE_allele_frequency.csv
    infile = '/NOBACKUP/scratch/xwu/grenet/merged_frequency/merged_hapFIRE_delta_p_normed.txt' 
    dt = data.table::fread(file = infile, skip = 1)
    return(dt)
}

prune_normdeltap <- function(dt) {
    # paste all columns together
    pastedt = apply(dt, 1, paste0, collapse = '')
    uniqueid = !duplicated(pastedt) # this need to be saved too 
    return(uniqueid)
}

subset_prunenormdeltap <- function(env_sites, dt) {
    # get real mergeid headers
    mergeids = colnames(read.csv("/NOBACKUP/scratch/xwu/grenet/merged_frequency/merged_hapFIRE_allele_frequency.csv", nrows = 1)) %>% 
        stringr::str_remove(.,'X') 
    myind = which(mergeids %in% env_sites$mergeid)
    stopifnot(all(mergeids[myind] == env_sites$mergeid)) # check sampling will not change orders
    # rows are LD pruned SNPs, columns are samples with mergeid falling in given generation (reflected in env_sites)
    outdt = dt[, ..myind] 
    print(dim(outdt))
    return(outdt)
}

# def variables --------
outdir = './data/lmm_pc/inputs/'
dir.create(outdir, recursive = T)

generations = 1:3

# load data --------
# load founder VCF PCs
PCs <- read.delim("/NOBACKUP/scratch/xwu/grenet/hapFIRE_updatedVCF/greneNet_final_v1.1_LDpruned_10PCs.txt",header = F,sep = " ")
# load and calculate scaled ecotype frequency shifts 
load('./data/local_adaptation/inputs/ecotypefreqs_deltaenv.rda')
ecop <- read.delim("/NOBACKUP/scratch/xwu/grenet/merged_frequency/merged_ecotype_frequency.txt")
ecop0 <- loadSomeRData(x = 'ecop0', file = )

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

# load merged, normalized allele frequency change
ndp = load_normdeltap()
pruneids = prune_normdeltap(dt = ndp)
table(pruneids)
# prune merged, normalized allele frequency change for completely the same rows
pndp = ndp[pruneids, ]
dim(pndp)
pndp[1:5,1:5]
save(pndp, file = paste0(outdir, 'pruned_merged_hapFIRE_delta_p_normed.rda'))

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
    # subset pruned, merged, normalized allele frequency shifts for this generation
    # no transpose performed, computationally too expensive
    deltap = subset_prunenormdeltap(env_sites = env_sites, dt = pndp)
    save(env_sites, pop_strc, deltap, file = paste0(outdir, 'envsites_popstrc_deltap_gen', mygen, '.rda'))
    return(invisible())
})

# output files --------
save(pruneids, file = paste0(outdir, 'LDpruneids_merged_hapFIRE_delta_p_normed.rda'))
rm(pndp, ndp) # remove large data files
save.image(paste0(outdir, 'step0_prepare_data.RData'))

# cleanup --------
date()
closeAllConnections()
