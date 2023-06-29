# Title: Prepare Fst matrix using grenedalf
# Author: Meixi Lin
# Date: Wed Jun 28 21:45:59 2023

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
mergep_bygen <- function(mydf, mygen, samp_gen) {
    mysamples = samp_gen[[mygen]]
    outdf = mydf[, ..mysamples]
    return(outdf)
}

# def variables --------
outdir = './data/fst/'
dir.create(outdir)

# load data --------
# the merged allele frequencies
load('./data-raw/snp_freq/merged_p_745_ldpruned.rda')
# the chromosome positions
load('./data-raw/snp_freq/seedmix_p0_231_ldpruned.rda')
# total flowers
load('./data-raw/snp_freq/snp_sampleinfo_745.rda')
# samples by generation
load('./data/lmm_loo/info/mergeid_bygen.rda')

# main --------
# split mergep_p by generation
mergep_pg <- lapply(1:3, mergep_bygen, mydf = mergep_p, samp_gen = samp_gen)

# split total flowers by generation
flowers_gen <- lapply(1:3, function(mygen){
    # get the snp samples
    df = snp_samples[snp_samples$generation == mygen, ]
    # get sample names are matching with mergep_pg
    if(!all(df$sample_name == colnames(mergep_pg[[mygen]]))) {
        stop('Sample names not matching')
    } else {
        write.table(df$total_flower_counts, 
                    file = paste0(outdir, 'nflowers_gen', mygen, '.txt'), 
                    quote = FALSE, row.names = FALSE, col.names = FALSE)
        return(invisible())
    }
})

# output grenedalf compliant files
mergep_pgdalf <- lapply(1:3, function(mygen) {
    df = mergep_pg[[mygen]]
    colnames(df) = paste0(colnames(df), '.freq')
    outdt = cbind(p0_p[,c('CHR','POS')], df)
    write.csv(outdt, file = paste0(outdir, 'mergep_gen', mygen, '.csv'), 
              quote = FALSE, row.names = FALSE)
    return(outdt)
})

# some checks --------
# TOFIX: some chr-pos sites should be invariant but does not show up now? 
lapply(mergep_pg, function(df){
    which(apply(df, 1, function(xx){all(xx == 0)}))
})


# cleanup --------
date()
closeAllConnections()



