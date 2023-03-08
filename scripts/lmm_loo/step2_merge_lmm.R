# Title: Merge the results of all the runs
# Author: Meixi Lin
# Date: Wed Mar  8 08:42:41 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses")
date()
sessionInfo()

library(data.table)
library(dplyr)
library(qqman)

# def functions --------
padj_lmres <- function(lmres, p_adjmethods = 'fdr') {
    # format lmres
    Pcols = lmres %>% dplyr::select(ends_with('_p'))
    adjPcols = apply(Pcols, 2, function(xx) {
        padj = p.adjust(xx, method = p_adjmethods)
    })
    colnames(adjPcols) = paste0(colnames(adjPcols), '_adj')
    # only looking at the LR_nofix_p so far
    outdt = cbind(lmres, adjPcols) %>%
        dplyr::mutate(pass_p = LR_nofix_p_adj < 0.05) %>%
        dplyr::arrange(LR_nofix_p_adj,desc(R2m))
    return(outdt)
}


# def variables --------
batches = 1:299
years = c('2018', '2019', '2020')
nsnps = 2948475

pval_loc = 'LR_nofix_p' # use the LRT output

outdir = 'data/lmm_loo/bio1/'

dir.create(outdir)

# iterate over years --------
for (myyear in years) {
    # load data
    lmresl <- lapply(batches, function(mybatch) {
        mylmresf = list.files(path = '/NOBACKUP/scratch/meixilin/grenenet/deltap_lmm', 
                              pattern = paste0('lmeout_', stringr::str_pad(mybatch, width = 3, side = 'left', pad = '0'), '_.+_', myyear),
                              full.names = TRUE)
        load(mylmresf)
        return(lmres)
    })
    lmresall = data.table::rbindlist(lmresl)
    # check nrow
    if (nrow(lmresall) != nsnps) {
        stop('Bad merging')
    }
    # adjust p value
    lmeresallp <- padj_lmres(lmresall)
    print(head(lmeresallp[, c('CHR', 'POS', 'BIC', 'R2m', 'beta', 'LR_nofix_p_adj')]))
    
    # output files
    save(lmeresallp, file = paste0(outdir, 'lmeout_bio1_all_', myyear, '.rda'))
    print(date())
}

# cleanup --------
date()
closeAllConnections()
