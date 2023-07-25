# Title: Load and interpret Fst results
# Author: Meixi Lin
# Date: Thu Jun 29 12:01:00 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses")
date()
sessionInfo()

library(heatmaply)
library(logger)
library(Matrix)

# def functions --------
load_fstm <- function(fstmf) {
    df = read.csv(file = fstmf)
    dfm = as.matrix(df[,-1])
    dimnames(dfm) = list(df$sample, df$sample)
    return(dfm)
}

find_na <- function(fstm) {
    allnav = apply(fstm, 1, function(xx) {all(is.na(xx) | xx == 0)})
    print(allnav[allnav])
    out = names(allnav[allnav])
    return(out)
}

plot_fst <- function(mygen, fstm, fsta) {
    pdf(file = paste0(plotdir, 'fst_histflower_gen', mygen, '.pdf'), height = 4, width = 8)
    par(mfrow = c(1,2))
    # histogram of fst
    hist(fstm)
    # correlations of fst and number of flowers
    myflowers = unlist(merged_samples_data[merged_samples_data$sample_name %in% names(fsta), 'total_flower_counts'])
    plot(x = myflowers, y = fsta, log = 'x')
    dev.off()
    return(invisible())
}

positive_semidefinite <- function(covfst) {
    eigenvalues <- eigen(covfst)$values
    return(all(eigenvalues >= 0))
}

full_rank <- function(covfst) {
    # Compute the rank of the matrix
    rank_A <- qr(covfst)$rank
    return(rank_A == min(dim(covfst)))
}

# get nearest positive semi-definite 
get_nearPD <- function(covfst, mygen) {
    # get nearPD
    covfst1 = nearPD(covfst, corr = TRUE, maxit = 1000)
    covfsti = covfst1$mat
    log_info('Transformed positive semidefinite: {positive_semidefinite(as.matrix(covfsti))}')
    # check if the fst is mostly preserved after transformation
    # for very high fst, the values were toned down
    png(filename = paste0(plotdir, 'nearPD_fst_gen', mygen, '.png'))
    plot(x = as.vector(as.matrix(covfst)), y = as.vector(as.matrix(covfsti)), 
         xlab = 'raw Fst', ylab = 'nearPD Fst')
    dev.off()
    return(covfsti)
}

# def variables --------
mygens = 1:3
plotdir = 'data/fst/plots/'
dir.create(plotdir)

# loop through each generations
for (mygen in mygens) {
    # load data --------
    fstm = load_fstm(paste0('data/fst/fst-matrix_gen', mygen, '.csv'))
    mergep_gen = data.table::fread(file = paste0('data/fst/mergep_gen', mygen, '.csv'))
    load('../metadata/data/merged_samples_data.rda')
    
    # some plotting --------
    # overall matrix
    heatmaply::heatmaply(fstm, Rowv = NA, Colv = NA, 
                         file = paste0(plotdir, 'fst_gen', mygen, '.html'))
    # summary info
    fsta = apply(fstm, 1, mean, na.rm = TRUE)
    plot_fst(mygen,fstm,fsta)
    
    log_info('Average Fst extreme values')
    print(tail(sort(fsta)))
    print(head(sort(fsta)))
    log_info('These samples had NaN')
    nasamps = find_na(fstm)
    log_info('Samples with high Fst')
    # not much pattern on what samples have high Fst
    print(merged_samples_data[merged_samples_data$sample_name %in% names(tail(sort(fsta))), ] %>% as.data.frame())
    log_info('Samples with negative Fst')
    # samples with negative Fst caused by the flower corrector
    print(merged_samples_data[merged_samples_data$sample_name %in% names(which(fsta < 0)),] %>% as.data.frame())
    log_info('Samples with Fst NaN')
    print(merged_samples_data[merged_samples_data$sample_name %in% nasamps, ] %>% as.data.frame())
    
    # output files --------
    notnasamp = setdiff(dimnames(fstm)[[1]], nasamps)
    fstmf = fstm[notnasamp, notnasamp]
    fstmf[fstmf < 0] = 0 # convert less than zero values to zero
    covfst = 1 - fstmf # convert to kinship like matrix. self-self is 1. 
    # check if positive semidefinite
    log_info('Covariance matrix full rank: {full_rank(covfst)}')
    log_info('Covariance matrix positive semidefinite: {positive_semidefinite(covfst)}')
    
    # construct covariance matrix to be positive definite
    covfst0 = Matrix::Matrix(covfst, sparse = TRUE)
    covfst = get_nearPD(covfst0, mygen)
    save(covfst, file = paste0('data/fst/covfst_gen', mygen, '.rda'))
}


# cleanup --------
date()
closeAllConnections()

