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

# def variables --------
mygen = 3

# load data --------
fstm = load_fstm(paste0('data/fst/fst-matrix_gen', mygen, '.csv'))
mergep_gen = data.table::fread(file = paste0('data/fst/mergep_gen', mygen, '.csv'))
load('data-raw/snp_freq/snp_sampleinfo_745.rda')

# main --------
heatmaply::heatmaply(fstm, Rowv = NA, Colv = NA)
# get average Fst per sample ========
fsta = apply(fstm, 1, mean, na.rm = TRUE)
tail(sort(fsta))
head(sort(fsta))
# not much pattern on what samples have high Fst
snp_samples[snp_samples$sample_name %in% names(tail(sort(fsta))), ]
hist(fsta)

# samples with very low Fst all had not a lot of flowers sampled
snp_samples[snp_samples$sample_name %in% names(which(fsta < 0)),]

# get na Fst ========
nasamps = find_na(fstm)
nasamps_mp = paste0(nasamps, '.freq')
# they had AF=0/1 compared to others
summary(mergep_gen[, ..nasamps_mp])
summary(mergep_gen[, 3:7])
# they all were sampled with total_flower_counts = 1 
all((snp_samples$sample_name %in% nasamps) == 
        (snp_samples$total_flower_counts == 1 & snp_samples$generation== 1))
snp_samples[snp_samples$sample_name %in% nasamps, ]

myids = snp_samples[snp_samples$total_flower_counts > 10 & snp_samples$generation== 1, 'sample_name']
fstm5 = fstm[myids, myids]



# rescale fst matrix for linear mixed model residual structure ========
notnasamp = setdiff(dimnames(fstm)[[1]], nasamps)
fstmf = fstm[notnasamp, notnasamp]
hist(fstmf)
summary(as.vector(fstmf))
heatmaply::heatmaply(fstmf, Rowv = NA, Colv = NA)

fstmfp = fstmf
sum(fstmf < 0) # 1708 had values less than zero
fstmfp[fstmfp < 0] = 0
outfst = 1 - fstmfp

# rescale does not change density
fstms = scales::rescale(fstm)

par(mfrow = c(3,1))
hist(fstmf)
hist(fstmfp)
plot(hist(as.vector(fstmf), na.rm = TRUE))
plot(hist(as.vector(fstmfp), na.rm = TRUE))
plot(density(as.vector(fstms), na.rm = TRUE))

hist(snp_samples[snp_samples$generation== 1, 'total_flower_counts'], breaks = 300, xlim = c(1,20))
outfst = 1-fstms
covfst = Matrix::Matrix(outfst, sparse = TRUE)
save(covfst, file = 'data/fst/covfst_gen3.rda')

# reproduce minus fst values ========
mergep1 = mergep_gen[,c('55_1_10.freq', '1_1_1.freq')]


# output files --------

# cleanup --------
date()
closeAllConnections()

# SCRATCH --------
# NA sample exploration


mergep_gen[1:5,..nasamps_mp]

# is there any problems


summary(fstm[, '55_1_10'])


hist(fstms)



