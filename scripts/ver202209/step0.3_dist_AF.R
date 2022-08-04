# Title: Caculate genetic distances from allele frequencies
# Author: Meixi Lin
# Date: Thu Aug  4 12:45:18 2022

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet")

library(poppr) # TODO: use the nei distance

date()
sessionInfo()

# def functions --------
calc_dist <- function(afdt, method = c('euclidean', 'nei'), filename, nsample = 285) {
    afdt = t(as.matrix(afdt)) # transpose to calculate by individuals
    # double check dimension
    if (dim(afdt)[1] != nsample) {
        stop('Wrong input dimension')
    }
    if (method == 'euclidean') {
        outdist = stats::dist(x = afdt, method = method)
    } else {if (method == 'nei') {
        outdist = poppr::nei.dist(x = afdt)
    } else {
        stop('Unspecified distance')
    }
    }
    saveRDS(object = outdist, file = filename)
    return(paste0(c(format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ' SUCCESS'))
}

# def variables --------


# load data --------
af285 = readRDS(file = './data/AF/ver202209/haplotype/AF284pSeed_0922.rds')

# there is missing data in two samples
miss285 = apply(af285, 2, function(xx) {any(is.na(xx))})
table(miss285, useNA = 'always')
which(miss285 == TRUE)

# find missing columns
misscol = af285[,'MLFH100220180130']
head(which(is.na(misscol)))
af285[392390:392410,'MLFH100220180130']

rm(miss285)
rm(misscol)

# remove sites with any missing data
miss285r = apply(af285, 1, function(xx) {any(is.na(xx))})
table(miss285r, useNA =  'always')
af285_nona = af285[!miss285r,] # 857365 sites
rm(af285)

# calculate meanaf for all data
meanaf285_nona = rowSums(af285_nona, na.rm = TRUE)/ncol(af285_nona)

png(file = 'plots/ver202209/AF/hist_meanaf285_nona.png', height = 4, width = 4, units = 'in')
hist(meanaf285_nona)
dev.off()

# main --------
# use different mafcutoffs
af285_nona_maf05 = af285_nona[meanaf285_nona > 0.05 & meanaf285_nona < 0.95, ] # 280092
dim(af285_nona_maf05)

af285_nona_maf10 = af285_nona[meanaf285_nona > 0.10 & meanaf285_nona < 0.90, ] # 180685
dim(af285_nona_maf10)

# output files --------
saveRDS(af285_nona, file = 'data/AF/ver202209/haplotype/cleaned/AF284pSeed_nona_0922.rds')
saveRDS(meanaf285_nona, file = 'data/AF/ver202209/haplotype/cleaned/meanAF284pSeed_nona_0922.rds')

# calculate distances and save the output
calc_dist(afdt = af285_nona, method = 'euclidean', filename = 'data/AF/ver202209/haplotype/distAF/distAF_maf00_euc_0922.rds')
calc_dist(afdt = af285_nona, method = 'nei', filename = 'data/AF/ver202209/haplotype/distAF/distAF_maf00_nei_0922.rds')
calc_dist(afdt = af285_nona_maf05, method = 'euclidean', filename = 'data/AF/ver202209/haplotype/distAF/distAF_maf05_euc_0922.rds')
calc_dist(afdt = af285_nona_maf05, method = 'nei', filename = 'data/AF/ver202209/haplotype/distAF/distAF_maf05_nei_0922.rds')
calc_dist(afdt = af285_nona_maf10, method = 'euclidean', filename = 'data/AF/ver202209/haplotype/distAF/distAF_maf10_euc_0922.rds')
calc_dist(afdt = af285_nona_maf10, method = 'nei', filename = 'data/AF/ver202209/haplotype/distAF/distAF_maf10_nei_0922.rds')

# cleanup --------
date()
closeAllConnections()

# SCRATCH --------
# afdt = af285_nona_maf10[sample(1:nrow(af285_nona_maf10), size = 10000),]
# # the original position had NaN
# 13899571,0.00535715
# 13899648,0.528198
# 13899930,0.0571844
# 13899954,0.0554402
# 13900058,nan
# 13900067,nan
# 13900105,nan
# 13900257,nan
# 13900263,nan
# 13900329,nan
# 13900340,nan
# 13900375,nan

