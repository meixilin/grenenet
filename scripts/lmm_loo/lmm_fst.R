# Title: Add Fst based covariance structure for accounting false positive
# Author: Meixi Lin
# Date: Thu Jun 29 22:33:42 2023
# Use 3rd generation only

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses")
date()
sessionInfo()

library(lme4)
library(lmerTest)
library(lme4qtl)
library(data.table)
library(dplyr)
require(Matrix)

# def functions --------
deltap_bygen <- function(mydf, mygen, samp_gen) {
    mysamples = samp_gen[[mygen]]
    outdf = mydf[, ..mysamples]
    return(outdf)
}

# def variables --------
mygen = 3
envvar = 'bio1'

# load data --------
# load LD pruned, scaled (divided by p(1-p)), change in allele freq SNPs
load('data-raw/snp_freq/merged_delta_p_745_ldpruned.rda')

# load the environmental information
load('data/lmm_loo/info/mergeid_bygen.rda')
load('data/lmm_loo/info/weather_bygen.rda')

# load the coordinates of chromosomes
load('data-raw/snp_freq/seedmix_p0_231_ldpruned.rda')

# load the fst matrix
load('data/fst/covfst_gen3.rda')

heatmaply(as.matrix(covfst), Rowv = NA, Colv = NA)

# main --------
mysamp = dimnames(covfst)[[1]]
metadty = metadt_gen[[mygen]] %>% 
    dplyr::filter(mergeid %in% mysamp)
deltapsy = deltap_p[,..mysamp]

# check sample orders
all(colnames(deltapsy) == metadty$mergeid)
all(colnames(deltapsy) == mysamp)

# get nearPD
covfst1 = nearPD(covfst)
covfsti = covfst1$mat

heatmaply(as.matrix(covfst), Rowv = NA, Colv = NA)

# get pvalues --------
# if set up things as having random slopes as well, then we get singular fit
pvector1 = rep(NA, times = nrow(deltap_p))
pvector2 = rep(NA, times = nrow(deltap_p))
for (ii in 1:nrow(deltap_p)) {
    mydata = cbind(unlist(deltapsy[ii,]), metadty[,c('mergeid', 'site', 'bio1', 'bio2', 'bio3', 'bio4')])
    colnames(mydata)[1] = 'deltap'
    m0 = lme4qtl::relmatLmer(deltap ~ (1|site) + (1|mergeid), data = mydata, relmat = list(mergeid = covfsti), REML = FALSE)
    m1 = lme4qtl::relmatLmer(deltap ~ bio1 + (1|site) + (1|mergeid), data = mydata, relmat = list(mergeid = covfsti), REML = FALSE)
    a1 = anova(m0,m1)
    m00 = lme4::lmer(deltap ~ (1|site), data = mydata, REML = FALSE)
    m01 = lme4::lmer(deltap ~ bio1 + (1|site), data = mydata, REML = FALSE)
    a01 = anova(m00, m01)
    pvector1[ii] = a1$`Pr(>Chisq)`[2]
    pvector2[ii] = a01$`Pr(>Chisq)`[2]
}

pdf(file = './data/lmm_loo/lmm_fst.pdf', width = 8, height = 4)
par(mfrow = c(1,2))
qqman::qq(pvector1, main = 'deltap ~ bio1 + (1|site) + (1|mergeid)')
qqman::qq(pvector2, main = 'deltap ~ bio1 + (1|site)')
dev.off()

# output files --------
save.image(file = './data/lmm_loo/lmm_fst.RData')

# cleanup --------
date()
closeAllConnections()
