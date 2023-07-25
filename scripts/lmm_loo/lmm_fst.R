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

# heatmaply(as.matrix(covfst), Rowv = NA, Colv = NA)

# main --------
# only include samples 
mysamp = dimnames(covfst)[[1]]
metadty = metadt_gen[[mygen]] %>% 
    dplyr::filter(mergeid %in% mysamp)
deltapsy = deltap_p[,..mysamp]

# check sample orders
all(colnames(deltapsy) == metadty$mergeid)
all(colnames(deltapsy) == mysamp)

# deltap_p = deltap_p[1:10,]

# get pvalues --------
pvector1 = rep(NA, times = nrow(deltap_p))
pvector2 = rep(NA, times = nrow(deltap_p))
pvector3 = rep(NA, times = nrow(deltap_p))
pvector4 = rep(NA, times = nrow(deltap_p))
for (ii in 1:nrow(deltap_p)) {
    mydata = cbind(unlist(deltapsy[ii,]), metadty[,c('mergeid', 'site', 'bio1', 'bio2', 'bio3', 'bio4')])
    colnames(mydata)[1] = 'deltap'
    # model 1
    m01 = lme4qtl::relmatLmer(deltap ~ (1|site), data = mydata, REML = FALSE)
    m11 = lme4qtl::relmatLmer(deltap ~ bio1 + (1|site), data = mydata, REML = FALSE)
    a1 = anova(m01, m11)
    # model 2
    m02 = lme4qtl::relmatLmer(deltap ~ (1|site) + (1|mergeid), data = mydata, REML = FALSE)
    m12 = lme4qtl::relmatLmer(deltap ~ bio1 + (1|site) + (1|mergeid), data = mydata, REML = FALSE)
    a2 = anova(m02,m12)
    # model 3
    m03 = lme4qtl::relmatLmer(deltap ~ (1|site) + (1|mergeid), data = mydata, relmat = list(mergeid = covfst), REML = FALSE)
    m13 = lme4qtl::relmatLmer(deltap ~ bio1 + (1|site) + (1|mergeid), data = mydata, relmat = list(mergeid = covfst), REML = FALSE)
    a3 = anova(m03,m13)
    
    # model 4
    m04 = lme4qtl::relmatLmer(deltap ~ bio2 + bio3 + bio4 + (1|site) + (1|mergeid), data = mydata, relmat = list(mergeid = covfst), REML = FALSE)
    m14 = lme4qtl::relmatLmer(deltap ~ bio1 + bio2 + bio3 + bio4 + (1|site) + (1|mergeid), data = mydata, relmat = list(mergeid = covfst), REML = FALSE)
    a4 = anova(m04,m14)  
    # only get p vectors for now
    pvector1[ii] = a1$`Pr(>Chisq)`[2]
    pvector2[ii] = a2$`Pr(>Chisq)`[2]
    pvector3[ii] = a3$`Pr(>Chisq)`[2]
    pvector4[ii] = a4$`Pr(>Chisq)`[2]
}

pdf(file = './data/lmm_loo/lmm_fst.pdf', width = 8, height = 8)
par(mfrow = c(2,2))
qqman::qq(pvector1, main = 'deltap ~ bio1 + (1|site)')
qqman::qq(pvector2, main = 'deltap ~ bio1 + (1|site) + (1|mergeid)')
qqman::qq(pvector3, main = 'deltap ~ bio1 + (1|site) + (1|mergeid) [covfst]')
qqman::qq(pvector4, main = 'deltap ~ bio1~4 + (1|site) + (1|mergeid) [covfst]')
dev.off()

# output files --------
save.image(file = './data/lmm_loo/lmm_fst.RData')

# cleanup --------
date()
closeAllConnections()
