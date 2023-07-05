# Title: Random slope for accounting false positive
# Author: Meixi Lin
# Date: Fri Jun 23 12:49:12 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses")
date()
sessionInfo()

library(lme4)
library(lmerTest)
library(data.table)
library(dplyr)

# def functions --------
deltap_bygen <- function(mydf, mygen, samp_gen) {
    mysamples = samp_gen[[mygen]]
    outdf = mydf[, ..mysamples]
    return(outdf)
}

# def variables --------

# load data --------
# load LD pruned, scaled (divided by p(1-p)), change in allele freq SNPs
load('data-raw/snp_freq/merged_delta_p_745_ldpruned.rda')

# load the environmental information
load('data/lmm_loo/info/mergeid_bygen.rda')
load('data/lmm_loo/info/weather_bygen.rda')

# load the coordinates of chromosomes
load('data-raw/snp_freq/seedmix_p0_231_ldpruned.rda')

# main --------
mygen = 3
envvar = 'bio1'
metadty = metadt_gen[[mygen]]
deltapsy = deltap_bygen(deltap_p, mygen, samp_gen)
# deltap_p = deltap_p[1:100,]

# get pvalues --------
# if set up things as having random slopes as well, then we get singular fit
pvector1 = rep(NA, times = nrow(deltap_p))
pvector2 = rep(NA, times = nrow(deltap_p))
pvector3 = rep(NA, times = nrow(deltap_p))
pvector4 = rep(NA, times = nrow(deltap_p))
for (ii in 1:nrow(deltap_p)) {
    mydata = cbind(unlist(deltapsy[ii,]), metadty[,c('mergeid', 'site', 'bio1', 'bio2', 'bio3', 'bio4')])
    colnames(mydata)[1] = 'deltap'
    # random intercept
    model1 = tryCatch({
        # nlme::lme(fixed = deltap ~ bio1, random = ~1|site, method = 'ML', data = mydata)
        lmerTest::lmer(deltap ~ bio1 + (1|site), data = mydata)
    }, error = function(cond) {
        print(cond)
        return(NULL)
    })
    # random slope and intercept
    model2 = tryCatch({
        # nlme::lme(fixed = deltap ~ bio1, random = ~bio1|site, method = 'ML', data = mydata, 
                  # control =list(msMaxIter = 1000, msMaxEval = 1000))
        #  (1 + bio1 || site) same as ((1|site) +(0+bio1|site))
        lmerTest::lmer(deltap ~ bio1 + (1 + bio1 || site), data = mydata)
        # lmerTest::lmer(deltap ~ bio1 + bio2 + bio3 + bio4 + (1|site), data = mydata, REML = FALSE)
    }, error = function(cond) {
        print(cond)
        return(NULL)
    })
    # random intercept but more variables 
    model3 = tryCatch({
        # nlme::lme(fixed = deltap ~ bio1, random = ~bio1|site, method = 'ML', data = mydata, 
        # control =list(msMaxIter = 1000, msMaxEval = 1000))
        #  (1 + bio1 || site) same as ((1|site) +(0+bio1|site))
        # lmerTest::lmer(deltap ~ bio1 + (1 + bio1 || site), data = mydata)
        lmerTest::lmer(deltap ~ bio1 + bio2 + bio3 + bio4 + (1|site), data = mydata)
    }, error = function(cond) {
        print(cond)
        return(NULL)
    })   
    # random intercept with per id random intercept 
    # this give the exact same pvector as model1. Seems like 1|mergeid was ignored
    model4 = tryCatch({
        # nlme::lme(fixed = deltap ~ bio1, random = ~bio1|site, method = 'ML', data = mydata, 
        # control =list(msMaxIter = 1000, msMaxEval = 1000))
        #  (1 + bio1 || site) same as ((1|site) +(0+bio1|site))
        # lmerTest::lmer(deltap ~ bio1 + (1 + bio1 || site), data = mydata)
        lmerTest::lmer(deltap ~ bio1 + (1|site) + (1|mergeid), data = mydata,
                       control=lmerControl(check.nobs.vs.nlev = "ignore",
                                           check.nobs.vs.rankZ = "ignore",
                                           check.nobs.vs.nRE="ignore"))
    }, error = function(cond) {
        print(cond)
        return(NULL)
    })   
    # 
    # pvector1[ii] = ifelse(is.null(model1), NA, summary(model1)$tTable['bio1', 'p-value'])
    # pvector2[ii] = ifelse(is.null(model2), NA, summary(model2)$tTable['bio1', 'p-value'])
    
    pvector1[ii] = ifelse(is.null(model1), NA, summary(model1)$coefficients['bio1', 'Pr(>|t|)'])
    pvector2[ii] = ifelse(is.null(model2), NA, summary(model2)$coefficients['bio1', 'Pr(>|t|)'])
    pvector3[ii] = ifelse(is.null(model3), NA, summary(model3)$coefficients['bio1', 'Pr(>|t|)'])
    pvector4[ii] = ifelse(is.null(model4), NA, summary(model4)$coefficients['bio1', 'Pr(>|t|)'])
}

pdf(file = './data/lmm_loo/lmm_randomslope.pdf', width = 8, height = 8)
par(mfrow = c(2,2))
qqman::qq(pvector1, main = 'deltap ~ bio1 + (1|site)')
qqman::qq(pvector2, main = 'deltap ~ bio1 + (1 + bio1 || site)')
qqman::qq(pvector3, main = 'deltap ~ bio1 + bio2 + bio3 + bio4 + (1|site)')
qqman::qq(pvector4, main = 'deltap ~ bio1 + (1|site) + (1|mergeid)')
dev.off()

# output files --------
save.image(file = './data/lmm_loo/lmm_randomslope.RData')

# cleanup --------
date()
closeAllConnections()
