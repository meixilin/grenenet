# Title: Run the linear mixed model
# Author: Meixi Lin
# Date: Mon Mar  6 14:23:46 2023
# ARCHIVED: March 2023 version

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses")

library(dplyr)
library(data.table)
library(nlme)
library(MuMIn)
library(logger)
date()
sessionInfo()

# def functions --------
deltap_bygen <- function(mydf, mygen, samp_gen) {
    mysamples = samp_gen[[mygen]]
    outdf = mydf[, ..mysamples]
    return(outdf)
}

get_LRT <- function(mydata, mymodel, fixed_effects) {
    # the linear mixed effect model without the fixed effect
    mymodel_nofix = tryCatch({
        nlme::lme(fixed = deltapp ~ 1, random = ~1|site, data = mydata, method = 'ML')
    }, error = function(cond) {
        return(NULL)
    })
    if (is.null(mymodel_nofix)) {
        outl = vector('list', length = 2)
    } else {
        # the linear model
        mymodel_lm = lm(formula = stats::as.formula(fixed_effects), data = mydata)
        lrt_nofix = anova(mymodel,mymodel_nofix)
        lrt_norand = anova(mymodel,mymodel_lm)
        outl = list(lrt_nofix, lrt_norand)
    }
    return(outl)
}

run_lme <- function(deltapp, metadty, modelvars) {
    # check if deltap == 0 
    if (all(deltapp == 0)) {
        outdt = rep(NA_real_, times = length(modelvars))
    } else {
        # assemble the data
        mydata = cbind(deltapp, metadty[,c('mergeid', 'site', envvar)])
        # the linear mixed effect model
        mymodel = tryCatch({
            nlme::lme(fixed = stats::as.formula(fixed_effects), random = ~1|site, 
                      data = mydata, method = 'ML')
        }, error = function(cond) {
            return(NULL)
        })
        if (is.null(mymodel)) {
            outdt = rep(NA_real_, times = length(modelvars))
        } else {
            lmesum = summary(mymodel)
            lmer2 = MuMIn::r.squaredGLMM(mymodel)
            # run LRT
            lrtl = get_LRT(mydata, mymodel, fixed_effects)
            lrt_nofix = lrtl[[1]]; lrt_norand = lrtl[[2]]
            if (is.null(lrt_nofix)) {
                # get relevant stats
                outdt = c(lmer2,
                          lmesum$tTable[envvar,"Value"],
                          lmesum$tTable[envvar,"p-value"],
                          lmesum$BIC,
                          rep(NA_real_, times = 4))
            } else {
                # get relevant stats
                outdt = c(lmer2,
                          lmesum$tTable[envvar,"Value"],
                          lmesum$tTable[envvar,"p-value"],
                          lmesum$BIC,
                          lrt_nofix[2, 'L.Ratio'],
                          lrt_nofix[2, 'p-value'],
                          lrt_norand[2, 'L.Ratio'],
                          lrt_norand[2, 'p-value'])
            }

        }
    }
    return(outdt)
}

format_lmres <- function(lmres0, modelvars, mygen) {
    colnames(lmres0) = modelvars
    lmres0[, 'gen' := mygen]
    return(lmres0)
}

# def variables --------
# envvar = as.character(args[2])
envvar = 'bio1' # only run bio1 for now
fixed_effects = paste0('deltapp ~ ', envvar)
print(fixed_effects)

# generations
gens = c(1,2,3)

# model variables collected
modelvars = c('R2m', 'R2c', 'beta', 'beta_p', 'BIC', 'LR_nofix','LR_nofix_p', 'LR_norand', 'LR_norand_p')

# load data --------
# load LD pruned, scaled (divided by p(1-p)), change in allele freq SNPs
load('data-raw/snp_freq/merged_delta_p_745_ldpruned.rda')

# load the environmental information
load('data/lmm_loo/info/mergeid_bygen.rda')
load('data/lmm_loo/info/weather_bygen.rda')

# load the coordinates of chromosomes
load('data-raw/snp_freq/seedmix_p0_231_ldpruned.rda')

# main --------
# loop through each generation
for (mygen in gens) {
    deltapsy = deltap_bygen(deltap_p, mygen, samp_gen)
    metadty = metadt_gen[[mygen]]
    # confirm that colnames is matching metadt
    if (!all(metadty$mergeid == colnames(deltapsy))) {
        stop('mergeid mismatch with weather data')
    }
    # iterate across genomic sites
    lmres0 = apply(deltapsy, 1, run_lme, metadty = metadty, modelvars = modelvars) %>%
        t() %>% data.table()
    lmres = format_lmres(lmres0,modelvars,mygen)
    log_success(paste0('Done linear mixed model ', fixed_effects, 
                       ' for gen = ', mygen, 
                       '. n = ', ncol(deltapsy)))
    
    # output the results by gen
    mylmm_file = paste0('data/lmm_loo/lmeout_', envvar, '_gen', mygen, '.rda')
    save(lmres, file = mylmm_file)
    log_info(paste0('Outputting to ', mylmm_file))
}

# output files --------

# cleanup --------
date()
closeAllConnections()
