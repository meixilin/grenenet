# Title: Run the linear mixed model (by batch)
# Author: Meixi Lin
# Date: Mon Mar  6 14:23:46 2023

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
get_deltapf <- function(mybatch) {
    mydeltapf = list.files(path = '/NOBACKUP/scratch/meixilin/grenenet/deltap_lmm', 
                           pattern = paste0('deltap_', stringr::str_pad(mybatch, width = 3, side = 'left', pad = '0'), '_'),
                           full.names = TRUE)
    
    if (length(mydeltapf) != 1) {
        stop('File not found')
    }
    log_info(paste0('Loading file: ', mydeltapf))
    return(mydeltapf)
}

deltap_byyear <- function(mydf, myyear, samp_year) {
    mysamples = samp_year[[myyear]]
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
            nlme::lme(fixed = stats::as.formula(fixed_effects), random = ~1|site, data = mydata, method = 'ML')
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

format_lmres <- function(lmres0, modelvars, mycoords, myyear) {
    colnames(lmres0) = modelvars
    # plot(lmres0$beta_p, lmres0$LR_nofix_p) LRT is pretty much perfect correlation with beta p-value
    # add chromosome info and year info
    lmres0 = cbind(mycoords, lmres0)
    lmres0[, 'year' := myyear]
    return(lmres0)
}

# def variables --------
args = commandArgs(trailingOnly=TRUE)
mybatch = as.integer(args[1])
# envvar = as.character(args[2])
envvar = 'bio1' # only run bio1 for now
fixed_effects = paste0('deltapp ~ ', envvar)
print(fixed_effects)

years = c('2018', '2019', '2020')

# model variables collected
modelvars = c('R2m', 'R2c', 'beta', 'beta_p', 'BIC', 'LR_nofix','LR_nofix_p', 'LR_norand', 'LR_norand_p')

# my file to load
mydeltap_file = get_deltapf(mybatch)
mylmm_prefix = gsub('.rda', '', gsub('deltap_lmm/deltap', 'deltap_lmm/lmeout', mydeltap_file))

# load data --------
load(mydeltap_file)

# load the batch information 
load('data/lmm_loo/info/mergeid_byyear.rda')
load('data/lmm_loo/info/weather_byyear.rda')
coordid = read.csv('data/lmm_loo/info/coord_splits.csv', row.names = 1)[mybatch,]

# load the coordinates of chromosomes
load('data/AF/seedmix_p0_231.rda')
mycoords = p0[coordid$start:coordid$end, c('CHR', 'POS')]
rm(p0)

# main --------
lmresl <- lapply(years, function(myyear){
    deltapsy = deltap_byyear(deltaps, myyear, samp_year)
    metadty = metadt_year[[myyear]]
    # confirm that colnames is matching metadt
    if (!all(metadty$mergeid == colnames(deltapsy))) {
        stop('mergeid mismatch with weather data')
    }
    # iterate across sites
    lmres0 = apply(deltapsy, 1, run_lme, metadty = metadty, modelvars = modelvars) %>%
        t() %>% data.table()
    lmres = format_lmres(lmres0,modelvars,mycoords,myyear)
    log_success(paste0('Done linear mixed model ', fixed_effects, ' for year = ', myyear))
    # output the results by year
    mylmm_file = paste0(mylmm_prefix, '_', myyear, '.rda')
    save(lmres, file = mylmm_file)
    log_info(paste0('Outputting to ', mylmm_file))
    return(invisible())
})

# output files --------

# cleanup --------
date()
closeAllConnections()
