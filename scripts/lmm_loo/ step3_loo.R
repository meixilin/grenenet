# Title: Leave one out test for lmeres
# Author: Meixi Lin
# Date: Wed Mar  8 22:32:50 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses")

library(data.table)
library(nlme)
library(caret)

date()
sessionInfo()

# def functions --------
deltap_byyear <- function(mydf, myyear, samp_year) {
    mysamples = samp_year[[myyear]]
    outdt = mydf[, ..mysamples]
    return(outdt)
}

find_similar_p0 <- function(p0, sourceid, backgroundid) {
    myp0 = p0[['p0']][sourceid]
    mychr = p0[['CHR']][sourceid]
    similarp0id = which((abs(p0[['p0']]-myp0)/myp0 < 1e-4) & p0[['CHR']] == mychr)
    candidate1 = intersect(similarp0id, backgroundid)
    if (length(candidate1) == 0) {
        stop('No background SNP available')
    }
    pairid = sample(candidate1, size = 1)
    # get not significant ones
    return(pairid)
}

subeset_by_id_year <- function(deltap, samp_year, id, myyear) {
    xx = deltap[id, ]
    xx_year = deltap_byyear(xx, myyear, samp_year)
    xx_year = unlist(xx_year)
    return(xx_year)
}

format_chrpos <- function(p0, ii) {
    mychr = p0[[ii, 'CHR']]; mypos = p0[[ii, 'POS']]
    out = paste(mychr, mypos, sep = '-')
    return(out)
}

find_deltap_pair <- function(lmeresallp, p0, deltap, backgroundid, samp_year, ii, myyear) {
    mychr = lmeresallp[[ii, 'CHR']]; mypos = lmeresallp[[ii, 'POS']]
    sourceid = which(p0[['CHR']] == mychr & p0[['POS']] == mypos)
    myp0 = p0[['p0']][sourceid]
    thisdeltap = subeset_by_id_year(deltap, samp_year, sourceid, myyear)
    pairid = find_similar_p0(p0, sourceid, backgroundid)
    pairdeltap = subeset_by_id_year(deltap, samp_year, pairid, myyear)
    thispos = format_chrpos(p0, sourceid)
    pairpos = format_chrpos(p0, pairid)
    return(list(thisdeltap, pairdeltap, thispos, pairpos))
}

get_mydata <- function(deltapp, metadty, envvar) {
    outdt = cbind(deltapp, metadty[,c('mergeid', 'site', envvar)])
    colnames(outdt)[1] = 'deltapp'
    return(outdt)
}

setup_traintest <- function(mydata) {
    mytraintest = lapply(base::split(mydata, mydata$site), function(xx) {
        testid = sample(1:nrow(xx), size = 1)
        testdt = xx[testid, ]
        traindt = xx[-testid, ]
        return(list(traindt, testdt))
    })
    mytrain = bind_rows(lapply(mytraintest, function(xx){xx[[1]]}))
    mytest = bind_rows(lapply(mytraintest, function(xx){xx[[2]]}))
    return(list(mytrain, mytest))
}

run_loo <- function(mydata, nrep = 200) {
    outdt = sapply(1:nrep, function(jj) {
        mytraintest = setup_traintest(mydata)
        mymodel = nlme::lme(fixed = stats::as.formula(fixed_effects), random = ~1|site, data = mytraintest[[1]], method = 'ML')
        mypre = predict(mymodel, newdata = mytraintest[[2]][,c('site', 'bio1')])
        # plot(mytraintest[[2]]$deltap, mypre)
        loostats = caret::postResample(pred = mypre, obs = mytraintest[[2]]$deltap)
        return(loostats)
    }) %>% t()
    return(outdt)
}

reformat_loores <- function(loores, mytype, mychrpos) {
    loores = loores %>%
        as.data.frame() %>%
        reshape2::melt(measure.vars = 1:3) %>% 
        dplyr::mutate(type = mytype,
                      pos = mychrpos) 
    return(loores)
}

run_loo_deltap <- function(deltapp, metadty, envvar, mytype, mychrpos) {
    # get mydata
    mydata = get_mydata(deltapp, metadty, envvar)
    loores = run_loo(mydata)
    loores = reformat_loores(loores, mytype, mychrpos)
    return(loores)
}

# def variables --------
myyear = '2018'
envvar = 'bio1'
fixed_effects = paste0('deltapp ~ ', envvar)

# load data --------
load(paste0('data/lmm_loo/bio1/lmeout_bio1_all_', myyear, '.rda'))
load('data/AF/merged_delta_p_776.rda')
load('data/AF/seedmix_p0_231.rda')
load('data/lmm_loo/info/mergeid_byyear.rda')
load('data/lmm_loo/info/weather_byyear.rda')

# set up background id
lmeresallp_bypos = data.table::copy(lmeresallp)
setorder(lmeresallp_bypos, CHR, POS)
all(p0[['CHR']]==lmeresallp_bypos[['CHR']], p0[['POS']]==lmeresallp_bypos[['POS']])
backgroundid = which(lmeresallp_bypos[['LR_nofix_p']] > 0.1)

# subset metadata 
metadty = metadt_year[[myyear]]

# main --------
lmeloo <- lapply(1:10, function(ii){
    # subset deltap
    deltapp = find_deltap_pair(lmeresallp, p0, deltap, backgroundid, samp_year, ii, myyear)
    print(deltapp[3:4])
    # significant results
    sigres = run_loo_deltap(deltapp[[1]], metadty, envvar, 'shift', deltapp[[3]]) 
    bakres = run_loo_deltap(deltapp[[2]], metadty, envvar, 'background', deltapp[[4]])
    outdt = rbind(sigres, bakres) %>% 
        dplyr::mutate(sigrank = ii)
    return(outdt)
})

lmeloodt = bind_rows(lmeloo)

# plot the comparisons
theme_set(cowplot::theme_cowplot())
pp <- ggplot(lmeloodt, aes(x = type, y = value)) +
    facet_grid(variable ~ sigrank, scales = 'free') + 
    geom_violin(alpha = 0.5, mapping = aes(fill = type, color = type)) +
    geom_boxplot(width = 0.2, mapping = aes(color = type)) +
    theme(legend.position = 'none')

# plot the second since it had the best performance
ii = 2
lmeresallp[ii, ]
deltapp = find_deltap_pair(lmeresallp, p0, deltap, backgroundid, samp_year, ii, myyear)
plotdt1 = get_mydata(deltapp[[1]], metadty, envvar) %>% 
    dplyr::mutate(chrpos = as.character(deltapp[[3]]))
plotdt2 = get_mydata(deltapp[[2]], metadty, envvar) %>% 
    dplyr::mutate(chrpos = deltapp[[4]])

pp2 <- ggplot(rbind(plotdt1, plotdt2), aes(x = bio1, ))


# output files --------
ggsave(filename = 'loo_bio1_2018.png', path = 'data/lmm_loo/bio1/', width = 10, height = 6)

# cleanup --------
date()
closeAllConnections()
