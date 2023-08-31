# Title: Calculate per ecotype strength of selection (Vs) 
# Author: Meixi Lin
# Date: Tue Aug 29 10:17:04 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses")

library(dplyr)

date()
sessionInfo()

# load data --------
load('./data/local_adaptation/inputs/ecotypefreqs_deltaenv.rda')

# def functions --------
# prepare data input (different from ecotype_Vs workflow)
format_lmtables <- function(ecop_l, ecop0, deltaenv, mysite, myenv, mygen) {
    tmpdeltaenv = deltaenv[deltaenv$site == mysite, c('ecotypeid', 'site', myenv)]
    tmpeco = ecop_l %>% 
        dplyr::filter(site == mysite, generation == mygen) %>% 
        dplyr::left_join(., ecop0, by = 'ecotypeid') %>% 
        dplyr::mutate(log_p1_p0 = log(frequency/p0)) %>% 
        dplyr::mutate(log_p1_p0 = ifelse(is.infinite(log_p1_p0), NA, log_p1_p0)) %>% # if frequency = 0 for this ecotype
        dplyr::left_join(.,  tmpdeltaenv, by = c('site', 'ecotypeid'))
    return(tmpeco)
}

null_lmsum <- function(tmpeco) {
    outdf = data.frame(matrix(ncol = 7))
    colnames(outdf) = c('site', 'envvar', 'generation', 'adj_r2', 'b_env', 'b_p', 'df')
    outdf$site <- unique(tmpeco$site)
    outdf$envvar <- colnames(tmpeco)[ncol(tmpeco)]
    outdf$generation <- unique(tmpeco$generation)
    return(outdf)
}

format_lmsum <- function(tmpeco, lmsum, myenv) {
    outdf = null_lmsum(tmpeco)
    # add info on linear model results
    outdf$adj_r2 <- lmsum$adj.r.squared
    outdf$b_env <- lmsum$coefficients[myenv, 'Estimate']
    outdf$b_p <- lmsum$coefficients[myenv, 'Pr(>|t|)']
    outdf$df <- lmsum$df[2]
    return(outdf)
}

run_lm <- function(tmpeco, myenv) {
    # if no environment associated with this ecotype
    if (all(is.na(tmpeco[, myenv]))) {
        lmres = null_lmsum(tmpeco)
    } else {
        # run linear regression
        mymodel = lm(formula = as.formula(paste0('log_p1_p0 ~ ', myenv)), data = tmpeco)
        lmres = format_lmsum(tmpeco, summary(mymodel), myenv)
    }
    return(lmres)
}


# def variables --------
envvars = colnames(deltaenv)[-c(1,2)]
generations = 1:3
sites = unique(deltaenv$site)

outdir = './data/local_adaptation/persite/'
dir.create(outdir, recursive = TRUE)
outdir1 = paste0(outdir, 'env_csv/')
dir.create(outdir1, recursive = TRUE)

# main --------
# run linear model for each ecotype, each bioclim var in each generation
# last parameter is `bioall` the euclidean distance on PCA scores
alldf = lapply(envvars, function(myenv) {
    envdf = lapply(generations, function(mygen) {
        lapply(1:length(sites), function(ii) {
            mysite = sites[ii]
            # print(paste(mysite, myenv, mygen))
            tmpeco = format_lmtables(ecop_l, ecop0, deltaenv, mysite, myenv, mygen)
            # some site/generation pairs do not have samples 
            if (nrow(tmpeco) == 0) {
                return(NULL)
            } else {
                lmout = run_lm(tmpeco, myenv)
                return(lmout)
            }
        }) %>% dplyr::bind_rows()
    }) %>% dplyr::bind_rows() 
    # add more information
    envdf = envdf %>%
        dplyr::mutate(tau = -1/b_env)
    write.csv(envdf, file = paste0(outdir1, 'la_lmres_',myenv,'_persite.csv'))
    return(envdf)
}) %>% dplyr::bind_rows()  

# output files --------
save(alldf, file = paste0(outdir, 'la_lmres_persite.rda'))
save.image(paste0(outdir, 'step2_site_tau.RData'))

# cleanup --------
date()
closeAllConnections()
