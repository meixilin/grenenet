# Title: Compare anova with linear models in performances
# Author: Meixi Lin
# Date: Tue Oct 11 11:18:33 2022

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet")
sink('./logs/DeltaAF_anova_lm_example_20221011.log')
library(data.table)
library(ggplot2)
library(dplyr)
library(cowplot)
library(qqman)
date()
sessionInfo()

theme_set(theme_cowplot())

# def functions --------
remove_badpos <- function(dt, badpos){
    print(dim(dt))
    finiteid = apply(dt, 1, function(xx){all(!is.finite(xx))})
    outdt = dt[!finiteid,]
    print(dim(outdt))
    if (length(intersect(badpos, rownames(outdt))) > 0) {
        stop('DT contains non-finite rows')
    }
    return(outdt)
}

lmsum_pval <- function(ff) {
    return(pf(ff[1],ff[2],ff[3],lower.tail = FALSE))
}

subset_deltadt <- function(deltadt, mycoords) {
    subdt = t(deltadt[mycoords,]) %>%
        as.data.frame(stringsAsFactor = FALSE) %>%
        tibble::rownames_to_column(var = 'siteday')
    forplot = reshape2::melt(data = subdt, id.vars = 'siteday',
                             variable.name = 'POS', value.name = 'deltaP')
    forplot = cbind(forplot, reshape2::colsplit(forplot$siteday, '_', names = c('SITE','DAY'))) %>%
        dplyr::left_join(., y = metadata[,c('SITE',envvar)], by = 'SITE') %>%
        dplyr::mutate(YEAR = substr(DAY, 1, 4),
                      POS = as.character(POS))
    return(forplot)
}

format_lmres <- function(lmres,genecoordsdt, p_adjmethods = 'bonferroni') {
    # format lmres
    Pcols = lmres %>% dplyr::select(starts_with('P_'))
    adjPcols = apply(Pcols, 2, function(xx) {
        padj = p.adjust(xx, method = p_adjmethods)
    })
    colnames(adjPcols) = paste0(colnames(adjPcols), '_adj')
    # only looking at the P_x_adj so far
    outdt = cbind(lmres, adjPcols) %>%
        tibble::rownames_to_column(var = 'POS') %>%
        dplyr::mutate(POS = as.integer(POS),
                      P_pass = P_x_adj < 0.001) %>%
        dplyr::arrange(desc(P_pass),desc(AdjR2))
    # find the gene coordinates
    outdt = outdt %>%
        dplyr::left_join(., genecoordsdt, by = 'POS')
    return(outdt)
}

gene_labeller <- function(forplot,lmres) {
    POSvalues <- unique(forplot$POS)
    if (length(POSvalues)!=2) {stop('Duplicated P_pass and POS')}
    dt <- lmres[lmres$POS %in% as.integer(POSvalues), ]
    POSlabs <- apply(dt, 1, function(xx){paste0('Chr1-', xx['POS'],'|Gene-',xx['GENES'])})
    names(POSlabs) <- !dt$P_pass
    return(POSlabs)
}

ggsaver <- function(pp, suffix, height, width) {
    ggsave(filename = paste0(prefix, envvar, suffix, '.pdf'), plot = pp, path = plotdir, height = height, width = width)
    ggsave(filename = paste0(prefix, envvar, suffix, '.png'), plot = pp, path = plotdir, height = height, width = width, units = 'in', dpi = 150)
}

# def variables --------
args = commandArgs(trailingOnly=TRUE)
envvar = as.character(args[1]) # environment variable to regress against
envvar = 'bio1'
# only use the scaled deltaP
deltadtpath = 'data/AF/ver202209/haplotype/DeltaP/scaled_delta_freq_uniq.rds'
prefix = 'scaled_deltaP_'

plotdir = paste0('plots/ver202209/AF/DeltaP_v3/', envvar, '/')
dir.create(plotdir, recursive = TRUE)

# these positions should not be input for linear models
badpos = c("8916","9094035","9527351","13914486","21327459","21368627")

# load data --------
# load deltadt
deltadt = readRDS(file = deltadtpath)
deltadt = remove_badpos(deltadt, badpos) # remove the six bad lines
deltadt[1:5,1:5]

# load the gene coordinates
genecoordsdt = readRDS(file = 'data/AF/ver202209/haplotype/DeltaP/genecoordsdt.rds')

# load metadata
metadata = readRDS(file = 'data/metadata/ver202209/metacensusclim_0922.rds')
# `day` should be `doy`, means days of the year
metadt = metadata[,c('SITE','DATE', 'year', 'day', envvar)] %>%
    dplyr::mutate(siteday = paste(SITE,DATE,sep='_')) %>%
                  # SITE = as.character(SITE)
    dplyr::distinct()

# double check that metadt is 1-1 match to colnames
all(metadt$siteday == colnames(deltadt))

# add a column where SITE are factors
metadt$SITE_f = factor(metadt$SITE)

# main --------
# moi's previous formula
# formula= s ~ site + plot %in% site + doy + LATITUDE

# test run a linear model
mydata = cbind(unlist(deltadt[1,]), metadt) # genomic site one
colnames(mydata)[1] = 'yy'

lmmodel = lm(formula = yy ~ SITE_f + year%in%SITE + day + bio1,
             data = mydata, na.action = na.exclude)
summary(lmmodel)
anova(lmmodel)

aovmodel = aov(formula = yy ~ SITE_f + year%in%SITE + day + bio1,
               data = mydata, na.action = na.exclude)
summary(aovmodel)

# an example on the importance of orders in anova models
aovmodel2 = aov(formula = yy ~ bio1 + SITE_f + year%in%SITE + day,
               data = mydata, na.action = na.exclude)
summary(aovmodel2)
aovmodel3 = aov(formula = yy ~ SITE_f + bio1 + year%in%SITE + day,
                data = mydata, na.action = na.exclude)
summary(aovmodel3)

# cleanup --------
date()
sink()
closeAllConnections()
