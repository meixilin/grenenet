# Title: Run replicability analyses using the leave-one-out method
# Author: Meixi Lin
# Date: Thu Oct 13 08:46:32 2022

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet")

library(ggplot2)
library(dplyr)
library(cowplot)
library(lme4)
library(lmerTest)
library(MuMIn)
date()
sessionInfo()

theme_set(theme_classic(base_size = 12))

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

ggsaver <- function(pp, suffix, height, width) {
    ggsave(filename = paste0(suffix, '.pdf'), plot = pp, path = plotdir, height = height, width = width)
    ggsave(filename = paste0(suffix, '.png'), plot = pp, path = plotdir, height = height, width = width, units = 'in', dpi = 150)
}

# def variables --------
# args = commandArgs(trailingOnly=TRUE)
# envvar = as.character(args[1]) # environment variable to regress against
envvar = 'bio1'

plotdir = paste0('plots/ver202209/Replicability_lmm/', envvar, '/')
dir.create(plotdir, recursive = TRUE)

# these positions should not be input for linear models
badpos = c("8916","9094035","9527351","13914486","21327459","21368627")

# load data --------
# load deltadt
deltadt = readRDS(file = './data/AF/ver202209/haplotype/scaled_delta_freq_1012.rds')
deltadt = remove_badpos(deltadt, badpos) # remove the six bad lines
deltadt[1:5,1:5]

# load the gene coordinates
genecoordsdt = readRDS(file = 'data/AF/ver202209/haplotype/DeltaP/genecoordsdt.rds')

# load metadata
metadata = readRDS(file = 'data/metadata/ver202209/metacensusclim_0922.rds')
# `day` should be `doy`, means days of the year
metadt = metadata[,c('id','SITE','DATE', 'PLOT', 'year', 'day', envvar)] %>%
    dplyr::mutate(SITE_f = as.character(SITE)) %>%
    dplyr::group_by(SITE) %>%
    dplyr::mutate(PLOT2 = ifelse(PLOT == max(PLOT), 'Rep2','Rep1')) %>%
    dplyr::ungroup()

# double check that metadt is 1-1 match to colnames
all(metadt$id == colnames(deltadt))

str(metadt)

# main --------
# moi's previous formula
# formula= s ~ site + plot %in% site + doy + LATITUDE

# test run a linear model
rowid = which(rownames(deltadt) == '7554251')
mydata = cbind(unlist(deltadt[rowid,]), metadt) # genomic site one
colnames(mydata)[1] = 'yy'

lmemodel = lmerTest::lmer(formula = yy ~ bio1 + (1|year) + (1|SITE_f/PLOT2),
             data = mydata, na.action = na.exclude)
ranef(lmemodel, drop = TRUE) # extract the random effect
summary(lmemodel)
anova(lmemodel)
r.squaredGLMM(lmemodel)

# plot fitted lmemodel
plot(x = mydata$bio1, y = fitted(lmemodel))
points(x = mydata$bio1, y = mydata$yy, col = 'blue')

# output files --------

# cleanup --------
date()
closeAllConnections()
