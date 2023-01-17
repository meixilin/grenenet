# Title: Run replicability analyses using the leave-one-out method
# Author: Meixi Lin
# Date: Thu Oct 13 08:46:32 2022
# Modification: Update to run the proper leave-one-out
# Date: Mon Dec 12 13:39:29 2022


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

ggsaver <- function(pp, suffix, height, width) {
    ggsave(filename = paste0(suffix, '.pdf'), plot = pp, path = plotdir, height = height, width = width)
    ggsave(filename = paste0(suffix, '.png'), plot = pp, path = plotdir, height = height, width = width, units = 'in', dpi = 150)
}

loadSomeRData <- function(x, file) {
    E <- new.env()
    load(file=file, envir=E)
    return(get(x, envir=E, inherits=F))
}


format_lmres <- function(lmres,genecoordsdt, p_adjmethods = 'bonferroni', p_cut = 0.05) {
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
                      P_pass = P_x_adj < p_cut) %>%
        dplyr::arrange(desc(P_pass),desc(Estimate_x))
    # find the gene coordinates
    outdt = outdt %>%
        dplyr::left_join(., genecoordsdt, by = 'POS')
    return(outdt)
}

run_lmm_D1 <- function(outsite) {
    tempmeta = metadt %>%
        dplyr::filter(SITE != outsite)
    biometa = metadt %>%
        dplyr::filter(SITE == outsite)
    tempyy = yy[names(yy) %in% tempmeta$id]

    # Account for experimental setup
    mydata = cbind(tempyy, tempmeta)
    lmemodel = suppressMessages(
        lmerTest::lmer(formula = tempyy ~ eval(as.name(envvar)) + (1|year) + (1|SITE_f/PLOT2),
                       data = mydata, na.action = na.exclude))
    # cannot keep year as a random effect because of singular fit
    # lmemodel = lmerTest::lmer(formula = yy ~ eval(as.name(envvar)) + year + (1|SITE_f/PLOT2),
    #                           data = mydata, na.action = na.exclude)
    lmesum = summary(lmemodel)
    lmepredict =
    lmer2 = MuMIn::r.squaredGLMM(lmemodel)
    # get relevant stats
    outdt = c(outsite,
              lmer2,
              lmesum$coefficients["eval(as.name(envvar))","Estimate"],
              log10(lmesum$coefficients["eval(as.name(envvar))","Pr(>|t|)"]),
              lmesum$AICtab)
    # keep the
    return(outdt)
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

# load the previous lmm result
lmres0.1 = loadSomeRData(x = 'lmres0', file = './plots/ver202209/LMM/bio1/part1/LMM_SSDAFbio1.RData')
lmres0.2 = loadSomeRData(x = 'lmres0', file = './plots/ver202209/LMM/bio1/part2/LMM_SSDAFbio1.RData')
lmres0 = rbind(lmres0.1,lmres0.2)

lmres = format_lmres(lmres0,genecoordsdt, p_adjmethods = 'fdr', p_cut = 0.1)
table(lmres$P_pass)
rm(lmres0.1,lmres0.2, lmres0)
lmresp = lmres %>%
    dplyr::filter(P_pass == TRUE) %>%
    dplyr::mutate(POS = as.character(POS),
                  log10P_x = log10(P_x))

# main --------
# for the top models, run the leave-one-out analyses
# run the lmm model with one site left out
deltadtp = deltadt[rownames(deltadt) %in% lmresp$POS, ]

for (ii in 1:nrow(lmresp)) {
    yy = deltadtp[lmresp$POS[ii],]

    # remove one site
    lmml = lapply(unique(metadt$SITE), run_lmm_D1)
    lmmD1 = lmml %>% dplyr::bind_rows()

    colnames(lmmD1) = c('outsite','R2_x','R2_all','Estimate_x','log10P_x','REML')

    forplot = lmmD1 %>%
        reshape2::melt(., id.vars = 'outsite')
    forplot2 = lmresp[ii, c(1, 2,3,4,10,6)] %>%
        reshape2::melt(., id.vars = 'POS')

    pp <- ggplot() +
        geom_violin(data = forplot, aes(x = variable, y = value)) +
        geom_point(data = forplot2, aes(x = variable, y = value), color = 'red', size = 2) +
        facet_wrap(. ~ variable, scales = 'free', nrow = 1)
    ggsave(filename = paste0('Replicability_lmm_POS', lmresp$POS[ii], '.png'),
           width = 8, height = 4, plot = pp, path = plotdir)
}


# output files --------

# cleanup --------
date()
closeAllConnections()
