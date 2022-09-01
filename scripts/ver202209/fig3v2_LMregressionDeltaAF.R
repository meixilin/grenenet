# Title: Run regression on the scaled deltaP files
# Plot it for presentation. Only run linear model.
# Update to include the possible confounding factors
# Author: Meixi Lin
# Date: Wed Aug 31 00:37:43 2022

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet")

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

# only use the scaled deltaP
deltadtpath = 'data/AF/ver202209/haplotype/DeltaP/scaled_delta_freq_uniq.rds'
prefix = 'scaled_deltaP_'

# envvar = 'bio1'
plotdir = paste0('plots/ver202209/AF/DeltaP_v2/', envvar, '/')
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

# main --------
# run a linear model on everything
lmres0 = apply(deltadt, 1, function(yy) {
    # Account for experimental setup
    mydata = cbind(yy, metadt)
    # lmmodel = lm(yy ~ eval(as.name(envvar)) + year + day %in% SITE + SITE, data = mydata, na.action = na.exclude)
    lmmodel = lm(yy ~ eval(as.name(envvar)) + year%in%SITE + day%in%SITE + SITE, data = mydata, na.action = na.exclude)
    lmsum = summary(lmmodel)
    # get overall p value
    outdt = c(lmsum$adj.r.squared,
              lmsum$coefficients["eval(as.name(envvar))","Estimate"],
              lmsum$coefficients["eval(as.name(envvar))","Pr(>|t|)"],
              lmsum_pval(lmsum$fstatistic))
    return(outdt)
}) %>% t() %>% data.frame()
dim(lmres0)

colnames(lmres0) = c('AdjR2','Estimate_x','P_x','P_model')
lmres = format_lmres(lmres0,genecoordsdt)
dim(lmres)
head(lmres)
tail(lmres)
# number of genomic positions that passed the 0.001 threshold and bonferroni correction
table(lmres$P_pass)

# tally the number of genes with significant result
generes = lmres %>%
    dplyr::filter(P_pass == TRUE) %>%
    dplyr::group_by(GENES) %>%
    dplyr::summarise(P_pass_count = n()) %>%
    dplyr::arrange(desc(P_pass_count))
head(as.data.frame(generes))

# plot the manhattan plots
pcutoff = 0.001/nrow(deltadt)
pp1 <- ggplot(data = lmres, aes(x = POS, y = -log10(P_x), color = P_pass)) +
    geom_point() +
    scale_color_manual(values = c('darkgray','red')) +
    geom_hline(yintercept = -log10(pcutoff), color = 'gray', linetype = 'dashed') +
    labs(x = 'Chr1 genomic position (bp)', y = '-log10(P)', title = envvar) +
    theme(legend.position = 'none')

# plot the example plot of allele frequency shift
# TODO: just picked the ones with the highest R^2 and the lowest R^2
forplot1 = subset_deltadt(deltadt, as.character(lmres[1,'POS'])) %>% dplyr::mutate(P_notpass = FALSE)
forplot2 = subset_deltadt(deltadt, as.character(lmres[nrow(lmres),'POS'])) %>% dplyr::mutate(P_notpass = TRUE)
forplot = bind_rows(forplot1,forplot2)
# Add Gene info in the facet
pp2 <- ggplot(data = forplot,
             mapping = aes(x = .data[[envvar]], y = deltaP, color = !P_notpass)) +
    geom_point(mapping = aes(shape = YEAR), alpha = 0.1) +
    scale_shape_manual(values = c(1,2)) +
    scale_color_manual(values = c('darkgray','red'), guide = 'none') +
    facet_wrap(. ~ P_notpass, nrow = 1, labeller = labeller(P_notpass = gene_labeller(forplot,lmres))) +
    labs(y = "\u0394p/(p0(1-p0))") +
    stat_smooth(method = 'lm', formula = y ~ x) +
    ggpmisc::stat_poly_eq() +
    theme(legend.position = 'top')

# output the figures
ggsaver(pp1, '_manhanttan',4,8)
options(warn = -1)
ggsaver(pp2, '_MaxMinR2',4,8) # pdf will have warnings
options(warn = 0)

# TODO: beware of the high P-value
png(filename = paste0(plotdir, prefix, envvar, '_QQplot.png'), width = 800, height = 400)
par(mfrow = c(1,2))
qqman::qq(lmres$P_x, xlim = c(0,6), ylim = c(0,12))
qqman::qq(lmres$P_x_adj, xlim = c(0,6), ylim = c(0,12))
dev.off()

# output files --------
rm(deltadt, forplot1, forplot2, genecoordsdt) # don't save the deltadt
write.csv(generes, file = paste0(plotdir, prefix, envvar, '_genes.csv'))
write.csv(lmres[lmres$P_pass, ], file = paste0(plotdir, prefix, envvar, '_lmres_Ppass.csv'))
save.image(file = paste0(plotdir, prefix, envvar, '.RData'))

# cleanup --------
date()
closeAllConnections()
