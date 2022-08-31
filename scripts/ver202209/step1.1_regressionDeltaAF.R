# Title: Run regression on all the deltadt/scalep files
# Author: Meixi Lin
# Date: Mon Aug 29 21:33:16 2022

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet")

library(data.table)
library(ggplot2)
library(dplyr)
library(cowplot)
date()
sessionInfo()

theme_set(theme_cowplot())

# def functions --------
plot_lm <- function(mycoords) {
    subdt = t(deltadt[mycoords,]) %>%
        as.data.frame(stringsAsFactor = FALSE) %>%
        tibble::rownames_to_column(var = 'siteday')
    forplot = reshape2::melt(data = subdt, id.vars = 'siteday', variable.factor = FALSE,
                               variable.name = 'POS', value.name = 'deltaP')
    forplot = cbind(forplot, reshape2::colsplit(forplot$siteday, '_', names = c('SITE','DAY'))) %>%
        dplyr::left_join(., y = metadata[,c('SITE',envvar)], by = 'SITE') %>%
        dplyr::mutate(YEAR = substr(DAY, 1, 4))

    pp <- ggplot(data = forplot,
                  mapping = aes(x = .data[[envvar]], y = deltaP)) +
        geom_point(mapping = aes(color = YEAR)) +
        ggpmisc::stat_poly_line() +
        ggpmisc::stat_poly_eq() +
        facet_wrap(. ~ POS, ncol = 4) +
        ggpubr::theme_pubr()
    return(pp)
}

plot_quadra <- function(mycoords) {
    subdt = t(deltadt[mycoords,]) %>%
        as.data.frame(stringsAsFactor = FALSE) %>%
        tibble::rownames_to_column(var = 'siteday')
    forplot = reshape2::melt(data = subdt, id.vars = 'siteday', variable.factor = FALSE,
                             variable.name = 'POS', value.name = 'deltaP')
    forplot = cbind(forplot, reshape2::colsplit(forplot$siteday, '_', names = c('SITE','DAY'))) %>%
        dplyr::left_join(., y = metadata[,c('SITE',envvar)], by = 'SITE') %>%
        dplyr::mutate(YEAR = substr(DAY, 1, 4))

    pp <- ggplot(data = forplot,
                 mapping = aes(x = .data[[envvar]], y = deltaP)) +
        geom_point(mapping = aes(color = YEAR)) +
        stat_smooth(method = "lm", formula = y ~ x + I(x^2)) +
        ggpmisc::stat_poly_eq(formula = y ~ x + I(x^2)) +
        facet_wrap(. ~ POS, ncol = 4) +
        ggpubr::theme_pubr()
    return(pp)
}

format_lmres <- function(lmres) {
    # format lmres
    Pcols = lmres %>% dplyr::select(starts_with('P_'))
    adjPcols = apply(Pcols, 2, function(xx) {
        padj = p.adjust(xx, method = 'bonferroni')
    })
    colnames(adjPcols) = paste0(colnames(adjPcols), '_adj')
    # only looking at the P_x_fdr so far
    outdt = cbind(lmres, adjPcols) %>%
        tibble::rownames_to_column(var = 'POS') %>%
        dplyr::mutate(POS = as.integer(POS),
                      P_pass = P_x_adj < 0.001) %>%
        dplyr::arrange(desc(P_pass),desc(AdjR2)) %>%
        tidyr::drop_na()
    return(outdt)
}

ggsaver <- function(pp, suffix, height, width) {
    ggsave(filename = paste0(prefix, envvar, suffix, '.pdf'), plot = pp, path = plotdir, height = height, width = width)
    ggsave(filename = paste0(prefix, envvar, suffix, '.png'), plot = pp, path = plotdir, height = height, width = width, units = 'in', dpi = 150)
}

get_prefix <- function(deltadtpath) {
    prefix = strsplit(deltadtpath,'/')[[1]]
    prefix = strsplit(prefix[length(prefix)],'_')[[1]]
    return(prefix[1])
}

# def variables --------
args = commandArgs(trailingOnly=TRUE)
deltadtpath = as.character(args[1])
envvar = as.character(args[2]) # environment variable to regress against
prefix = get_prefix(deltadtpath)
print(prefix)
# deltadtpath = 'data/AF/ver202209/haplotype/DeltaP/scaled_delta_freq_uniq.rds'
# envvar = 'bio1'
plotdir = paste0('plots/ver202209/AF/DeltaP/', envvar, '/')
dir.create(plotdir, recursive = TRUE)

# load data --------
# load deltadt
deltadt = readRDS(file = deltadtpath)
deltadt[1:5,1:5]
table(apply(deltadt, 1, function(xx){any(is.na(xx))})) # some CHR POS are NA

# load the gene coordinates
genecoords = readRDS(file = 'data/AF/ver202209/haplotype/DeltaP/genecoords.rds')

# load metadata
metadata = readRDS(file = 'data/metadata/ver202209/metacensusclim_0922.rds')
metadt = metadata[,c('SITE','DATE', envvar)] %>%
    dplyr::mutate(siteday = paste(SITE,DATE,sep='_')) %>%
    dplyr::distinct()

# main --------
# run a linear model on everything
lmres0 = apply(deltadt, 1, function(yy) {
    outdt = tryCatch({
        lmsum = summary(lm(yy ~ metadt[,envvar], na.action=na.exclude))
        outdt = c(lmsum$adj.r.squared, lmsum$coefficients[2,"Estimate"], lmsum$coefficients[2,"Pr(>|t|)"])
        return(outdt)
        }, error = function(cond){
            outdt = c(NA_real_,NA_real_,NA_real_)
            return(outdt)
        })
    return(outdt)
}) %>% t() %>% data.frame()
colnames(lmres0) = c('AdjR2','Estimate_x','P_x')
lmres = format_lmres(lmres0)
head(lmres)
table(lmres$P_pass) # number of genomic positions that passed the 0.001 threshold and bonferroni correction

# plot the manhattan plots
pp1 <- ggplot(data = lmres, aes(x = POS, y = -log10(P_x_adj), color = P_pass)) +
    geom_point() +
    scale_color_manual(values = c('darkgray','red')) +
    geom_hline(yintercept = -log10(0.001), color = 'gray', linetype = 'dashed') +
    labs(x = 'Chr1 genomic position (bp)', y = '-log10(P_adj)', title = envvar) +
    theme(legend.position = 'none')

# plot the example


# run a quadratic model
lmrsq2 = apply(deltadt, 1, function(yy) {
    rsq2 = tryCatch({
        data = data.frame(yy, metadt[,envvar], metadt[,envvar]^2)
        colnames(data) = c('yy','xx','xx2')
        rsq2 = summary(lm(yy ~ xx + xx2, data = data, na.action=na.exclude))$r.squared
    },
        error = function(cond){return(NA)})
    return(rsq2)
})

png(paste0(plotdir, 'HIST_lmrsq2.png'), height = 500, width = 500)
hist(lmrsq2, xlab = paste0('lmrsq: deltaP ~ ', envvar, ' + ', envvar, '^2'))
dev.off()

# print the highest and lowest values
length(lmrsq2)
lmrsq2 = sort(lmrsq2, na.last = NA)
length(lmrsq2)
head(lmrsq2)
tail(lmrsq2)

# plot the picked ones for linear model ========
pp1 = plot_lm(mycoords = names(head(lmrsq, n = 12)))
ggsaver(pp1,'_minLMx',7,10)
pp2 = plot_lm(mycoords = names(tail(lmrsq, n = 12)))
ggsaver(pp2,'_maxLMx',7,10)

# plot the picked ones for quadratic model ========
pp3 = plot_quadra(mycoords = names(head(lmrsq2, n = 12)))
ggsaver(pp3,'_minLMx2',7,10)
pp4 = plot_quadra(mycoords = names(tail(lmrsq2, n = 12)))
ggsaver(pp4,'_maxLMx2',7,10)


# output files --------
rm(deltadt) # don't save the deltadt
save.image(file = paste0(plotdir, envvar, '.RData'))

# cleanup --------
date()
closeAllConnections()
