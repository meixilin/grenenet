# Title: Run regression on all the deltadt files
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

date()
sessionInfo()

# def functions --------

# def variables --------

# load data --------
# get the coordinates of the allele frequency (should be only the chr1)
afcoords = readRDS(file = 'data/AF/ver202209/haplotype/SNPcoords_922975.rds')
deltadt = readRDS(file = 'data/AF/ver202209/haplotype/DeltaP/delta_freq.rds')
# genes = readRDS(file = 'data/TAIR10/Chr1_genes_gff.rds')

# # subset by the genes
# afcoordsi = as.integer(afcoords)
# genecoords = apply(genes, 1, function(xx) {
#     start = as.integer(xx['start'])
#     end = as.integer(xx['end'])
#     outids = afcoords[afcoordsi >= start & afcoordsi <= end]
#     return(outids)
# })
# names(genecoords) = genes$ID
# saveRDS(genecoords, file = 'data/AF/ver202209/haplotype/DeltaP/genecoords.rds')
genecoords, file = 'data/AF/ver202209/haplotype/DeltaP/genecoords.rds')
# load metadata
metadata = readRDS(file = 'data/metadata/ver202209/metacensusclim_0922.rds')

# main --------
# sample
colnames(deltadt) = c('siteday',afcoords)

# run a linear model on everything
deltadt1 = deltadt[,2:100]

out1 = apply(deltadt[,2:10000], 2, function(yy) {
    rsq = summary(lm(unlist(yy) ~ metadt$LATITUDE))$r.squared
    return(rsq)
})
hist(out1)
max(out1,na.rm=TRUE)

names(out) = colnames(deltadt1)

out = apply(deltadt[,2:10000], 2, function(yy) {
    data = data.frame(yy, metadt$LATITUDE, metadt$LATITUDE^2)
    colnames(data) = c('yy','xx','xx2')
    rsq2 = summary(lm(yy ~ xx + xx2, data = data))$r.squared
    return(rsq2)
})
hist(out)

mycoord = c('siteday',names(which(out==max(out, na.rm = TRUE))))
forplot = data.table::melt(data = deltadt[,..mycoords], id.vars = 'siteday', variable.factor = FALSE,
                           variable.name = 'POS', value.name = 'deltaP')
forplot = cbind(forplot, reshape2::colsplit(forplot$siteday, '_', names = c('SITE','DAY'))) %>%
    dplyr::left_join(., y = metadata[,c('SITE','LATITUDE','bio1')], by = 'SITE') %>%
    dplyr::mutate(YEAR = substr(DAY, 1, 4),
                  POS = as.integer(POS))

pp <- ggplot(data = forplot,
             mapping = aes(x = LATITUDE, y = deltaP)) +
    geom_point(mapping = aes(color = YEAR)) +
    ggpmisc::stat_poly_line() +
    ggpmisc::stat_poly_eq() +
    labs(title = mycoord[2], subtitle = 'AT1G01390: PIP5K') +
    # geom_violin() +
    ggpubr::theme_pubr()

names(out) = colnames(deltadt1)






















outdt = dplyr::bind_rows(out, .id = 'POS')

# get the metadatadt
metadt = metadata %>%
    select(SITE,DATE,starts_with('bio'),'LATITUDE') %>%
    dplyr::mutate(siteday = paste(SITE,DATE,sep='_')) %>%
    dplyr::distinct()
all.equal(deltadt$siteday, metadt$siteday)

# sample the first 10
mycoords = c('siteday',genecoords[[mygene]][1:10])
forplot = data.table::melt(data = deltadt[,..mycoords], id.vars = 'siteday', variable.factor = FALSE,
                           variable.name = 'POS', value.name = 'deltaP')
# plot by latitude
forplot = cbind(forplot, reshape2::colsplit(forplot$siteday, '_', names = c('SITE','DAY'))) %>%
    dplyr::left_join(., y = metadata[,c('SITE','LATITUDE','bio1')], by = 'SITE') %>%
    dplyr::mutate(YEAR = substr(DAY, 1, 4),
                  POS = as.integer(POS))

pp <- ggplot(data = forplot[forplot$YEAR == '2018', ],
             mapping = aes(x = LATITUDE, y = deltaP, color = SITE)) +
    geom_point() +
    facet_wrap(. ~ POS) +
    theme_bw() +
    theme(legend.position = 'none')

outdt = scalep[mycoords[-1],] %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = 'POS')
colnames(outdt)[-1] = deltadt$siteday
forplot = data.table::melt(data = outdt, id.vars = 'POS', variable.factor = FALSE,
                           variable.name = 'siteday', value.name = 'scaleP')
forplot = cbind(forplot, reshape2::colsplit(forplot$siteday, '_', names = c('SITE','DAY'))) %>%
    dplyr::left_join(., y = metadata[,c('SITE','LATITUDE','bio1')], by = 'SITE') %>%
    dplyr::mutate(YEAR = substr(DAY, 1, 4),
                  POS = as.integer(POS))
forplot2 = forplot %>%
    dplyr::group_by(siteday,bio1, YEAR) %>%
    dplyr::summarise(meanP = median(scaleP),
                     nsample = n())

# plot the regions of allele frequency changes
pp <- ggplot(data = forplot[forplot$YEAR == '2019', ],
             mapping = aes(x = LATITUDE, y = scaleP,group = SITE)) +
    # geom_point() +
    geom_violin() +
    theme_bw() +
    theme(legend.position = 'none')

pp <- ggplot() +
    geom_point(data = forplot2, mapping = aes(x = bio1, y = abs(meanP), color = YEAR)) +
    geom_line() +
    theme_bw()
ggsave(filename = 'test.png', plot = pp, height = 4, width = 10)

# output files --------
saveRDS(afcoords, file = 'data/AF/ver202209/haplotype/SNPcoords_922975.RDS')

# cleanup --------
date()
closeAllConnections()
