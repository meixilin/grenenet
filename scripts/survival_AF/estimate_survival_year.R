# Title: Estimate the year of survival for all the samples by the fastq info
# Author: Meixi Lin
# Date: Wed Oct 12 21:08:12 2022

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet")
date()
sessionInfo()

library(ggplot2)
library(dplyr)
# def functions --------

# def variables --------

# load data --------
load(file = '../../ath_evo/grenephase1/data/fastq_info.rda')
# load metadata
metadata = readRDS(file = 'data/metadata/ver202209/metacensusclim_0922.rds')

# main --------
sampledt = data.frame(SITE = substr(fastq_info$sampleid,5,6),
                      PLOT = substr(fastq_info$sampleid,7,8),
                      DATE = substr(fastq_info$sampleid,9,16),
                      year = substr(fastq_info$sampleid,9,12))
# survival dt by all samples
survivaldt = sampledt %>%
    dplyr::group_by(SITE, PLOT) %>%
    dplyr::summarise(nyear = n_distinct(year))

# survival dt by the largest
survivalsite = survivaldt %>%
    dplyr::group_by(SITE) %>%
    dplyr::summarise(nyear = max(nyear))

# plot
ggplot(data = survivaldt, aes(x = SITE, fill = factor(nyear), group = nyear)) +
    geom_bar(position="stack", stat="count") +
    labs(y = 'Number of Plots', fill = 'Year has fastq\n(Survival yr)')

# plot on the map
survivalsite$SITE = as.integer(survivalsite$SITE)
mapdt = dplyr::left_join(metadata[,c('SITE','LONGITUDE','LATITUDE')], survivalsite, by = 'SITE') %>%
    dplyr::distinct()
n_distinct(metadata$SITE)

library(maps)
theme_set(cowplot::theme_cowplot())
wrd <- map_data("world")

pp <- ggplot() +
    geom_polygon(data = wrd, aes(x=long, y = lat, group = group), fill = 'gray80') +
    geom_point(data = mapdt, size = 2, aes(x = LONGITUDE, y = LATITUDE, color = as.factor(nyear))) +
    coord_cartesian(xlim = range(mapdt$LONGITUDE), ylim = range(mapdt$LATITUDE)) +
    scale_color_manual(values = RColorBrewer::brewer.pal(n = 11, name = 'Spectral')[c(1,4,8,11)]) +
    theme(legend.position = 'right')

# output files --------
write.csv(survivaldt, file = './data/metadata/ver202209/survival_byplot.csv')
write.csv(survivalsite, file = './data/metadata/ver202209/survival_bysite.csv')

# cleanup --------
date()
closeAllConnections()
