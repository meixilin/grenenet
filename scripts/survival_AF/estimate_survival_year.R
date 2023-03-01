# Title: Estimate the year of survival for all the samples by the fastq info
# Author: Meixi Lin
# Date: Wed Oct 12 21:08:12 2022
# Modification: Rerun the data
# Date: Thu Feb  9 10:02:03 2023


# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses")
date()
sessionInfo()

library(ggplot2)
library(dplyr)
library(diffdf)

library(maps)
theme_set(cowplot::theme_cowplot())
wrd <- map_data("world")

# def functions --------

# def variables --------

# load data --------
load(file = '../metadata/data/fastq_info.rda')
# load metadata
load(file = '../metadata/data/samples_data.rda')
load(file = '../metadata/data/locations_data.rda')

# main --------
# survival dt by all samples (filter usesample or no, the results are the same)
survivaldt = samples_data %>%
    dplyr::filter(usesample) %>%
    dplyr::group_by(site, plot) %>%
    dplyr::summarise(nyear = n_distinct(year))

# survival dt by the largest (only keep three years)
survivalsite = survivaldt %>%
    dplyr::group_by(site) %>%
    dplyr::summarise(nyear = max(nyear)) %>%
    dplyr::mutate(nyear = ifelse(nyear > 3, 3, nyear))

# plot
ggplot(data = survivaldt, aes(x = site, fill = factor(nyear), group = nyear)) +
    geom_bar(position="stack", stat="count") +
    labs(y = 'Number of Plots', fill = 'Year has fastq\n(Survival yr)')

# plot on the map
mapdt = dplyr::full_join(locations_data[,c('site','longitude','latitude', 'survivalyear')], survivalsite, by = 'site') %>%
    dplyr::mutate(survivalyear = ifelse(is.na(survivalyear), nyear, survivalyear))

pp <- ggplot() +
    geom_polygon(data = wrd, aes(x=long, y = lat, group = group), fill = 'gray80') +
    geom_point(data = mapdt, size = 2, aes(x = longitude, y = latitude, color = as.factor(survivalyear))) +
    coord_cartesian(xlim = range(mapdt$longitude), ylim = range(mapdt$latitude)) +
    scale_color_manual(values = RColorBrewer::brewer.pal(n = 4, name = 'Spectral')) +
    theme(legend.position = 'right')

# output files --------
dir.create('./data/metadata/survival_years/', recursive = TRUE)
write.csv(mapdt, file = './data/metadata/survival_years/survival_bysite.csv', quote = TRUE, row.names = FALSE)

# cleanup --------
date()
closeAllConnections()
