# Title: 
# Author: Meixi Lin
# Date: Thu Feb  9 09:59:49 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses/data/manuscript/")

dir.create('plot_worldclim')

source('/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses/scripts/manuscript/config_manuscript.R')

library(dplyr)
library(ggplot2)

date()
sessionInfo()

# def functions --------

# def variables --------

# load data --------
load('../../../metadata/data/worldclim_sitesdata.rda')
load('../../../metadata/data/worldclim_ecotypesdata.rda')
load('../../../metadata/data/ecotypes_data.rda')
survivaldt = read.csv('../metadata/survival_years/survival_bysite.csv')

# main --------
forplot1 = dplyr::left_join(worldclim_sitesdata, survivaldt, by = 'site')
hull_forplot1 <- forplot1 %>%
    group_by(survivalyear) %>%
    slice(chull(bio12, bio1))

hull_ecotypes <- worldclim_ecotypesdata %>%
    dplyr::select(bio12, bio1) %>% 
    dplyr::mutate(ecogroup = 'ecotype') %>%
    tidyr::drop_na() %>%
    group_by(ecogroup) %>%
    slice(chull(bio12, bio1))
 
pp <- ggplot() +
    geom_polygon(data = hull_ecotypes, 
                 mapping = aes(x = bio12, y = bio1), 
                 alpha = 0.2, fill = 'gray') +
    geom_point(data = worldclim_ecotypesdata, mapping = aes(x = bio12, y = bio1), pch = 3) +
    geom_polygon(data = hull_forplot1, 
                 mapping = aes(x = bio12, y = bio1, fill = as.factor(survivalyear)), 
                 alpha = 0.2) +
    geom_text(data = forplot1, mapping = aes(x = bio12, y = bio1, color = as.factor(survivalyear), label = site), size = 3) +
    labs(x = 'Annual precipitation (mm)', y = 'Annual temperature (ºC)') +
    scale_color_manual(values = RColorBrewer::brewer.pal(n = 4, name = 'Spectral')) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(n = 4, name = 'Spectral')) +
    theme(legend.position = 'none')

# output files --------
ggsave(filename = 'worldclim.pdf', width = 5, height = 5, path = 'plot_worldclim/')
save.image(file ='plot_worldclim/plotdata.rda')

# cleanup --------
date()
closeAllConnections()
