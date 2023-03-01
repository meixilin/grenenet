# Title: Plot number of flowers by site temperature
# Author: Meixi Lin
# Date: Tue Feb 28 13:53:11 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses")
date()
sessionInfo()

library(ggplot2)
library(dplyr)

source('./scripts/manuscript/config_manuscript.R')

# def functions --------

# def variables --------

# load data --------
load(file = '../metadata/data/samples_data.rda')
# only keep data with `usesample == TRUE`
samples_data = samples_data %>% 
    dplyr::filter(usesample)

load(file = '../metadata/data/worldclim_sitesdata.rda')

# sort the locations by annual temperature
samples_data = dplyr::left_join(samples_data, worldclim_sitesdata, by = 'site')

# main --------
# get total number of flowers
plotdt = samples_data %>%
    dplyr::group_by(site, plot, generation, bio1) %>% 
    dplyr::summarise(totalflowers = sum(flowerscollected)) %>% 
    dplyr::arrange(bio1) %>% 
    dplyr::mutate(site_bio1 = paste(site, format(bio1, digits = 3), sep = ': '))
# order the locations by annual temperature
plotdt$site_bio1 = factor(plotdt$site_bio1, levels = unique(plotdt$site_bio1))

pp1 <- ggplot(data = plotdt, aes(x = generation, y = totalflowers, color = plot, group = plot)) + 
    geom_line() +
    geom_point() + 
    facet_wrap(. ~ site_bio1, ncol = 6)

# use examples with at least two years of data ========
sites_remove = c(23, 26, 33, 37, 57, 60)

plotdt2 = plotdt %>% 
    dplyr::filter(!(site %in% sites_remove)) %>% 
    dplyr::group_by(site, generation, bio1) %>% 
    dplyr::summarise(meanflowers = mean(totalflowers),
                     minflowers = mean(totalflowers) - sd(totalflowers),
                     maxflowers = mean(totalflowers) + sd(totalflowers)) %>% 
    dplyr::mutate(tempcat = case_when(bio1 < 8 ~ '[5,8)',
                                      bio1 >= 8 & bio1 < 10 ~ '[8,10)',
                                      bio1 >= 10 & bio1 < 15 ~ '[10,15)',
                                      TRUE ~ '[15,20)'))
plotdt2$tempcat = factor(plotdt2$tempcat, c('[5,8)', '[8,10)', '[10,15)', '[15,20)'))

pp2 <- ggplot(data = plotdt2, aes(x = generation, y = meanflowers, ymin = minflowers, ymax = maxflowers, 
                                 color = bio1, group = site)) + 
    # geom_errorbar(width = 0.2) +
    geom_line(linetype = 'dashed') +
    scale_color_viridis_c() +
    geom_point() + 
    facet_wrap(. ~ tempcat, nrow = 1) +
    scale_x_continuous(breaks = c(1, 2, 3))

# output files --------
outdir = 'data/manuscript/plot_nflowers_site/'
dir.create(outdir)

ggsave(filename = 'nflowers_siteALL.pdf', width = 8, height = 8, plot = pp1, path = outdir)
ggsave(filename = 'nflowers_site_summary.pdf', width = 8, height = 4, plot = pp2, path = outdir)

save.image(file = paste0(outdir, 'plotdata.rda'))

# cleanup --------
date()
closeAllConnections()
