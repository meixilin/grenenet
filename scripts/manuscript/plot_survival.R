# Title: 
# Author: Meixi Lin
# Date: Thu Feb  9 10:23:37 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses/data/manuscript/")
dir.create('plot_survival')

source('/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses/scripts/manuscript/config_manuscript.R')

library(rgdal)
library(raster)
library(ggplot2)

date()
sessionInfo()

# def functions --------

# def variables --------

# load data --------
# Defines the x axes required
x_lines <- seq(-120,180, by = 60)
data("wrld_simpl", package = "maptools")

# load survival by site
survivaldt = read.csv('../metadata/survival_years/survival_bysite.csv')

# load ecotypes
load('../../../metadata/data/ecotypes_data.rda')

# main --------
pp <- ggplot() +
    geom_polygon(data = wrld_simpl, aes(x = long, y = lat, group = group), fill = "gray", colour = "gray") +
    geom_point(data = ecotypes_data, mapping = aes(x = longitude, y = latitude), pch = 3, color = 'black') +
    geom_point(data = survivaldt, mapping = aes(x = longitude, y = latitude, color = as.factor(survivalyear)), size = 3) +
    # Convert to polar coordinates
    coord_map("ortho", orientation = c(90, 0, 0)) +
    scale_y_continuous(breaks = seq(0, 90, by = 30), labels = NULL) +
    scale_color_manual(values = RColorBrewer::brewer.pal(n = 4, name = 'Spectral')) +
    # Removes Axes and labels
    scale_x_continuous(breaks = NULL) +
    labs(x = '', y = '',color = '') + 
    theme(panel.background = element_blank(),
          panel.grid.major = element_line(size = 0.25, linetype = 'dashed',
                                          colour = "black"),
          axis.ticks=element_blank(),
          axis.line = element_blank())

# output files --------
ggsave(filename = 'survival.pdf', width = 5, height = 5, path = 'plot_survival/')
save.image(file ='plot_survival/plotdata.rda')

# cleanup --------
date()
closeAllConnections()

