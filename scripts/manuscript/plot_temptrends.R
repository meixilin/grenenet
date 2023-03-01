# Title: 
# Author: Meixi Lin
# Date: Thu Feb  9 10:23:37 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses/data/manuscript/")
dir.create('plot_temptrends')

source('/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses/scripts/manuscript/config_manuscript.R')

library(ggplot2)

date()
sessionInfo()

# def functions --------

# def variables --------

# load data --------
# load weatherstation info
load('../../../metadata/data/weatherstation_info.rda')
load('../../../metadata/data/weatherstation_data.rda')
load('../../../metadata/data/locations_data.rda')
weatherstation_data = weatherstation_data %>% 
    dplyr::mutate(weatherdate = as.POSIXct(weatherdate, format="%Y-%m-%d"))

ibuttons = read.csv(file = '../../../metadata/data/temp_ib_data.csv', row.names = 1) %>% 
    dplyr::mutate(datetime_right = as.POSIXct(datetime_right, format="%Y-%m-%d %H")) 

year3 = weatherstation_data %>% dplyr::filter(site == 46) # this is germany data. distance to site = 46 km
year1 = weatherstation_data %>% dplyr::filter(site == 26) # this is israel data. distance to site = 43 km
forplot = rbind(year3, year1)
year3ib = ibuttons %>% dplyr::filter(site == 46)
year1ib = ibuttons %>% dplyr::filter(site == 26)


# main --------
pp1 <- ggplot() +
    geom_line(data = year3, mapping = aes(x = weatherdate, y = temp), color = '#2B83BA') + 
    geom_line(data = year1, mapping = aes(x = weatherdate, y = temp), color = '#FDAE61') + 
    geom_line(data = year3ib, mapping = aes(x = datetime_right, y = value), color = '#2B83BA', alpha = 0.2) + 
    geom_line(data = year1ib, mapping = aes(x = datetime_right, y = value), color = '#FDAE61', alpha = 0.2) + 
    scale_x_datetime(date_breaks = '1 year') + 
    labs(x = 'Year', y = 'Annual temperature (ºC)')

pp2 <- ggplot() +
    geom_col(data = forplot, mapping = aes(x = weatherdate, y = prcp, fill = as.factor(site))) + 
    facet_wrap(. ~ site, ncol = 1) +
    scale_x_datetime(date_breaks = '1 year') +
    scale_fill_manual(values = c('#FDAE61','#2B83BA')) + 
    theme(legend.position = 'none') + 
    labs(x = 'Year', y = 'Annual precipitation (mm)')
    
# output files --------
ggsave(filename = 'temperature.pdf', width = 6, height = 3, plot = pp1, path = 'plot_temptrends/')
ggsave(filename = 'precipitation.pdf', width = 6, height = 4, plot = pp2, path = 'plot_temptrends/')
save.image(file ='plot_temptrends/plotdata.rda')

# cleanup --------
date()
closeAllConnections()

