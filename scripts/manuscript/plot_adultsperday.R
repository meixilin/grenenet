# Title: Plot Figure 1C with updated data
# Author: Meixi Lin
# Date: Thu Jan 26 16:11:01 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/central/groups/carnegie_poc/meixilin/grenenet/analyses/data/manuscript/")
dir.create('plot_adultsperday')

source('/central/groups/carnegie_poc/meixilin/grenenet/analyses/scripts/manuscript/config_manuscript.R')

library(stringi)
library(dplyr)
library(ggplot2)

# def functions --------

# def variables --------

# load data --------
load(paste0(metadatadir, 'samples_data.rda'))

# 2339 samples
samples_data = samples_data %>% 
    dplyr::filter(usesample) %>% 
    dplyr::filter(year != 2021) 
samples_date = as.Date(samples_data$date, '%Y%m%d')

# main --------
# convert to the format week of the month and day of the week
# DayOfWeek: 1-based, 1 denotes Sunday
# GMT: enforce to 00:00:00 for the date time
dayofyear = stringi::stri_datetime_fields(samples_date, tz = 'GMT')$DayOfYear

# merge the day of year by the number of flowers sampled
plotdf = cbind(samples_data, dayofyear)

pp <- ggplot(data = plotdf, aes(x = dayofyear, y = flowerscollected, color = flowerscollected)) +
    geom_point() + 
    scale_color_gradientn(colours = RColorBrewer::brewer.pal(9,"Greens")[3:9]) +
    # scale_y_continuous(expand = c(0,0)) + 
    labs(x = 'calendar day in a year')
    
# output files --------
ggsave(filename = 'adultsperday.pdf', width = 10, height = 2, path = 'plot_adultsperday/')
save.image(file ='plot_adultsperday/plotdata.rda')

# cleanup --------
date()
closeAllConnections()

