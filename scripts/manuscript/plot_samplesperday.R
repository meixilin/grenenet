# Title: Plot Figure 1B with updated data
# Author: Meixi Lin
# Date: Thu Jan 26 13:12:28 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/central/groups/carnegie_poc/meixilin/grenenet/analyses/data/manuscript/")
dir.create('plot_samplesperday')

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
    dplyr::filter(usesample) 
samples_date = as.Date(samples_data$date, '%Y%m%d')

# main --------
# convert to the format week of the month and day of the week
# DayOfWeek: 1-based, 1 denotes Sunday
# GMT: enforce to 00:00:00 for the date time
samples_datedf = stringi::stri_datetime_fields(samples_date, tz = 'GMT') %>% 
    dplyr::group_by(Year, Month, Day, WeekOfMonth, DayOfWeek) %>% 
    dplyr::summarize(nsample = n()) %>%
    dplyr::filter(Year != 2021) 
samples_datedf$DayOfWeek = factor(samples_datedf$DayOfWeek, levels = rev(c(2,3,4,5,6,7,1)))
samples_datedf$Month = factor(samples_datedf$Month, 
                              labels = c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'))
pp <- ggplot(data = samples_datedf, aes(x = WeekOfMonth, y = DayOfWeek, fill = nsample)) +
    geom_tile() + 
    facet_grid(Year ~ Month) +
    scale_fill_gradientn(colours = RColorBrewer::brewer.pal(9,"Greens")[3:9]) + 
    scale_y_discrete(labels = rev(c('M','T','W','T','F','S','S'))) + 
    theme(strip.background = element_blank(), 
          axis.line=element_blank(),
          axis.ticks=element_blank(),
          axis.text.x = element_blank(),
          axis.title = element_blank(),
          strip.text = element_text(size = 15)) 
    
# output files --------
ggsave(filename = 'samplesperday.pdf', width = 10, height = 5, path = 'plot_samplesperday')
save.image(file ='plot_samplesperday/plotdata.rda')

# cleanup --------
date()
closeAllConnections()

