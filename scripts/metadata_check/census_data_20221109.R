# Title:
# Author: Meixi Lin
# Date: Wed Nov  9 16:30:22 2022

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
load('/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/ath_evo/grenephase1/data/locations_data.rda')

load(file = '/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/ath_evo/grenephase1/data/recordssorted.rda')
samples_dt = read.csv(file = '/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/ath_evo/grenephase1/data/samples_test.csv') %>%
    dplyr::mutate(year = substr(DATE, 1, 4)) %>%
    dplyr::mutate(site = stringr::str_pad(SITE, 2, side = 'left', pad = '0'))
census_dt = read.csv(file = '/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/ath_evo/grenephase1/data/census_test.csv') %>%
    dplyr::mutate(year = substr(date, 1, 4)) %>%
    dplyr::mutate(site = stringr::str_pad(site, 2, side = 'left', pad = '0'))
ibuttons_dt = read.csv(file = '/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/ath_evo/grenephase1/data/ibuttons_info.csv')
str(census_dt)
census_summary = census_dt %>%
    dplyr::group_by(site, plot, year) %>%
    dplyr::summarise(n_dates = n())

census_tally2 = census_dt %>%
    dplyr::group_by(site, plot, year) %>%
    dplyr::summarise(n_dates = n() > 0) %>%
    dplyr::group_by(site, plot) %>%
    dplyr::summarise(n_years = sum(n_dates))
census_tally2 =dplyr::full_join(census_tally2, sites)

pp <- ggplot(data = census_tally2, mapping = aes(x = plot, y = site, fill = factor(n_years))) +
    geom_tile() +
    scale_fill_brewer(palette = 'Greens') +
    labs(fill = '# of years') +
    theme_bw(base_size = 14)

ggsave(filename = 'census_byyear.png', plot = pp, path = './plots/SCRATCH/', width = 5, height = 5)

ibuttons_tally2 = ibuttons_dt %>%
    dplyr::mutate(year = substr(date, 1, 4)) %>%
    dplyr::mutate(site = stringr::str_pad(site, 2, side = 'left', pad = '0')) %>%
    dplyr::group_by(site, year) %>%
    dplyr::summarise(n_dates = n())

pp <- ggplot(data = ibuttons_tally2, mapping = aes(x = year, y = site, fill = factor(n_dates))) +
    geom_tile() +
    scale_fill_brewer(palette = 'Greens') +
    labs(fill = '# of entries') +
    theme_bw(base_size = 14)

ggsave(filename = 'ibuttons_byyear.png', plot = pp, path = './plots/SCRATCH/', width = 5, height = 5)


# use the sample or census data


# main --------

# output files --------

# cleanup --------
date()
closeAllConnections()
