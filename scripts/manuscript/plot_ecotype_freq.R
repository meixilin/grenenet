# Title: plot_ecotype_freq geom_area
# Author: Meixi Lin
# Date: Mon Feb  5 08:32:00 2024

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses")

library(dplyr)
library(ggplot2)
library(randomcoloR)

date()
sessionInfo()

# def functions --------
source('scripts/manuscript/config_manuscript.R')

ggsaver <- function(myname, mypp, mywidth, myheight, mypath) {
    ggsave(filename = paste0(myname, '.pdf'), plot = mypp, width = mywidth, height = mywidth, path = mypath)
    ggsave(filename = paste0(myname, '.png'), plot = mypp, width = mywidth, height = mywidth, path = mypath, dpi = 72)
}

add_p0 <- function(df, myp0) {
    # fill in all the rows. if there is 12 plots x 231 ecotypes = 2772 rows for p0df
    p0df = dplyr::full_join(df %>% dplyr::distinct(ecotypeid, site, plot), myp0, by = 'ecotypeid') %>% 
        dplyr::mutate(generation = 0,
                      sample_name = paste(site, generation, plot, sep = '_')) %>%
        dplyr::rename(frequency = p0)
    # p0df should be n times 231
    stopifnot(nrow(p0df) %% 231 == 0)
    outdf = dplyr::bind_rows(df, p0df)
    return(outdf)
}

site_plot <- function(mysite) {
    xx = ecop_l231 %>% dplyr::filter(site == mysite)
    xx = add_p0(xx, ecop0_231)
    pp <-  ggplot(xx, 
                  aes(x = as.character(generation), y = frequency, fill = as.factor(ecotypeid), group = as.factor(ecotypeid))) + 
        geom_area() + 
        facet_wrap(. ~ plot) +
        labs(x = 'generation', title = paste0('site', mysite)) + 
        scale_x_discrete(expand = c(0,0)) + 
        scale_fill_manual(values = mycolors) + 
        theme(legend.position = 'none')
    ggsaver(paste0('ecotype_freq_site', mysite), pp, 8, 8, paste0(plotdir, 'persite/'))
    
    # print top ecotypes
    xx %>% 
        dplyr::filter(generation == 3) %>% 
        dplyr::group_by(plot) %>% 
        dplyr::top_n(., 3, frequency) %>%
        as.data.frame()
    return(pp)
}

# def variables --------
plotdir = 'data/manuscript/plot_ecotype_freq/'
dir.create(plotdir)
dir.create(paste0(plotdir, 'persite/'))

# load data --------
# the experiments ecotype frequencies
rdata = loadSomeRData(c('ecop_l231', 'ecop0_231'), './data/local_adaptation/inputs/step0_prepare_data.RData')
ecop_l231 = rdata$ecop_l231; ecop0_231 = rdata$ecop0_231
# set color scheme (for now this is not reproducible because of randomcoloR but the color is saved as mycolors in one run)
mycolors = randomcoloR::randomColor(count = 231); names(mycolors) = ecop0_231$ecotypeid

# main --------
# individual ecotype plots
pplist <- lapply(unique(ecop_l231$site), site_plot)

# plot allsite summary
ecop_mean = ecop_l231 %>% 
    dplyr::group_by(site, generation, ecotypeid) %>% 
    dplyr::summarise(frequency = mean(frequency),
                     nplot = n(),
                     .groups = 'drop') %>% 
    dplyr::mutate(plot = 'all',
                  sample_name = paste(site, generation, plot, sep = '_'))
ecop_meanp0 = add_p0(ecop_mean, ecop0_231)

# check that the above mean sums up to 1 each generation at each site
ecop_mean_checker = ecop_mean %>% 
    dplyr::group_by(site, generation) %>% 
    dplyr::summarise(sumfreq = sum(frequency))
summary(ecop_mean_checker$sumfreq) # PASSED! All 1. 

# now plot allsite
ppall <- ggplot(ecop_meanp0, 
                aes(x = as.character(generation), y = frequency, fill = as.factor(ecotypeid), group = as.factor(ecotypeid))) + 
    geom_area() + 
    facet_wrap(. ~ site) +
    labs(x = 'generation', title = 'allsite averaged plots') + 
    scale_x_discrete(expand = c(0,0)) + 
    scale_fill_manual(values = mycolors) + 
    theme(legend.position = 'none')
ggsaver('ecotype_freq_allsite', ppall, 8, 8, plotdir)

# print the site-generation pairs that had less than 5 plots averaged
ecop_mean %>% 
    dplyr::distinct(site, generation, nplot) %>% 
    dplyr::filter(nplot < 6) %>% 
    dplyr::arrange(nplot, site, generation) %>% 
    as.data.frame()

# output files --------
save.image(file = paste0(plotdir, 'plotdata.rda'))

# cleanup --------
date()
closeAllConnections()

# Rscript --vanilla scripts/manuscript/plot_ecotype_freq.R &> data/manuscript/plot_ecotype_freq/plot_ecotype_freq.log

