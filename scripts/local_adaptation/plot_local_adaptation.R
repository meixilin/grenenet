# Title: Plot results from local adaptation analyses
# Author: Meixi Lin
# Date: Wed Aug 30 21:21:11 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses")

library(ggplot2)
library(dplyr)

date()
sessionInfo()

theme_set(cowplot::theme_cowplot())

# function --------
loadSomeRData <- function(x, file) {
    E <- new.env()
    load(file=file, envir=E)
    return(get(x, envir=E, inherits=F))
}

# load data --------
# note that ecop_l in this dataset have 4 ecotypes removed due to lack of worldclim data. 
# therefore it does not sum to 1 per sample
load('./data/local_adaptation/inputs/ecotypefreqs_deltaenv.rda')
load('../metadata/data/worldclim_ecotypesdata.rda')
load('../metadata/data/worldclim_sitesdata.rda')
load('../metadata/data/ecotypes_data.rda')
load('../metadata/data/locations_data.rda')
load('../metadata/data/merged_samples_data.rda')
# per ecotype adaptation
res_eco = loadSomeRData('alldf', './data/local_adaptation/perecotype/la_lmmres_perecotype.rda')
# per site adaptation
res_site = loadSomeRData('alldf', './data/local_adaptation/persite/la_lmres_persite.rda')

data("wrld_simpl", package = "maptools")
# 
# # modify data --------
# # calculate deltaecop
# ecodp = dplyr::left_join(ecop_l, ecop0, by = 'ecotypeid') %>% 
#     dplyr::mutate(dp = frequency - p0,
#                   p1_p0 = log10(frequency/p0))
# 
# hist(ecodp$dp)

# A: local adaptation with every ecotypes and every site --------
# bio11: average temperature at coldest quarter. 
# [reason:] average R2 is among the highest and the highest bio variable. 
# other high R2 observations were also temperature data at winter, e.g. tmax02, tmax01, etc. 
plotadf = dplyr::left_join(ecop_l, deltaenv_raw[,c('ecotypeid','site', 'bio11')], by = c('ecotypeid', 'site')) %>% 
    dplyr::left_join(., worldclim_ecotypesdata[,c('ecotypeid', 'bio11')], by = 'ecotypeid', suffix = c('.diff', '.eco')) %>% 
    dplyr::sample_frac() 
ggplot(data = plotadf, aes(x = bio11.diff, y = frequency, color = bio11.eco)) + 
    geom_hline(yintercept = 1/231, linetype = 'dashed', color = 'gray') + 
    geom_point(alpha = 0.8) + 
    facet_wrap(generation ~ ., ncol = 1) +
    scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 9, 'Greens'))

# B: map of local adaptation using bioall at ecotype level --------
plotbdf = res_eco %>% 
    dplyr::filter(envvar == 'bioall') %>% 
    dplyr::left_join(., ecotypes_data[,c('ecotypeid', 'longitude','latitude')], by = 'ecotypeid') %>% 
    dplyr::arrange(desc(b_env))

ggplot() +
    geom_polygon(data = wrld_simpl, aes(x = long, y = lat, group = group), fill = "lightgray", colour = "lightgray") +
    geom_point(data = plotbdf, mapping = aes(x = longitude, y = latitude, color = -b_env, size = -log10(b_p)), alpha = 0.6) + 
    coord_map(projection = 'ortho', xlim = range(plotbdf$longitude), ylim = range(plotbdf$latitude), orientation = c(90,0,0)) + 
    scale_size_continuous(range = c(0.1,5)) + 
    viridis::scale_color_viridis() + 
    facet_wrap(generation ~ ., ncol = 1) 

# C: map of local adaptation using bioall at experimental site level --------
# projection = 'ortho', 
plotcdf = res_site %>% 
    dplyr::filter(envvar == 'bioall') %>% 
    dplyr::left_join(., locations_data[,c('site', 'longitude','latitude')], by = 'site') %>% 
    dplyr::arrange(desc(b_env))

ggplot() +
    geom_polygon(data = wrld_simpl, aes(x = long, y = lat, group = group), fill = "lightgray", colour = "lightgray") +
    geom_point(data = plotcdf, mapping = aes(x = longitude, y = latitude, color = -b_env, size = -log10(b_p)), alpha = 0.6) + 
    coord_map(projection = 'ortho', xlim = range(plotcdf$longitude), ylim = range(plotcdf$latitude), orientation = c(90,0,0)) + 
    scale_size_continuous(range = c(0.1,5)) + 
    viridis::scale_color_viridis() + 
    facet_wrap(generation ~ ., ncol = 1) 


# D: correlations of tau and population size --------
# proxy: number of unique ecotypes in this sample
# any ecotypes with frequency less than maximum flowers collected at that sample are removed.
# load the ecotypes without filtering
ecop = read.delim('/NOBACKUP/scratch/xwu/grenet/merged_frequency/merged_ecotype_frequency.txt')
xx = ecop_l %>% 
    dplyr::group_by(sample_name) %>% 
    dplyr::summarise(sumf = sum(frequency, na.rm = T))

xx = ecop_l %>% 
    dplyr::filter(sample_name == '52_3_10')

necodf = ecop_l %>% 
    dplyr::left_join(worldclim_sitesdata[,c('site', 'bio11')], by = 'site') %>%
    dplyr::left_join(merged_samples_data[,c('sample_name', 'total_flower_counts')], by = 'sample_name') %>% 
    dplyr::filter(frequency > 1/total_flower_counts) %>%
    dplyr::group_by(sample_name, site, generation, plot, bio11) %>% 
    dplyr::summarise(n_eco = n(), .groups = 'drop') %>% 
    dplyr::arrange(bio11) 
necodf$site = factor(necodf$site, levels = unique(necodf$site))

ggplot(necodf, aes(x = bio11, y = n_eco, color = bio11)) + 
    geom_point() +
    facet_wrap(generation ~ ., ncol = 1) +
    scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 9, 'Greens')[3:9])



plotddf = res_site %>% 
    dplyr::left_join(merged_samples_data, by =c('site', 'generation')) %>% 
    dplyr::select(site, generation, plot, envvar, total_flower_counts, b_env) %>% 
    dplyr::filter(envvar == 'bio11')

ggplot(plotddf, aes(x = total_flower_counts, y = b_env, color = as.factor(site))) + 
    geom_point() +
    facet_wrap(generation ~ ., ncol = 1) 


# supp: R2 values in the models --------
# summarize results 
# which environmental variables have the highest average R2m
sum_eco = res_eco %>% 
    dplyr::group_by(envvar) %>% 
    dplyr::summarise(ave_R2 = mean(R2m, na.rm = TRUE),
                     med_R2 = median(R2m, na.rm = TRUE),
                     n_R215 = sum(R2m > 0.15, na.rm = TRUE)) %>% 
    dplyr::arrange(desc(ave_R2))
sum_eco

# which environmental variables have the highest average R2 per site
sum_site = res_site %>% 
    dplyr::group_by(envvar) %>% 
    dplyr::summarise(ave_R2 = mean(adj_r2, na.rm = TRUE),
                     med_R2 = median(adj_r2, na.rm = TRUE),
                     n_R215 = sum(adj_r2 > 0.15, na.rm = TRUE)) %>% 
    dplyr::arrange(desc(ave_R2))
sum_site

plots1df = res_eco %>% 
    dplyr::mutate(envtype = ifelse(envvar %in% c('latitude', 'longitude'), 'loc', gsub("\\d", "", envvar)),
                  envvar = factor(envvar, levels = sum_eco$envvar)) %>% 
    dplyr::arrange(envvar)

ggplot(data = plots1df, mapping = aes(x = envvar, y = R2m, fill = envtype)) + 
    geom_boxplot() + 
    scale_fill_brewer(palette = 'Set3') + 
    theme(axis.text.x = element_text(angle = 90))

plots2df = res_site %>% 
    dplyr::mutate(envtype = ifelse(envvar %in% c('latitude', 'longitude'), 'loc', gsub("\\d", "", envvar)),
                  envvar = factor(envvar, levels = sum_site$envvar)) %>% 
    dplyr::arrange(envvar)

ggplot(data = plots2df, mapping = aes(x = envvar, y = adj_r2, fill = envtype)) + 
    geom_boxplot() + 
    scale_fill_brewer(palette = 'Set3') + 
    theme(axis.text.x = element_text(angle = 90))

# supp: correlations of tau and population size --------
# proxy 1: number of flowers in merged samples (NO CORRELATIONS)
plotddf = res_site %>% 
    dplyr::left_join(merged_samples_data, by =c('site', 'generation')) %>% 
    dplyr::select(site, generation, plot, envvar, total_flower_counts, b_env) %>% 
    dplyr::filter(envvar == 'bio11')

ggplot(plotddf, aes(x = total_flower_counts, y = b_env, color = as.factor(site))) + 
    geom_point() +
    scale_x_log10() + 
    facet_wrap(generation ~ ., ncol = 1) 


# output files --------

# cleanup --------
date()
closeAllConnections()
    