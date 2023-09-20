# Title: Plot results from local adaptation analyses
# Author: Meixi Lin
# Date: Wed Aug 30 21:21:11 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses")

library(cowplot)
library(ggplot2)
library(dplyr)

date()
sessionInfo()

theme_set(cowplot::theme_cowplot(font_size = 8))

# function --------
loadSomeRData <- function(x, file) {
    E <- new.env()
    load(file=file, envir=E)
    return(get(x, envir=E, inherits=F))
}

save_figure <- function(plotdf, pp, myname, mywid, myhei) {
    save(plotdf, pp, file = paste0(outdir, myname, '.rda'))
    ggsave(filename = paste0(outdir, myname, '.pdf'), plot = pp, useDingbats = TRUE,
           width = mywid, height = myhei)
    return(invisible())
}

# define variables --------
outdir = './data/local_adaptation/plots/'
dir.create(outdir)

# environmental variables to plot
# bio11: average temperature at coldest quarter for persite analysis. 
# [reason:] average R2 is among the highest and it is a bio variable. 
# other high R2 observations were also temperature data at winter, e.g. tmax02, tmax01, etc. 
myvar = 'bio11'
# when differences for specific biovariables are not needed. `bioall` was used as a proxy for multidimensional environmental space distance. 
allvar = 'bioall' 

# plotting settings
myscale_size = scale_size_continuous(limits = c(-0.002,0.5), range = c(0.1,4), breaks = c(0,0.2,0.4))


# load data --------
# NOTE: ecop_l in this dataset have 4 ecotypes removed due to lack of worldclim data. 
load('./data/local_adaptation/inputs/ecotypefreqs_deltaenv.rda')
ecop_l231 = loadSomeRData('ecop_l231', './data/local_adaptation/inputs/step0_prepare_data.RData')
deltaenv_raw = loadSomeRData('deltaenv_raw', './data/local_adaptation/inputs/step0_prepare_data.RData')

# metadata
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


# supp: R2 values in the models --------
# summarize results 
# which environmental variables have the highest average R2m at ecotype levels
sum_eco = res_eco %>% 
    dplyr::group_by(envvar) %>% 
    dplyr::summarise(ave_R2 = mean(R2m, na.rm = TRUE),
                     med_R2 = median(R2m, na.rm = TRUE),
                     n_R215 = sum(R2m > 0.15, na.rm = TRUE)) %>% 
    dplyr::arrange(desc(ave_R2))
sum_eco

# which environmental variables have the highest average R2 at site levels
sum_site = res_site %>% 
    dplyr::group_by(envvar) %>% 
    dplyr::summarise(ave_R2 = mean(adj_r2, na.rm = TRUE),
                     med_R2 = median(adj_r2, na.rm = TRUE),
                     n_R215 = sum(adj_r2 > 0.15, na.rm = TRUE)) %>% 
    dplyr::arrange(desc(ave_R2))
sum_site

# A: local adaptation with every ecotypes and every site --------
# dplyr::sample_frac so the overlay of points are not forced sorted by ecotypeid
plotadf = dplyr::left_join(ecop_l, deltaenv_raw[,c('ecotypeid','site', myvar)], by = c('ecotypeid', 'site')) %>% 
    dplyr::left_join(., worldclim_ecotypesdata[,c('ecotypeid', myvar)], by = 'ecotypeid', suffix = c('.diff', '.eco')) %>% 
    dplyr::sample_frac(size = 1) 
ppA <- ggplot(data = plotadf, aes(x = eval(as.name(paste0(myvar, '.diff'))), y = frequency, color = eval(as.name(paste0(myvar, '.eco'))))) + 
    geom_hline(yintercept = 1/231, linetype = 'dashed', color = 'lightgray') + 
    geom_point(alpha = 0.6) + 
    facet_wrap(generation ~ ., ncol = 1) +
    scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 9, 'Greens')) +
    labs(color = 'bio11',
         x = 'climate distance (bio11)',
         y = 'ecotype frequency') +
    theme(legend.position = 'top',
          legend.title = element_blank(),
          legend.key.width = unit(0.7, 'cm'))

# B: map of local adaptation using bioall at ecotype level --------
plotbdf = res_eco %>% 
    dplyr::filter(envvar == allvar) %>% 
    dplyr::left_join(., ecotypes_data[,c('ecotypeid', 'longitude','latitude')], by = 'ecotypeid') %>% 
    dplyr::arrange(desc(b_env))

ppB <- ggplot() +
    geom_polygon(data = wrld_simpl, aes(x = long, y = lat, group = group), fill = "lightgray", colour = "lightgray") +
    geom_point(data = plotbdf, mapping = aes(x = longitude, y = latitude, color = -b_env, size = R2m), alpha = 0.6) + 
    coord_map(projection = 'ortho', xlim = range(plotbdf$longitude), ylim = range(plotbdf$latitude), orientation = c(90,0,0)) + 
    myscale_size + 
    guides(size = 'none') + 
    viridis::scale_color_viridis() + 
    facet_wrap(generation ~ ., ncol = 1) +
    theme(axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = 'top',
          legend.title = element_blank(),
          legend.key.width = unit(0.7, 'cm')) + 
    labs(color = '1/Vs',
         x = 'ecotype locations')

# C: map of local adaptation using bioall at experimental site level --------
plotcdf = res_site %>% 
    dplyr::filter(envvar == allvar) %>% 
    dplyr::left_join(., locations_data[,c('site', 'longitude','latitude')], by = 'site') %>% 
    dplyr::arrange(desc(b_env))
# summary(plotcdf)

ppC <- ggplot() +
    geom_polygon(data = wrld_simpl, aes(x = long, y = lat, group = group), fill = "lightgray", colour = "lightgray") +
    geom_point(data = plotcdf, mapping = aes(x = longitude, y = latitude, color = -b_env, size = adj_r2), alpha = 0.6) + 
    coord_map(projection = 'ortho', xlim = range(plotcdf$longitude), ylim = range(plotcdf$latitude), orientation = c(90,0,0)) + 
    myscale_size + 
    guides(size = 'none') + 
    viridis::scale_color_viridis() + 
    facet_wrap(generation ~ ., ncol = 1) +
    theme(axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = 'top',
          legend.title = element_blank(),
          legend.key.width = unit(0.7, 'cm')) + 
    labs(color = '1/Vs',
         x = 'experiment site locations')

# D: correlations of Vs and climate --------
plotddf = res_eco %>% 
    dplyr::filter(envvar == myvar) %>% 
    dplyr::left_join(worldclim_ecotypesdata[,c('ecotypeid', myvar)], by = 'ecotypeid')

ppD <- ggplot(data = plotddf, aes(x = eval(as.name(myvar)), y = -b_env, size = R2m, color = as.factor(generation))) + 
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'lightgray') + 
    geom_point(alpha = 0.6) +
    myscale_size +
    labs(color = 'generation',
         size = bquote('R'^2),
         x = 'ecotype bio11',
         y = 'selection strength (bio11)') + 
    theme(legend.position = 'none')

# E: correlations of tau and climate --------
plotedf = res_site %>% 
    dplyr::filter(envvar == myvar) %>% 
    dplyr::left_join(worldclim_sitesdata[,c('site', myvar)], by = 'site')

ppE <- ggplot(data = plotedf, aes(x = eval(as.name(myvar)), y = -b_env, size = adj_r2, color = as.factor(generation))) + 
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'lightgray') + 
    geom_point(alpha = 0.6) +
    myscale_size +
    labs(color = 'generation',
         size = bquote('R'^2),
         x = 'experiment site bio11',
         y = 'selection strength (bio11)') + 
    theme(legend.position = 'none')

# F: correlations of tau and population size --------
# proxy: number of unique ecotypes in this sample
# use the ecotypes without filtering (ecop_l231). Make a small cutoff so ecotypes <= total_flowers
necodf = ecop_l231 %>% 
    dplyr::left_join(merged_samples_data[,c('sample_name', 'total_flower_counts')], by = 'sample_name') %>% 
    dplyr::mutate(freqcut = ifelse(1/total_flower_counts < 0.02, 1/total_flower_counts, 0.02),  
                  freq_s = frequency > freqcut,
                  freq_0 = frequency > 0) %>%
    dplyr::group_by(sample_name, site, generation, plot, total_flower_counts) %>% 
    dplyr::summarise(neco = sum(freq_s),
                     neco0 = sum(freq_0), 
                     .groups = 'drop') 
# plot(necodf$total_flower_counts, necodf$neco)

# now calculate the number of ecotypes on average 
ave_neco = necodf %>% 
    dplyr::group_by(site, generation) %>% 
    dplyr::summarise(ave_neco = mean(neco),
                     ave_nflower = mean(total_flower_counts))
# plot(ave_neco$ave_nflowers, ave_neco$ave_neco)

# plot correlations with b_env
plotfdf = res_site %>% 
    dplyr::left_join(ave_neco, by =c('site', 'generation')) %>% 
    dplyr::filter(envvar == allvar)

# ave_nflower has less correlations. use ave_neco here
ppF <- ggplot(data = plotfdf, aes(x = ave_neco, y = -b_env, size = adj_r2, color = as.factor(generation))) + 
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'lightgray') + 
    geom_point(alpha = 0.6) +
    myscale_size +
    # scale_x_log10() +
    labs(color = 'generation',
         size = bquote('R'^2),
         x = 'number of unique ecotypes',
         y = 'selection strength (all)') + 
    theme(legend.position = 'none')

# assemble final figure --------
ptop <- plot_grid(ppA, ppB, ppC, labels = "AUTO", align = 'h', nrow = 1)
pbot <- plot_grid(ppD, ppE, ppF, labels = c('D','E','F'), align = 'hv', nrow = 1)
botleg <- ggpubr::as_ggplot(ggpubr::get_legend(ppF + theme(legend.position = 'right')))

ggsave(filename = paste0(outdir, 'ABC_local_adaptation.pdf'), width = 180, height = 160, units = 'mm', plot = ptop)
ggsave(filename = paste0(outdir, 'DEF_local_adaptation.pdf'), width = 180, height = 60, units = 'mm', plot = pbot)
ggsave(filename = paste0(outdir, 'DEFleg_local_adaptation.pdf'), width = 20, height = 60, units = 'mm', plot = botleg)

# output files --------
save(plotadf, plotbdf, plotcdf, plotddf, plotedf, plotfdf, file = paste0(outdir, 'plotdf.rda'))
save.image(paste0(outdir, 'plot_local_adaptation.RData'))

# cleanup --------
date()
closeAllConnections()



