# Title: Run an NMDS on allele frequencies
# Author: Meixi Lin
# Date: Thu Aug  4 00:20:37 2022

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet")

date()
library(vegan)
library(ggplot2)
library(dplyr)
sessionInfo()

# def functions --------
plotter_ord <- function(dt,xx,yy,cc,legpos = 'none', colpalette) {
    pp <- ggplot(dt, aes(x = .data[[xx]], y = .data[[yy]], color = .data[[cc]])) +
        geom_point() +
        geom_point(data = dt[dt$id == 'seed_mix',], shape = 9, size = 5, color = 'black') +
        scale_color_manual(values = colpalette) +
        ggpubr::theme_pubr() +
        theme(legend.position = legpos)
    return(pp)
}

plot_ord <- function(forplot, cc, colpalette, addhull = FALSE) {
    pp1 <- plotter_ord(dt = forplot, xx = 'MDS1', yy = 'MDS2',
                       cc = cc, colpalette = colpalette)
    pp2 <- plotter_ord(dt = forplot, xx = 'MDS1', yy = 'MDS3',
                       cc = cc, colpalette = colpalette)
    pp3 <- plotter_ord(dt = forplot, xx = 'MDS3', yy = 'MDS2',
                       cc = cc, colpalette = colpalette)
    leg <- ggpubr::get_legend(pp1 + theme(legend.position = 'bottom'))
    if (addhull) {
        pp1 <- add_hull(pp1)
        pp2 <- add_hull(pp2)
        pp3 <- add_hull(pp3)
    }
    pp <- ggpubr::ggarrange(pp1, pp3, pp2, leg, nrow = 2, ncol = 2)
    return(pp)
}

add_hull <- function(pp) {
    ppbuild <- ggplot_build(pp)
    # get color scale
    mycolors = unique(ppbuild$data[[1]]$colour)
    names(mycolors) = unique(pp$data[,pp$labels$colour])
    hull_forplot <- pp$data %>%
        group_by(.data[[pp$labels$colour]]) %>%
        slice(chull(.data[[pp$labels$x]], .data[[pp$labels$y]]))
    pph <- pp +
        geom_polygon(data = hull_forplot, alpha = 0.1, aes(fill = .data[[pp$labels$colour]])) +
        scale_fill_manual(values = mycolors)
    return(pph)
}

format_mds <- function(myMDS,allmeta0922) {
    forplot = myMDS$point %>%
        data.frame() %>%
        tibble::rownames_to_column(var = 'id') %>%
        dplyr::left_join(.,
                         y = allmeta0922[,c('id','SITE','year', 'country')], by = 'id') %>%
        dplyr::mutate(SITE = as.character(SITE))
    forplot[forplot$id == 'seed_mix', c('SITE','year','country')] = c('seed_mix', 'start', 'seed_mix')
    return(forplot)
}

ggsaver <- function(filename, plot, height, width) {
    ggsave(filename = filename, plot = plot, path = plotdir,
           height = height, width = width, dpi = 72)
    return(invisible())
}

# def variables --------
maflist = c('00','05','10')
methodlist = c('euc','nei')
mafmethodl = expand.grid(maflist,methodlist, stringsAsFactors = FALSE)
colnames(mafmethodl) = c('maf','method')
mafmethodn = unlist(apply(mafmethodl, 1, function(xx) {paste0('maf',xx['maf'], '_', xx['method'])}))


plotdir = './plots/ver202209/AF/NMDS/'
dir.create(plotdir)
dir.create('data/AF/ver202209/haplotype/NMDS')

source('./scripts/ver202209/config.R')

# load data --------



# metadata ========
allmeta0922 = readRDS(file = 'data/metadata/ver202209/metacensusclim_0922.rds')
# append the country location
allmeta0922$country <- maps::map.where(database="world", allmeta0922$LONGITUDE, allmeta0922$LATITUDE)

# main --------
# loop through all six combinations
mdslist <- apply(mafmethodl, 1, function(xx) {
    maf = xx['maf'];method = xx['method']
    # load distance
    distfile = paste0('./data/AF/ver202209/haplotype/distAF/distAF_maf', maf, '_', method, '_0922.rds')
    mydist = readRDS(distfile)
    # run mds (Species scores not available expected given we only give dist)
    myMDS <- vegan::metaMDS(mydist, k = 3, trymax = 10)
    print(paste(distfile, 'MDS summary:'))
    print(myMDS)
    return(myMDS)
})
names(mdslist) = mafmethodn

# color the NMDS by SITE, year, SITE/year
nmds_site <- mapply(function(myMDS,mdsid) {
    # format forplot
    forplot <- format_mds(myMDS, allmeta0922)
    # run plotting command
    pp1 <- plot_ord(forplot, cc = 'SITE', colpalette = site_cols)
    ggsaver(filename = paste0('nmds_site_',mdsid, '.png'), plot = pp1, height = 7, width = 8)
    # add hulls
    pp2 <- plot_ord(forplot, cc = 'SITE', colpalette = site_cols, addhull = TRUE)
    ggsaver(filename = paste0('nmdsH_site_',mdsid, '.png'), plot = pp2, height = 7, width = 8)
    return(list(pp1,pp2))
}, myMDS = mdslist, mdsid = names(mdslist))

# color the NMDS by year
nmds_year <- mapply(function(myMDS,mdsid) {
    # format forplot
    forplot <- format_mds(myMDS, allmeta0922)
    # run plotting command
    pp1 <- plot_ord(forplot, cc = 'year', colpalette = year_cols)
    ggsaver(filename = paste0('nmds_year_',mdsid, '.png'), plot = pp1, height = 7, width = 8)
    # add hulls
    pp2 <- plot_ord(forplot, cc = 'year', colpalette = year_cols, addhull = TRUE)
    ggsaver(filename = paste0('nmdsH_year_',mdsid, '.png'), plot = pp2, height = 7, width = 8)
    return(list(pp1,pp2))
}, myMDS = mdslist, mdsid = names(mdslist))

# color by sites in 2018
nmds_site2018 <- mapply(function(myMDS,mdsid) {
    # format forplot
    forplot <- format_mds(myMDS, allmeta0922) %>%
        dplyr::filter(year != '2019')
    # run plotting command
    pp1 <- plot_ord(forplot, cc = 'SITE', colpalette = site_cols)
    ggsaver(filename = paste0('nmds_site2018_',mdsid, '.png'), plot = pp1, height = 7, width = 8)
    # add hulls
    pp2 <- plot_ord(forplot, cc = 'SITE', colpalette = site_cols, addhull = TRUE)
    ggsaver(filename = paste0('nmdsH_site2018_',mdsid, '.png'), plot = pp2, height = 7, width = 8)
    return(list(pp1,pp2))
}, myMDS = mdslist, mdsid = names(mdslist))

# color by sites in 2019
nmds_site2019 <- mapply(function(myMDS,mdsid) {
    # format forplot
    forplot <- format_mds(myMDS, allmeta0922) %>%
        dplyr::filter(year != '2018')
    # run plotting command
    pp1 <- plot_ord(forplot, cc = 'SITE', colpalette = site_cols)
    ggsaver(filename = paste0('nmds_site2019_',mdsid, '.png'), plot = pp1, height = 7, width = 8)
    # add hulls
    pp2 <- plot_ord(forplot, cc = 'SITE', colpalette = site_cols, addhull = TRUE)
    ggsaver(filename = paste0('nmdsH_site2019_',mdsid, '.png'), plot = pp2, height = 7, width = 8)
    return(list(pp1,pp2))
}, myMDS = mdslist, mdsid = names(mdslist))

# plot the NMDS by geographical clusters (countries)
nmds_country <- mapply(function(myMDS,mdsid) {
    # format forplot
    forplot <- format_mds(myMDS, allmeta0922)
    # run plotting command
    pp1 <- plot_ord(forplot, cc = 'country', colpalette = country_cols)
    ggsaver(filename = paste0('nmds_country_',mdsid, '.png'), plot = pp1, height = 7, width = 8)
    # add hulls
    pp2 <- plot_ord(forplot, cc = 'country', colpalette = country_cols, addhull = TRUE)
    ggsaver(filename = paste0('nmdsH_country_',mdsid, '.png'), plot = pp2, height = 7, width = 8)
    return(list(pp1,pp2))
}, myMDS = mdslist, mdsid = names(mdslist))

# output files --------
# save the mds object
saveRDS(mdslist, file = 'data/AF/ver202209/haplotype/NMDS/mdslist.rds')
# save rdata
save.image(file = './rdata/ver202209/fig2.1_nmds_af.RData')

# cleanup --------
date()
closeAllConnections()
