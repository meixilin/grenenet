# Title: Run replicability analyses using the leave-one-out method
# Author: Meixi Lin
# Date: Thu Oct 13 08:46:32 2022

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet")

library(ggplot2)
library(dplyr)
library(cowplot)
date()
sessionInfo()

theme_set(theme_classic(base_size = 12))

# def functions --------
loadSomeRData <- function(x, file) {
    E <- new.env()
    load(file=file, envir=E)
    return(get(x, envir=E, inherits=F))
}

ggsaver <- function(pp, suffix, height, width) {
    ggsave(filename = paste0(suffix, '.pdf'), plot = pp, path = plotdir, height = height, width = width)
    ggsave(filename = paste0(suffix, '.png'), plot = pp, path = plotdir, height = height, width = width, units = 'in', dpi = 150)
}
# def variables --------
# args = commandArgs(trailingOnly=TRUE)
# envvar = as.character(args[1]) # environment variable to regress against
envvar = 'bio1'

plotdir = paste0('plots/ver202209/Replicability/', envvar, '/')
dir.create(plotdir, recursive = TRUE)

# these positions should not be input for linear models
badpos = c("8916","9094035","9527351","13914486","21327459","21368627")

# load data --------
deltadt = readRDS(file = './data/AF/ver202209/haplotype/scaled_delta_freq_1012.rds')
metadata = readRDS(file = './data/metadata/ver202209/metacensusclim_0922.rds')
deltadt = deltadt[-which(rownames(deltadt) %in% badpos),]
dim(deltadt)

# load the linear model output (significant alleles)
lmres = read.csv(file = './plots/ver202209/AF/DeltaP_v2/bio1/scaled_deltaP_bio1_lmres_Ppass.csv')

# main --------
# plot the changes by time and replicates
pplist <- lapply(1:nrow(lmres), function(xx){
    thislm = lmres[xx, ]
    scalep = deltadt[as.character(thislm$POS),]
    all(metadata$id == colnames(deltadt))
    myheader = paste0('CHR1-',thislm$POS, ' AdjR2=',round(thislm$AdjR2, digits = 3),
                      ' PxAdj=',round(thislm$P_x_adj, digits = 3))
    plotdt = cbind(scalep, metadata[,c(envvar, 'SITE', 'PLOT', 'year', 'DATE')]) %>%
        dplyr::mutate(DATE = as.Date(DATE, format = '%Y%m%d'),
                      PLOT = as.character(PLOT)) %>%
        dplyr::arrange(eval(as.name(envvar)))
    plotdt$SITE = paste0(plotdt$SITE, '-',envvar,':',plotdt[,envvar])
    plotdt$SITE = factor(plotdt$SITE, levels = unique(plotdt$SITE))

    # plot the allele frequency shift
    pp <- ggplot(data = plotdt, aes(x = DATE, y = scalep, color = PLOT, group = PLOT)) +
        geom_point() +
        geom_line() +
        geom_hline(yintercept = 0, linetype = 'dashed', color = 'gray') +
        facet_wrap(. ~ SITE) +
        labs(title = myheader) +
        theme(axis.text.x = element_text(angle = 90))
    ggsaver(pp=pp,suffix = paste0('Trajectory_', envvar,'-', thislm$POS), height = 8, width = 9)

    # what if we average the year values
    plotdt2 = plotdt %>%
        dplyr::group_by(SITE,year,PLOT) %>%
        dplyr::summarise(scalepA = mean(scalep),
                         nsample = n(),
                         scalepSD = sd(scalep),
                         .groups = 'drop') %>%
        dplyr::group_by(SITE) %>%
        dplyr::mutate(PLOT = ifelse(as.integer(PLOT) == max(as.integer(PLOT)), 'Rep2','Rep1'))
    startdt = data.frame(SITE = rep(unique(plotdt2$SITE),each = 2),
                         PLOT = rep(c('Rep1','Rep2'), times = length(unique(plotdt2$SITE))))
    startdt$year = '2017'; startdt$nsample=1; startdt$scalepA = 0; startdt$scalepSD = NA
    plotdt2 = bind_rows(startdt,plotdt2)
    pp2 <- ggplot(data = plotdt2, aes(x = year, y = scalepA, color = PLOT, group = PLOT)) +
        geom_hline(yintercept = 0, linetype = 'dashed', color = 'gray') +
        geom_pointrange(alpha = 0.8, aes(ymin=scalepA-scalepSD, ymax=scalepA+scalepSD))+
        geom_line(alpha = 0.8) +
        facet_wrap(. ~ SITE) +
        labs(title = myheader) +
        theme(axis.text.x = element_text(angle = 90))
    ggsaver(pp=pp2,suffix = paste0('AveTrajectory_', envvar,'-', thislm$POS), height = 8, width = 9)

    # test for consistency
    testdt = data.frame(SITE = rep(unique(plotdt2$SITE),each = 4),
                        PLOT = rep(c('Rep1','Rep1','Rep2','Rep2'), times = length(unique(plotdt2$SITE))),
                        year = rep(c('2018','2019','2018','2019'), times = length(unique(plotdt2$SITE))),
                        Dscalep = NA)
    for (ii in 1:nrow(testdt)) {
        if (testdt[ii,'year'] == 2018) {
            testdt[ii,'Dscalep'] = plotdt2[which(plotdt2$SITE == testdt[ii,'SITE'] &
                                           plotdt2$PLOT == testdt[ii,'PLOT'] &
                                           plotdt2$year == '2018'),'scalepA']
        } else {
            id = which(plotdt2$SITE == testdt[ii,'SITE'] &
                           plotdt2$PLOT == testdt[ii,'PLOT'] &
                           plotdt2$year == '2019')
            if(length(id) == 0) {
                testdt[ii,'Dscalep'] = NA
            } else {
            testdt[ii,'Dscalep'] = plotdt2[id,'scalepA'] -
                                   plotdt2[which(plotdt2$SITE == testdt[ii,'SITE'] &
                                                      plotdt2$PLOT == testdt[ii,'PLOT'] &
                                                      plotdt2$year == '2018'),'scalepA']
            }
        }
    }

    # pass testdt
    testdtout = testdt %>%
        dplyr::mutate(Dsign = Dscalep > 0) %>%
        dplyr::group_by(SITE, year) %>%
        dplyr::mutate(allsign = all(Dsign == TRUE))
    print(thislm$POS)
    print(table(testdtout$allsign, useNA = 'always'))
    return(list(pp,pp2))
})

# output files --------

# cleanup --------
date()
closeAllConnections()
