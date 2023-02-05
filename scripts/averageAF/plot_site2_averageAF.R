# Title: Plot the tested weighted average in site 2 (exploratory)
# Author: Meixi Lin
# Date: Fri Jan 27 10:11:57 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/central/groups/carnegie_poc/meixilin/grenenet/analyses")

library(data.table)
library(reshape2)

# def functions --------
plot_pairs <- function(mysite, mysitedt) {
    png(filename = paste0(plotdir, 'site',mysite, '_pairs_year2018_chr1_10kSNPs.png'), width = 500, height = 500)
    pairs(mysitedt[sample(1:nrow(mysitedt), 10000),-1], pch = '.')
    dev.off()
}

get_name <- function(xx) {
    nn = stringr::str_split(xx, '/')[[1]]
    nn = stringr::str_split(nn[length(nn)], '\\.')[[1]][1]
    return(nn)
}

load_chr1_plots <- function(myspy) {
    myfiles = list.files(path = '/central/groups/carnegie_poc/lczech/grenephase1/hafpipe-231/frequencies/', 
                         pattern = paste0('MLFH',myspy, '.+.1.afSite$'), full.names = TRUE)
    mydtl = lapply(myfiles, function(xx) {
        dt = data.table::fread(file = xx)[,2]
        return(dt)
    })
    outdt = do.call(cbind, mydtl)
    colnames(outdt) = unname(sapply(myfiles, get_name))
    return(outdt)
}

plot_hist <- function(mydt, mysp, nreadsdt, nflowersdt) {
    png(filename = paste0(plotdir, 'site',mysp, '_hist_year2018_chr1_10kSNPs.png'), width = 500, height = 500)
    pp <- par(mfrow = c(2,2))
    mydt = data.table::setDF(mydt)
    for (ii in 1:ncol(mydt)) {
        hist(mydt[,ii], main = colnames(mydt)[ii])
        text('topright')
    }
    par(pp)
    dev.off()
}

# def variables --------
plotdir = './data/averageAF/plots/'
dir.create(plotdir)

# load data --------
site2_2018 = data.table::fread('./data/averageAF/FH_site2_year2018_chr1_averageAF.csv')
site10_2018 = data.table::fread('./data/averageAF/FH_site10_year2018_chr1_averageAF.csv')
site13_2018 = data.table::fread('./data/averageAF/FH_site13_year2018_chr1_averageAF.csv')

# load original data
site10.1_2018 = load_chr1_plots('10012018')
site10.4_2018 = load_chr1_plots('10042018')

# main --------
plot_pairs(2, site2_2018)
plot_pairs(10, site10_2018)
plot_pairs(13, site13_2018)

# since site10 plot1 - plot3 had good consistencies, check with site2 plot1 - plot3
all.equal(site2_2018[,1], site10_2018[,1])
diffsite_2018 = cbind(site2_2018[,1:4], site10_2018[,2:4])
plot_pairs('2vs10', diffsite_2018)

# check with site10, year2018, within plot1 and within plot4
plot_pairs('2plot1', cbind(site10_2018[,1], site10.1_2018))
plot_hist(site10.1_2018, '1001')

plot_pairs('2plot4', cbind(site10_2018[,1], site10.4_2018))
plot_hist(site10.4_2018, '1004')

# check coverage and number of flowers
load('../metadata/data/samples_data.rda')
dd = samples_data %>%
    dplyr::filter(year == 2018, site == 10, plot %in% c(1,4)) %>%
    dplyr::select(sampleid,flowerscollected)

load('data/fastq_check/fastqc_summary.rda')
fastqcdt %>% 
    dplyr::filter(sampleid %in% dd$sampleid) %>% 
    dplyr::distinct(sampleid,.keep_all = TRUE) %>% 
    dplyr::select(sampleid,coverage,nreads)

# output files --------

# cleanup --------
date()
closeAllConnections()