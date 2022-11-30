# Title: Plot fastq info
# Author: Meixi Lin
# Date: Wed Nov  9 14:10:33 2022

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet")
date()
sessionInfo()

library(ggplot2)
library(dplyr)

theme_set(cowplot::theme_cowplot())

# def functions --------

# def variables --------

# load data --------
samplenames = read.delim(file = './data/metadata/fastq_info/samplenames.txt', header = FALSE)
lucasdt = read.delim(file = './data/metadata/fastq_info/seqtable_cal.tsv')
load(file = '../../ath_evo/grenephase1/data/fastq_info.rda')
setdiff(fastq_info$sampleid, lucasdt$sample)
# MLFH570120190821
setdiff(lucasdt$sample, fastq_info$sampleid)
fastq_info[fastq_info$sampleid == 'MLFH570120190821',]
# sample names check
setdiff(samplenames$V1, fastq_info$sampleid)
setdiff(fastq_info$sampleid, samplenames$V1)
# fastq names check
setdiff()

# main --------
# heatmap of files
sample_list = data.frame(site = substr(fastq_info$sampleid,5,6),
                         plot = as.integer(substr(fastq_info$sampleid,7,8)),
                         date = substr(fastq_info$sampleid,9,16),
                         year = substr(fastq_info$sampleid,9,12),
                         sampleid = fastq_info$sampleid)

sites = data.frame(site = unique(sample_list$site))
sample_tally = sample_list %>%
    dplyr::group_by(site, plot, year) %>%
    dplyr::summarise(n_dates = n())

sample_tally2 = sample_list %>%
    dplyr::group_by(site, plot, year) %>%
    dplyr::summarise(n_dates = n() > 0) %>%
    dplyr::group_by(site, plot) %>%
    dplyr::summarise(n_years = sum(n_dates))

pp <- ggplot(data = sample_tally, mapping = aes(x = plot, y = site, fill = n_dates)) +
    facet_wrap(. ~ year, nrow = 1) +
    geom_tile() +
    scale_fill_brewer(pal = 'Reds')

ggsave(filename = 'fastq_byyear.png', plot = pp, path = './plots/SCRATCH/', width = 10, height = 4.8)
table(sample_tally2$n_years)

pp <- ggplot(data = sample_tally2, mapping = aes(x = plot, y = site, fill = factor(n_years))) +
    geom_tile() +
    scale_fill_brewer(palette = 'Greens') +
    labs(fill = '# of years') +
    theme_bw(base_size = 14)

ggsave(filename = 'fastq_byyear.png', plot = pp, path = './plots/SCRATCH/', width = 5, height = 5)

# output files --------

# cleanup --------
date()
closeAllConnections()
