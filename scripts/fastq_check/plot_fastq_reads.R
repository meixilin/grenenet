# Title: Plot fastq reads
# Author: Meixi Lin
# Date: Thu Jan 26 10:53:36 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/central/groups/carnegie_poc/meixilin/grenenet/analyses/")

library(ggplot2)
library(stringr)
library(dplyr)

# def functions --------

# def variables --------
arag = 135*1e+6

# load data --------
fastqcdt = read.delim(file = 'data/fastq_check/fastqc_totalseq.txt', header = FALSE) %>% 
    # dplyr::filter(str_detect(V1, 'R1')) %>%
    dplyr::mutate(V3 = V2*2*150/arag)
fastqcs = reshape2::colsplit(fastqcdt$V1, pattern = '-', names = c('sampleid','seqid','trim','readid'))
fastqcdt = cbind(fastqcdt, fastqcs)
colnames(fastqcdt)[1:3] = c('fileid','nreads','coverage')

load('../metadata/data/samples_data.rda')
plotdt = dplyr::left_join(fastqcdt, samples_data, by = 'sampleid')

# main --------

# output files --------
png(filename = 'data/fastq_check/fastqc_coverage.png', width = 500, height = 500)
hist(fastqcdt$coverage)
dev.off()

png(filename = 'data/fastq_check/fastqc_coverage_nflowers.png', width = 500, height = 500)
plot(plotdt$flowerscollected, plotdt$nreads,xlim = c(0,20))
dev.off()

sum(fastqcdt$coverage < 2) # 10% of all the data

# cleanup --------
date()
closeAllConnections()


