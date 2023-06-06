# Title: Load ecotypes frequencies
# Author: Meixi Lin
# Date: Thu Jun  1 14:58:11 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses")
date()
sessionInfo()

library(ggplot2)
library(dplyr)

# def functions --------
read_ecotypes <- function(filename, samplename) {
    df = read.delim(file = filename, header = FALSE)
    colnames(df) = c('ecotypeid', samplename)
    df = df %>% 
        dplyr::arrange(ecotypeid)
    return(df)
}

my_weighted_mean <- function(eco_freq, nflowers) {
    out_freq = stats::weighted.mean(eco_freq, nflowers)
    return(out_freq)
}

calc_weighted_mean <- function(ecop, merged_samples_data, samples_data) {
    # set up output data frame, row is ecotype by id, column is merged sample_name
    outdf = apply(merged_samples_data, 1, function(thissample) {
        sampleidlist = strsplit(thissample['sampleidlist'], ';')[[1]]
        subecop = as.matrix(ecop[, sampleidlist])
        subsamples_data = samples_data[samples_data$sampleid %in% sampleidlist, ]
        # confirm that the order is correct
        if (!(all(colnames(subecop) == subsamples_data$sampleid))) {
            stop('samples_data and ecop not matching')
        }
        out_freq = apply(subecop, 1, my_weighted_mean, nflowers = subsamples_data$flowerscollected)
        return(out_freq)
    })
    dimnames(outdf)[[2]] = paste0('X', merged_samples_data$sample_name) 
    
    return(outdf)
}

# def variables --------
ecop0f = '/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/ath_evo/grenephase1-data/ecotype_frequency/avg_founder_ecotype_frequency.txt'
ecopf = '/NOBACKUP/scratch/xwu/grenet/hapFE/ecotype_frequency/'
ecodeltapf = '/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/ath_evo/grenephase1-data/ecotype_frequency/delta_ecotype_frequency.csv'

# load data --------
load('../metadata/data/merged_samples_data.rda')
load('../metadata/data/ecotypes_data.rda')
load('../metadata/data/samples_data.rda')

# load starting ecotypes --------
ecop0 = read_ecotypes(filename = ecop0f, samplename = 'p0') 
str(ecop0)

# confirm that the ecotypeid are the same Yes. 
all(ecop0$ecotypeid == ecotypes_data$ecotypeid) 

# load ecotype frequency per sample --------
ecopfs = list.files(path = ecopf)
samps = sapply(ecopfs, function(xx) {sub("^(.*?)_.*", "\\1", xx)})
ecopl = lapply(ecopfs, function(xx) {
    df = read_ecotypes(filename = paste0(ecopf, xx), samplename = 'temp')
    return(df$temp)
}) 
ecop = cbind(ecop0$ecotypeid, dplyr::bind_cols(ecopl))
colnames(ecop) = c('ecotypeid', samps)
str(ecop)

# check that the samples are the same as the samples in the metadata
# yes, all samples were included in this case
setdiff(samps, samples_data$sampleid) 
setdiff(samples_data$sampleid, samps)

# merge ecotype frequency by the number of flowers sampled -------
ecomergedp = calc_weighted_mean(ecop, merged_samples_data, samples_data)

# load ecotype frequency change --------
ecodeltap = read.csv(file = ecodeltapf)
# the order for ecodeltapf is sorted based on ecop0f's original orders

# reorder the ecodeltap sort the ecotypeid
ecop0r = read.delim(file = ecop0f, header = FALSE)
ecodeltap = cbind(ecop0r[, 'V1'], ecodeltap)
colnames(ecodeltap)[1] = 'ecotypeid'
ecodeltap = ecodeltap %>% dplyr::arrange(ecotypeid)
str(ecodeltap)

# do a sanity check on the frequency changes
all(sort(dimnames(ecodeltap)[[2]]) == sort(dimnames(ecomergedp)[[2]]))
# reorder the ecodeltap by site, generation, plot
ecodeltap = ecodeltap[, dimnames(ecomergedp)[[2]]]

# calculate deltap from ecomergedp
# SUCCESS! Reproduced Xing's results. 
ecodeltap1 = ecomergedp - ecop0$p0
all.equal(as.matrix(ecodeltap[,-1]), ecodeltap1)

# reshape the merged ecotype frequency to site, generation and plots --------
ecomergedp = cbind(ecop0[, 'ecotypeid'], as.data.frame(ecomergedp))
colnames(ecomergedp)[1] = 'ecotypeid'

ecomergedp_long = ecomergedp %>% 
    reshape2::melt(id.var = 'ecotypeid', value.name = 'frequency', variable.name = 'sample_name') %>% 
    dplyr::mutate(sample_name = stringr::str_remove(sample_name, '^X')) %>% 
    dplyr::left_join(., y = merged_samples_data[, c('sample_name', 'site', 'generation', 'plot')], by = 'sample_name')

# output files --------
save(ecop0, file = './data-raw/ecotype_freq/ecotype_p0.rds')
save(ecop, file = './data-raw/ecotype_freq/ecotype_freq.rds')
save(ecodeltap, file = './data-raw/ecotype_freq/merged_ecotype_deltap.rds')
save(ecomergedp_long, file = './data-raw/ecotype_freq/merged_ecotype_freq_long.rds')

# plot the starting seeds --------
ecop0pair = cbind(ecop0, ecotypes_data$seedsperplot/sum(ecotypes_data$seedsperplot))
colnames(ecop0pair) = c('ecotypeid', 'hapfire', 'weight')
summary(ecop0pair)
ecop0pair %>% slice_max(hapfire, n = 5)
ecop0pair %>% slice_min(hapfire, n = 5)
ecop0pair %>% slice_max(weight, n = 5)
ecop0pair %>% slice_min(weight, n = 5)

ecop0pair = ecop0pair %>% reshape2::melt(id.var = 'ecotypeid')
pp = ggplot(data = ecop0pair, aes(x = value, fill = variable)) +
    geom_histogram(bins = 100) +
    labs(x = 'ecotype frequency', y = '# of ecotypes')

ggsave(filename = './data-raw/ecotype_freq/ecotype_p0_hapfire.pdf',plot = pp, height = 4, width = 4)

# cleanup --------
date()
closeAllConnections()

