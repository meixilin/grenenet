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

find_plots <- function(samples_data, mysite, mygen) {
    # already filtered for the usesample
    outdf = samples_data %>% 
        dplyr::filter(site == mysite, generation == mygen)
    outplots = unique(outdf$plot)
    return(outplots)
}

find_samples <- function(mysite, myyear) {
    def find_plots(df: pd.DataFrame,ss: int,yy: int):
        outdf = df.query('site == @ss and year == @yy and usesample == True')
        outplots = outdf['plot'].unique().tolist() # cannot access without quoting
        return outplots
}

weight_ecotypes <- function(mysite, myplot, myyear) {
    
}

# def variables --------
ecop0f = '/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/ath_evo/grenephase1-data/ecotype_frequency/avg_founder_ecotype_frequency.txt'
ecopf = '/NOBACKUP/scratch/xwu/grenet/hapFE/ecotype_frequency/'
ecodeltapf = '/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/ath_evo/grenephase1-data/ecotype_frequency/delta_ecotype_frequency.csv'

# load data --------
load('../metadata/data/ecotypes_data.rda')
load('../metadata/data/samples_data.rda')
load('./data-raw/snp_freq/snp_sampleinfo_745.rda')

samples_data = samples_data %>% dplyr::filter(usesample)
allsites = unique(samples_data$site) # 31 sites

# load starting ecotypes --------
ecop0 = read_ecotypes(filename = ecop0f, samplename = 'p0') 

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


# load ecotype frequency change --------
ecodeltap = read.csv(file = ecodeltapf)

# output files --------
save(ecop0, file = './data-raw/ecotype_freq/ecotype_p0.rds')
save(ecop, file = './data-raw/ecotype_freq/ecotype_freq.rds')
save(ecodeltap, file = './data-raw/ecotype_freq/merged_ecotype_deltap.rds')

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
