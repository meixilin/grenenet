# Title: Load allele frequency data and convert to an RDS
# Author: Meixi Lin
# Date: Thu Aug  4 00:21:08 2022
# Rscript --vanilla ./scripts/ver202209/step0.2_load_AF.R &> ./logs/ver202209/step0.2_load_AF.log

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet")

date()
sessionInfo()

# def functions --------
# load hafpipe derived sites
# TODO: quick and dirty solution
load_af <- function(filelist,samplelist) {
    dtlist <- plyr::laply(filelist, function(xx) {read.csv(file = xx, colClasses = c('NULL', 'numeric'))},
                          .progress = "none")
    dt <- dplyr::bind_cols(dtlist)
    colnames(dt) = samplelist
    return(dt)
}

# def variables --------
indir = './data-raw/AF/ver202209/haplotype/'
outdir = './data/AF/ver202209/haplotype/'

# load metadata
meta0922 <- readRDS(file = './data/metadata/ver202209/metadata_0922.rds')

# list of samples to load
sample_filelist <- paste0(indir, paste0('samples_frequency/', meta0922$id, '-1.sorted.bam.1.afSite'))

# load data --------
afdt = load_af(filelist = sample_filelist, samplelist = meta0922$id)
afdt
class(afdt)
# convert back to a plain data frame
afdt = as.data.frame(afdt)
class(afdt)

# load the additional seed mix
seedmixdt = read.delim(file = paste0(indir, 'seed_mix_chr1_frequency.txt'))
colnames(seedmixdt) = c('pos','seed_mix')

# main --------
# append the position info for afdt in the rownames
rownames(afdt) = seedmixdt$pos

# make another allele frq dt with both info
afdts = cbind(afdt,seedmixdt$seed_mix)
colnames(afdts)[length(colnames(afdts))] = 'seed_mix'
tail(colnames(afdts))

# output files --------
saveRDS(afdt, file = paste0(outdir, 'AF284_0922.rds'))
saveRDS(afdts, file = paste0(outdir, 'AF284pSeed_0922.rds'))
saveRDS(seedmixdt, file = paste0(outdir, 'AFSeed_0922.rds'))

# cleanup --------
date()
closeAllConnections()
