# Title: Load allele frequency and other data, convert to an RDS
# Author: Meixi Lin
# Date: Wed May 31 14:06:04 PDT 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses/")

date()
sessionInfo()

library(data.table)
setDTthreads(6)

# def functions --------

# def variables --------
# source data directory # 745 binned sites
deltapf = '/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/ath_evo/grenephase1-data/snp_frequency/haf-pipe/merged_delta_p.csv'
p0f = '/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/ath_evo/grenephase1-data/snp_frequency/haf-pipe/p0_genomewide_frequency.txt'
# ld prune directory
ldf = '/NOBACKUP/scratch/xwu/grenet/genomic_evolution/plink.prune.in'
sampf = '/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/ath_evo/grenephase1-data/meta_table/merged_sample_table.csv'

# load deltap --------
deltap = data.table::fread(input = deltapf)

# look at NA occurrences (only 3096 rows out of 3363913 rows)
allnarows = apply(deltap, 1, function(xx){all(is.na(xx))})
anynarows = apply(deltap, 1, function(xx){any(is.na(xx))})
table(allnarows)
table(anynarows)
all.equal(allnarows, anynarows)

# remove Inf values 
deltap[, names(deltap) := lapply(.SD, function(xx) replace(xx, is.infinite(xx), NA))]
allnarows = apply(deltap, 1, function(xx){all(is.na(xx))})
all.equal(allnarows, anynarows)
rm(anynarows)

# load p0 --------
p0 = data.table::fread(input = p0f) 
setnames(p0, c('CHR', 'POS', 'p0'))
p0[, CHR_POS := paste(CHR, POS, sep = '_')]
setcolorder(p0, c('CHR_POS', 'CHR', 'POS', 'p0'))

# The deltap NAs (after correcting for Inf) are the 0/1 AF (invariant sites) in the starting mix
table(p0[allnarows, 'p0']) 
invarrows = (p0$p0 == 0 | p0$p0 == 1)
all.equal(allnarows, invarrows)

# load LD pruned locations --------
# 24337 locations
prunedloc = read.delim(file = ldf, header = FALSE) 
# which rowid is the prunedloc 24273 rows
prunedrow = which(p0$CHR_POS %in% prunedloc$V1)
# which are not included, likely due to hafpipe artifacts
setdiff(prunedloc$V1, p0$CHR_POS)

# prune the deltap
deltap_p = deltap[prunedrow, ]
p0_p = p0[prunedrow, ]

# load the sample info --------
snp_samples = read.csv(file = sampf)
summary(snp_samples)
save(snp_samples, file = 'data-raw/snp_freq/snp_sampleinfo_745.rda')

# output files --------
# 3363913 rows x 745 columns. rows are genomic sites. columns are samples. 
save(deltap, file = 'data-raw/snp_freq/merged_delta_p_745.rda') 
# 3363913 rows x 3 columns. columns 1 and 2 are CHR and POS. column 3 is the p0 from 231 founders. 
save(p0, file = 'data-raw/snp_freq/seedmix_p0_231.rda')

# 24273 rows x 745 columns. rows are genomic sites. columns are samples. 
str(deltap_p)
save(deltap_p, file = 'data-raw/snp_freq/merged_delta_p_745_ldpruned.rda') 
# 24273 rows x 3 columns. columns 1 and 2 are CHR and POS. column 3 is the p0 from 231 founders. 
str(p0_p)
save(p0_p, file = 'data-raw/snp_freq/seedmix_p0_231_ldpruned.rda')

# plot some basics --------
# plot histogram of p0
png(filename = 'data-raw/snp_freq/seedmix_p0_231_hist.png')
hist(p0$p0)
dev.off()

png(filename = 'data-raw/snp_freq/seedmix_p0_231_hist_LDpruned.png')
hist(p0_p$p0)
dev.off()

# cleanup --------
date()
closeAllConnections()

# SCRATCH --------
# this time it's only 745 columns because site 57 is now defined as generations 1 to 3
# load('data/AF/merged_delta_p_776.rda')
# n776 = colnames(deltap)
# > setdiff(n776,n745)
# [1] "57_4_1"  "57_4_3"  "57_4_4"  "57_4_6"  "57_4_7"  "57_4_10" "57_4_12" "57_4_13"
# [9] "57_4_15" "57_4_16" "57_4_18" "57_6_1"  "57_6_3"  "57_6_4"  "57_6_6"  "57_6_7" 
# [17] "57_6_10" "57_6_12" "57_6_13" "57_6_16" "57_5_3"  "57_5_4"  "57_5_6"  "57_5_7" 
# [25] "57_5_9"  "57_5_10" "57_5_12" "57_5_13" "57_5_15" "57_5_16" "57_5_18"
