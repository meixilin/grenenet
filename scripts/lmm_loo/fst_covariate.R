# Title: 
# Author: Meixi Lin
# Date: Tue Jun 20 15:45:39 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses")
date()
sessionInfo()

# def functions --------

# def variables --------

# load data --------
load('./data-raw/snp_freq/snp_sampleinfo_745.rda')

# main --------
colnames(mergep) = paste0(colnames(mergep), '.freq')
mergep$chr = '1'
mergep$pos = seq(1,1000)

mergep = mergep[,c(747,746,1:745)]

# output files --------
write.csv(mergep, file = './data/mergep.csv', quote = FALSE, row.names = FALSE)
write.table(snp_samples$total_flower_counts, file = './data/flowers.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)
mergep = mergep[,-c(1,2)]

which(apply(mergep, 1, function(xx){all(xx == 0)}))

# grenedalf --------
grenedalf fst --window-type=genome --method=unbiased-nei --frequency-table-path 'mergep.csv' --pool-sizes 'flowers.txt'

# cleanup --------
date()
closeAllConnections()
