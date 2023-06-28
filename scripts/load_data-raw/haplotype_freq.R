# Title: Load the haplotype frequency data 
# Author: Meixi Lin
# Date: Wed Jun 14 14:41:55 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses")
date()
sessionInfo()

library(data.table)
setDTthreads(6)

# def functions --------
get_founder_freq <- function(haplop0d) {
    founder <- as.data.frame(matrix(nrow = 39313,ncol=8))
    for(i in 1:8){
        name <- paste(haplop0d, "s",i,"_clustered_haplotype_frequency.txt",sep = "")
        tmp <- read.delim(name,header=F)
        founder[,i]<- tmp$V2
    }
    founder_headers = reshape2::colsplit(tmp$V1, pattern = '_haplotype_cluster', 
                                         names = c('CHR_RANGE', 'CLUSTER'))
    p0 <- apply(founder,1,mean)
    founder_freq <- cbind(founder_headers, p0)
    return(founder_freq)
}

# def variables --------
haplopf = '/NOBACKUP/scratch/xwu/grenet/selection_signals/haplotype/merged_haplotype_frequency.csv'

haplop0d = '/NOBACKUP/scratch/xwu/grenet/hapFE/seedmix/'
    
# load data --------
haplop = data.table::fread(haplopf)

haplop0 = get_founder_freq(haplop0d)

# main --------
haplop[1:5, 1:5]
haplop0[1:5, ]

# output files --------
# 39313 rows x 745 columns. rows are genomic sites. columns are samples. 
save(haplop, file = 'data-raw/haplotype_freq/merged_haplotype_p_745.rda') 
save(haplop0, file = 'data-raw/haplotype_freq/seedmix_haplotype_p0_231.rda') 

# cleanup --------
date()
closeAllConnections()

# SCRATCH ---------
# sanity check the haplotype id (TODO)
