# Title: Test why there are over representation for certain variant types
# Author: Meixi Lin
# Date: Thu Mar  2 14:00:52 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses")

library(data.table)
setDTthreads(threads = 6)

date()
sessionInfo()

# def functions --------

# def variables --------

# load data --------
load('data/AF/seedmix_p0_231.rda')
p0 = p0[CHR == 1, ] # only chromosome 1

# load the merged allele frequency for the first chromosome
testp = data.table::fread(input = '/NOBACKUP/scratch/xwu/grenet/haf-pipe/hafpipe_231/merged_allele_frequency.csv', 
                          nThread = 6, nrows = 894424) 
testp[1:5,1:5]

# main --------
# start from p0 inference
p0[, 'id'] = 1:nrow(p0)
setorder(p0, p0, CHR)

# split by p0
p0l = split(p0, by = 'p0')

# remove the ones that have only one value
p0len = sapply(p0l, nrow)
p0tb = sort(table(p0len), decreasing = T)
head(p0tb)
write.csv(p0tb, file = '/NOBACKUP/scratch/meixilin/deltap_qc/p0_tallyCHR1.csv')

# get the index that have the duplicated values
p0l = p0l[unname(p0len > 1)]

# get the p0$id that are duplicated by group
p0id = lapply(p0l, function(xx) {
    return(xx$id)
})

head(p0id)

# look into the testp and see if it's the same 
testpunique = sapply(p0id, function(ii) {
    tempdt = testp[ii,]
    temphash = apply(tempdt, 1, paste, collapse = '')
    mylen = length(unique(temphash)) 
    if (mylen > 1) {
        write.csv(p0[id %in% ii, ], file = paste0('/NOBACKUP/scratch/meixilin/deltap_qc/p0_', ii[1], '.csv'))
        write.csv(tempdt, file = paste0('/NOBACKUP/scratch/meixilin/deltap_qc/merged_p_', ii[1], '.csv'))
        print(head(p0[id %in% ii, ]))
        print(head(tempdt[, 1:5]))
    }
    return(mylen)
})

testptb = sort(table(testpunique), decreasing = T)
head(testptb)
write.csv(testptb, file = '/NOBACKUP/scratch/meixilin/deltap_qc/testp_tallyCHR1.csv')

# if you have the same starting allele frequency, then most of the time, 
# all the sites are having the same allele frequency. 
# these sites are not right next to each other, they are completely the same because of the HAFpipe setup

# explain why this happened --------
library(Rlab)

# five founder haplotypes, 10 SNPs
set.seed(7)
snptable = matrix(nrow = 5, ncol = 10, data = rbern(n = 50, prob  = 0.5))
heatmap(snptable, Rowv = NA, Colv = NA)

# twenty sample, five haplotype frequency
set.seed(7)
samplehap = matrix(nrow = 20, ncol = 5, data = runif(n = 20*5))
heatmap(samplehap, Rowv = NA, Colv = NA)

# twenty sample, 10 SNPs (output files)
samplesnps = samplehap %*% snptable
# in SNP positions 1 and 7, 4 and 6 the allele frequency is 100% matching across all samples. 
# because the SNPxhaplotype table is the same across haplotypes. 

# if you have the same founder vcf columns, then you are going to get the same rows in the end if you are in the same haplotype window, in other words, the SNPs are bined in the 100 bp window. 

# cleanup --------
date()
closeAllConnections()
