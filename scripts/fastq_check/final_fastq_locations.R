# Title: Final check on the fastq locations
# Author: Meixi Lin
# Date: Fri Jan 27 13:19:38 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/central/groups/carnegie_poc/meixilin/grenenet/analyses")

# def functions --------

# def variables --------

# load data --------
seqtable = read.delim(file = '/central/groups/carnegie_poc/lczech/grenephase1/seqtable_cal.tsv')
load('../metadata/data/fastq_info.rda')

# main --------
fastq_gg = fastq_info %>% 
    dplyr::mutate(fq1 = str_replace(r1srafolder,
                                    '/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/',
                                    '/central/groups/carnegie_poc/Moilab/'),
                  fq2 = str_replace(r2srafolder,
                                    '/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/',
                                    '/central/groups/carnegie_poc/Moilab/')) %>%
    dplyr::select(sampleid,unit,fq1,fq2,platform) %>%
    dplyr::mutate(platform = '-') %>%
    dplyr::rename(sample = sampleid) 
fastq_g = fastq_gg %>%
    dplyr::filter(sample != 'MLFH570120190821')

adddt = fastq_gg %>%
    dplyr::filter(sample == 'MLFH570120190821')

all(fastq_g$sample == seqtable$sample)
all(fastq_g$fq1 == seqtable$fq1)
all(fastq_g$fq2 == seqtable$fq2)
all(fastq_g$unit == seqtable$unit)

# the only difference is the MLFH570120190821 file


# output files --------
write.table(adddt, file = 'data-raw/fastq_check/seqtable_add.tsv', quote = F, row.names = F, sep ='\t')
# cleanup --------
date()
closeAllConnections()

