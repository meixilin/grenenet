# Title:
# Author: Meixi Lin
# Date: Thu Aug 25 13:48:27 2022

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet")

library(rtracklayer)
date()
sessionInfo()

# def functions --------

# def variables --------
my_tags <- c("ID", "Name", "geneID")

outdir = './data/TAIR10/'
dir.create(outdir)
# dir.create(plotdir)

# load data --------
gff = readGFF(filepath = '/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/arabidopsisthaliana_reference/TAIR10_GFF3_genes_transposons.gff',
              tags = my_tags)

# main --------
# keep only the genes and only things on chr1
table(gff$type, useNA = 'always')
table(gff$seqid, useNA = 'always')
genes = gff[gff$type == 'gene' & gff$seqid == 'Chr1',]

# check that there is no duplicate genes
any(duplicated(genes$ID)) # YES!

# output files --------
saveRDS(genes, file = paste0(outdir, 'Chr1_genes_gff.rds'))

# cleanup --------
date()
closeAllConnections()
