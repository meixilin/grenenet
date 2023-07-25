# Title: Load GFF file, generate a dictionary of gene locations
# Author: Meixi Lin
# Date: Tue Jul 18 23:31:21 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses")

library(rtracklayer)
library(GenomicRanges)
date()
sessionInfo()

# def functions --------

# def variables --------
my_tags <- c("ID", "Name", "geneID")

outdir = './data-raw/TAIR10/'
dir.create(outdir)
# dir.create(plotdir)

# load data --------
gff = readGFF(filepath = '/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/arabidopsisthaliana_reference/TAIR10_GFF3_genes_transposons.gff',
              tags = my_tags)

# main --------
# keep only the genes 
table(gff$type, useNA = 'always')
table(gff$seqid, useNA = 'always')
genes = gff[gff$type == 'gene' & 
            gff$seqid %in% c('Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5'),]

# check that there is no duplicate genes
any(duplicated(genes$ID)) # No duplicates

# generate genomic ranges --------
genesGR = GenomicRanges::GRanges(seqnames = genes$seqid,
                                 ranges = IRanges(start = genes$start, end = genes$end),
                                 strand = genes$strand, 
                                 ID = genes$ID)
genesGR

# get a chromosome cumsum table --------
chrdf = gff[gff$type == 'chromosome', ]
chrdf$CHR = c(1,2,3,4,5,NA,NA)
chrdf$sumstart = c(1, cumsum(chrdf$end) + chrdf$start)[1:7]

# output files --------
saveRDS(genes, file = paste0(outdir, 'TAIR10_genes_gff.rds'))
saveRDS(genesGR, file = paste0(outdir, 'TAIR10_genesGR_gff.rds'))
write.csv(chrdf, file = paste0(outdir, 'TAIR10_chr_len.csv'))

# cleanup --------
date()
closeAllConnections()


