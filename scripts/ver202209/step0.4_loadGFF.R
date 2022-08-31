# Title: Load GFF file, generate a dictionary of gene locations
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
# make a dataframe for gene coordinates
format_genecoords <- function(genecoords) {
    outdtl = vector('list',length=length(genecoords))
    for (ii in seq_along(genecoords)) {
        if (length(genecoords[[ii]]) == 0) {
            next
        } else {
            outdtl[[ii]] = data.frame(gene = names(genecoords)[ii],
                                      POS = genecoords[[ii]])
        }
    }
    outdt = dplyr::bind_rows(outdtl) %>%
        dplyr::mutate(POS = as.integer(POS))
    return(outdt)
}

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

# generate gene coords --------
# See step 1.0_loadDeltaAF.R

# append the genecoords --------
format_genecoords <- function(genecoords) {
    outdtl = vector('list',length=length(genecoords))
    for (ii in seq_along(genecoords)) {
        if (length(genecoords[[ii]]) == 0) {
            next
        } else {
            outdtl[[ii]] = data.frame(gene = names(genecoords)[ii],
                                      POS = genecoords[[ii]])
        }
    }
    outdt = dplyr::bind_rows(outdtl) %>%
        dplyr::mutate(POS = as.integer(POS))
    # if there are multiple genes per location, bin them
    cleandt = outdt %>%
        dplyr::group_by(POS) %>%
        dplyr::summarise(GENES = paste(gene, collapse = ';'))
    if (any(duplicated(cleandt$POS))) {
        stop('Duplicated POS')
    }
    return(cleandt)
}

genecoords = readRDS(file = 'data/AF/ver202209/haplotype/DeltaP/genecoords.rds')
genecoordsdt = format_genecoords(genecoords)

# output files --------
saveRDS(genes, file = paste0(outdir, 'Chr1_genes_gff.rds'))
saveRDS(genecoordsdt, file = 'data/AF/ver202209/haplotype/DeltaP/genecoordsdt.rds')

# cleanup --------
date()
closeAllConnections()


