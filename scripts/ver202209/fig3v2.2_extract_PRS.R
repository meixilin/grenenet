# Title: Pull out the PRS scores
# Author: Meixi Lin
# Date: Thu Sep  1 13:52:06 2022

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet")
date()
sessionInfo()

# def functions --------
# TODO: can be optimized
read_gemma <- function(filepath, pos) {
    con <- file(filepath, 'r')
    ii = 1
    # make an empty list
    dtlines <- vector('list', length = length(pos)+1)
    while(TRUE) {
        line <- readLines(con = con, n = 1)
        if (ii == 1) {
            dtlines[[ii]] = line
            ii = ii+1
        } else{
            if (as.integer(strsplit(line, '\t')[[1]][3]) %in% pos) {
                dtlines[[ii]] = line
                ii = ii+1
            }
        }
        # stop when we are at the second chromosome
        if (substr(line,1,1) == '2') {
            break
        }
    }
    close(con)

    # convert the dtlines into delimited files
    outdt = do.call(cbind, (lapply(dtlines, function(xx){strsplit(xx,'\t')[[1]]}))) %>%
        as.data.frame() %>%
        tibble::column_to_rownames(var = 'V1') %>%
        t() %>%
        as.data.frame() %>%
        dplyr::mutate(beta = as.numeric(beta), ps = as.integer(ps))
    return(outdt)
}

loadSomeRData <- function(x, file) {
    E <- new.env()
    load(file=file, envir=E)
    return(get(x, envir=E, inherits=F))
}


# def variables --------
envvar = 'bio1'

# load data --------
# load the significant sites from the bio1
lmres = read.csv(file = paste0('plots/ver202209/AF/DeltaP_v2/', envvar, 'scaled_deltaP_', envvar, '_lmres_Ppass.csv'), row.names = 1)

# TODO: GEMMA LMM might not be the best.
# # Tried to use the BSLMM output
# gemma_out = read.delim(file = '~/safedata/natvar/phenotypes/Exposito-Alonso_Nature_2019_PID_31462776/1001/bslmm_rSurvival_fruit_mlp/output/Exposito-Alonso_Nature_2019_PID_31462776_bslmm_rSurvival_fruit_mlp.assoc.txt.param.txt')
# intersect(gemma_out$ps, lmres$POS) # No intersection of sites

# load the GEMMA LMM output from the natvar 2019 experiments
gemma_out = read_gemma(filepath = '~/safedata/natvar/phenotypes/Exposito-Alonso_Nature_2019_PID_31462776/1001/rSurvival_fruit_mlp/output/Exposito-Alonso_Nature_2019_PID_31462776.lmm.assoc.txt', pos = lmres$POS)

# get the scaled deltaP
deltadt = readRDS(file = 'data/AF/ver202209/haplotype/DeltaP/scaled_delta_freq_uniq.rds')
deltadt = deltadt[rownames(deltadt) %in% as.character(lmres$POS),]

# get the environment
metadt = loadSomeRData(x = 'metadt', file = paste0('plots/ver202209/AF/DeltaP_v2/', envvar, '/scaled_deltaP_', envvar, '.RData'))

# main --------
# calculate PRS for all the sites
prsdt = lapply(1:nrow(deltadt), function(ii) {
    deltap = deltadt[ii,]
    pos = rownames(deltadt)[ii]
    if (any(gemma_out$ps == as.integer(pos))) {
        beta = gemma_out[gemma_out$ps == as.integer(pos), 'beta_all']
        prs = beta*deltap
    } else{
        prs = deltap; prs[1,] = NA
    }
})
prsdt = dplyr::bind_rows(prsdt)

# get colSums for each site
site_prs = colSums(prsdt)

# join dt
dt = dplyr::left_join(lmres, gemma_out, by = c('POS' = 'ps'))

all(metadt$siteday == names(site_prs))

# plot the site env with the site_prs
forplot = data.frame(env = metadt[,envvar], PRS = $)
pp <-
plot(metadt$bio1, site_prs) # not correlated?

#


# output files --------
save.image()

# cleanup --------
date()
closeAllConnections()
