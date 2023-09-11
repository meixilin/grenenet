# Title: Run linear mixed models in parallel
# Author: Meixi Lin
# Date: Thu Sep  7 12:02:10 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses")

library(dplyr)
library(data.table)
library(doParallel)
cl <- makeCluster(20)
registerDoParallel(cl)

date()
sessionInfo()

# def functions --------
# column headers are defined in env_sites
# deltaps: scaled and merged change in allele frequnecy [dp/(p0(1-p0))]
# total rows: 3235480 no batching performed
load_deltaps <- function(env_sites) {
    # this file had a header V1,V2,V3, so you need to skip = 1 to read the first SNP (but fread also automatically find the header)
    infile = '/NOBACKUP/scratch/xwu/grenet/merged_frequency/merged_hapFIRE_delta_p_normed.txt' 
    dt = data.table::fread(file = infile, nThread = 20, skip = 1)
    # get real mergeid headers
    mergeids = colnames(read.csv("/NOBACKUP/scratch/xwu/grenet/merged_frequency/merged_hapFIRE_allele_frequency.csv", nrows = 1)) %>% 
        stringr::str_remove(.,'X') 
    myind = which(mergeids %in% env_sites$mergeid)
    stopifnot(all(mergeids[myind] == env_sites$mergeid)) # check sampling will not change orders
    # rows are SNPs, columns are samples with mergeid falling in given generation (reflected in env_sites)
    outdt = as.matrix(dt[, ..myind])
    dimnames(outdt)[[2]] = NULL
    outdt = t(outdt)
    print(dim(outdt))
    return(outdt)
}

# assemble lmm objects
prep_lmm <- function(yy, env_sites, envvar, pop_strc) {
    mydata = cbind(yy, env_sites[,c('mergeid', 'site', envvar)], pop_strc)
    return(mydata)
}

# get lmm model results
format_lmm <- function(mymodel, envvar) {
    lmesum = summary(mymodel)
    lmer2 = MuMIn::r.squaredGLMM(mymodel) # contains two rows
    outdt = c(lmer2,
              lmesum$tTable[envvar, "Value"],
              lmesum$tTable[envvar, "p-value"],
              lmesum$BIC)
    return(outdt)
}

# def variables --------
# not batching, run full genome in one run
args = commandArgs(trailingOnly=TRUE)
envvar = as.character(args[1])
mygen = as.character(args[2])
# envvar = 'bio1'
# mygen = '2'

myfm0 = as.formula(paste0('yy ~ ', envvar))
myfm1 = as.formula(paste0('yy ~ ', paste(c(paste0('PC',1:10), envvar), collapse = ' + ')))
myfm2 = as.formula(paste0('yy ~ ', paste(c(paste0('PC',1:5), envvar), collapse = ' + ')))
print(myfm0); print(myfm1); print(myfm2)

outdir = '/NOBACKUP/scratch/meixilin/grenenet/lmm_pc/testpc/'
dir.create(outdir, recursive = TRUE)

# load data --------
load(paste0('./data/lmm_pc/inputs/envsites_popstrc_gen', mygen, '.rda'))
deltaps = load_deltaps(env_sites)
date()

# main --------
# run the linear mixed model in a parallel for loop
lmeres <- foreach(yy = deltaps, .combine='rbind') %dopar% {
    .GlobalEnv$myfm0 <- myfm0
    .GlobalEnv$myfm1 <- myfm1 # fix a global env bug
    .GlobalEnv$myfm2 <- myfm2
    mydata = prep_lmm(yy, env_sites, envvar, pop_strc) 
    model0 = nlme::lme(fixed = myfm0, random = ~ 1|site, data = mydata) # old model without pop strc
    model1 = nlme::lme(fixed = myfm1, random = ~ 1|site, data = mydata) # 10 popstr PCs
    model2 = nlme::lme(fixed = myfm2, random = ~ 1|site, data = mydata) # 5 popstr PCs
    unlist(lapply(list(model0, model1, model2), format_lmm, envvar = envvar))
}

dimnames(lmeres)[[2]] = c('R2m.0', 'R2c.0', 'beta.0', 'beta_p.0', 'BIC.0',
                          'R2m.1', 'R2c.1', 'beta.1', 'beta_p.1', 'BIC.1',
                          'R2m.2', 'R2c.2', 'beta.2', 'beta_p.2', 'BIC.2')
print(dim(lmeres))

# preliminary plotting
png(filename = paste0(outdir, 'QQlme_', envvar, '_gen', mygen, '.png'), width = 9, height = 4, res = 300, units = 'in')
par(mfrow = c(1,3))
qqman::qq(lmeres[,4], main = paste0('yy ~ ', envvar))
qqman::qq(lmeres[,9], main = paste0('yy ~ PC1:10 + ', envvar))
qqman::qq(lmeres[,14], main = paste0('yy ~ PC1:5 + ', envvar))
dev.off()

# output files --------
save(lmeres, file = paste0(outdir, 'lmeres_', envvar, '_gen', mygen, '.rda'))

# cleanup --------
stopCluster(cl)
date()
closeAllConnections()

# ARCHIVE --------
# batched reading files
# load_deltaps <- function(batchid, batchsize, env_sites) {
#     start = (batchid-1)*batchsize + 1
#     # this file had a header V1,V2,V3, so you need to skip = 1 to read the first SNP
#     infile = '/NOBACKUP/scratch/xwu/grenet/merged_frequency/merged_hapFIRE_delta_p_normed.txt' 
#     dt = data.table::fread(file = infile, skip = start, nrows = batchsize)
#     # get real mergeid headers
#     mergeids = colnames(read.csv("/NOBACKUP/scratch/xwu/grenet/merged_frequency/merged_hapFIRE_allele_frequency.csv", nrows = 1)) %>% 
#         stringr::str_remove(.,'X') 
#     myind = which(mergeids %in% env_sites$mergeid)
#     stopifnot(length(myind) == nrow(env_sites))
#     # rows are SNPs, columns are samples with mergeid falling in given generation (reflected in env_sites)
#     outdt = as.matrix(dt[, ..myind])
#     dimnames(outdt)[[2]] = NULL
#     outdt = t(outdt)
#     return(outdt)
# }



