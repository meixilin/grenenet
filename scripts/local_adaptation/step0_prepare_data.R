# Title: Prepare data for ecotype cline analyses
# Author: Meixi Lin
# Date: Mon Aug 28 18:57:35 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses")

library(dplyr)
date()
sessionInfo()

# def functions --------
# wrapper to read and calculate starting ecotype frequencies
read_ecop0 <- function(indir) {
    founder <- as.data.frame(matrix(nrow = 231,ncol=8))
    for(i in 1:8){
        name <- paste(indir, "s",i,"_ecotype_frequency.txt",sep = "")
        tmp <- read.delim(name,header=F)
        founder[,i]<- tmp$V2
    }
    founder_freq <- apply(founder,1,mean)
    outdf = data.frame(ecotypeid = tmp$V1,
                       p0 = founder_freq)
    return(outdf)
}

euclidean_denv <- function(row, ecotypedf) {
    return((row - ecotypedf)^2)
}

raw_denv <- function(row, ecotypedf) {
    return((row - ecotypedf))
}

# calculate distance in the forms of (x-x0)^2
calc_env_dist <- function(ecotypedf, sitedf, myfun) {
    sitedist = apply(sitedf[,-1], 1, myfun, ecotypedf = ecotypedf) %>% t()
    sitedist = data.frame(sitedist) 
    sitedist$site = sitedf[,1] 
    return(sitedist)
}

# pca on bio1 to bio19
pca_bio <- function(env_ecotypes, env_sites) {
    # extract only bio1 to bio19 variables
    bio_ecotypes <- env_ecotypes %>%  
        dplyr::mutate(ecotypeid = paste0('eco', ecotypeid)) %>% 
        tibble::column_to_rownames(var = 'ecotypeid') %>% 
        dplyr::select(starts_with('bio'))
    bio_sites <- env_sites %>% 
        dplyr::mutate(site = paste0('site', site)) %>% 
        tibble::column_to_rownames(var = 'site') %>% 
        dplyr::select(starts_with('bio'))
    # run PCA
    pcabio <- prcomp(rbind(bio_ecotypes, bio_sites), scale. = TRUE)
    scoresbio <- pcabio$x
    return(scoresbio)
}

calc_pca_dist <- function(scoresbio, neco, nsite, deltaenv) {
    # calculate euclidean distance^2 on PCA space
    bioall <- sapply(1:neco, function(ii) {
        sapply(1:nsite, function(jj) {
            scorespair <- scoresbio[c(ii, jj+neco), ]
            tmpdist <- as.numeric(dist(scorespair, method = 'euclidean'))
            return(tmpdist^2)
        }) 
    }) 
    
    colnames(bioall) = as.integer(stringr::str_remove(dimnames(scoresbio)[[1]][1:neco], 'eco'))
    rownames(bioall) = as.integer(stringr::str_remove(dimnames(scoresbio)[[1]][(neco+1):(neco+nsite)], 'site'))
    
    deltaenv_a = bioall %>% 
        reshape2::melt()
    
    # sanity check
    if (!(all(deltaenv$site == deltaenv_a$Var1) & all(deltaenv$ecotypeid == deltaenv_a$Var2))) {
        stop('Mismatching ecotype and site names')
    }
    return(deltaenv_a$value)
}

# def variables --------
outdir = './data/local_adaptation/inputs/'
dir.create(outdir, recursive = TRUE)

# load data --------
# ecotype data
ecop <- read.delim('/NOBACKUP/scratch/xwu/grenet/merged_frequency/merged_ecotype_frequency.txt') # md5: 4272ec578c1c45511f45d9d92916fdb2
summary(colSums(ecop))
ecop0 <- read_ecop0('/NOBACKUP/scratch/xwu/grenet/hapFIRE_frequencies/seed_mix/')
# add one column on the left for ecop ecotypeid
ecop <- cbind(ecop0[,'ecotypeid'], ecop)
colnames(ecop)[1]= 'ecotypeid'
# reorder the ecop to match the metadata
ecop = ecop %>% dplyr::arrange(ecotypeid)
ecop0 = ecop0 %>% dplyr::arrange(ecotypeid)

# load and reformat climate data
load('../metadata/data/worldclim_ecotypesdata.rda')
load('../metadata/data/worldclim_sitesdata.rda')
load('../metadata/data/ecotypes_data.rda')
load('../metadata/data/locations_data.rda')
all(colnames(worldclim_ecotypesdata)[-1] == colnames(worldclim_sitesdata)[-1])

# add longitude and latitude for this study
# TODO no reasonable altitude data available for ecotypes
env_ecotypes = dplyr::left_join(worldclim_ecotypesdata, 
                                ecotypes_data[,c('ecotypeid', 'longitude', 'latitude')], by = 'ecotypeid')
env_sites = dplyr::left_join(worldclim_sitesdata, 
                             locations_data[,c('site', 'longitude', 'latitude')], by = 'site')

# load samples data
load('../metadata/data/merged_samples_data.rda')

# only keep the 31 sites that have any sequence data
env_sites = env_sites %>%
    dplyr::filter(site %in% unique(merged_samples_data$site))

# only analyze the ecotypes that had environmental data (not too close to the oceans)
env_ecotypes = env_ecotypes %>% 
    tidyr::drop_na()

# this will remove the ecotypes that had missing data
# keep the original ecop and ecop0 in RData
# TODO: if wanted, we can obtain worldclim for the ecotypes using a nearest neighbor approach
ecop_231 = ecop
ecop0_231 = ecop0

ecop = ecop %>% 
    dplyr::filter(ecotypeid %in% env_ecotypes$ecotypeid)
ecop0 = ecop0 %>% 
    dplyr::filter(ecotypeid %in% env_ecotypes$ecotypeid)

# check that the ordering is correct
stopifnot(all(env_ecotypes$ecotypeid == ecop$ecotypeid))
stopifnot(all(ecop0$ecotypeid == ecop$ecotypeid))

# main --------
# calculate distance of ecotypes and experimental sites ========
# for each ecotype and each bioclim var, we calculate the euclidean distance to each site
deltaenv = lapply(1:nrow(env_ecotypes), function(ii) {
    tmp = unlist(env_ecotypes[ii,-1]) # env vector for this ecotype
    out = calc_env_dist(ecotypedf = tmp, sitedf = env_sites, myfun = euclidean_denv) %>% 
        dplyr::mutate(ecotypeid = env_ecotypes[ii,'ecotypeid']) %>% 
        dplyr::relocate(ecotypeid, site)
    return(out)
}) %>% dplyr::bind_rows()

# add another deltaenv calculation for plotting
deltaenv_raw = lapply(1:nrow(env_ecotypes), function(ii) {
    tmp = unlist(env_ecotypes[ii,-1]) # env vector for this ecotype
    out = calc_env_dist(ecotypedf = tmp, sitedf = env_sites, myfun = raw_denv) %>% 
        dplyr::mutate(ecotypeid = env_ecotypes[ii,'ecotypeid']) %>% 
        dplyr::relocate(ecotypeid, site)
    return(out)
}) %>% dplyr::bind_rows()

# reshape the ecotype frequency data to a long format ========
ecop_l = ecop %>% 
    reshape2::melt(id.var = 'ecotypeid', value.name = 'frequency', variable.name = 'sample_name') %>% 
    dplyr::mutate(sample_name = stringr::str_remove(sample_name, '^X')) %>% 
    dplyr::left_join(., y = merged_samples_data[, c('sample_name', 'site', 'generation', 'plot')], by = 'sample_name')

# keep this long form in the step0_prepare_data.RData
ecop_l231 = ecop_231 %>% 
    reshape2::melt(id.var = 'ecotypeid', value.name = 'frequency', variable.name = 'sample_name') %>% 
    dplyr::mutate(sample_name = stringr::str_remove(sample_name, '^X')) %>% 
    dplyr::left_join(., y = merged_samples_data[, c('sample_name', 'site', 'generation', 'plot')], by = 'sample_name')

dim(ecop)
dim(ecop_231)

# calculate environmental distances as a whole ========
scoresbio <- pca_bio(env_ecotypes, env_sites)
neco <- sum(grepl('eco', dimnames(scoresbio)[[1]]))
nsite <- dim(scoresbio)[1] - neco

# calculate euclidean distance^2 on PCA space
bioall <- calc_pca_dist(scoresbio, neco, nsite, deltaenv)
deltaenv$bioall = bioall

hist(deltaenv$bioall) # most distances are small 

# TODO: can also use the Mahalanobis Distance

# output files --------
save(ecop0, ecop_l, deltaenv, file = paste0(outdir, 'ecotypefreqs_deltaenv.rda'))
save.image(file = paste0(outdir, 'step0_prepare_data.RData'))

# cleanup --------
date()
closeAllConnections()

