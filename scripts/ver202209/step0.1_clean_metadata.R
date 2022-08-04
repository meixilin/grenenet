# Title: Load and clean metadata needed for september plotting
# Author: Meixi Lin
# Date: Wed Aug  3 12:06:43 2022
# Rscript --vanilla ./scripts/ver202209/step0.1_clean_metadata.R &> ./logs/ver202209/step0.1_clean_metadata.log

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet")

library(ggplot2)
library(dplyr)

sessionInfo()

# def functions --------
read_samples <- function(filename) {
    temp <- readLines(con =  filename)
    outchr <- sapply(temp, function(xx) {base::strsplit(xx, split = " ")[[1]][1]})
    outchr <- unname(outchr)
    return(outchr)
}

match_records <- function(id) {
    # split before the date
    id0 = base::substr(id,1,8)
    # find all values with id0 in samples.fam
    sampleid0 = samples[stringr::str_starts(samples, id0)]
    # find all values with id0 in recordssorted
    recordsid0 = recordssorted[stringr::str_starts(recordssorted$id, id0), 'id']
    # set differences
    diffid_s = setdiff(sampleid0, recordsid0)
    diffid_r = setdiff(recordsid0, sampleid0)
    # make sure there is only one mismatch
    if (any(diffid_s != id)) {
        # print(id)
        outdt = data.frame(id = id,
                           id_s.newid = paste(diffid_s, collapse = ' '),
                           id_r.diff_dates = paste(diffid_r, collapse = ' '),
                           resolved = FALSE)
    } else {
        # find similar dates
        id_date = as.Date(base::substr(id,9,16), format = '%Y%m%d')
        diffid_r_date = as.Date(base::substr(diffid_r,9,16), format = '%Y%m%d')
        diffdates = abs(as.integer(difftime(time1= id_date, time2 = diffid_r_date, units = 'days')))
        names(diffdates) = diffid_r
        # return the one with the minimal differences in dates
        outdt = data.frame(id = id,
                           id_s.newid = names(diffdates)[diffdates == min(diffdates)],
                           id_r.diff_dates = as.character(min(diffdates)),
                           resolved = TRUE)
    }
    return(outdt)
}

# update metadata to include FH2501 2019 data to some extend
update_metadata <- function(metadata) {
    outmetadata = metadata[is.na(metadata$CODES),] %>%
        dplyr::mutate(CODES = substr(samples_fix, 3,4),
                      SITE = substr(samples_fix, 5,6),
                      PLOT = as.integer(substr(samples_fix, 7,8)),
                      DATE = substr(samples_fix, 9, 16),
                      SAMPLE_ID = paste(CODES, SITE, PLOT, DATE, sep = '-'),
                      D = as.Date(DATE, format = '%Y%m%d'),
                      year = substr(samples_fix, 9, 12),
                      startyear = as.Date(paste0(year, '0101'), format = '%Y%m%d'),
                      doy = as.integer(difftime(D,startyear)))
    # return the metadata
    metadata[is.na(metadata$CODES),] = outmetadata
    return(metadata)
}

date_fixer <- function(DATE) {
    outDATE = sapply(DATE, function(xx) {
        if (is.na(xx)) {yy = NA_character_}
        else {if (!is.na(as.Date(xx, format = '%m-%d-%y', optional = TRUE))) {
            yy = format(as.Date(xx, format = '%m-%d-%y', optional = TRUE), '%Y%m%d')
        } else {if (!is.na(as.integer(xx))) {
            yy = format(as.Date(as.character(as.integer(xx)), format = '%Y%m%d', optional = TRUE),'%Y%m%d')
        }
            else {stop(xx); yy = NULL}}}
    })
    outDATE = unname(outDATE)
    return(outDATE)
}

total_plant_number_optional_fixer <- function(var) {
    var[var == '>100'] = '100.0'
    var[var == 'more than 20'] = '20.0'
    var[stringr::str_detect(var, pattern = "total number of plants")] = NA
    var[var == 'many but probably not A. thaliana'] = NA
    # check ending
    end = unname(sapply(var[!is.na(var)], function(xx) {substr(xx, start = nchar(xx)-1,stop = nchar(xx))}))
    if (any(end != '.0')) {
        stop('Wrong TOTAL_PLANT_NUMBER_OPTIONAL format')
    }
    var = as.integer(var)
    return(var)
}

output_data <- function(object,outdir, filename) {
    # write csv file, quote everything
    write.csv(x = object, file = paste0(outdir, filename, '.csv'))
    # write RDS file
    saveRDS(object = object, file = paste0(outdir, filename, '.rds'))
    return(NULL)
}

# def variables --------
outdir = './data/metadata/ver202209/'

# load data --------
# load ver202209 284 samples
# (note the vars are not well separated due to elevation formats)
samples = read_samples(filename = './data-raw/metadata/ver202209/preliminary.fam')
# double check the samples are not duplicated
any(duplicated(samples))
# match the samples with their metadata (two tier)
# samples -- recordssorted$id; recordssorted$SITE -- sites.clim$SITE_CODE
# load(file = 'data-raw/metadata/ver202209/grenephase1-data/seqtable.rda')
load(file = 'data-raw/metadata/ver202209/grenephase1-data/recordssorted.rda')
load(file = 'data-raw/metadata/ver202209/grenephase1-data/sites.clim.rda')
load(file = 'data-raw/metadata/ver202209/grenephase1-data/census.rda')

# get a list of samples in the haplotype-AF samples
af_samples = list.files(path ='data-raw/AF/ver202209/haplotype/samples_frequency/')
af_samples = unname(sapply(af_samples, function(xx) {strsplit(x = xx, split = '-')[[1]][1]}))

# make sure the samples are the same (YES!)
all.equal(sort(samples), sort(af_samples))

# main --------
# samples -- recordssorted$id ========
# match samples with metadata recordssorted
metadata0 = dplyr::left_join(data.frame(samples), y = recordssorted, by = c('samples' = 'id'))
# some of the recordssorted did not include the samples in plink files
mismatch = metadata0[is.na(metadata0$CODES),'samples']
mismatch

# attempt to fix the mismatches
fixed_mismatch <- dplyr::bind_rows(lapply(mismatch, match_records))
fixed_mismatch

fixer_mismatch <- fixed_mismatch[fixed_mismatch$resolved == TRUE, 'id']
names(fixer_mismatch) <- fixed_mismatch[fixed_mismatch$resolved == TRUE, 'id_s.newid']
fixer_mismatch
samples_fix = samples
samples_fix[sapply(fixer_mismatch, function(xx) {unname(which(samples == xx))})] = names(fixer_mismatch)

# get metadata again
metadata1 = dplyr::left_join(data.frame(samples, samples_fix), y = recordssorted, by = c('samples_fix' = 'id'))

# fill in the FH2501/FH2502 for 2019 information
metadata = update_metadata(metadata = metadata1)
# rename some columns
colnames(metadata)[colnames(metadata) %in% c('samples', 'samples_fix', 'doy')] = c('id', 'id_fix', 'day')
meta0922 = metadata %>%
    dplyr::mutate(NUMBER_FLOWERS_COLLECTED = as.numeric(NUMBER_FLOWERS_COLLECTED),
                  SITE = as.integer(SITE),
                  PLOT = as.integer(PLOT)) %>%
    dplyr::arrange(id_fix)
str(meta0922)

table(meta0922[,c('SITE','PLOT','year')])
table(meta0922$PLOT, useNA = 'always')

# recordssorted$SITE -- sites.clim$SITE_CODE ========
# 30 out of 45 sites were included, so only subset the sites that were included with climate data
clim0922 = sites.clim[sites.clim$SITE_CODE %in% metadata$SITE, ] %>%
    dplyr::select(SITE_CODE, LATITUDE, LONGITUDE,
                  starts_with('bio'), starts_with('prec'), starts_with('tmin'), starts_with('tmax')) %>%
    dplyr::rename(SITE = SITE_CODE)

str(clim0922)

# census data ========
# fix census data in general first (so we can join by the RECORD_ID/SAMPLE_ID)
census_edit = census[,1:11] %>%
    dplyr::filter(!is.na(SITE), !is.na(PLOT))
any(is.na(census_edit$RECORD_ID))
# rename census_edit
colnames(census_edit)
colnames(census_edit) = c('column_label','SITE','PLOT','DATE','RECORD_ID','DIAGONAL_PLANT_NUMBER',
                          'OFF_DIAGONAL_PLANT_NUMBER','TOTAL_PLANT_NUMBER_OPTIONAL',
                          'MEAN_FRUITS_PER_PLANT_OPTIONAL','SD_FRUITS_PER_PLANT_OPTIONAL', 'COMMENTS')
# basic info
table(census_edit$SITE,useNA = 'always')
table(census_edit$PLOT,useNA = 'always')
table(census_edit$DATE,useNA = 'always')
table(census_edit$DIAGONAL_PLANT_NUMBER,useNA = 'always')
table(census_edit$TOTAL_PLANT_NUMBER_OPTIONAL,useNA = 'always')
table(census_edit$MEAN_FRUITS_PER_PLANT_OPTIONAL,useNA = 'always')
table(census_edit$SD_FRUITS_PER_PLANT_OPTIONAL,useNA = 'always')
# table(census_edit$COMMENTS,useNA = 'always')

# fix some formatting
census_edit = census_edit %>%
    dplyr::mutate(across(.cols = DIAGONAL_PLANT_NUMBER:SD_FRUITS_PER_PLANT_OPTIONAL, .fns = ~ na_if(.x, '?'))) %>%
    dplyr::mutate(across(.cols = DIAGONAL_PLANT_NUMBER:SD_FRUITS_PER_PLANT_OPTIONAL, .fns = ~ na_if(.x, 'na'))) %>%
    dplyr::mutate(across(.cols = DIAGONAL_PLANT_NUMBER:SD_FRUITS_PER_PLANT_OPTIONAL, .fns = ~ na_if(.x, 'NA'))) %>%
    dplyr::mutate(DIAGONAL_PLANT_NUMBER = as.integer(DIAGONAL_PLANT_NUMBER),
                  OFF_DIAGONAL_PLANT_NUMBER = as.integer(OFF_DIAGONAL_PLANT_NUMBER),
                  across(.cols = ends_with('FRUITS_PER_PLANT_OPTIONAL'), .fns = as.numeric))

# fix the date (loop through to avoid NA issues)
census_edit$DATE = date_fixer(census_edit$DATE)

# fix the TOTAL_PLANT_NUMBER_OPTIONAL
census_edit$TOTAL_PLANT_NUMBER_OPTIONAL = total_plant_number_optional_fixer(census_edit$TOTAL_PLANT_NUMBER_OPTIONAL)

census_edit = census_edit %>%
    dplyr::mutate(SAMPLE_ID = paste('FH', SITE, PLOT, DATE, sep = '-')) %>%
    dplyr::arrange(SITE, PLOT, DATE)

# print final data structure
str(census_edit)

# now generate the census0922 data by left_join with metadata
census0922 = dplyr::left_join(x = meta0922, y = census_edit %>% dplyr::select(-SITE, -PLOT, -DATE), by = 'SAMPLE_ID')
table(is.na(census0922$RECORD_ID)) # Only 35 contained the census data from the same day
# SCRATCH
# View(census0922[,c('DIAGONAL_PLANT_NUMBER1','DIAGONAL_PLANT_NUMBER')])
# View(census0922[,c('DATE','DATE1')])

# generate a large metadata file with metadata + census + climate
allmeta0922 = dplyr::left_join(x = census0922, y = clim0922, by = 'SITE')
str(allmeta0922)

# output files --------
output_data(object = meta0922, outdir = outdir, filename = 'metadata_0922')
output_data(object = clim0922, outdir = outdir, filename = 'climate_0922')
output_data(object = census_edit, outdir = outdir, filename = 'census_edit')
output_data(object = census0922, outdir = outdir, filename = 'census_0922')
output_data(object = allmeta0922, outdir = outdir, filename = 'metacensusclim_0922')

save.image(file = './rdata/ver202209/step0.1_clean_metadata.RData')

# cleanup --------
date()
closeAllConnections()
