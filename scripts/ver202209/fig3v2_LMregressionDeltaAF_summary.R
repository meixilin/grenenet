# Title: Summarize the LM regression results
# Author: Meixi Lin
# Date: Thu Sep  1 13:30:07 2022

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet")
date()
sessionInfo()

# def functions --------

# def variables --------

# load data --------
envvars = c("LATITUDE","LONGITUDE","bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19","prec1","prec2","prec3","prec4","prec5","prec6","prec7","prec8","prec9","prec10","prec11","prec12","tmin1","tmin2","tmin3","tmin4","tmin5","tmin6","tmin7","tmin8","tmin9","tmin10","tmin11","tmin12","tmax1","tmax2","tmax3","tmax4","tmax5","tmax6","tmax7","tmax8","tmax9","tmax10","tmax11","tmax12")

lmreslist <- lapply(envvars, function(env) {
    dt = read.csv(file = paste0('./plots/ver202209/AF/DeltaP_v2/', env, '/scaled_deltaP_', env, '_lmres_Ppass.csv'), row.names = 1) %>%
        dplyr::summarise(nPpass = n(),
                         maxR2 = max(AdjR2, na.rm = TRUE),
                         meanR2 = mean(AdjR2, na.rm = TRUE)) %>%
        dplyr::mutate(env = env)
})

# main --------
lmresdt = dplyr::bind_rows(lmreslist) %>%
    dplyr::arrange(desc(nPpass))

head(lmresdt)

# output files --------
write.csv(lmresdt, file = './plots/ver202209/AF/DeltaP_v2/summary_lm.csv')

# cleanup --------
date()
closeAllConnections()
