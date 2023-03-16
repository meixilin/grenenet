# Title: Plot LMM for bio1
# Author: Meixi Lin
# Date: Wed Mar  8 09:32:30 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses")

library(ggplot2)
library(dplyr)

date()
sessionInfo()

# def functions --------
source('scripts/manuscript/config_manuscript.R')

plot_manhattan <- function(dt, outfile, plotdir) {
    pp1 <- ggplot(data = dt, aes(x = POS, y = -log10(LR_nofix_p), color = as.factor(pass_p2), size = as.factor(pass_p2))) +
        geom_point() +
        scale_size_manual(values = c(1, 0.2, 0.1)) + 
        scale_color_manual(values = c('green','darkgreen','darkgray')) +
        facet_wrap(. ~ CHR, ncol = 1, scales = 'free_x') + 
        geom_hline(yintercept = -log10(pcutoff), color = 'darkgray', linetype = 'dashed') +
        labs(x = 'genomic position (bp)', y = '-log10(P)') +
        theme(legend.position = 'none')
    ggsave(filename =  outfile, plot = pp1, path = plotdir, height = 8, width = 8)
    return(pp1)
}

plot_qqman <- function(dt, outprefix, plotdir) {
    png(filename = paste0(plotdir, outprefix, '_LR_nofix_p.png'), width = 400, height = 400)
    qqman::qq(dt[['LR_nofix_p']], xlim = c(0,9), ylim = c(0,9))
    dev.off()
    png(filename = paste0(plotdir, outprefix, '_beta_p.png'), width = 400, height = 400)
    qqman::qq(dt[['beta_p']], xlim = c(0,9), ylim = c(0,9))
    dev.off()
    return(invisible())
}

# def variables --------
pcutoff = 5e-8
plotdir = 'data/manuscript/plot_lmm_bio1/'
dir.create(plotdir)

years = c('2018', '2019', '2020')

# load data --------
lmeresl <- lapply(years, function(myyear) {
    load(paste0('data/lmm_loo/bio1/lmeout_bio1_all_', myyear, '.rda'))
    lmeresallp$pass_p2 = ifelse(lmeresallp$LR_nofix_p < pcutoff, 1, ifelse(lmeresallp$pass_p == TRUE, 2, 3))
    print(table(lmeresallp$pass_p2))
    return(lmeresallp)
})

# 1       2       3 
# 3    6253 2939798 
# 
# 1       2       3 
# 18   10947 2934880 
# 
# 1       2       3 
# 5    3248 2942714 

# main --------
for (ii in 1:length(years)) {
    myyear = years[ii]
    mydt = lmeresl[[ii]]
    pp = plot_manhattan(mydt, outfile = paste0('lmeout_bio1_all_', myyear, '_manhattan.png'), plotdir = plotdir)
    qq = plot_qqman(mydt, outprefix = paste0('lmeout_bio1_all_', myyear, '_QQplot'), plotdir = plotdir)
}

# output files --------


# cleanup --------
date()
closeAllConnections()
