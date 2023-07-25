# Title: Plot LMM for bio1
# Author: Meixi Lin
# Date: Wed Mar  8 09:32:30 2023
# Modification: Update for new LMM setup
# Date: Tue Jul 25 12:45:46 2023
# P-value vector derived from: lmm_randomslope.R script output pvector1
# generation == 3
# lmerTest::lmer(deltap ~ bio1 + (1|site), data = mydata)

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses")

library(ggplot2)
library(dplyr)
library(lmerTest)
library(qqman)

date()
sessionInfo()

# def functions --------
source('scripts/manuscript/config_manuscript.R')

format_plotdt <- function(pv) {
    # chromosome positions
    load('data-raw/snp_freq/seedmix_p0_231_ldpruned.rda')
    # adjust p-values
    pvc = p.adjust(pv, method = 'bonferroni')
    plotdt = cbind(p0_p[,c('CHR', 'POS')], pv, pvc) 
    plotdt$id = 1:nrow(plotdt)
    plotdt = dplyr::left_join(plotdt, chrdf[,c('CHR', 'sumstart')], by = 'CHR') %>% 
        dplyr::mutate(psig = pvc < 0.05,
                      sumpos = POS + sumstart - 1)
    data.table::setDF(plotdt)
    return(plotdt)
}

select_snps <- function(deltapsy, metadty, plotdt, id) {
    mydata = cbind(unlist(deltapsy[id,]), metadty[,c('mergeid', 'site', 'bio1')])
    colnames(mydata)[1] = 'deltap'
    mydata$id = id
    mydata$psig = plotdt[plotdt$id == id, 'psig']
    return(mydata)
}

extract_regline <- function(mymodel) {
    ss = summary(mymodel)
    output = c(aa = ss$coefficients['(Intercept)', 'Estimate'],
               bb = ss$coefficients['bio1', 'Estimate'])
    return(output)
}

plot_example <- function(mydf, mymodel, mycol, mytext) {
    regline = extract_regline(mymodel)
    pp <- ggplot(data = mydf,
           mapping = aes(x = bio1, y = deltap)) +
        geom_abline(intercept = regline['aa'], slope = regline['bb'], 
                    color = mycol, linetype = 'dashed') + 
        annotate('text', label = mytext, size = 4,
                 x = 8.2, y = 30, vjust = 1, hjust = 0) +
        geom_point(color = mycol) +
        labs(x = 'Annual temperature (ºC)', y = expression(Delta*p/p[0]*(1-p[0]))) +
        ylim(-2.3,30)
    return(pp)
}

# def variables --------
pcutoff = 5e-8
plotdir = './data/manuscript/plot_lmm_bio1/'
dir.create(plotdir)

mygen = 3 # third generation

# load data --------
lmedata <- loadSomeRData(c('pvector1', 'deltapsy', 'metadty'), './data/lmm_loo/lmm_randomslope.RData')
pvector <- lmedata[[1]]; deltapsy <- lmedata[[2]]; metadty <- lmedata[[3]]
# add the data for chromosome lengths 
chrdf = read.csv('data-raw/TAIR10/TAIR10_chr_len.csv', row.names = 1)[1:5,]
plotdt <- format_plotdt(pvector)

# main --------
# manhattan plot
pdf(file = paste0(plotdir, 'lmm_gen3_bio1_site.pdf'), width = 8, height = 4)
qqman::manhattan(plotdt, chr = 'CHR', bp = 'POS', p = 'pv', snp = 'id', highlight = plotdt[plotdt$psig, 'id'], 
                 genomewideline = FALSE, suggestiveline = FALSE)
dev.off()

# qq plot
pdf(file = paste0(plotdir, 'lmm_gen3_bio1_qq.pdf'), width = 5, height = 5)
qqman::qq(plotdt$pv)
dev.off()

# add example clines 
sig_snp = select_snps(deltapsy, metadty, plotdt, id = plotdt[plotdt$pv == min(plotdt$pv),'id']) 
nosig_snp = select_snps(deltapsy, metadty, plotdt, id = plotdt[plotdt$pv == max(plotdt$pv),'id'])

# run the model again to get the regression lines
sig_model = lmerTest::lmer(deltap ~ bio1 + (1|site), data = sig_snp) # this model had singular fit though
nosig_model = lmerTest::lmer(deltap ~ bio1 + (1|site), data = nosig_snp)

sig_pp <- plot_example(sig_snp,sig_model,'#238b45','y = -4.24 + 0.33x')
nosig_pp <- plot_example(nosig_snp,nosig_model,'darkgray','y = 0.24 - 0.00x')

ggsave(filename = 'lmm_gen3_bio1_sig_chr2-872619.pdf', path = plotdir, height = 5, width = 5, plot = sig_pp)
ggsave(filename = 'lmm_gen3_bio1_nosig_chr5-14564240.pdf', path = plotdir, height = 5, width = 5, plot = nosig_pp)

# output files --------
save.image(file = paste0(plotdir, 'lmm_plotdata.RData'))

# cleanup --------
date()
closeAllConnections()



