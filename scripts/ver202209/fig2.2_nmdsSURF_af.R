# Title: Ordisurf and envfit on the NMDS plane
# Author: Meixi Lin
# Date: Fri Aug  5 16:31:56 2022

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet")

library(dplyr)
library(vegan)
date()
sessionInfo()

# def functions --------
pnger <- function(filename,height,width) {
    png(filename = paste0(plotdir, filename), width = width, height = height, units = 'in', res = 150)
}

# def variables --------
plotdir = './plots/ver202209/AF/NMDS_surf_maf05euc/'
dir.create(plotdir)

# variables to test against
vars = paste0('bio',1:19)

source('./scripts/ver202209/config.R')

# load data --------
mdslist = readRDS(file = 'data/AF/ver202209/haplotype/NMDS/mdslist.rds')

# use the MAF = 0.05, distance measure = euclidean
mymds = mdslist$maf05_euc
allmeta0922 = readRDS(file = 'data/metadata/ver202209/metacensusclim_0922.rds')

# sanity check
all.equal(rownames(mymds$points)[-285],allmeta0922$id)

# main --------
# run envfit
envdt = rbind(allmeta0922[,vars], rep(NA, times = length(vars)))
envfitres <- vegan::envfit(ord = mymds, env = envdt, na.rm = TRUE, permutations = 2499)
envfitres
mycolors = site_cols[c(as.character(allmeta0922$SITE), 'seed_mix')]
pnger(filename = 'envfit_NMDS.png', height = 6, width = 6)
plot(x = mymds$points[,1], y = mymds$points[,2], col = mycolors, xlab = 'MDS1', ylab = 'MDS2')
plot(envfitres)
dev.off()

# run ordisurf
ordisurflist <- lapply(vars, function(xx) {
    out <- vegan::ordisurf(x = mymds, y = envdt[,xx], plot = FALSE)
})
names(ordisurflist) = vars

# plot every surf result
ordisurfplot <- mapply(function(surfob, myvar) {
    mycolors = site_cols[c(as.character(allmeta0922$SITE), 'seed_mix')]
    print(myvar)
    surfsum = summary(surfob)
    print(surfsum)
    r2 = surfsum$r.sq
    reml = unname(surfsum$sp.criterion)
    mytitle = paste0(myvar, ' R2=', round(r2, digits = 3),
                     ' REML=',round(reml, digits = 3))
    pnger(filename = paste0('ordisurf_',myvar,'.png'), height = 6, width = 6)
    plot(surfob$model$x1, surfob$model$x2, xlab = 'MDS1', ylab = 'MDS2', main = mytitle, col = mycolors)
    plot(surfob, main = myvar, add = TRUE, col = 'gray20')
    dev.off()
    output = c(r2,reml)
    return(output)
}, surfob = ordisurflist, myvar = names(ordisurflist))
ordisurfsummary = data.frame(t(ordisurfplot))
colnames(ordisurfsummary) = c('R2', 'REML')
ordisurfsummary %>%
    dplyr::arrange(desc(R2))

# output files --------
saveRDS(envfitres, file = './data/AF/ver202209/haplotype/NMDS/envfitres.rds')
write.csv(ordisurfplot, file = './data/AF/ver202209/haplotype/NMDS/ordisurf_res.csv')
save.image(file = './rdata/ver202209/fig2.2_nmdsSURF_af.RData')

# cleanup --------
date()
closeAllConnections()

