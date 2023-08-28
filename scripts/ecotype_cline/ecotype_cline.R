# Title: Ecotype clustering with deviation from the source environments
# Based on local-adaptation-stabilizing-selection.Rmd
# Author: Meixi Lin
# Date: Fri Jun  2 09:49:43 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses")
date()
sessionInfo()

library(ggplot2)
library(dplyr)
library(reshape2)

# def functions --------
calc_env_dist <- function(ecotypedf, sitedf) {
    sitedist = apply(sitedf[,-1], 1, function(row) (row - ecotypedf)^2) %>% t()
    sitedist = data.frame(sitedist) 
    sitedist$site = sitedf[,1] 
    return(sitedist)
}

format_lmtables <- function(ecop0, deltaenv, myecotype, myenv, mygen){
    myp0 = ecop0[ecop0$ecotypeid == myecotype, 'p0']
    tmpdeltaenv = deltaenv[deltaenv$ecotypeid == myecotype, c('ecotypeid', 'site', myenv)]
    tmpeco = ecomergedp_long %>% 
        dplyr::filter(ecotypeid == myecotype, generation == mygen) %>% 
        dplyr::mutate(log_p1_p0 = log(frequency/myp0)) %>% 
        dplyr::mutate(log_p1_p0 = ifelse(is.infinite(log_p1_p0), NA, log_p1_p0)) %>% 
        dplyr::left_join(.,  tmpdeltaenv, by = c('site', 'ecotypeid'))
    return(tmpeco)
}

format_lmsumm <- function(tmpeco) {
    outdf = data.frame(matrix(ncol = 7))
    colnames(outdf) = c('ecotypeid', 'envvar', 'generation', 'adj_r2', 'b_env', 'b_p', 'df')
    outdf$ecotypeid <- unique(tmpeco$ecotypeid)
    outdf$envvar <- colnames(tmpeco)[ncol(tmpeco)]
    outdf$generation <- unique(tmpeco$generation)
    return(outdf)
}

format_lmsum <- function(tmpeco, lmsum) {
    outdf = data.frame(matrix(ncol = 7))
    colnames(outdf) = c('ecotypeid', 'envvar', 'generation', 'adj_r2', 'b_env', 'b_p', 'df')
    outdf$adj_r2 <- lmsum$adj.r.squared
    outdf$b_env <- lmsum$coefficients[myenv, 'Estimate']
    outdf$b_p <- lmsum$coefficients[myenv, 'Pr(>|t|)']
    outdf$df <- lmsum$df[2]
    outdf$ecotypeid <- unique(tmpeco$ecotypeid)
    outdf$envvar <- colnames(tmpeco)[ncol(tmpeco)]
    outdf$generation <- unique(tmpeco$generation)
    return(outdf)
}

run_lm <- function(tmpeco, myenv) {
    if (all(is.na(tmpeco[, myenv]))) {
        lmres = format_lmsumm(tmpeco)
    } else {
        mymodel = lm(formula = as.formula(paste0('log_p1_p0 ~ ', myenv)), data = tmpeco)
        lmres = format_lmsum(tmpeco, summary(mymodel))
    }
    return(lmres)
}

plot_vs <- function(alldf, myenv) {
    mydf = alldf[alldf$envvar == myenv,]
    pp <- ggplot(mydf, aes(x = as.factor(ecotypeid), y = b_env, 
                           color = pass_both, shape = as.factor(generation))) + 
        geom_point() + 
        scale_color_manual(values = c('lightgray', 'darkgreen')) + 
        labs(x = 'ecotypeid', y = '-1/Vs', shape = 'generation', title = myenv, color = 'R2 > 0.15 & P < 0.05') +  
        theme(axis.text.x = element_text(angle = 90, size = 4),
              legend.position = 'top')
    ggsave(filename = paste0('vs_', myenv, '.pdf'), plot = pp, path = plotdir, height = 4, width = 8)
    return(invisible())
}

plot_map_vs <- function(plotdf, myenv, mygen) {
    # put this on a map
    mydf = plotdf[plotdf$generation == mygen & plotdf$envvar == myenv, ]
    data("wrld_simpl", package = "maptools")
    pp <- ggplot() +
        geom_polygon(data = wrld_simpl, aes(x = long, y = lat, group = group), fill = "lightgray", colour = "gray") +
        geom_point(data = mydf, 
                   mapping = aes(x = longitude, y = latitude, color = (-1/b_env), shape = pass_both),
                   size = 2) + 
        # Convert to polar coordinates
        # coord_map("ortho", orientation = c(90, 0, 0)) +
        # scale_y_continuous(breaks = seq(0, 90, by = 30), labels = NULL) +
        scale_colour_gradient(low = 'white', high = 'darkgreen', trans = 'pseudo_log') +
        # Removes Axes and labels
        # scale_x_continuous(breaks = NULL) +
        coord_cartesian(xlim = c(-24, 89), ylim = c(15,64)) + 
        labs(x = '', y = '',color = 'Vs', shape = 'P < 0.05', title = paste0(myenv, 'gen', mygen)) + 
        theme(panel.background = element_blank(),
              panel.grid.major = element_line(size = 0.25, linetype = 'dashed',
                                              colour = "black"),
              axis.ticks=element_blank(),
              axis.line = element_blank())
    ggsave(filename = paste0('mapvs_', myenv, 'gen', mygen, '.pdf'), plot = pp, path = plotdir, height = 4, width = 8)
    return(pp)
}

plot_bestfit <- function(alldf, myenv, mygen) {
    mydf = alldf[alldf$envvar == myenv & alldf$generation == mygen,]
    besteco = mydf %>% dplyr::slice_max(order_by = adj_r2)
    worsteco = mydf %>% dplyr::slice_min(order_by = adj_r2)
    tmpeco1 = format_lmtables(ecop0, deltaenv, besteco$ecotypeid, myenv, mygen)
    tmpeco2 = format_lmtables(ecop0, deltaenv, worsteco$ecotypeid, myenv, mygen)
    plotdf = rbind(tmpeco1, tmpeco2)
    pp <- ggplot(plotdf, aes(x = !!as.name(myenv), y = log_p1_p0, color = as.factor(ecotypeid))) + 
        geom_point() + 
        stat_smooth(method = "lm", formula = y ~ x, geom = "smooth") + 
        ggpubr::stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~"))) + 
        facet_wrap(. ~ ecotypeid, scales = 'free_x') +
        labs(x = paste0('(x - x0)^2 ', myenv), y = 'log(p1/p0)', color = 'ecotypeid', title = paste0) + 
        theme(legend.position = 'top')
}

# def variables --------
outdir = './data/ecotype_cline/'
plotdir = paste0(outdir, 'plots/')
dir.create(outdir)
dir.create(plotdir)

theme_set(theme_bw())

# load data --------
load('../metadata/data/worldclim_ecotypesdata.rda')
load('../metadata/data/worldclim_sitesdata.rda')
load('../metadata/data/ecotypes_data.rda')
load('../metadata/data/locations_data.rda')
all(colnames(worldclim_ecotypesdata)[-1] == colnames(worldclim_sitesdata)[-1])

load('./data-raw/ecotype_freq/merged_ecotype_freq_long.rds')
load('./data-raw/ecotype_freq/ecotype_p0.rds')

# add longitude and latitude
worldclim_ecotypesdata = dplyr::left_join(worldclim_ecotypesdata, 
                                          ecotypes_data[,c('ecotypeid', 'longitude', 'latitude')], by = 'ecotypeid')
worldclim_sitesdata = dplyr::left_join(worldclim_sitesdata, 
                                       locations_data[,c('site', 'longitude', 'latitude')], by = 'site')

envvars = colnames(worldclim_ecotypesdata)[-1]
generations = unique(ecomergedp_long$generation)

# main --------
# for each ecotype and each bioclim var, we calculate the euclidean distance to each site
deltaenv = lapply(1:nrow(worldclim_ecotypesdata), function(ii) {
    tmp = unlist(worldclim_ecotypesdata[ii,-1])
    out = calc_env_dist(ecotypedf = tmp, sitedf = worldclim_sitesdata) %>% 
        dplyr::mutate(ecotypeid = worldclim_ecotypesdata[ii,'ecotypeid']) %>% 
        dplyr::relocate(ecotypeid, site)
    return(out)
}) %>% dplyr::bind_rows()

# run linear model for each ecotype, each bioclim var 
alldf = data.frame()
# envvars = envvars[c(1:19, 104:105)]
envvars = 'bio1'
for (myenv in envvars) {
    for (mygen in generations) {
        lmdf = lapply(1:nrow(ecop0), function(ii) {
            # print(ii)
            myecotype = ecop0[ii, 'ecotypeid']
            tmpeco = format_lmtables(ecop0, deltaenv, myecotype, myenv, mygen)
            lmout = run_lm(tmpeco, myenv)
            return(lmout)
        }) %>% dplyr::bind_rows()
        alldf = rbind(alldf, lmdf)
    }
}

# add more corrections
alldf = alldf %>% 
    dplyr::mutate(pass_r2 = adj_r2 > 0.15,
                  pass_bp = b_p < 0.05/nrow(alldf),
                  pass_both = pass_r2 & pass_bp,
                  vs = -1/b_env,
                  vs_pos = vs > 0) 
head(alldf)

lmressum = alldf %>% 
    dplyr::group_by(envvar, generation) %>% 
    dplyr::summarise(n_passr2 = sum(pass_r2, na.rm = T),

                     max_r2 = max(adj_r2, na.rm = T)) 

alldf %>% 
    dplyr::group_by(envvar, generation) %>% 
    dplyr::summarise(n_vs_pos = sum(vs_pos, na.rm = T),
                     total_ecotypes = n(),
                     prp_vs_pos = n_vs_pos/n()) 

summary(alldf$vs)

plotdf = dplyr::left_join(alldf, ecotypes_data, by = 'ecotypeid')

# plotting --------
for (myenv in envvars) {
    plot_vs(alldf, myenv)
}

# put things on a map
for (myenv in envvars) {
    for (mygen in generations) {
        plot_map_vs(plotdf, myenv, mygen)
    }
}

# output files --------
save.image(file = paste0(outdir, 'ecotype_cline_bio1.RData'))

# cleanup --------
date()
closeAllConnections()
