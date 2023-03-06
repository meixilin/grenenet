# Title: Set up and test leave-one-out linear mixed model
# Author: Meixi Lin
# Date: Sun Mar  5 20:03:34 2023

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

setwd("/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses")

library(dplyr)
library(nlme)
library(Metrics)

date()
sessionInfo()

# def functions --------

# def variables --------
outdir = 'data/lmm_loo/test/'
dir.create(outdir, recursive = TRUE)

# load data --------
# # load and only select the first 30000 sites
# load('data/AF/merged_delta_p_776.rda')
# deltapt = deltap[1:30000, ]
# save(deltapt,file = paste0(outdir, 'merged_delta_p_776_30k.rda'))
# load('data/AF/seedmix_p0_231.rda')
# p0t = p0[1:30000, ]
# save(p0t,file = paste0(outdir, 'seedmix_p0_231_30k.rda'))
# rm(deltap, p0)

# load the data
load('data/lmm_loo/test/merged_delta_p_776_30k.rda')
load('data/lmm_loo/test/seedmix_p0_231_30k.rda')

# set up sample metadata
metadt = reshape2::colsplit(colnames(deltapt), '_', names = c('site', 'generation', 'plot')) %>% 
    dplyr::group_by(site, generation) %>%
    dplyr::mutate(nplot = n())

# load climate data
load('../metadata/data/weatherstation_bioclim.rda')

# which is the first year, and have more than one plot
firstyearsample = which(metadt$generation == 1 & metadt$nplot > 1)

# set up the first year relavant deltap
deltapt1 = deltapt[,..firstyearsample]

# set up first year metadata
metadt1 = metadt %>% 
    dplyr::filter(generation == 1, nplot > 1) %>% 
    dplyr::left_join(., weatherstation_bioclim[weatherstation_bioclim$year == 2018, c('site', 'year', 'bio1')], by = 'site')

length(unique(metadt1$site)) # 30 sites

# main --------
for (ii in 1:3000) {
    mydata = data.frame(deltap = unname(unlist(deltapt1[ii,])), 
                        id = colnames(deltapt1),
                        bio1 = metadt1$bio1,
                        site = metadt1$site)
    if (all(mydata$deltap == 0)) {
        next
    }
    # set up train and test
    mytraintest = lapply(base::split(mydata, mydata$site), function(xx) {
        testid = sample(1:nrow(xx), size = 1)
        testdt = xx[testid, ]
        traindt = xx[-testid, ]
        return(list(traindt, testdt))
    })
    mytrain = bind_rows(lapply(mytraintest, function(xx){xx[[1]]}))
    mytest = bind_rows(lapply(mytraintest, function(xx){xx[[2]]}))
    
    mymodel = nlme::lme(fixed = deltap ~ bio1, random = ~1|site, data = mytrain, method = 'ML')
    mymodel1 = nlme::lme(fixed = deltap ~ 1, random = ~1|site, data = mytrain, method = 'ML')
    mymodel0 = lm(deltap ~ bio1, data = mytrain)
    anova(mymodel, mymodel0)
    anova(mymodel, mymodel1)
    anova(mymodel0, mymodel)
    
    
    summary(mymodelx)
    summary(mymodel)
    mysum = summary(mymodel)
    myr2 = MuMIn::r.squaredGLMM(mymodel)
    if(mysum$tTable[2, 'p-value'] < 0.001 & myr2[1, 'R2c'] > 0.2) {
        break
    }
    
    # get the significance of the model using leave-one-out
    library(caret)
    
    lmm <- list(library = 'nlme',
                type = 'Regression',
                )
    
    #specify the cross-validation method
    ctrl <- trainControl(method = "LOOCV")
    
    #fit a regression model and use LOOCV to evaluate performance
    model <- train(deltap ~ bio1, data = mytrain, method = "lm", trControl = ctrl)
    model
    # make prediction
    mymodel2 = nlme::lme(fixed = deltap ~ bio1, random = ~ 1, data = mytrain)
    
    testpre = predict(mymodel,newdata = mytest[,c('site', 'bio1')])
    table(sign(testpre) == sign(mytest$deltap))
    View(cbind(testpre, mytest$deltap))
    plot(mytest$deltap, testpre)
    summary(lm(testpre ~ mytest$deltap))
    
    # calculate mse
    caret::postResample(pred = testpre, obs = mytest$deltap)
    sqrt(mean((testpre - mytest$deltap)^2)) # hand calculate RMSE
    summary(lm(testpre ~ mytest$deltap)) # hand calculate R2
    caret::RMSE(testpre, )
}



# output files --------

# cleanup --------
date()
closeAllConnections()
