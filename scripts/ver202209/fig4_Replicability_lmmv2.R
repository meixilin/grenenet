# load the data
deltadt = readRDS(file = './data/AF/ver202209/haplotype/scaled_delta_freq_1012.rds')
yy = deltadt['17546696',]
metadt = metadt[,c('id','bio1', 'year','SITE_f', 'PLOT2')]

# run lmm
tempmeta = metadt %>%
    dplyr::filter(SITE_f != '1')
biometa = metadt %>%
    dplyr::filter(SITE_f == '1')
tempyy = yy[names(yy) %in% tempmeta$id]

# Account for experimental setup
mydata = cbind(tempyy, tempmeta)
lmemodel = lme4::lmer(formula = tempyy ~ bio1 + (1|year) + (1|SITE_f/PLOT2),
                      data = mydata, na.action = na.exclude)

yypredict = predict(lmemodel)

predyy = predict(lmemodel, newdata = biometa, re.form = ~ (1|year))
names(predyy) = biometa$id
all(names(predyy) == bioyy$id)
bioyy = data.frame(yy[!(names(yy) %in% tempmeta$id)]) %>%
    tibble::rownames_to_column(var = 'id')
colnames(bioyy) = 'bioyy'
valdt = cbind(bioyy, predyy, biometa)

ggplot(data = valdt, mapping = aes(x = year, y = bioyy)) +
    geom_point()


lmemodel = suppressMessages(
    lmerTest::lmer(formula = tempyy ~ eval(as.name(envvar)) + )
# cannot keep year as a random effect because of singular fit
# lmemodel = lmerTest::lmer(formula = yy ~ eval(as.name(envvar)) + year + (1|SITE_f/PLOT2),
#                           data = mydata, na.action = na.exclude)
lmesum = summary(lmemodel)
lmepredict =
    lmer2 = MuMIn::r.squaredGLMM(lmemodel)
# get relevant stats
outdt = c(outsite,
          lmer2,
          lmesum$coefficients["eval(as.name(envvar))","Estimate"],
          log10(lmesum$coefficients["eval(as.name(envvar))","Pr(>|t|)"]),
          lmesum$AICtab)
# keep the
return(outdt)
