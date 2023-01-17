# asked by xing wu
# date: 2022/12/13

vcfheader = read.table(file = '~/1001g.txt', header = FALSE)

ecotypesvcf = unname(unlist(vcfheader[,10:ncol(vcfheader)]))
ecotypesvcf = as.integer(ecotypesvcf)
# 224 ecotypes

# load metadata

load(file = '../../ath_evo/grenephase1/data/ecotypes_data.rda')
str(ecotypes_data)
inter = intersect(ecotypes_data$ecotypeid, ecotypesvcf)

ecotypes_dt = read.csv(file = '../../ath_evo/grenephase1/data/ecotypes_data.csv')
inter2 = intersect(ecotypes_dt$ecotypeid, ecotypesvcf)

inter3 = intersect(ecotypes_dt$ecotypeid, ecotypesvcf)

setdiff(inter3, inter2)
setdiff(inter2, inter3)
