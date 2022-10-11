setwd('./scripts/metadata_check/')

require(dplyr)
# 376 files in release 02/05 was duplicated
dup2 = readLines(con = 'duplicated_reads_release02.csv')
dup5 = readLines(con = 'duplicated_reads_release05.csv')
any(duplicated(dup2))
any(duplicated(dup5))

dupr2 = dup2 %>%
    reshape2::colsplit(., pattern =  '/', names = paste0('V',1:14))
dupr5 = dup5 %>%
    reshape2::colsplit(., pattern =  '/', names = paste0('V',1:14))
any(duplicated(dupr2$V14))
any(duplicated(dupr5$V14))

md5r2 = read.table(file = 'duplicated_reads_release02.md5sum.txt')
md5r5 = read.table(file = 'duplicated_reads_release05.md5sum.txt')
md5r2$filename = reshape2::colsplit(md5r2$V2, pattern =  '/', names = paste0('V',1:14))[,14]
md5r5$filename = reshape2::colsplit(md5r5$V2, pattern =  '/', names = paste0('V',1:14))[,14]

all(md5r2$filename==md5r5$filename)
# TRUE
all(md5r2$V1==md5r5$V1)
# TRUE
