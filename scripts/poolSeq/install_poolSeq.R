# R (>= 3.3.1)
# data.table (>= 1.9.4)
# foreach (>= 1.4.2)
# stringi (>= 0.4-1)
# matrixStats (>= 0.14.2)
# Rcpp (>= 1.0.3)

library(data.table)
library(foreach)
library(stringi)
library(matrixStats)
library(Rcpp)

sessionInfo()
# [1] Rcpp_1.0.7         matrixStats_0.56.0 stringi_1.4.6      foreach_1.5.0      
# data.table_1.14.2 

# cd /home/mlin/Software/poolSeq_R
# wget https://github.com/ThomasTaus/poolSeq/archive/refs/tags/v0.3.5.tar.gz
install.packages("/home/mlin/Software/poolSeq_R/poolSeq_v0.3.5.tar.gz", repos=NULL, type="source")

library(poolSeq)


example(wf.traj)
