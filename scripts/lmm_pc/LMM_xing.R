## LMM model to test delta p and climate association

library(data.table)
library(lme4)
library(afex)
library(nlme)

setwd("/NOBACKUP/scratch/xwu/grenet/LMM_snp_climate/")

meta <- read.delim("../merged_frequency/merged_sample_table.csv",sep=",")

snp_p0 <- read.delim("../merged_frequency/average_seedmix_p0.txt",header=F)

delta_snp_p <- fread("../merged_frequency/merged_hapFIRE_delta_p.txt",nThread = 10)

delta_snp_p_norm <- as.matrix(delta_snp_p) / snp_p0$V3 / (1 - snp_p0$V3)

fwrite(delta_snp_p_norm,"../merged_frequency/merged_hapFIRE_delta_p_normed.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T,nThread = 10)

delta_snp_p_norm = fread("../merged_frequency/merged_hapFIRE_delta_p_normed.txt", nThread = 10, skip = 10000, nrows = 100)

delta_ecotype_p <- read.delim("../merged_frequency/delta_ecotype_freq.txt")

PC <- read.delim("../hapFIRE_updatedVCF/greneNet_final_v1.1_LDpruned_10PCs.txt",header = F,sep = " ")

bioclim <- read.delim("../metadata/worldclim_sitesdata.csv",sep=",")



### parse out the generation 3 data #### 

gen3_index <- which(meta$generation==3 & meta$total_flower_counts > 5)

meta_gen3 <- meta[gen3_index,]

bio1 <- c()
bio19 <- c()

for (i in 1:nrow(meta_gen3)){
  bio1 <- c(bio1,bioclim$bio1[which(meta_gen3$site[i] == bioclim$site)])
  bio19 <- c(bio19,bioclim$bio19[which(meta_gen3$site[i] == bioclim$site)])
}

meta_gen3$bio1 <- bio1
meta_gen3$bio19 <- bio19

delta_snp_p_norm_gen3 <- delta_snp_p_norm[,..gen3_index]

non_inform_index <- apply(delta_snp_p_norm_gen3,1,function(x) sum(x<0.01) > 145)

sum(non_inform_index)



### TEST LMM with randomly sampled SNPs 

set.seed(1)
random_snps <- sample(1:3235480,100000)

delta_snp_p_norm_gen3_random <- delta_snp_p_norm_gen3[random_snps,]

delta_snp_p_norm_gen3_random = delta_snp_p_norm_gen3 
delta_snp_p_norm_gen3_random = as.matrix(delta_snp_p_norm_gen3_random)

'''
the mixed model 

delta_snp_p ~ (1|site) + PC %*% delta_ecotype_p + bio1 

'''
p_value <- c()

pop_strc <- t(as.matrix(delta_ecotype_p[,gen3_index])) %*% as.matrix(PC)

for(i in 1:100){
  model_0 <- lmer(delta_snp_p_norm_gen3_random[i,] ~ pop_strc + (1|as.character(meta_gen3$site)), method = 'ML')
  model_1 <-  lmer(delta_snp_p_norm_gen3_random[i,] ~ pop_strc + (1|as.character(meta_gen3$site)) + meta_gen3$bio1, method = 'ML')
  anova_results <- anova(model_0,model_1)
  p_value <- c(p_value,anova_results$`Pr(>Chisq)`[2])
}

p_value_lmm <- c()
### LmerTest
for(i in 1:1000){
  lmm <- lmer(delta_snp_p_norm_gen3_random[i,] ~ pop_strc[,1:5] + (1|as.character(meta_gen3$site)) + scale(meta_gen3$bio19))
  results <- summary(lmm)
  p_value_lmm <- c(p_value_lmm,results$coefficients[7,5])
}
qq(p_value_lmm)

p_value_lmm2 <- c()
### LmerTest
for(i in 1:1000){
  lmm <- lmer(delta_snp_p_norm_gen3_random[i,] ~ (1|as.character(meta_gen3$site)) + meta_gen3$bio1)
  results <- summary(lmm)
  p_value_lmm2 <- c(p_value_lmm2,results$coefficients[2,5])
}
qq(p_value_lmm2)


