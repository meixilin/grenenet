### deltaP trend and within site replicability 

setwd("/NOBACKUP/scratch/xwu/grenet/deltaP_trend/")

library(data.table)


meta <- read.delim("../haf-pipe/hafpipe_231/merged_sample_table.csv",sep=",",stringsAsFactors = F)
meta$flowerless5 <- meta$total_flower_counts <= 5
delta_p <- fread("../haf-pipe/hafpipe_231/merged_delta_p.csv",sep=",",nThread = 10)

p0 <- c()
for (ch in 1:5){
  print(ch)
  tmp <- fread(paste("../seed_mix/haf-pipe/231_vcf/chr",ch,"_p0.txt",sep=""))
  p0 <- c(p0,tmp$V2)
}


### read the snp positions and LD pruned snps

snps <- read.delim("../genomic_evolution/chromosome_locations.txt",header=F)
ld_pruned <- read.delim("../genomic_evolution/plink.prune.in",header=F)

ld_pruned_index <- which(snps %in% ld_pruned$V1)

### remove the samples with <=5 flowers 
delta_p_matrix <- as.matrix(delta_p)

ld_delta_p <- delta_p_matrix[ld_pruned_index,which(meta$flowerless5==F)]

meta_flower5 <- meta[meta$flowerless5==F,]
dim(ld_delta_p)

write.table(ld_delta_p,"ld_pruned_deltaP_flower5.txt",row.names = F,quote = F,append = F,sep="\t")

#### 


ld_delta_p <- fread("ld_pruned_deltaP_flower5.txt",sep="\t")
ld_delta_p <- as.data.frame(ld_delta_p)

scale <- p0$p0[ld_pruned_index] * (1-p0$p0[ld_pruned_index])
norm_delta_p <-ld_delta_p / scale


sum(meta_flower5$generation==1)

unique_sites <- unique(meta_flower5$site[meta_flower5$generation==1])
snp_correlation <- c()
site <- c()

for (i in 1:length(unique_sites)){
  index <- which(meta_flower5$site==unique_sites[i] & meta_flower5$generation==1)
  if (length(index) > 1){
    data <- norm_delta_p[,index]
    m <- cor(data)
    m <- m[upper.tri(m)]
    snp_correlation <- c(snp_correlation,m)
    site <- c(site,rep(unique_sites[i],length(m)))
  }
}
print(length(site)) #1408

unique_sites <- unique(meta_flower5$site[meta_flower5$generation==2])
for (i in 1:length(unique_sites)){
  index <- which(meta_flower5$site==unique_sites[i] & meta_flower5$generation==2)
  if (length(index) > 1){
    data <- norm_delta_p[,index]
    m <- cor(data)
    m <- m[upper.tri(m)]
    snp_correlation <- c(snp_correlation,m)
    site <- c(site,rep(unique_sites[i],length(m)))
  }
}
print(length(site)) #2252

unique_sites <- unique(meta_flower5$site[meta_flower5$generation==3])
for (i in 1:length(unique_sites)){
  index <- which(meta_flower5$site==unique_sites[i] & meta_flower5$generation==3)
  if (length(index) > 1){
    data <- norm_delta_p[,index]
    m <- cor(data)
    m <- m[upper.tri(m)]
    snp_correlation <- c(snp_correlation,m)
    site <- c(site,rep(unique_sites[i],length(m)))
  }
}
print(length(site)) #2953
gen <- c(rep(1,1408),rep(2,844),rep(3,680))

data <- as.data.frame(matrix(nrow=2932,ncol=3))
colnames(data) <- c('generation','site','cor')
data$generation <- as.factor(gen)
data$site <- as.factor(site)
data$cor <- snp_correlation
ggplot(data,aes(x=reorder(site,cor),y=cor,fill=generation))+
  geom_boxplot(outlier.color = "red",outlier.size = 0.5,position = position_dodge(width = 0.7,preserve = "single"))+
  theme_classic()
ggsave("site_snp_correlation.pdf",device = "pdf",units = "in",height = 6,width = 10)

