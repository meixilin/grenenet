# R version 4.2.1
# library(lme4)
# library(devtools)
# devtools::install_git('https://github.com/perpdgo/lme4GS/',subdir='pkg_src/lme4GS')
# install.packages('BGLR')

library(lme4GS)

library(BGLR)
library(lme4GS)

#Example 1, wheat 
data(wheat)
X<-wheat.X

Z<-scale(X,center=TRUE,scale=TRUE)
# tcrossprod(x) is formally equivalent to, but faster than, the call x %*% t(x),
G<-tcrossprod(Z)/ncol(Z)
A<-wheat.A
rownames(G)<-colnames(G)<-rownames(A)
y<-as.vector(wheat.Y)
X[1:5,1:5]
G[1:5,1:5] # pair-wise sample matrix
A[1:5,1:5] # pair-wise sample matrix

data=data.frame(y=y,m_id=rep(rownames(G),4),a_id=rep(rownames(A),4))
all(data$m_id == data$a_id)

out<-lme4GS::lmerUvcov(y ~ (1|m_id) + (1|a_id), data=data,
               Uvcov=list(m_id=list(K=G),a_id=list(K=A)))

summary(out)

plot(y,predict(out))