library(doParallel)
registerDoParallel(cores = 4)
library(glmnet)
library(ITRSelect)
devtools::load_all()

ds=read.csv('~/Desktop/stard.csv')
ds=ds[!(is.na(ds$S2)),]
ds=ds[!(is.na(ds$S1)),]
ds=ds[!(is.na(ds$Y1)),]
ds=ds[!ds$S2==Inf,]
ds=ds[!ds$S2==-Inf,]

ind=!(is.na(ds$Y2))
ds2=ds[ind,]

ds$Y1[ind]=NA

## noise variables
p=5;n=1027
set.seed(3333)
noise1 <- matrix(rbinom(n = n*p,1,prob = runif(1,0.3,0.8)),nrow = 1027)

index=(1: n*p)[rbinom(n*p,1,0.2)]
noise2=noise1
noise2[index]=abs(noise2[index]-1)

X2=cbind(ds2$A1,ds2$Q2,ds2$S2,ds2$P2,noise2[ind,])
A2=ds2$A2;Y2=ds2$Y2

m=glm(A2~ds2$P2,binomial)
w2=abs(A2-fitted(m))

X1=cbind(ds$Q1,ds$S1,ds$P1,noise1)
A1=ds$A1;Y1=ds$Y1

m=glm(A1~ds$P1,binomial)
w1=abs(A1-fitted(m))
