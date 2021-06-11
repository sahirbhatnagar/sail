library(doParallel)
registerDoParallel(cores = 4)
library(glmnet)
library(ITRSelect)
devtools::load_all()
k=1
l1=l2=l3=l4=s1=s2=s3=s4=lr1=lr2=lr3=lr4=sr1=sr2=sr3=sr4=ar1=ar2=ar3=ar4=a1=a2=a3=a4=
  al1=al2=al3=al4=as1=as2=as3=as4=list()

while (k<=10) {
  ds=g_data(1000,10);X=ds[[1]];Y=ds[[2]];w=ds[[3]];A=ds[[4]]
  ## 1. both wrong
  coeff=abs(coef(lm(Y~cbind(X,A,A*X)))[-1])
  pfac=1/coeff
  pfac=pfac/sum(pfac)*21
  m=cv.glmnet(x=cbind(X,A,A*X),y=Y,nfolds = 3,parallel = T)
  lasso=lasso_refit=coef(m,s='lambda.min')
  index <- (1:21)[lasso[2:(22)]!=0]
  lasso_refit[index+1]=coef(lm(Y~cbind(X,A,A*X)[,index]))[-1]
  lasso=lasso[12:22];lasso_refit=lasso_refit[12:22]
  m=cv.glmnet(x=cbind(X,A,A*X),y=Y,nfolds = 4,parallel = T, penalty.factor=pfac)
  l1[[k]]=lasso;lr1[[k]]=lasso_refit
  al1[[k]]=coef(m,s='lambda.min')[12:22]

  ## pal
  m=PAL(Y~X|A,refit = F,penalty = 'LASSO',pi1.est=.5)
  a1[[k]]=m$beta1.est
  m=PAL(Y~X|A,refit = T,penalty = 'LASSO',pi1.est=.5)
  ar1[[k]]=m$beta1.est

  ##sail
  pfac=1/coeff
  pfac[12:21]=coeff[1:10]*coeff[11]/coeff[12:21]
  pfac=c(pfac[11],pfac[-11])
  pfac=pfac/sum(pfac)*21
  m=cv.sail(y=Y,e=A,x=X,nfolds = 3,parallel = T,basis=function(i) i)
  sail=coef(m,s='lambda.min');sail_refit=sail
  index <- (1:21)[sail[2:(22)]!=0]
  sail_refit[index+1]=coef(lm(Y~cbind(X,A,A*X)[,index]))[-1]
  sail=sail[12:22];sail_refit=sail_refit[12:22]
  m=cv.sail(y=Y,e=A,x=X,nfolds = 4,penalty.factor=pfac,
            parallel = T,basis=function(i) i)
  s1[[k]]=sail;sr1[[k]]=sail_refit; as1[[k]]=coef(m,s='lambda.min')[12:22]

  ## 2.w correct; free incorrect

  ## pal
  m=PAL(Y~X|A,refit = F,penalty = 'LASSO')
  a2[[k]]=m$beta1.est
  m=PAL(Y~X|A,refit = T,penalty = 'LASSO')
  ar2[[k]]=m$beta1.est

  ##sail

  coeff=abs(coef(lm(Y~cbind(X,A,A*X),weights=w))[-1])
  pfac=1/coeff
  pfac[12:21]=coeff[1:10]*coeff[11]/coeff[12:21]
  pfac=c(pfac[11],pfac[-11])
  pfac=pfac/sum(pfac)*21

  m=cv.sail(y=Y,e=A,x=X,nfolds = 4,weights = w,
            parallel = T,basis=function(i) i)
  sail=coef(m,s='lambda.min');sail_refit=sail
  index <- (1:21)[sail[2:(22)]!=0]
  sail_refit[index+1]=coef(lm(Y~cbind(X,A,A*X)[,index],weights = w))[-1]
  sail=sail[12:22];sail_refit=sail_refit[12:22]
  m=cv.sail(y=Y,e=A,x=X,nfolds = 4,penalty.factor=pfac,weights = w,
            parallel = T,basis=function(i) i)
  s2[[k]]=sail;sr2[[k]]=sail_refit;as2[[k]]=coef(m,s='lambda.min')[12:22]

  ## 3  w incorrect; free correct

  coeff=abs(coef(lm(Y~cbind(exp(X[,1])+log(abs(X[,1])),X,A,A*cbind(exp(X[,1])+log(abs(X[,1])),X))))[-1])
  pfac=1/coeff
  pfac=pfac/sum(pfac)*23

  m=cv.glmnet(x=cbind(exp(X[,1])+log(abs(X[,1])),X,A,A*cbind(exp(X[,1])+log(abs(X[,1])),X)),
              parallel = T,y=Y,nfolds = 4)

  lasso=lasso_refit=coef(m,s='lambda.min')
  index <- (1:23)[lasso[2:(24)]!=0]
  lasso_refit[index+1]=coef(lm(Y~cbind(exp(X[,1])+log(abs(X[,1])),X,A,A*cbind(exp(X[,1])+log(abs(X[,1])),X))[,index]))[-1]
  lasso=lasso[13:24];lasso_refit=lasso_refit[13:24]
  m=cv.glmnet(x=cbind(exp(X[,1])+log(abs(X[,1])),X,A,A*cbind(exp(X[,1])+log(abs(X[,1])),X)),penalty.factor=pfac,
              parallel = T,y=Y,nfolds = 4)
  l3[[k]]=lasso;lr3[[k]]=lasso_refit;al3[[k]]=coef(m,s='lambda.min')[13:24]

  ##pal
  m=PAL(Y~cbind(exp(X[,1])+log(abs(X[,1])),X)|A,refit = F,penalty = 'LASSO',pi1.est=.5)
  a3[[k]]=m$beta1.est
  m=PAL(Y~cbind(exp(X[,1])+log(abs(X[,1])),X)|A,refit = T,penalty = 'LASSO',pi1.est=.5)
  ar3[[k]]=m$beta1.est

  ##sail
  coeff=abs(coef(lm(Y~cbind(exp(X[,1])+log(abs(X[,1])),X,A,A*cbind(exp(X[,1])+log(abs(X[,1])),X))))[-1])
  pfac=1/coeff
  pfac[13:23]=coeff[1:11]*coeff[12]/coeff[13:23]
  pfac=c(pfac[12],pfac[-12])
  pfac=pfac/sum(pfac)*23

  m=cv.sail(y=Y,e=A,x=cbind(exp(X[,1])+log(abs(X[,1])),X),nfolds = 4,
            parallel = T,basis=function(i) i)
  sail=coef(m,s='lambda.min');sail_refit=sail
  index <- (1:23)[sail[2:(24)]!=0]
  sail_refit[index+1]=coef(lm(Y~cbind(exp(X[,1])+log(abs(X[,1])),X,A,A*cbind(exp(X[,1])+log(abs(X[,1])),X))[,index]))[-1]
  sail=sail[13:24];sail_refit=sail_refit[13:24]
  m=cv.sail(y=Y,e=A,x=cbind(exp(X[,1])+log(abs(X[,1])),X),nfolds = 4,penalty.factor=pfac,
            parallel = T,basis=function(i) i)
  s3[[k]]=sail;sr3[[k]]=sail_refit; as3[[k]]=coef(m,s='lambda.min')[13:24]

  ## 4. both correct

  ##sail
  coeff=abs(coef(lm(Y~cbind(exp(X[,1])+log(abs(X[,1])),X,A,A*cbind(exp(X[,1])+log(abs(X[,1])),X)),weights = w))[-1])
  pfac=1/coeff
  pfac[13:23]=coeff[1:11]*coeff[12]/coeff[13:23]
  pfac=c(pfac[12],pfac[-12])
  pfac=pfac/sum(pfac)*23

  m=cv.sail(y=Y,e=A,x=cbind(exp(X[,1])+log(abs(X[,1])),X),nfolds = 4,weights = w,
            parallel = T,basis=function(i) i)
  sail=coef(m,s='lambda.min');sail_refit=sail
  index <- (1:23)[sail[2:(24)]!=0]
  sail_refit[index+1]=coef(lm(Y~cbind(exp(X[,1])+log(abs(X[,1])),X,A,A*cbind(exp(X[,1])+log(abs(X[,1])),X))[,index],weights = w))[-1]
  sail=sail[13:24];sail_refit=sail_refit[13:24]
  m=cv.sail(y=Y,e=A,x=cbind(exp(X[,1])+log(abs(X[,1])),X),nfolds = 4,weights = w,penalty.factor=pfac,
            parallel = T,basis=function(i) i)
  s4[[k]]=sail;sr4[[k]]=sail_refit; as4[[k]]=coef(m,s='lambda.min')[13:24]

  ## pal
  m=PAL(Y~cbind(exp(X[,1])+log(abs(X[,1])),X)|A,refit = F,penalty = 'LASSO')
  a4[[k]]=m$beta1.est
  m=PAL(Y~cbind(exp(X[,1])+log(abs(X[,1])),X)|A,refit = T,penalty = 'LASSO')
  ar4[[k]]=m$beta1.est

  k=k+1
}

save.image("100.RData")






