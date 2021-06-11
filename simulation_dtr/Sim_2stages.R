library(doParallel)
library(ITRSelect)
registerDoParallel(cores = 4)
library(glmnet)
devtools::load_all()
k=1
l2=l1=s2=s1=a2=a1=lr2=lr1=sr2=sr1=ar2=ar1=al2=al1=as2=as1=list()

expit=function(x){exp(x)/(1+exp(x))}
q1=0.5*(expit(2)+expit(0));q2=0.5*(expit(-2)+expit(0))
q3=0.5*(expit(2)-expit(0));q4=0.5*(expit(0)-expit(-2))
f1=-.5;f2=2.5
-1+(q1-q2)*(2-f2)

while (k<=100) {
  ds=g2(50000,10);Y=ds[[1]];X1=ds[[2]];X2=ds[[3]];A1=ds[[4]];A2=ds[[5]]

  ## lasso and refit
  m=cv.glmnet(x=cbind(X2,A2,A2*X2),y=Y,nfolds = 4,parallel = T)
  lasso2=coef(m,s='lambda.min');lasso2_refit=lasso2
  index <- (1:21)[lasso2[2:(22)]!=0]
  lasso2_refit[index+1]=coef(lm(Y~cbind(X2,A2,A2*X2)[,index]))[-1]
  lasso2=lasso2[12:22];lasso2_refit=lasso2_refit[12:22]

  aopt=as.numeric(cbind(1,X2)%*%lasso2>0)
  yopt=Y+(aopt-A2)* cbind(1,X2)%*%lasso2

  m=cv.glmnet(x=cbind(X1,A1,A1*X1),y=yopt,nfolds = 4,parallel = T)
  lasso1=coef(m,s='lambda.min')[12:22]

  aopt=as.numeric(cbind(1,X2)%*%lasso2_refit>0)
  yopt=Y+(aopt-A2)* cbind(1,X2)%*%lasso2_refit
  m=cv.glmnet(x=cbind(X1,A1,A1*X1),y=yopt,nfolds = 4,parallel = T)
  lasso1_refit=coef(m,s='lambda.min')
  index <- (1:21)[lasso1_refit[2:(22)]!=0]
  lasso1_refit[index+1]=coef(lm(yopt~cbind(X1,A1,A1*X1)[,index]))[-1]
  lasso1_refit=lasso1_refit[12:22]

  ## ad lasso
  pfac2=abs(1/coef(lm(Y~cbind(X2,A2,A2*X2)))[-1])
  pfac2=pfac2/sum(pfac2)*21

  m=cv.glmnet(x=cbind(X2,A2,A2*X2),y=Y,nfolds = 4,parallel = T,penalty.factor=pfac2)
  adlasso2=coef(m,s='lambda.min')[12:22]

  aopt=as.numeric(cbind(1,X2)%*%adlasso2>0)
  yopt=Y+(aopt-A2)* cbind(1,X2)%*%adlasso2

  pfac=abs(1/coef(lm(yopt~cbind(X1,A1,A1*X1)))[-1])
  pfac=pfac/sum(pfac)*21
  m=cv.glmnet(x=cbind(X1,A1,A1*X1),y=yopt,nfolds = 4,parallel = T,penalty.factor=pfac)
  al2[[k]]=adlasso2
  al1[[k]]=coef(m,s='lambda.min')[12:22]


  ### sail
  m=cv.sail(y=Y,e=A2,x=X2,nfolds = 4,
             parallel = T,basis=function(i) i)

  sail2=coef(m,s='lambda.min');sail2_refit=sail2
  index <- (1:21)[sail2[2:(22)]!=0]
  sail2_refit[index+1]=coef(lm(Y~cbind(X2,A2,A2*X2)[,index]))[-1]
  sail2=sail2[12:22];sail2_refit=sail2_refit[12:22]

  aopt=as.numeric(cbind(1,X2)%*%sail2>0)
  yopt=Y+(aopt-A2)* cbind(1,X2)%*%sail2

  m=cv.sail(y=yopt,e=A1,x=X1,nfolds = 4,
             parallel = T,basis=function(i) i)
  sail1=coef(m,s='lambda.min')[12:22]

  aopt=as.numeric(cbind(1,X2)%*%sail2_refit>0)
  yopt=Y+(aopt-A2)* cbind(1,X2)%*%sail2_refit

  m=cv.sail(y=yopt,e=A1,x=X1,nfolds = 4,
            parallel = T,basis=function(i) i)
  sail1_refit=coef(m,s='lambda.min')
  index <- (1:21)[sail1_refit[2:(22)]!=0]
  sail1_refit[index+1]=coef(lm(yopt~cbind(X1,A1,A1*X1)[,index]))[-1]
  sail1_refit=sail1_refit[12:22]

  ## ad sail

  pfac2=c(pfac2[11],pfac2[-11])
  m=cv.sail(y=Y,e=A2,x=X2,nfolds = 4,penalty.factor=pfac2,
            parallel = T,basis=function(i) i)

  adsail2=coef(m,s='lambda.min')[12:22]

  aopt=as.numeric(cbind(1,X2)%*%adsail2>0)
  yopt=Y+(aopt-A2)* cbind(1,X2)%*%adsail2

  pfac=abs(1/coef(lm(yopt~cbind(X1,A1,A1*X1)))[-1])
  pfac=pfac/sum(pfac)*21
  pfac=c(pfac[11],pfac[-11])
  m=cv.sail(y=yopt,e=A1,x=X1,nfolds = 4,penalty.factor=pfac,
            parallel = T,basis=function(i) i)
  as2[[k]]=adsail2
  as1[[k]]=coef(m,s='lambda.min')[12:22]

  ## a-learning
  m=PAL(Y~X2|A2,refit = F,penalty = 'LASSO',pi1.est=.5)
  pal2=m$beta1.est
  aopt=as.numeric(cbind(1,X2)%*%pal2>0)
  yopt=Y+(aopt-A2)* cbind(1,X2)%*%pal2
  m=PAL(yopt~X1|A1,refit = F,penalty = 'LASSO',pi1.est=.5)
  pal1=m$beta1.est

  m=PAL(Y~X2|A2,penalty = 'LASSO',pi1.est=.5)
  pal2_refit=m$beta1.est
  aopt=as.numeric(cbind(1,X2)%*%pal2_refit>0)
  yopt=Y+(aopt-A2)* cbind(1,X2)%*%pal2_refit
  m=PAL(yopt~X1|A1,penalty = 'LASSO',pi1.est=.5)
  pal1_refit=m$beta1.est

  l2[[k]]=lasso2;l1[[k]]=lasso1
  s2[[k]]=sail2;s1[[k]]=sail1
  a2[[k]]=pal2;a1[[k]]=pal1

  lr2[[k]]=lasso2_refit;lr1[[k]]=lasso1_refit
  sr2[[k]]=sail2_refit;sr1[[k]]=sail1_refit
  ar2[[k]]=pal2_refit;ar1[[k]]=pal1_refit
  k=k+1
}

for (i in 1:400) {
  l2[[i]]=res[[i]][[1]]
  l1[[i]]=res[[i]][[2]]
  s2[[i]]=res[[i]][[3]]
  s1[[i]]=res[[i]][[4]]
  a2[[i]]=res[[i]][[5]]
  a1[[i]]=res[[i]][[6]]
  lr2[[i]]=res[[i]][[7]]
  lr1[[i]]=res[[i]][[8]]
  sr2[[i]]=res[[i]][[9]]
  sr1[[i]]=res[[i]][[10]]
  ar2[[i]]=res[[i]][[11]]
  ar1[[i]]=res[[i]][[12]]
  al2[[i]]=res[[i]][[13]]
  al1[[i]]=res[[i]][[14]]
  as2[[i]]=res[[i]][[15]]
  as1[[i]]=res[[i]][[16]]
}


lasso2=do.call(cbind,l2);lasso1=do.call(cbind,l1)
sail2=do.call(cbind,s2);sail1=do.call(cbind,s1)
pal2=do.call(cbind,a2);pal1=do.call(cbind,a1)

lasso2_refit=do.call(cbind,lr2);lasso1_refit=do.call(cbind,lr1)
sail2_refit=do.call(cbind,sr2);sail1_refit=do.call(cbind,sr1)
pal2_refit=do.call(cbind,ar2);pal1_refit=do.call(cbind,ar1)


rowMeans(lasso2);rowMeans(sail2);rowMeans(pal2)
rowMeans(lasso2_refit);rowMeans(sail2_refit);rowMeans(pal2_refit)

rowMeans(lasso1);rowMeans(sail1);rowMeans(pal1)
rowMeans(lasso1_refit);rowMeans(sail1_refit);rowMeans(pal1_refit)

# Selection Rate

lasso_sel2=apply(lasso2!=0, 1, mean);lasso_sel1=apply(lasso1!=0, 1, mean)
sail_sel2=apply(sail2!=0, 1, mean);sail_sel1=apply(sail1!=0, 1, mean)
pal_sel2=apply(pal2!=0, 1, mean);pal_sel1=apply(pal1!=0, 1, mean)

lasso_sel2_refit=apply(lasso2_refit!=0, 1, mean);lasso_sel1_refit=apply(lasso1_refit!=0, 1, mean)
sail_sel2_refit=apply(sail2_refit!=0, 1, mean);sail_sel1_refit=apply(sail1_refit!=0, 1, mean)
pal_sel2_refit=apply(pal2_refit!=0, 1, mean);pal_sel1_refit=apply(pal1_refit!=0, 1, mean)


sel_rate=cbind(sail_sel1,lasso_sel1,pal_sel1,sail_sel1_refit,lasso_sel1_refit,
               pal_sel1_refit, sail_sel2, lasso_sel2,pal_sel2,sail_sel2_refit,lasso_sel2_refit,
               pal_sel2_refit)

colnames(sel_rate)=c('lasso (stage 1)', 'relaxed lasso (stage 1)',
                     'sail (stage 1)','relaxed sail (stage 1)', 'PAL (stage 1)',
                     'relaxed PAL (stage 1)',
                     'lasso (stage 2)', 'relaxed lasso (stage 2)',
                     'sail (stage 2)','relaxed sail (stage 2)', 'PAL (stage 2)',
                     'relaxed PAL (stage 2)')

rownames(sel_rate)=c('A','X1','Noise1','Noise2','Noise3','Noise4','Noise5','Noise6','Noise7','Noise8',
                    'Noise9')
sel_rate=round(sel_rate*100,2)

## Error Rate
set.seed(999);ds_test=g2(10000,10);Xtest_1=ds_test[[2]];Xtest_2=ds_test[[3]]

opt1=-1.190399+1.5*Xtest_1[,1]>0;opt2=1-1.5*Xtest_2[,1]>0

lasso_opt1=apply(lasso1, 2, function(i) cbind(1,Xtest_1)%*%i)>0
lasso_opt2=apply(lasso2, 2, function(i) cbind(1,Xtest_2)%*%i)>0

lasso_opt1_refit=apply(lasso1_refit, 2, function(i) cbind(1,Xtest_1)%*%i)>0
lasso_opt2_refit=apply(lasso2_refit, 2, function(i) cbind(1,Xtest_2)%*%i)>0

sail_opt1=apply(sail1, 2, function(i) cbind(1,Xtest_1)%*%i)>0
sail_opt2=apply(sail2, 2, function(i) cbind(1,Xtest_2)%*%i)>0

sail_opt1_refit=apply(sail1_refit, 2, function(i) cbind(1,Xtest_1)%*%i)>0
sail_opt2_refit=apply(sail2_refit, 2, function(i) cbind(1,Xtest_2)%*%i)>0

pal_opt1=apply(pal1, 2, function(i) cbind(1,Xtest_1)%*%i)>0
pal_opt2=apply(pal2, 2, function(i) cbind(1,Xtest_2)%*%i)>0

pal_opt1_refit=apply(pal1_refit, 2, function(i) cbind(1,Xtest_1)%*%i)>0
pal_opt2_refit=apply(pal2_refit, 2, function(i) cbind(1,Xtest_2)%*%i)>0

### error!!

lasso_error1=1-mean(colMeans(apply(lasso_opt1, 2, function(i) i==opt1)))
lasso_error2=1-mean(colMeans(apply(lasso_opt2, 2, function(i) i==opt2)))

sail_error1=1-mean(colMeans(apply(sail_opt1, 2, function(i) i==opt1)))
sail_error2=1-mean(colMeans(apply(sail_opt2, 2, function(i) i==opt2)))

pal_error1=1-mean(colMeans(apply(pal_opt1, 2, function(i) i==opt1)))
pal_error2=1-mean(colMeans(apply(pal_opt2, 2, function(i) i==opt2)))

## refit one

lasso_error1_refit=1-mean(colMeans(apply(lasso_opt1_refit, 2, function(i) i==opt1)))
lasso_error2_refit=1-mean(colMeans(apply(lasso_opt2_refit, 2, function(i) i==opt2)))

sail_error1_refit=1-mean(colMeans(apply(sail_opt1_refit, 2, function(i) i==opt1)))
sail_error2_refit=1-mean(colMeans(apply(sail_opt2_refit, 2, function(i) i==opt2)))

pal_error1_refit=1-mean(colMeans(apply(pal_opt1_refit, 2, function(i) i==opt1)))
pal_error2_refit=1-mean(colMeans(apply(pal_opt2_refit, 2, function(i) i==opt2)))

err_rate=cbind(sail_error1,lasso_error1,pal_error1,sail_error1_refit,lasso_error1_refit,
               pal_error1_refit,sail_error2,lasso_error2,pal_error2,
               sail_error2_refit,lasso_error2_refit,pal_error2_refit)


## total error
lopt1=apply(lasso_opt1, 2, function(i) i==opt1)==T
lopt2=apply(lasso_opt2, 2, function(i) i==opt2)==T
k=0
for (i in 1:length(lopt1)) {
  if (lopt1[i]==F|lopt2[i]==F) {
    k=k+1
  }
}

lasso_error=k/length(lopt1)

sopt1=apply(sail_opt1, 2, function(i) i==opt1)==T
sopt2=apply(sail_opt2, 2, function(i) i==opt2)==T
k=0
for (i in 1:length(lopt1)) {
  if (sopt1[i]==F|sopt2[i]==F) {
    k=k+1
  }
}
sail_error=k/length(lopt1)

aopt1=apply(pal_opt1, 2, function(i) i==opt1)==T
aopt2=apply(pal_opt2, 2, function(i) i==opt2)==T
k=0
for (i in 1:length(lopt1)) {
  if (aopt1[i]==F|aopt2[i]==F) {
    k=k+1
  }
}

pal_error=k/length(lopt1)

lopt1=apply(lasso_opt1_refit, 2, function(i) i==opt1)==T
lopt2=apply(lasso_opt2_refit, 2, function(i) i==opt2)==T
k=0
for (i in 1:length(lopt1)) {
  if (lopt1[i]==F|lopt2[i]==F) {
    k=k+1
  }
}

lasso_error_refit=k/length(lopt1)

sopt1=apply(sail_opt1_refit, 2, function(i) i==opt1)==T
sopt2=apply(sail_opt2_refit, 2, function(i) i==opt2)==T
k=0
for (i in 1:length(lopt1)) {
  if (sopt1[i]==F|sopt2[i]==F) {
    k=k+1
  }
}
sail_error_refit=k/length(lopt1)

aopt1=apply(pal_opt1_refit, 2, function(i) i==opt1)==T
aopt2=apply(pal_opt2_refit, 2, function(i) i==opt2)==T
k=0
for (i in 1:length(lopt1)) {
  if (aopt1[i]==F|aopt2[i]==F) {
    k=k+1
  }
}

pal_error_refit=k/length(lopt1)


error_t=c(sail_error,lasso_error,pal_error,
          sail_error_refit,lasso_error_refit,pal_error_refit)

colnames(err_rate)=c('lasso (stage 1)', 'relaxed lasso (stage 1)',
                     'sail (stage 1)','relaxed sail (stage 1)', 'PAL (stage 1)',
                     'relaxed PAL (stage 1)',
                     'lasso (stage 2)', 'relaxed lasso (stage 2)',
                     'sail (stage 2)','relaxed sail (stage 2)', 'PAL (stage 2)',
                     'relaxed PAL (stage 2)')


err_rate=round(err_rate*100,1)

rownames(err_rate)=c('Error Rate')

err_rate=round(err_rate,2)

# Value Function
mu=-1+2*X1[,1]+X2[,1]-A1+1.5*X1[,1]*A1+A2-1.5*X2[,1]*A2


Value_T=mean(-1+Xtest_2[,1]+2*Xtest_1[,1]-opt1+1.5*Xtest_1[,1]*opt1+opt2-1.5*Xtest_2[,1]*opt2)

v_0=mean(-1+Xtest_2[,1]+2*Xtest_1[,1])

v_1=mean(-1+Xtest_2[,1]+2*Xtest_1[,1]-1+1.5*Xtest_1[,1]*1+1-1.5*Xtest_2[,1]*1)

value_lasso=mean(-1+Xtest_2[,1]+2*Xtest_1[,1]-lasso_opt1+1.5*Xtest_1[,1]*lasso_opt1+lasso_opt2-1.5*Xtest_2[,1]*lasso_opt2)

value_sail=mean(-1+Xtest_2[,1]+2*Xtest_1[,1]-sail_opt1+1.5*Xtest_1[,1]*sail_opt1+sail_opt2-1.5*Xtest_2[,1]*sail_opt2)

value_pal=mean(-1+Xtest_2[,1]+2*Xtest_1[,1]-pal_opt1+1.5*Xtest_1[,1]*pal_opt1+pal_opt2-1.5*Xtest_2[,1]*pal_opt2)


## refit

value_lasso_refit=mean(-1+Xtest_2[,1]+2*Xtest_1[,1]-lasso_opt1_refit
                       +1.5*Xtest_1[,1]*lasso_opt1_refit+lasso_opt2_refit-1.5*Xtest_2[,1]*lasso_opt2_refit)


value_sail_refit=mean(-1+Xtest_2[,1]+2*Xtest_1[,1]-sail_opt1_refit
                      +1.5*Xtest_1[,1]*sail_opt1_refit+sail_opt2_refit-1.5*Xtest_2[,1]*sail_opt2_refit)


value_pal_refit=mean(-1+Xtest_2[,1]+2*Xtest_1[,1]-pal_opt1_refit
                     +1.5*Xtest_1[,1]*pal_opt1_refit+pal_opt2_refit-1.5*Xtest_2[,1]*pal_opt2_refit)


value_function=rbind(value_sail,value_lasso,value_pal,
                     value_sail_refit,value_lasso_refit,
                     value_pal_refit)

rownames(value_function)=c('True Value Function', 'Estimated Value Function (lasso)',
                           'Estimated Value Function (relaxed lasso)',
                           'Estimated Value Function (sail)', 'Estimated Value Function (relaxed sail)',
                           'Estimated Value Function (PAL)','Estimated Value Function (relaxed PAL)')

colnames(value_function)=c('Value Function')

value_function=round(value_function,3)

save.image("~/Desktop/weights/multiple.RData")
load("~/Desktop/weights/multiple.RData")

## dwols setting

library(doParallel)
registerDoParallel(cores = 4)
library(glmnet)
devtools::load_all()
k=1
l2=l1=s2=s1=a2=a1=lr2=lr1=sr2=sr1=ar2=ar1=al1=al2=as1=as2=list()
while (k<=100) {
  ds=g3(8000,10);X1=ds[[1]];X2=ds[[2]];Y=ds[[3]];w1=ds[[4]];w2=ds[[5]];A1=ds[[6]];A2=ds[[7]]
  pfac2=coef(lm(Y~cbind(X2,A2,A2*X2),weights = w2))[-1]
  pfac2=abs(1/pfac2)
  pfac2=pfac2/sum(pfac2) *21

  m=cv.glmnet(x=cbind(X2,A2,A2*X2),y=Y,nfolds = 4,parallel = T,weights = w2,
              penalty.factor=pfac2)
  lasso2=coef(m,s='lambda.min');lasso2_refit=lasso2
  index <- (1:21)[lasso2[2:(22)]!=0]
  lasso2_refit[index+1]=coef(lm(Y~cbind(X2,A2,A2*X2)[,index],weights = w2))[-1]
  lasso2=lasso2[12:22];lasso2_refit=lasso2_refit[12:22]

  aopt=as.numeric(cbind(1,X2)%*%lasso2>0)
  yopt=Y+(aopt-A2)* cbind(1,X2)%*%lasso2
  pfac=abs(1/coef(lm(yopt~cbind(X1,A1,A1*X1),weights = w1))[-1])
  pfac=pfac/sum(pfac)*21

  m=cv.glmnet(x=cbind(X1,A1,A1*X1),y=yopt,nfolds = 4,parallel = T,weights = w1,penalty.factor=pfac)
  lasso1=coef(m,s='lambda.min')[12:22]

  aopt=as.numeric(cbind(1,X2)%*%lasso2_refit>0)
  yopt=Y+(aopt-A2)* cbind(1,X2)%*%lasso2_refit
  pfac=abs(1/coef(lm(yopt~cbind(X1,A1,A1*X1),weights = w1))[-1])
  pfac=pfac/sum(pfac)*21

  m=cv.glmnet(x=cbind(X1,A1,A1*X1),y=yopt,nfolds = 4,parallel = T,weights = w1,penalty.factor=pfac)
  lasso1_refit=coef(m,s='lambda.min')
  index <- (1:21)[lasso1_refit[2:(22)]!=0]
  lasso1_refit[index+1]=coef(lm(yopt~cbind(X1,A1,A1*X1)[,index],weights = w1))[-1]
  lasso1_refit=lasso1_refit[12:22]

  ### sail
  pfac2=c(pfac2[11],pfac2[-11])
  m=cv.sail(y=Y,e=A2,x=X2,nfolds = 4,weights = w2,
            parallel = T, basis=function(i) i)

  sail2=coef(m,s='lambda.min');sail2_refit=sail2
  index <- (1:21)[sail2[2:(22)]!=0]
  sail2_refit[index+1]=coef(lm(Y~cbind(X2,A2,A2*X2)[,index],weights = w2))[-1]
  sail2=sail2[12:22];sail2_refit=sail2_refit[12:22]

  aopt=as.numeric(cbind(1,X2)%*%sail2>0)
  yopt=Y+(aopt-A2)* cbind(1,X2)%*%sail2
  pfac=abs(1/coef(lm(yopt~cbind(X1,A1,A1*X1),weights = w1))[-1])
  pfac=pfac/sum(pfac)*21
  pfac=c(pfac[11],pfac[-11])

  m=cv.sail(y=yopt,e=A1,x=X1,nfolds = 4,weights = w1,
            parallel = T,basis=function(i) i)
  sail1=coef(m,s='lambda.min')[12:22]

  aopt=as.numeric(cbind(1,X2)%*%sail2_refit>0)
  yopt=Y+(aopt-A2)* cbind(1,X2)%*%sail2_refit
  pfac=abs(1/coef(lm(yopt~cbind(X1,A1,A1*X1),weights = w1))[-1])
  pfac=pfac/sum(pfac)*21
  pfac=c(pfac[11],pfac[-11])

  m=cv.sail(y=yopt,e=A1,x=X1,nfolds = 4,weights = w1,penalty.factor=pfac,
            parallel = T,basis=function(i) i)
  sail1_refit=coef(m,s='lambda.min')
  index <- (1:21)[sail1_refit[2:(22)]!=0]
  sail1_refit[index+1]=coef(lm(yopt~cbind(X1,A1,A1*X1)[,index],weights = w1))[-1]
  sail1_refit=sail1_refit[12:22]

  ## a-learning
  m=PAL(Y~X2|A2,refit = F,penalty = 'LASSO')
  pal2=m$beta1.est
  aopt=as.numeric(cbind(1,X2)%*%pal2>0)
  yopt=Y+(aopt-A2)* cbind(1,X2)%*%pal2
  m=PAL(yopt~X1|A1,refit = F,penalty = 'LASSO')
  pal1=m$beta1.est

  m=PAL(Y~X2|A2,penalty = 'LASSO')
  pal2_refit=m$beta1.est
  aopt=as.numeric(cbind(1,X2)%*%pal2_refit>0)
  yopt=Y+(aopt-A2)* cbind(1,X2)%*%pal2_refit
  m=PAL(yopt~X1|A1,penalty = 'LASSO')
  pal1_refit=m$beta1.est

  l2[[k]]=lasso2;l1[[k]]=lasso1
  s2[[k]]=sail2;s1[[k]]=sail1
  a2[[k]]=pal2;a1[[k]]=pal1

  lr2[[k]]=lasso2_refit;lr1[[k]]=lasso1_refit
  sr2[[k]]=sail2_refit;sr1[[k]]=sail1_refit
  ar2[[k]]=pal2_refit;ar1[[k]]=pal1_refit
  k=k+1
}


save.image('dwols.RData')
lasso2=do.call(cbind,l2);lasso1=do.call(cbind,l1)
sail2=do.call(cbind,s2);sail1=do.call(cbind,s1)
pal2=do.call(cbind,a2);pal1=do.call(cbind,a1)

lasso2_refit=do.call(cbind,lr2);lasso1_refit=do.call(cbind,lr1)
sail2_refit=do.call(cbind,sr2);sail1_refit=do.call(cbind,sr1)
pal2_refit=do.call(cbind,ar2);pal1_refit=do.call(cbind,ar1)


rowMeans(lasso2);rowMeans(sail2);rowMeans(pal2)
rowMeans(lasso2_refit);rowMeans(sail2_refit);rowMeans(pal2_refit)

rowMeans(lasso1);rowMeans(sail1);rowMeans(pal1)
rowMeans(lasso1_refit);rowMeans(sail1_refit);rowMeans(pal1_refit)


# Selection Rate

lasso_sel2=apply(lasso2!=0, 1, mean);lasso_sel1=apply(lasso1!=0, 1, mean)
sail_sel2=apply(sail2!=0, 1, mean);sail_sel1=apply(sail1!=0, 1, mean)
pal_sel2=apply(pal2!=0, 1, mean);pal_sel1=apply(pal1!=0, 1, mean)

lasso_sel2_refit=apply(lasso2_refit!=0, 1, mean);lasso_sel1_refit=apply(lasso1_refit!=0, 1, mean)
sail_sel2_refit=apply(sail2_refit!=0, 1, mean);sail_sel1_refit=apply(sail1_refit!=0, 1, mean)
pal_sel2_refit=apply(pal2_refit!=0, 1, mean);pal_sel1_refit=apply(pal1_refit!=0, 1, mean)


sel_rate=cbind(lasso_sel1,lasso_sel1_refit,sail_sel1,sail_sel1_refit,
               pal_sel1,pal_sel1_refit,  lasso_sel2,lasso_sel2_refit,sail_sel2,
               sail_sel2_refit,pal_sel2,pal_sel2_refit)

colnames(sel_rate)=c('lasso (stage 1)', 'relaxed lasso (stage 1)',
                     'sail (stage 1)','relaxed sail (stage 1)', 'PAL (stage 1)',
                     'relaxed PAL (stage 1)',
                     'lasso (stage 2)', 'relaxed lasso (stage 2)',
                     'sail (stage 2)','relaxed sail (stage 2)', 'PAL (stage 2)',
                     'relaxed PAL (stage 2)')

rownames(sel_rate)=c('A','X1','Noise1','Noise2','Noise3','Noise4','Noise5','Noise6','Noise7','Noise8',
                     'Noise9')

sel_rate=round(sel_rate*100,0)

## Error Rate
set.seed(999);ds_test=g3(10000,10);Xtest_1=ds_test[[1]];Xtest_2=ds_test[[2]]

opt1=0.8-2*Xtest_1[,1]>0;opt2=1-1.5*Xtest_2[,1]>0

lasso_opt1=apply(lasso1, 2, function(i) cbind(1,Xtest_1)%*%i)>0
lasso_opt2=apply(lasso2, 2, function(i) cbind(1,Xtest_2)%*%i)>0

lasso_opt1_refit=apply(lasso1_refit, 2, function(i) cbind(1,Xtest_1)%*%i)>0
lasso_opt2_refit=apply(lasso2_refit, 2, function(i) cbind(1,Xtest_2)%*%i)>0

sail_opt1=apply(sail1, 2, function(i) cbind(1,Xtest_1)%*%i)>0
sail_opt2=apply(sail2, 2, function(i) cbind(1,Xtest_2)%*%i)>0

sail_opt1_refit=apply(sail1_refit, 2, function(i) cbind(1,Xtest_1)%*%i)>0
sail_opt2_refit=apply(sail2_refit, 2, function(i) cbind(1,Xtest_2)%*%i)>0

pal_opt1=apply(pal1, 2, function(i) cbind(1,Xtest_1)%*%i)>0
pal_opt2=apply(pal2, 2, function(i) cbind(1,Xtest_2)%*%i)>0

pal_opt1_refit=apply(pal1_refit, 2, function(i) cbind(1,Xtest_1)%*%i)>0
pal_opt2_refit=apply(pal2_refit, 2, function(i) cbind(1,Xtest_2)%*%i)>0

### error!!

lasso_error1=1-mean(colMeans(apply(lasso_opt1, 2, function(i) i==opt1)))
lasso_error2=1-mean(colMeans(apply(lasso_opt2, 2, function(i) i==opt2)))

sail_error1=1-mean(colMeans(apply(sail_opt1, 2, function(i) i==opt1)))
sail_error2=1-mean(colMeans(apply(sail_opt2, 2, function(i) i==opt2)))

pal_error1=1-mean(colMeans(apply(pal_opt1, 2, function(i) i==opt1)))
pal_error2=1-mean(colMeans(apply(pal_opt2, 2, function(i) i==opt2)))

## refit one

lasso_error1_refit=1-mean(colMeans(apply(lasso_opt1_refit, 2, function(i) i==opt1)))
lasso_error2_refit=1-mean(colMeans(apply(lasso_opt2_refit, 2, function(i) i==opt2)))

sail_error1_refit=1-mean(colMeans(apply(sail_opt1_refit, 2, function(i) i==opt1)))
sail_error2_refit=1-mean(colMeans(apply(sail_opt2_refit, 2, function(i) i==opt2)))

pal_error1_refit=1-mean(colMeans(apply(pal_opt1_refit, 2, function(i) i==opt1)))
pal_error2_refit=1-mean(colMeans(apply(pal_opt2_refit, 2, function(i) i==opt2)))



err_rate=cbind(lasso_error1,lasso_error1_refit,sail_error1,sail_error1_refit,
               pal_error1,pal_error1_refit,lasso_error2,lasso_error2_refit,
               sail_error2,sail_error2_refit,pal_error2,pal_error2_refit)


## total error
lopt1=apply(lasso_opt1, 2, function(i) i==opt1)==T
lopt2=apply(lasso_opt2, 2, function(i) i==opt2)==T
k=0
for (i in 1:length(lopt1)) {
  if (lopt1[i]==F|lopt2[i]==F) {
    k=k+1
  }
}

lasso_error=k/length(lopt1)

sopt1=apply(sail_opt1, 2, function(i) i==opt1)==T
sopt2=apply(sail_opt2, 2, function(i) i==opt2)==T
k=0
for (i in 1:length(lopt1)) {
  if (sopt1[i]==F|sopt2[i]==F) {
    k=k+1
  }
}
sail_error=k/length(lopt1)

aopt1=apply(pal_opt1, 2, function(i) i==opt1)==T
aopt2=apply(pal_opt2, 2, function(i) i==opt2)==T
k=0
for (i in 1:length(lopt1)) {
  if (aopt1[i]==F|aopt2[i]==F) {
    k=k+1
  }
}

pal_error=k/length(lopt1)

lopt1=apply(lasso_opt1_refit, 2, function(i) i==opt1)==T
lopt2=apply(lasso_opt2_refit, 2, function(i) i==opt2)==T
k=0
for (i in 1:length(lopt1)) {
  if (lopt1[i]==F|lopt2[i]==F) {
    k=k+1
  }
}

lasso_error_refit=k/length(lopt1)

sopt1=apply(sail_opt1_refit, 2, function(i) i==opt1)==T
sopt2=apply(sail_opt2_refit, 2, function(i) i==opt2)==T
k=0
for (i in 1:length(lopt1)) {
  if (sopt1[i]==F|sopt2[i]==F) {
    k=k+1
  }
}
sail_error_refit=k/length(lopt1)

aopt1=apply(pal_opt1_refit, 2, function(i) i==opt1)==T
aopt2=apply(pal_opt2_refit, 2, function(i) i==opt2)==T
k=0
for (i in 1:length(lopt1)) {
  if (aopt1[i]==F|aopt2[i]==F) {
    k=k+1
  }
}

pal_error_refit=k/length(lopt1)


error_t=c(sail_error,lasso_error,pal_error,
          sail_error_refit,lasso_error_refit,pal_error_refit)

colnames(err_rate)=c('lasso (stage 1)', 'relaxed lasso (stage 1)',
                     'sail (stage 1)','relaxed sail (stage 1)', 'PAL (stage 1)',
                     'relaxed PAL (stage 1)',
                     'lasso (stage 2)', 'relaxed lasso (stage 2)',
                     'sail (stage 2)','relaxed sail (stage 2)', 'PAL (stage 2)',
                     'relaxed PAL (stage 2)')

rownames(err_rate)=c('Error Rate')

err_rate=round(err_rate*100,2)

# Value Function

Value_T=mean(.5+2*(Xtest_1[,1]+2*Xtest_1[,2]))

v_0=mean(.5+2*(Xtest_1[,1]+2*Xtest_1[,2])
              +(0-opt1)*(.8-2*Xtest_1[,1])+
                (0-opt2)*(1-1.5*Xtest_2[,1]))


v_1=mean(.5+2*(Xtest_1[,1]+2*Xtest_1[,2])
         +(1-opt1)*(.8-2*Xtest_1[,1])+
           (1-opt2)*(1-1.5*Xtest_2[,1]))



value_lasso=mean(colMeans(.5+2*(Xtest_1[,1]+2*Xtest_1[,2])
                          +(lasso_opt1-opt1)*(.8-2*Xtest_1[,1])+
                            (lasso_opt2-opt2)*(1-1.5*Xtest_2[,1])))

value_sail=mean(colMeans(.5+2*(Xtest_1[,1]+2*Xtest_1[,2])
                         +(sail_opt1-opt1)*(.8-2*Xtest_1[,1])+
                           (sail_opt2-opt2)*(1-1.5*Xtest_2[,1])))

value_pal=mean(colMeans(.5+2*(Xtest_1[,1]+2*Xtest_1[,2])
                         +(pal_opt1-opt1)*(.8-2*Xtest_1[,1])+
                           (pal_opt2-opt2)*(1-1.5*Xtest_2[,1])))


value_lasso_refit=mean(colMeans(.5+2*(Xtest_1[,1]+2*Xtest_1[,2])
                          +(lasso_opt1_refit-opt1)*(.8-2*Xtest_1[,1])+
                            (lasso_opt2_refit-opt2)*(1-1.5*Xtest_2[,1])))

value_sail_refit=mean(colMeans(.5+2*(Xtest_1[,1]+2*Xtest_1[,2])
                         +(sail_opt1_refit-opt1)*(.8-2*Xtest_1[,1])+
                           (sail_opt2_refit-opt2)*(1-1.5*Xtest_2[,1])))

value_pal_refit=mean(colMeans(.5+2*(Xtest_1[,1]+2*Xtest_1[,2])
                        +(pal_opt1_refit-opt1)*(.8-2*Xtest_1[,1])+
                          (pal_opt2_refit-opt2)*(1-1.5*Xtest_2[,1])))

value_function=cbind(value_sail,value_lasso,value_pal,
                     value_sail_refit,value_lasso_refit,
                     value_pal_refit)

error_t=100*error_t

value_function=rbind(value_function,error_t)

colnames(value_function)=c('pdWOLS','Q-Learning (lasso)','PAL',
                           'refitted pdWOLS','Q-Learning (refitted lasso)',
                           'refitted PAL')

rownames(value_function)=c('Value Function','Total Error')

value_function=round(value_function,2)

save.image("~/Desktop/weights/dwols.RData")
load("~/Desktop/weights/dwols.RData")

