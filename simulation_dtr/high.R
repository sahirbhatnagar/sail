library(doParallel)
library(ITRSelect)
registerDoParallel(cores = 4)
library(glmnet)
devtools::load_all()
k=1
l=s=a=list()
lpfac=c(rep(1,300),0,rep(1,300))
lpfac=lpfac/sum(lpfac) *601
spfac=c(lpfac[301],lpfac[-301])

while (k<=100) {
  ds=gen_high(200,300);Y=ds[[2]];X=ds[[1]];A=ds[[3]]
  m=cv.glmnet(x=cbind(X,A,A*X),y=Y,nfolds = 4,parallel = T, penalty.factor=lpfac)
  lasso=coef(m,s='lambda.min')
  lasso_refit=lasso
  index <- (1:601)[lasso[2:602]!=0]
  lasso_refit[index+1]=coef(lm(Y~cbind(X,A,A*X)[,index]))[-1]
  lasso=lasso_refit[302:602]

  ### sail
  m=cv.sail(y=Y,e=A,x=X,nfolds = 4,penalty.factor=spfac,
            parallel = T,basis=function(i) i)
  sail=coef(m,s='lambda.min')
  sail_refit=sail
  index <- (1:601)[sail[2:602]!=0]
  sail_refit[index+1]=coef(lm(Y~cbind(X,A,A*X)[,index]))[-1]
  sail=sail_refit[302:602]

  ## a-learning
  m=PAL(Y~X|A,refit = T,penalty = 'LASSO',pi1.est=.5)
  pal=m$beta1.est
  l[[k]]=lasso; s[[k]]=sail; a[[k]]=pal
  k=k+1
}

save.image("~/Desktop/weights/high.RData")
load("~/Desktop/weights/high.RData")

lasso=do.call(cbind,l);sail=do.call(cbind,s);pal=do.call(cbind,a)

rowMeans(lasso);rowMeans(sail);rowMeans(pal)
d=cbind(rowMeans(lasso),rowMeans(sail),rowMeans(pal))
colnames(d)=c('lasso','sail','pal')

# Selection Rate
true_sel=c(T,T,T,T,rep(F,297))

lasso_sel=apply(lasso!=0, 2, function (i) i==true_sel)
sail_sel=apply(sail!=0, 2, function (i) i==true_sel)
pal_sel=apply(pal!=0, 2, function (i) i==true_sel)

fn_lasso=1-mean(apply(lasso_sel[1:4,],2, mean))
fn_sail=1-mean(apply(sail_sel[1:4,],2, mean))
fn_pal=1-mean(apply(pal_sel[1:4,],2, mean))

fn=c(fn_lasso,fn_sail,fn_pal)

fp_lasso=1-mean(apply(lasso_sel[-(1:4),],2, mean))
fp_sail=1-mean(apply(sail_sel[-(1:4),],2, mean))
fp_pal=1-mean(apply(pal_sel[-(1:4),],2, mean))

fp=c(fp_lasso,fp_sail,fp_pal)



## Error Rate
set.seed(999);ds_test=gen_high(10000,300);Xtest=ds_test[[1]]
psi=c(1,-1.5,-1.5,1)
opt=cbind(1,Xtest[,1:3])%*%psi>0

lasso_opt=apply(lasso, 2, function(i) cbind(1,Xtest)%*%i)>0

sail_opt=apply(sail, 2, function(i) cbind(1,Xtest)%*%i)>0

pal_opt=apply(pal, 2, function(i) cbind(1,Xtest)%*%i)>0


### error!!

lasso_error=1-mean(colMeans(apply(lasso_opt, 2, function(i) i==opt)))

sail_error=1-mean(colMeans(apply(sail_opt, 2, function(i) i==opt)))

pal_error=1-mean(colMeans(apply(pal_opt, 2, function(i) i==opt)))

err_rate=cbind(lasso_error,sail_error,pal_error)

colnames(err_rate)=c('lasso','sail','pal')

rownames(err_rate)=c('Error Rate')

err_rate=round(err_rate,2)

# Value Function
tfree=3+3*exp(Xtest[,1])+3*Xtest[,2]+3*Xtest[,3]+3*Xtest[,4]

Value_T=mean(tfree+opt*(cbind(1,Xtest[,1:3])%*%psi))

value_lasso=mean(tfree+lasso_opt*drop((cbind(1,Xtest[,1:3])%*%psi)))
value_sail=mean(tfree+sail_opt*drop((cbind(1,Xtest[,1:3])%*%psi)))
value_pal=mean(tfree+pal_opt*drop((cbind(1,Xtest[,1:3])%*%psi)))

value_function=c(value_lasso,
                     value_sail,value_pal)

rownames(value_function)=c('True Value Function', 'Estimated Value Function (lasso)',
                           'Estimated Value Function (sail)','Estimated Value Function (PAL)')

colnames(value_function)=c('Value Function')

value_function=round(value_function,3)

save.image("~/Desktop/weights/high.RData")

vr=value_function/Value_T

res=rbind(fn,fp,err_rate,vr)

colnames(res)=c('lasso','pdWOLS','PAL')
rownames(res)=c('False Negative','False Positive','Error Rate','Value Ratio')

res=round(res*100,2)

kable(res,  'latex', booktabs = T,
      caption = "Variable Selection Results (%)")%>%
  kable_styling(latex_options = "HOLD_position",position = 'center')

### plot

psi_0=c(lasso[1,],sail[1,],pal[1,])
psi_1=c(lasso[2,],sail[2,],pal[2,])
psi_2=c(lasso[3,],sail[3,],pal[3,])
psi_3=c(lasso[4,],sail[4,],pal[4,])

Estimates=c(psi_0,psi_1,psi_2,psi_3)
methods=c(rep('lasso',100),rep('pdWOLS',100),rep('PAL',100))
methods=rep(methods,4)

par=rep(c('psi0','psi1','psi2','psi3'),each=300)

ds=data.frame(Estimates,methods,par)
library(latex2exp)
levels(ds$par)=c(TeX('$\\psi_{0}$'),TeX('$\\psi_{1}$'),
                 TeX('$\\psi_{2}$'),TeX('$\\psi_{3}$'))

true=rep(c(1,-1.5,-1.5,1),each=300)

ind=rep(c(1,2,3,4),each=300)

ds=cbind(ds,true,ind)

library(ggplot2)
ggplot(ds, aes(x=ind, y=Estimates, fill=methods)) +
  geom_boxplot()+
  facet_grid(vars(par),scales='free',labeller = label_parsed)+
  geom_hline(aes(yintercept=true),linetype="dashed")+
  xlab('Blip Parameters')+
  ylab('Estimates')








