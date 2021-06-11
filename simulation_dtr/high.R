library(doParallel)
library(ITRSelect)
registerDoParallel(cores = 4)
library(glmnet)
devtools::load_all()

ds=gen_high(200,400);Y=ds[[2]];X=ds[[1]];A=ds[[3]]
m=cv.glmnet(x=cbind(X,A,A*X),y=Y,nfolds = 4)
lasso=coef(m,s='lambda.min')
l=lasso[402:802]
lasso_refit=lasso
index <- (1:801)[lasso[2:802]!=0]
lasso_refit[index+1]=coef(lm(Y~cbind(X,A,A*X)[,index]))[-1]
lr=lasso_refit[402:802]

### sail
m=cv.sail(y=Y,e=A,x=X,nfolds = 4,parallel=T,
          basis=function(i) i)
sail=coef(m,s='lambda.min')
s=sail[402:802]
sail_refit=sail
index <- (1:801)[sail[2:802]!=0]
sail_refit[index+1]=coef(lm(Y~cbind(X,A,A*X)[,index]))[-1]
sr=sail_refit[402:802]

load("~/Desktop/weights/high.RData")
l=lr=s=sr=a=ar=list()
for (i in 1:400) {
  l[[i]]=res[[i]][[1]]
  lr[[i]]=res[[i]][[2]]
  s[[i]]=res[[i]][[3]]
  sr[[i]]=res[[i]][[4]]
  a[[i]]=res[[i]][[5]]
  ar[[i]]=res[[i]][[6]]
}
rm(res);rm(lpfac);rm(spfac)
save.image("~/Desktop/weights/high.RData")

lasso=do.call(cbind,l);sail=do.call(cbind,s);pal=do.call(cbind,a)
lasso_refit=do.call(cbind,lr);sail_refit=do.call(cbind,sr);pal_refit=do.call(cbind,ar)

d=cbind(rowMeans(lasso),rowMeans(sail),rowMeans(pal),
        rowMeans(lasso_refit),rowMeans(sail_refit),rowMeans(pal_refit))

colnames(d)=c('lasso','sail','pal','relaxed lasso','relaxed sail','relaxed pal')

# Selection Rate
true_sel=c(T,T,rep(F,399))

lasso_sel=apply(lasso!=0, 2, function (i) i==true_sel)
sail_sel=apply(sail!=0, 2, function (i) i==true_sel)
pal_sel=apply(pal!=0, 2, function (i) i==true_sel)

fn_lasso=1-mean(apply(lasso_sel[1:2,],2, mean))
fn_sail=1-mean(apply(sail_sel[1:2,],2, mean))
fn_pal=1-mean(apply(pal_sel[1:2,],2, mean))

fn=c(fn_lasso,fn_sail,fn_pal)

fp_lasso=1-mean(apply(lasso_sel[-(1:2),],2, mean))
fp_sail=1-mean(apply(sail_sel[-(1:2),],2, mean))
fp_pal=1-mean(apply(pal_sel[-(1:2),],2, mean))

fp=c(fp_lasso,fp_sail,fp_pal)

## Error Rate
set.seed(999);ds_test=gen_high(10000,400);Xtest=ds_test[[1]]
psi=c(1,-1.5)
opt=cbind(1,Xtest[,1])%*%psi>0

lasso_opt=apply(lasso, 2, function(i) cbind(1,Xtest)%*%i)>0

sail_opt=apply(sail, 2, function(i) cbind(1,Xtest)%*%i)>0

pal_opt=apply(pal, 2, function(i) cbind(1,Xtest)%*%i)>0

lasso_opt_refit=apply(lasso_refit, 2, function(i) cbind(1,Xtest)%*%i)>0

sail_opt_refit=apply(sail_refit, 2, function(i) cbind(1,Xtest)%*%i)>0

pal_opt_refit=apply(pal_refit, 2, function(i) cbind(1,Xtest)%*%i)>0


### error!!

lasso_error=1-mean(colMeans(apply(lasso_opt, 2, function(i) i==opt)))

sail_error=1-mean(colMeans(apply(sail_opt, 2, function(i) i==opt)))

pal_error=1-mean(colMeans(apply(pal_opt, 2, function(i) i==opt)))

lasso_error_refit=1-mean(colMeans(apply(lasso_opt_refit, 2, function(i) i==opt)))

sail_error_refit=1-mean(colMeans(apply(sail_opt_refit, 2, function(i) i==opt)))

pal_error_refit=1-mean(colMeans(apply(pal_opt_refit, 2, function(i) i==opt)))

err_rate=cbind(lasso_error,sail_error,pal_error,
               lasso_error_refit,sail_error_refit,pal_error_refit)

colnames(err_rate)=c('Q-Learning (lasso)','pdWOLS','PAL',
                     'Q-Learning (relaxed lasso)','relaxed pdWOLS','relaxed PAL')

rownames(err_rate)=c('Error Rate')


# Value Function
tfree=1-2*Xtest[,1]+2*sin(pi*Xtest[,1])-exp(Xtest[,2])-2*Xtest[,2]

Value_T=mean(tfree+opt*(cbind(1,Xtest[,1])%*%psi))
Value_0=mean(tfree+0*(cbind(1,Xtest[,1])%*%psi))
Value_1=mean(tfree+1*(cbind(1,Xtest[,1])%*%psi))

value_lasso=mean(tfree+lasso_opt*drop((cbind(1,Xtest[,1])%*%psi)))
value_sail=mean(tfree+sail_opt*drop((cbind(1,Xtest[,1])%*%psi)))
value_pal=mean(tfree+pal_opt*drop((cbind(1,Xtest[,1])%*%psi)))

value_lasso_refit=mean(tfree+lasso_opt_refit*drop((cbind(1,Xtest[,1])%*%psi)))
value_sail_refit=mean(tfree+sail_opt_refit*drop((cbind(1,Xtest[,1])%*%psi)))
value_pal_refit=mean(tfree+pal_opt_refit*drop((cbind(1,Xtest[,1])%*%psi)))

value_function=c(value_lasso,value_sail,value_pal,
                 value_lasso_refit,value_sail_refit,value_pal_refit)


save.image("~/Desktop/weights/high.RData")

vr=value_function/Value_T
value_function=round(value_function,2)
res=rbind(fn,fp,err_rate,vr,value_function)
res=t(res)

colnames(res)=c('False Negative','False Positive','Error Rate','VR','Value Function')

res[,1:4]=round(res[,1:4]*100,2)


library(kableExtra)
kable(res,  'latex', booktabs = T,
      caption = "Variable Selection Results")%>%
  kable_styling(latex_options = "HOLD_position",position = 'center')

### plot

psi_0=c(sail[1,],lasso[1,],pal[1,],
        sail_refit[1,],lasso_refit[1,],pal_refit[1,])
psi_1=c(sail[2,],lasso[2,],pal[2,],
        sail_refit[2,],lasso_refit[2,],pal_refit[2,])


m=400
Estimates=c(psi_0,psi_1)
methods=factor(c(rep('pdWOLS',m),rep('QL',m),rep('PAL',m),
                 rep('refitted pdWOLS',m),rep('refitted QL',m),rep('refitted PAL',m)),
               levels=c('pdWOLS','QL','PAL','refitted pdWOLS',
                        'refitted QL','refitted PAL'))

methods=rep(methods,2)

par=rep(c('psi0','psi1'),each=400*6)

ds=data.frame(Estimates,methods,par)
library(latex2exp)
levels(ds$par)=c(TeX('$\\psi_{0}$'),TeX('$\\psi_{1}$'))

true=rep(c(1,-1.5),each=2400)

ind=rep(c('1','2'),each=2400)

ds=cbind(ds,true,ind)

ccols <- c(RColorBrewer::brewer.pal(9, "Blues")[c(4,5,6)],
           RColorBrewer::brewer.pal(9, "RdPu")[c(4,5,6)])


library(ggplot2)
ggplot(ds, aes(x=ind, y=Estimates, fill=methods)) +
  geom_boxplot()+
  cowplot::theme_cowplot()+
  facet_grid(vars(par),scales='free',labeller = label_parsed)+
  geom_hline(aes(yintercept=true),linetype="dashed")+
  scale_fill_manual(values=ccols)+
  xlab('Blip Parameters')+
  ylab('Estimates')


### new plot
methods=factor(c(rep('pdWOLS',m),rep('QL',m),rep('PAL',m),
                 rep('refitted pdWOLS',m),rep('refitted QL',m),rep('refitted PAL',m)),
               levels=c('pdWOLS','QL','PAL','refitted pdWOLS',
                        'refitted QL','refitted PAL'))

ds=data.frame(psi_0,psi_1,methods,par)

library(latex2exp)

ccols <- c(RColorBrewer::brewer.pal(9, "Blues")[c(4,5,6)],
           RColorBrewer::brewer.pal(9, "RdPu")[c(4,5,6)])

library(ggplot2)
### no ggplot
library(latex2exp)
par(mfrow = c(1, 2))
boxplot(sail[1,],lasso[1,],pal[1,],
         sail_refit[1,],lasso_refit[1,],pal_refit[1,],ylab=TeX('$\\psi_{0}$'),
        names=c('pdWOLS','QL','PAL','RpdWOLS',
                'RQL','RPAL'))
abline(h = 1,lty = 2)

boxplot(sail[2,],lasso[2,],pal[2,],
        sail_refit[2,],lasso_refit[2,],pal_refit[2,],ylab=TeX('$\\psi_{1}$'),
        names=c('pdWOLS','QL','PAL','RpdWOLS',
                'RQL','RPAL'))

abline(h = -1.5,lty = 2)





