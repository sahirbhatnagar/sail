## Analyze the results

load("~/Desktop/weights/2000.RData")

l2=l1; lr2=lr1;al2=al1;l4=l3; lr4=lr3;al4=al3
lasso1=do.call(cbind,l1);sail1=do.call(cbind,s1);pal1=do.call(cbind,a1)
lasso2=do.call(cbind,l2);sail2=do.call(cbind,s2);pal2=do.call(cbind,a2)
lasso3=do.call(cbind,l3);sail3=do.call(cbind,s3);pal3=do.call(cbind,a3)
lasso4=do.call(cbind,l4);sail4=do.call(cbind,s4);pal4=do.call(cbind,a4)

adlasso1=do.call(cbind,al1);adsail1=do.call(cbind,as1)
adlasso2=do.call(cbind,al2);adsail2=do.call(cbind,as2)
adlasso3=do.call(cbind,al3);adsail3=do.call(cbind,as3)
adlasso4=do.call(cbind,al4);adsail4=do.call(cbind,as4)

lasso1_refit=do.call(cbind,lr1);sail1_refit=do.call(cbind,sr1);pal1_refit=do.call(cbind,ar1)
lasso2_refit=do.call(cbind,lr2);sail2_refit=do.call(cbind,sr2);pal2_refit=do.call(cbind,ar2)
lasso3_refit=do.call(cbind,lr3);sail3_refit=do.call(cbind,sr3);pal3_refit=do.call(cbind,ar3)
lasso4_refit=do.call(cbind,lr4);sail4_refit=do.call(cbind,sr4);pal4_refit=do.call(cbind,ar4)

cbind(rowMeans(lasso2),rowMeans(sail2),rowMeans(pal2),rowMeans(pal2_refit),rowMeans(adlasso2),rowMeans(adsail2))

cbind(rowMeans(lasso1),rowMeans(sail1),rowMeans(pal1),rowMeans(pal1_refit),rowMeans(adlasso1),rowMeans(adsail1))

cbind(rowMeans(lasso3),rowMeans(sail3),rowMeans(pal3),rowMeans(pal3_refit),rowMeans(adlasso3),rowMeans(adsail3))

cbind(rowMeans(lasso4),rowMeans(sail4),rowMeans(pal4),rowMeans(pal4_refit),rowMeans(adlasso4),rowMeans(adsail4))


psi_0=c(lasso1[1,],lasso2[1,],lasso3[1,],lasso4[1,],
        sail1[1,],sail2[1,],sail3[1,],sail4[1,],
        pal1[1,],pal2[1,],pal3[1,],pal4[1,],
        lasso1_refit[1,],lasso2_refit[1,],lasso3_refit[1,],lasso4_refit[1,],
        sail1_refit[1,],sail2_refit[1,],sail3_refit[1,],sail4_refit[1,],
        pal1_refit[1,],pal2_refit[1,],pal3_refit[1,],pal4_refit[1,])

psi_1=c(lasso1[2,],lasso2[2,],lasso3[3,],lasso4[3,],
        sail1[2,],sail2[2,],sail3[3,],sail4[3,],
        pal1[2,],pal2[2,],pal3[3,],pal4[3,],
        lasso1_refit[2,],lasso2_refit[2,],lasso3_refit[3,],lasso4_refit[3,],
        sail1_refit[2,],sail2_refit[2,],sail3_refit[3,],sail4_refit[3,],
        pal1_refit[2,],pal2_refit[2,],pal3_refit[3,],pal4_refit[3,])

psi_2=c(lasso1[3,],lasso2[3,],lasso3[4,],lasso4[4,],
        sail1[3,],sail2[3,],sail3[4,],sail4[4,],
        pal1[3,],pal2[3,],pal3[4,],pal4[4,],
        lasso1_refit[3,],lasso2_refit[3,],lasso3_refit[4,],lasso4_refit[4,],
        sail1_refit[3,],sail2_refit[3,],sail3_refit[4,],sail4_refit[4,],
        pal1_refit[3,],pal2_refit[3,],pal3_refit[4,],pal4_refit[4,])

psi_3=c(lasso1[4,],lasso2[4,],lasso3[5,],lasso4[5,],
        sail1[4,],sail2[4,],sail3[5,],sail4[5,],
        pal1[4,],pal2[4,],pal3[5,],pal4[5,],
        lasso1_refit[4,],lasso2_refit[4,],lasso3_refit[5,],lasso4_refit[5,],
        sail1_refit[4,],sail2_refit[4,],sail3_refit[5,],sail4_refit[5,],
        pal1_refit[4,],pal2_refit[4,],pal3_refit[5,],pal4_refit[5,])

estimates=c(psi_0,psi_1,psi_2,psi_3)

m=200
methods=c(rep('lasso',4*m),rep('pdWOLS',4*m),rep('PAL',4*m),
          rep('relaxed lasso',4*m),rep('relaxed pdWOLS',4*m),rep('relaxed PAL',4*m))
methods=rep(methods,4)


scene=rep(c('Scenario 1','Scenario 2','Scenario 3','Scenario 4'),each=m)
scene=rep(scene,4*6)


par=rep(c('psi0','psi1','psi2','psi3'),each=m*24)

pdata1=data.frame(estimates,methods,scene,par)

# write.csv(pdata1, file = "~/Desktop/weights/p200.csv")
# write.csv(pdata1, file = "~/Desktop/weights/p800.csv")
# write.csv(pdata1, file = "~/Desktop/weights/p2000.csv")

d1=read.csv('~/Desktop/weights/p200.csv')
d2=read.csv('~/Desktop/weights/p800.csv')
d3=read.csv('~/Desktop/weights/p2000.csv')

ds=rbind(d1,d2,d3)
sample=rep(c(200,800,2000),each=19200)
ds=cbind(ds,sample)
ds=ds[,-1]
library(latex2exp)

levels(ds$par)=c(TeX('$\\psi_0$'),TeX('$\\psi_1$'),
                 TeX('$\\psi_2$'),TeX('$\\psi_3$'))

true=rep(c(1,-1.5,-1.5,1),each=m*24)
true=rep(true,3)
ds=cbind(ds,true)

ds$sample=factor(ds$sample)
levels(ds$sample)=c(TeX('$n=200$'),TeX('$n=800$'),
                    TeX('$n=2000$'))


library(ggplot2)
ggplot(ds, aes(x=scene, y=estimates, fill=methods)) +
  geom_boxplot()+
  facet_grid(par~sample,scales='free',labeller = label_parsed)+
  geom_hline(aes(yintercept=true),linetype="dashed")+
  xlab('Scenarios')+
  ylab('Estimates')


## selection results
load("~/Desktop/weights/2000.RData")

lasso1=do.call(cbind,l1);sail1=do.call(cbind,s1);pal1=do.call(cbind,a1)
lasso2=do.call(cbind,l2);sail2=do.call(cbind,s2);pal2=do.call(cbind,a2)
lasso3=do.call(cbind,l3);sail3=do.call(cbind,s3);pal3=do.call(cbind,a3)
lasso4=do.call(cbind,l4);sail4=do.call(cbind,s4);pal4=do.call(cbind,a4)


lasso1_refit=do.call(cbind,lr1);sail1_refit=do.call(cbind,sr1);pal1_refit=do.call(cbind,ar1)
lasso2_refit=do.call(cbind,lr2);sail2_refit=do.call(cbind,sr2);pal2_refit=do.call(cbind,ar2)
lasso3_refit=do.call(cbind,lr3);sail3_refit=do.call(cbind,sr3);pal3_refit=do.call(cbind,ar3)
lasso4_refit=do.call(cbind,lr4);sail4_refit=do.call(cbind,sr4);pal4_refit=do.call(cbind,ar4)

### Selection Rate
lasso_sel1=apply(lasso1!=0, 1, mean);lasso_sel2=apply(lasso2!=0, 1, mean)
lasso_sel3=apply(lasso3!=0, 1, mean);lasso_sel4=apply(lasso4!=0, 1, mean)
lasso_sel1_refit=apply(lasso1_refit!=0, 1, mean);lasso_sel2_refit=apply(lasso2_refit!=0, 1, mean)
lasso_sel3_refit=apply(lasso3_refit!=0, 1, mean);lasso_sel4_refit=apply(lasso4_refit!=0, 1, mean)

sail_sel1=apply(sail1!=0, 1, mean);sail_sel2=apply(sail2!=0, 1, mean)
sail_sel3=apply(sail3!=0, 1, mean);sail_sel4=apply(sail4!=0, 1, mean)
sail_sel1_refit=apply(sail1_refit!=0, 1, mean);sail_sel2_refit=apply(sail2_refit!=0, 1, mean)
sail_sel3_refit=apply(sail3_refit!=0, 1, mean);sail_sel4_refit=apply(sail4_refit!=0, 1, mean)

pal_sel1=apply(pal1!=0, 1, mean);pal_sel2=apply(pal2!=0, 1, mean)
pal_sel3=apply(pal3!=0, 1, mean);pal_sel4=apply(pal4!=0, 1, mean)
pal_sel1_refit=apply(pal1_refit!=0, 1, mean);pal_sel2_refit=apply(pal2_refit!=0, 1, mean)
pal_sel3_refit=apply(pal3_refit!=0, 1, mean);pal_sel4_refit=apply(pal4_refit!=0, 1, mean)

select_rate=cbind(lasso_sel1,pal_sel1,sail_sel1,lasso_sel2,pal_sel2,sail_sel2)

select_rate=rbind(select_rate[1,],c(NA,NA,NA,NA,NA,NA),select_rate[-1,])

select_rate=cbind(select_rate,lasso_sel3,pal_sel3,sail_sel3,lasso_sel4,pal_sel4,sail_sel4)

library(dplyr)
library(kableExtra)
colnames(select_rate)=rep(c('lasso','sail'),4)
rownames(select_rate)=c('psi_0','Noise','psi_1','psi_2','psi_3',
                        'Noise','Noise','Noise','Noise','Noise','Noise','Noise')
kable(select_rate,  'latex', booktabs = T,
      caption = "Variable Selection Rate of the Blip Parameters (Including the Noise Variables)")%>%
  add_header_above(c('','Scenario 1'=2, 'Scenario 2'=2,'Scenario 3'=2, 'Scenario 4'=2))%>%
  kable_styling(latex_options = "HOLD_position",position = 'center')


## Error Rate

set.seed(999);ds_test=g_data(10000,10);Xtest=ds_test[[1]];rm(ds_test);true_psi=c(1,-1.5)
opt=as.numeric(cbind(1,Xtest[,1])%*%true_psi>=0)
table(opt)

lasso1_opt=apply(lasso1, 2, function(i) cbind(1,Xtest)%*%i)>0
lasso2_opt=apply(lasso2, 2, function(i) cbind(1,Xtest)%*%i)>0
lasso3_opt=apply(lasso3, 2, function(i) cbind(1,exp(Xtest[,1]),Xtest)%*%i)>0
lasso4_opt=apply(lasso4, 2, function(i) cbind(1,exp(Xtest[,1]),Xtest)%*%i)>0

lasso1_opt_refit=apply(lasso1_refit, 2, function(i) cbind(1,Xtest)%*%i)>0
lasso2_opt_refit=apply(lasso2_refit, 2, function(i) cbind(1,Xtest)%*%i)>0
lasso3_opt_refit=apply(lasso3_refit, 2, function(i) cbind(1,exp(Xtest[,1]),Xtest)%*%i)>0
lasso4_opt_refit=apply(lasso4_refit, 2, function(i) cbind(1,exp(Xtest[,1]),Xtest)%*%i)>0

sail1_opt=apply(sail1, 2, function(i) cbind(1,Xtest)%*%i)>0
sail2_opt=apply(sail2, 2, function(i) cbind(1,Xtest)%*%i)>0
sail3_opt=apply(sail3, 2, function(i) cbind(1,exp(Xtest[,1]),Xtest)%*%i)>0
sail4_opt=apply(sail4, 2, function(i) cbind(1,exp(Xtest[,1]),Xtest)%*%i)>0

sail1_opt_refit=apply(sail1_refit, 2, function(i) cbind(1,Xtest)%*%i)>0
sail2_opt_refit=apply(sail2_refit, 2, function(i) cbind(1,Xtest)%*%i)>0
sail3_opt_refit=apply(sail3_refit, 2, function(i) cbind(1,exp(Xtest[,1]),Xtest)%*%i)>0
sail4_opt_refit=apply(sail4_refit, 2, function(i) cbind(1,exp(Xtest[,1]),Xtest)%*%i)>0

pal1_opt=apply(pal1, 2, function(i) cbind(1,Xtest)%*%i)>0
pal2_opt=apply(pal2, 2, function(i) cbind(1,Xtest)%*%i)>0
pal3_opt=apply(pal3, 2, function(i) cbind(1,exp(Xtest[,1]),Xtest)%*%i)>0
pal4_opt=apply(pal4, 2, function(i) cbind(1,exp(Xtest[,1]),Xtest)%*%i)>0

pal1_opt_refit=apply(pal1_refit, 2, function(i) cbind(1,Xtest)%*%i)>0
pal2_opt_refit=apply(pal2_refit, 2, function(i) cbind(1,Xtest)%*%i)>0
pal3_opt_refit=apply(pal3_refit, 2, function(i) cbind(1,exp(Xtest[,1]),Xtest)%*%i)>0
pal4_opt_refit=apply(pal4_refit, 2, function(i) cbind(1,exp(Xtest[,1]),Xtest)%*%i)>0


adsail1_opt=apply(adsail1, 2, function(i) cbind(1,Xtest)%*%i)>0
adsail2_opt=apply(adsail2, 2, function(i) cbind(1,Xtest)%*%i)>0
adsail3_opt=apply(adsail3, 2, function(i) cbind(1,exp(Xtest[,1]),Xtest)%*%i)>0
adsail4_opt=apply(adsail4, 2, function(i) cbind(1,exp(Xtest[,1]),Xtest)%*%i)>0

adlasso1_opt=apply(adlasso1, 2, function(i) cbind(1,Xtest)%*%i)>0
adlasso2_opt=apply(adlasso2, 2, function(i) cbind(1,Xtest)%*%i)>0
adlasso3_opt=apply(adlasso3, 2, function(i) cbind(1,exp(Xtest[,1]),Xtest)%*%i)>0
adlasso4_opt=apply(adlasso4, 2, function(i) cbind(1,exp(Xtest[,1]),Xtest)%*%i)>0

## error

lasso_error1=1-mean(colMeans(apply(lasso1_opt, 2, function(i) i==opt)))
lasso_error2=1-mean(colMeans(apply(lasso2_opt, 2, function(i) i==opt)))
lasso_error3=1-mean(colMeans(apply(lasso3_opt, 2, function(i) i==opt)))
lasso_error4=1-mean(colMeans(apply(lasso4_opt, 2, function(i) i==opt)))

sail_error1=1-mean(colMeans(apply(sail1_opt, 2, function(i) i==opt)))
sail_error2=1-mean(colMeans(apply(sail2_opt, 2, function(i) i==opt)))
sail_error3=1-mean(colMeans(apply(sail3_opt, 2, function(i) i==opt)))
sail_error4=1-mean(colMeans(apply(sail4_opt, 2, function(i) i==opt)))

pal_error1=1-mean(colMeans(apply(pal1_opt, 2, function(i) i==opt)))
pal_error2=1-mean(colMeans(apply(pal2_opt, 2, function(i) i==opt)))
pal_error3=1-mean(colMeans(apply(pal3_opt, 2, function(i) i==opt)))
pal_error4=1-mean(colMeans(apply(pal4_opt, 2, function(i) i==opt)))

## adaptive

## refit one

lasso_error1_refit=1-mean(colMeans(apply(lasso_opt1_refit, 2, function(i) i==opt1)))
lasso_error2_refit=1-mean(colMeans(apply(lasso_opt2_refit, 2, function(i) i==opt2)))

sail_error1_refit=1-mean(colMeans(apply(sail_opt1_refit, 2, function(i) i==opt1)))
sail_error2_refit=1-mean(colMeans(apply(sail_opt2_refit, 2, function(i) i==opt2)))

pal_error1_refit=1-mean(colMeans(apply(pal_opt1_refit, 2, function(i) i==opt1)))
pal_error2_refit=1-mean(colMeans(apply(pal_opt2_refit, 2, function(i) i==opt2)))


## Value Function

set.seed(999);ds_test=g_data(10000,10);Xtest=ds_test[[1]];rm(ds_test);true_psi=c(1,-1.5)
opt=as.numeric(cbind(1,Xtest[,1])%*%true_psi>=0)
table(opt)
tfree=1-2*(exp(Xtest[,1])+log(abs(Xtest[,1])))+2*Xtest[,2]

v_true=mean(tfree+opt*(cbind(1,Xtest[,1])%*%true_psi))
v_0=mean(tfree+0*(cbind(1,Xtest[,1])%*%true_psi))
v_1=mean(tfree+(cbind(1,Xtest[,1])%*%true_psi))

v_lasso1=mean(tfree+lasso1_opt*drop((cbind(1,Xtest[,1])%*%true_psi)))
v_lasso2=mean(tfree+lasso2_opt*drop((cbind(1,Xtest[,1])%*%true_psi)))
v_lasso3=mean(tfree+lasso3_opt*drop((cbind(1,Xtest[,1])%*%true_psi)))
v_lasso4=mean(tfree+lasso4_opt*drop((cbind(1,Xtest[,1])%*%true_psi)))

v_sail1=mean(tfree+sail1_opt*drop((cbind(1,Xtest[,1])%*%true_psi)))
v_sail2=mean(tfree+sail2_opt*drop((cbind(1,Xtest[,1])%*%true_psi)))
v_sail3=mean(tfree+sail3_opt*drop((cbind(1,Xtest[,1])%*%true_psi)))
v_sail4=mean(tfree+sail4_opt*drop((cbind(1,Xtest[,1])%*%true_psi)))

v_pal1=mean(tfree+pal1_opt*drop((cbind(1,Xtest[,1])%*%true_psi)))
v_pal2=mean(tfree+pal2_opt*drop((cbind(1,Xtest[,1])%*%true_psi)))
v_pal3=mean(tfree+pal3_opt*drop((cbind(1,Xtest[,1])%*%true_psi)))
v_pal4=mean(tfree+pal4_opt*drop((cbind(1,Xtest[,1])%*%true_psi)))


scene=c('Scenerio 1','Scenerio 1','Scenerio 2','Scenerio 2',
        'Scenerio 3','Scenerio 3','Scenerio 4','Scenerio 4')
method=rep(c('lasso','pdWOLS'),4)
ds=data.frame(mean_value,scene,method)

ggplot(ds, aes(fill=method, y=mean_value, x=scene)) +
  geom_bar(position="dodge", stat="identity")+
  ylab('Value')+
  xlab('Scenarios')


library(kableExtra)
ds=rbind(mean_error,mean_value)
rownames(ds)=c('Error rate','Value')
colnames(ds)=rep(c('lasso','pdWOLS'),4)
kable(round(ds,2),  'latex', booktabs = T,
      caption = "Error Rate and Value Function of the Estimated Rules")%>%
  add_footnote(c('The value function under the true rule is 9.37'),notation = 'symbol')%>%
  add_header_above(c('','Scenario 1'=2, 'Scenario 2'=2,'Scenario 3'=2, 'Scenario 4'=2))%>%
  kable_styling(latex_options = "HOLD_position",position = 'center')




