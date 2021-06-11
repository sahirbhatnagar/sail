## Analyze the results

l1=l2=l3=l4=s1=s2=s3=s4=lr1=lr2=lr3=lr4=sr1=sr2=sr3=sr4=ar1=ar2=ar3=ar4=a1=a2=a3=a4=
  al1=al2=al3=al4=as1=as2=as3=as4=list()

for (i in 1:400) {
  l1[[i]]=res[[i]][[1]]
  lr1[[i]]=res[[i]][[2]]
  al1[[i]]=res[[i]][[3]]
  a1[[i]]=res[[i]][[4]]

  ar1[[i]]=res[[i]][[5]]

  s1[[i]]=res[[i]][[6]]

  sr1[[i]]=res[[i]][[7]]
  as1[[i]]=res[[i]][[8]]

  a2[[i]]=res[[i]][[9]]

  ar2[[i]]=res[[i]][[10]]

  s2[[i]]=res[[i]][[11]]

  sr2[[i]]=res[[i]][[12]]
  as2[[i]]=res[[i]][[13]]

  l3[[i]]=res[[i]][[14]]

  lr3[[i]]=res[[i]][[15]]

  al3[[i]]=res[[i]][[16]]

  a3[[i]]=res[[i]][[17]]

  ar3[[i]]=res[[i]][[18]]

  s3[[i]]=res[[i]][[19]]

  sr3[[i]]=res[[i]][[20]]

  as3[[i]]=res[[i]][[21]]

  a4[[i]]=res[[i]][[22]]

  ar4[[i]]=res[[i]][[23]]

  s4[[i]]=res[[i]][[24]]
  sr4[[i]]=res[[i]][[25]]
  as4[[i]]=res[[i]][[26]]
}


l2=l1; lr2=lr1;al2=al1;l4=l3; lr4=lr3;al4=al3

lasso1=do.call(cbind,l1);sail1=do.call(cbind,s1);pal1=do.call(cbind,a1)
lasso2=do.call(cbind,l2);sail2=do.call(cbind,s2);pal2=do.call(cbind,a2)
lasso3=do.call(cbind,l3);sail3=do.call(cbind,s3);pal3=do.call(cbind,a3)
lasso4=do.call(cbind,l4);sail4=do.call(cbind,s4);pal4=do.call(cbind,a4)
#
# adlasso1=do.call(cbind,al1);adsail1=do.call(cbind,as1)
# adlasso2=do.call(cbind,al2);adsail2=do.call(cbind,as2)
# adlasso3=do.call(cbind,al3);adsail3=do.call(cbind,as3)
# adlasso4=do.call(cbind,al4);adsail4=do.call(cbind,as4)

lasso1_refit=do.call(cbind,lr1);sail1_refit=do.call(cbind,sr1);pal1_refit=do.call(cbind,ar1)
lasso2_refit=do.call(cbind,lr2);sail2_refit=do.call(cbind,sr2);pal2_refit=do.call(cbind,ar2)
lasso3_refit=do.call(cbind,lr3);sail3_refit=do.call(cbind,sr3);pal3_refit=do.call(cbind,ar3)
lasso4_refit=do.call(cbind,lr4);sail4_refit=do.call(cbind,sr4);pal4_refit=do.call(cbind,ar4)

psi_0=c(sail1[1,],sail2[1,],sail3[1,],sail4[1,],
        lasso1[1,],lasso2[1,],lasso3[1,],lasso4[1,],
        pal1[1,],pal2[1,],pal3[1,],pal4[1,],
        sail1_refit[1,],sail2_refit[1,],sail3_refit[1,],sail4_refit[1,],
        lasso1_refit[1,],lasso2_refit[1,],lasso3_refit[1,],lasso4_refit[1,],

        pal1_refit[1,],pal2_refit[1,],pal3_refit[1,],pal4_refit[1,])

psi_1=c(sail1[2,],sail2[2,],sail3[3,],sail4[3,],
        lasso1[2,],lasso2[2,],lasso3[3,],lasso4[3,],

        pal1[2,],pal2[2,],pal3[3,],pal4[3,],
        sail1_refit[2,],sail2_refit[2,],sail3_refit[3,],sail4_refit[3,],
        lasso1_refit[2,],lasso2_refit[2,],lasso3_refit[3,],lasso4_refit[3,],

        pal1_refit[2,],pal2_refit[2,],pal3_refit[3,],pal4_refit[3,])


estimates=c(psi_0,psi_1)

m=400
methods=factor(c(rep('pdWOLS',4*m),rep('QL',4*m),rep('PAL',4*m),
          rep('RpdWOLS',4*m),rep('RQL',4*m),rep('RPAL',4*m)),
          levels=c('pdWOLS','QL','PAL','RpdWOLS',
                   'RQL','RPAL'))

methods=rep(methods,2)

scene=rep(c('1','2','3','4'),each=m)
scene=rep(scene,2*6)

par=rep(c('psi0','psi1'),each=m*24)

pdata1=data.frame(estimates,methods,scene,par)
pdata2=data.frame(estimates,methods,scene,par)
pdata3=data.frame(estimates,methods,scene,par)


# write.csv(pdata1, file = "~/Desktop/weights/p100.csv")
# write.csv(pdata1, file = "~/Desktop/weights/p500.csv")
# write.csv(pdata1, file = "~/Desktop/weights/p2000.csv")

# d1=read.csv('~/Desktop/weights/p100.csv')
# d2=read.csv('~/Desktop/weights/p500.csv')
# d3=read.csv('~/Desktop/weights/p2000.csv')
#
# levels(d1$methods)

ds=rbind(pdata1,pdata2,pdata3)
sample=rep(c(100,500,2000),each=19200)
ds=cbind(ds,sample)
library(latex2exp)

levels(ds$par)=c(TeX('$\\psi_0$'),TeX('$\\psi_1$'))
m=400
true=rep(c(1,-1.5),each=m*24)
true=rep(true,3)
ds=cbind(ds,true)

ds$sample=factor(ds$sample)
levels(ds$sample)=c(TeX('$n=100$'),TeX('$n=500$'),
                    TeX('$n=2000$'))

library(ggplot2)

ccols <- c(RColorBrewer::brewer.pal(9, "Blues")[c(4,5,6)],
           RColorBrewer::brewer.pal(9, "RdPu")[c(4,5,6)])

ggplot(ds, aes(x=scene, y=estimates, fill=methods)) +
  geom_boxplot()+
  cowplot::theme_cowplot() +
  cowplot::panel_border()+
  facet_grid(par~sample,scales='free',labeller = label_parsed)+
  geom_hline(aes(yintercept=true),linetype="dashed")+
  scale_fill_manual(values=ccols)+
  xlab('Scenarios')+
  ylab('Estimates')


## selection results

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

adlasso_sel1=apply(adlasso1!=0, 1, mean);adlasso_sel2=apply(adlasso2!=0, 1, mean)
adlasso_sel3=apply(adlasso3!=0, 1, mean);adlasso_sel4=apply(adlasso4!=0, 1, mean)

adsail_sel1=apply(adsail1!=0, 1, mean);sail_sel2=apply(adsail2!=0, 1, mean)
adsail_sel3=apply(adsail3!=0, 1, mean);sail_sel4=apply(adsail4!=0, 1, mean)

select_rate=cbind(sail_sel1,lasso_sel1,pal_sel1,sail_sel2,lasso_sel2,pal_sel2)

select_rate=rbind(select_rate[1,],c(NA,NA,NA,NA,NA,NA),select_rate[-1,])

select_rate=cbind(select_rate,sail_sel3,lasso_sel3,pal_sel3,sail_sel4,lasso_sel4,pal_sel4)



library(dplyr)
library(kableExtra)
colnames(select_rate)=rep(c('pdWOLS','Q-Learning (lasso)','PAL'),4)
rownames(select_rate)=c('psi_0','Noise','psi_1','Noise','Noise',
                        'Noise','Noise','Noise','Noise','Noise','Noise','Noise')
kable(round(select_rate,2)*100,  'latex', booktabs = T,
      caption = "Variable Selection Rate of the Blip Parameters (Including the Noise Variables)")%>%
  add_header_above(c('','Scenario 1'=3, 'Scenario 2'=3,'Scenario 3'=3, 'Scenario 4'=3))%>%
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


# adsail1_opt=apply(adsail1, 2, function(i) cbind(1,Xtest)%*%i)>0
# adsail2_opt=apply(adsail2, 2, function(i) cbind(1,Xtest)%*%i)>0
# adsail3_opt=apply(adsail3, 2, function(i) cbind(1,exp(Xtest[,1]),Xtest)%*%i)>0
# adsail4_opt=apply(adsail4, 2, function(i) cbind(1,exp(Xtest[,1]),Xtest)%*%i)>0
#
# adlasso1_opt=apply(adlasso1, 2, function(i) cbind(1,Xtest)%*%i)>0
# adlasso2_opt=apply(adlasso2, 2, function(i) cbind(1,Xtest)%*%i)>0
# adlasso3_opt=apply(adlasso3, 2, function(i) cbind(1,exp(Xtest[,1]),Xtest)%*%i)>0
# adlasso4_opt=apply(adlasso4, 2, function(i) cbind(1,exp(Xtest[,1]),Xtest)%*%i)>0

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

lasso_error1_refit=1-mean(colMeans(apply(lasso1_opt_refit, 2, function(i) i==opt)))
lasso_error2_refit=1-mean(colMeans(apply(lasso2_opt_refit, 2, function(i) i==opt)))
lasso_error3_refit=1-mean(colMeans(apply(lasso3_opt_refit, 2, function(i) i==opt)))
lasso_error4_refit=1-mean(colMeans(apply(lasso4_opt_refit, 2, function(i) i==opt)))

sail_error1_refit=1-mean(colMeans(apply(sail1_opt_refit, 2, function(i) i==opt)))
sail_error2_refit=1-mean(colMeans(apply(sail2_opt_refit, 2, function(i) i==opt)))
sail_error3_refit=1-mean(colMeans(apply(sail3_opt_refit, 2, function(i) i==opt)))
sail_error4_refit=1-mean(colMeans(apply(sail4_opt_refit, 2, function(i) i==opt)))

pal_error1_refit=1-mean(colMeans(apply(pal1_opt_refit, 2, function(i) i==opt)))
pal_error2_refit=1-mean(colMeans(apply(pal2_opt_refit, 2, function(i) i==opt)))
pal_error3_refit=1-mean(colMeans(apply(pal3_opt_refit, 2, function(i) i==opt)))
pal_error4_refit=1-mean(colMeans(apply(pal4_opt_refit, 2, function(i) i==opt)))

error=rbind(sail_error1,lasso_error1,pal_error1,sail_error1_refit,lasso_error1_refit,pal_error1_refit,
            sail_error2,lasso_error2,pal_error2,sail_error2_refit,lasso_error2_refit,pal_error2_refit,
            sail_error3,lasso_error3,pal_error3,sail_error3_refit,lasso_error3_refit,pal_error3_refit,
            sail_error4,lasso_error4,pal_error4,sail_error4_refit,lasso_error4_refit,pal_error4_refit)



## Value Function

opt=as.numeric(cbind(1,Xtest[,1])%*%true_psi>=0)
tfree=.5-2*Xtest[,1]-0.6*exp(Xtest[,1])-2*Xtest[,2]

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

v_lasso1_refit=mean(tfree+lasso1_opt_refit*drop((cbind(1,Xtest[,1])%*%true_psi)))
v_lasso2_refit=mean(tfree+lasso2_opt_refit*drop((cbind(1,Xtest[,1])%*%true_psi)))
v_lasso3_refit=mean(tfree+lasso3_opt_refit*drop((cbind(1,Xtest[,1])%*%true_psi)))
v_lasso4_refit=mean(tfree+lasso4_opt_refit*drop((cbind(1,Xtest[,1])%*%true_psi)))

v_sail1_refit=mean(tfree+sail1_opt_refit*drop((cbind(1,Xtest[,1])%*%true_psi)))
v_sail2_refit=mean(tfree+sail2_opt_refit*drop((cbind(1,Xtest[,1])%*%true_psi)))
v_sail3_refit=mean(tfree+sail3_opt_refit*drop((cbind(1,Xtest[,1])%*%true_psi)))
v_sail4_refit=mean(tfree+sail4_opt_refit*drop((cbind(1,Xtest[,1])%*%true_psi)))

v_pal1_refit=mean(tfree+pal1_opt_refit*drop((cbind(1,Xtest[,1])%*%true_psi)))
v_pal2_refit=mean(tfree+pal2_opt_refit*drop((cbind(1,Xtest[,1])%*%true_psi)))
v_pal3_refit=mean(tfree+pal3_opt_refit*drop((cbind(1,Xtest[,1])%*%true_psi)))
v_pal4_refit=mean(tfree+pal4_opt_refit*drop((cbind(1,Xtest[,1])%*%true_psi)))

value=c(v_sail1,v_lasso1,v_pal1,v_sail1_refit,v_lasso1_refit,v_pal1_refit,
        v_sail2,v_lasso2,v_pal2,v_sail2_refit,v_lasso2_refit,v_pal2_refit,
        v_sail3,v_lasso3,v_pal3,v_sail3_refit,v_lasso3_refit,v_pal3_refit,
        v_sail4,v_lasso4,v_pal4,v_sail4_refit,v_lasso4_refit,v_pal4_refit)

res=cbind(error,value/v_true,value)

res[,1:2]=round(res[,1:2]*100,2)
res[,3]=round(res[,3],3)


scene=rep(c('1','2','3','4'),each=6)

methods=rep(c('pdWOLS','Q-L (lasso)','PAL',
              'refitted pdWOLS','Q-L (refitted lasso)','refitted PAL'),4)

res=cbind(methods,res)

rownames(res)=scene


library(kableExtra)

kable(res,  'latex', booktabs = T,align = "c",
      caption = "Error Rate and Value Function of the Estimated Rules")%>%
  column_spec(1:2) %>%
  collapse_rows(columns = 1, latex_hline = "major", valign = "middle")
  kable_styling(latex_options = "HOLD_position",position = 'center')



