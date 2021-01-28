##remove n=2000
ds=rbind(pdata1,pdata2)
sample=rep(c(100,500),each=19200)
ds=cbind(ds,sample)
library(latex2exp)

levels(ds$par)=c(TeX('$\\psi_0$'),TeX('$\\psi_1$'))
m=400
true=rep(c(1,-1.5),each=m*24)
true=rep(true,2)
ds=cbind(ds,true)

ds$sample=factor(ds$sample)
levels(ds$sample)=c(TeX('$n=100$'),TeX('$n=500$'))

library(ggplot2)

ccols <- c(RColorBrewer::brewer.pal(9, "Blues")[c(4,5,6)],
           RColorBrewer::brewer.pal(9, "RdPu")[c(4,5,6)])

ggplot(ds, aes(x=scene, y=estimates, fill=methods)) +
  geom_boxplot()+
  facet_wrap(par~sample,scales='free',labeller = label_parsed)+
  geom_hline(aes(yintercept=true),linetype="dashed")+
  scale_fill_manual(values=ccols)+
  xlab('Scenarios')+
  ylab('Estimates')


## remove scenario 1

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

lasso2=do.call(cbind,l2);sail2=do.call(cbind,s2);pal2=do.call(cbind,a2)
lasso3=do.call(cbind,l3);sail3=do.call(cbind,s3);pal3=do.call(cbind,a3)
lasso4=do.call(cbind,l4);sail4=do.call(cbind,s4);pal4=do.call(cbind,a4)

lasso2_refit=do.call(cbind,lr2);sail2_refit=do.call(cbind,sr2);pal2_refit=do.call(cbind,ar2)
lasso3_refit=do.call(cbind,lr3);sail3_refit=do.call(cbind,sr3);pal3_refit=do.call(cbind,ar3)
lasso4_refit=do.call(cbind,lr4);sail4_refit=do.call(cbind,sr4);pal4_refit=do.call(cbind,ar4)

psi_0=c(sail2[1,],sail3[1,],sail4[1,],
        lasso2[1,],lasso3[1,],lasso4[1,],
        pal2[1,],pal3[1,],pal4[1,],
        sail2_refit[1,],sail3_refit[1,],sail4_refit[1,],
        lasso2_refit[1,],lasso3_refit[1,],lasso4_refit[1,],

        pal2_refit[1,],pal3_refit[1,],pal4_refit[1,])

psi_1=c(sail2[2,],sail3[3,],sail4[3,],
        lasso2[2,],lasso3[3,],lasso4[3,],

        pal2[2,],pal3[3,],pal4[3,],
        sail2_refit[2,],sail3_refit[3,],sail4_refit[3,],
        lasso2_refit[2,],lasso3_refit[3,],lasso4_refit[3,],

        pal2_refit[2,],pal3_refit[3,],pal4_refit[3,])


estimates=c(psi_0,psi_1)

m=400
methods=factor(c(rep('pdWOLS',3*m),rep('QL',3*m),rep('PAL',3*m),
                 rep('RpdWOLS',3*m),rep('RQL',3*m),rep('RPAL',3*m)),
               levels=c('pdWOLS','QL','PAL','RpdWOLS',
                        'RQL','RPAL'))

methods=rep(methods,2)

scene=rep(c('2','3','4'),each=m)
scene=rep(scene,2*6)

par=rep(c('psi0','psi1'),each=m*18)

pdata1=data.frame(estimates,methods,scene,par)
pdata2=data.frame(estimates,methods,scene,par)
pdata3=data.frame(estimates,methods,scene,par)

ds=rbind(pdata1,pdata2,pdata3)
sample=rep(c(100,500,2000),each=14400)
ds=cbind(ds,sample)
library(latex2exp)

levels(ds$par)=c(TeX('$\\psi_0$'),TeX('$\\psi_1$'))
m=400
true=rep(c(1,-1.5),each=m*18)
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
  cowplot::theme_cowplot()+
  facet_grid(par~sample,scales='free',labeller = label_parsed)+
  geom_hline(aes(yintercept=true),linetype="dashed")+
  scale_fill_manual(values=ccols)+
  xlab('Scenarios')+
  ylab('Estimates')

### remove 2000 and scenario 1

ds=rbind(pdata1,pdata2)
sample=rep(c(100,500),each=14400)
ds=cbind(ds,sample)
library(latex2exp)

levels(ds$par)=c(TeX('$\\psi_0$'),TeX('$\\psi_1$'))
m=400
true=rep(c(1,-1.5),each=m*18)
true=rep(true,2)
ds=cbind(ds,true)

ds$sample=factor(ds$sample)
levels(ds$sample)=c(TeX('$n=100$'),TeX('$n=500$'))

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





