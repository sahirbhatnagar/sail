load("~/Desktop/weights/dwols.RData")
load("~/Desktop/weights/ql.RData")

psi_01=c(sail1[1,],lasso1[1,],pal1[1,],
         sail1_refit[1,],lasso1_refit[1,],pal1_refit[1,])

psi_11=c(sail1[2,],lasso1[2,],pal1[2,],
         sail1_refit[2,],lasso1_refit[2,],pal1_refit[2,])

psi_02=c(sail2[1,],lasso2[1,],pal2[1,],
         sail2_refit[1,],lasso2_refit[1,],pal2_refit[1,])
psi_12=c(sail2[2,],lasso2[2,],pal2[2,],
         sail2_refit[2,],lasso2_refit[2,],pal2_refit[2,])

Estimates=c(psi_01,psi_11,psi_02,psi_12)
m=400
methods=factor(c(rep('pdWOLS',m),rep('QL',m),rep('PAL',m),
                 rep('refitted pdWOLS',m),rep('refitted QL',m),rep('refitted PAL',m)),
               levels=c('pdWOLS','QL','PAL','refitted pdWOLS',
                        'refitted QL','refitted PAL'))
methods=rep(methods,4)

stage=rep(c('Stage 1','Stage 2'),each=12*m)

par=rep(c('psi01','psi11','psi02','psi12'),each=6*m)

ds=data.frame(Estimates,methods,stage,par)
library(latex2exp)
levels(ds$par)=c(TeX('$\\psi_{01}$'),TeX('$\\psi_{02}$'),
                 TeX('$\\psi_{11}$'),TeX('$\\psi_{12}$'))

# true=rep(c(0.8,-2,1,-1.5),each=6*m)
# true=rep(c(-1.190399,1.5,1,-1.5),each=6*m)

ds=cbind(ds,true)

library(ggplot2)

ccols <- c(RColorBrewer::brewer.pal(9, "Blues")[c(4,5,6)],
           RColorBrewer::brewer.pal(9, "RdPu")[c(4,5,6)])

ggplot(ds, aes(x=stage,y=Estimates, fill=methods)) +
  geom_boxplot()+
  facet_grid(vars(par),scales='free',labeller = label_parsed)+
  geom_hline(aes(yintercept=true),linetype="dashed")+
  scale_fill_manual(values=ccols)+
  xlab('Stage')+
  ylab('Estimates')

## new plot
library(latex2exp)
par(mfrow = c(2, 2))
boxplot(sail1[1,],lasso1[1,],pal1[1,],
        sail1_refit[1,],lasso1_refit[1,],pal1_refit[1,],ylab=TeX('$\\psi_{01}$'),cex.axis=1.1,
        names=c('pdWOLS','QL','PAL','RpdWOLS',
                'RQL','RPAL'))
abline(h = -1.195,lty = 2)

boxplot(sail1[2,],lasso1[2,],pal1[2,],
        sail1_refit[2,],lasso1_refit[2,],pal1_refit[2,],ylab=TeX('$\\psi_{11}$'),cex.axis=1.1,
        names=c('pdWOLS','QL','PAL','RpdWOLS',
                'RQL','RPAL'))
abline(h = 1.5,lty = 2)

boxplot(sail2[1,],lasso2[1,],pal2[1,],
        sail2_refit[1,],lasso2_refit[1,],pal2_refit[1,],ylab=TeX('$\\psi_{02}$'),cex.axis=1.1,
        names=c('pdWOLS','QL','PAL','RpdWOLS',
                'RQL','RPAL'))
abline(h = 1,lty = 2)

boxplot(sail2[2,],lasso2[2,],pal2[2,],
        sail2_refit[2,],lasso2_refit[2,],pal2_refit[2,],ylab=TeX('$\\psi_{12}$'),cex.axis=1.1,
        names=c('pdWOLS','QL','PAL','RpdWOLS',
                'RQL','RPAL'))
abline(h = -1.5,lty = 2)




library(kableExtra)
library(dplyr)


sel_rate=sel_rate[,1:9]

kable(sel_rate,  'latex', booktabs = T,
      caption = "Variable Selection Rate (%)")%>%
  add_header_above(c('','Stage 1'=6, 'Stage 2'=3))%>%
  kable_styling(latex_options = "HOLD_position",position = 'center')

err_rate=cbind(err_rate[1:6],err_rate[7:12])
rownames(err_rate)=c('pdWOLS','Q-Learning (lasso)','PAL',
                     'refitted pdWOLS','Q-Learning (refitted lasso)',
                     'refitted PAL')

colnames(err_rate)=c('ER (Stage 1)','ER (Stage 2)')
err_rate=t(err_rate)

res=rbind(err_rate,value_function)
res=rbind(res[4,],res[1:3,])
rownames(res)[1]='Total Error'
res[1:3,]=round(res[1:3,],2)
res=t(res)

kable(res,  'latex', booktabs = T,
      caption = "Value Ratio and Error Rate (%)")%>%
  kable_styling(latex_options = "HOLD_position",position = 'center')


err_rate

