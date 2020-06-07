load("~/Desktop/weights/dwols.RData")
load("~/Desktop/weights/multiple.RData")

psi_01=c(lasso1[1,],lasso1_refit[1,],sail1[1,],sail1_refit[1,],pal1[1,],pal1_refit[1,])
psi_11=c(lasso1[2,],lasso1_refit[2,],sail1[2,],sail1_refit[2,],pal1[2,],pal1_refit[2,])

psi_02=c(lasso2[1,],lasso2_refit[1,],sail2[1,],sail2_refit[1,],pal2[1,],pal2_refit[1,])
psi_12=c(lasso2[2,],lasso2_refit[2,],sail2[2,],sail2_refit[2,],pal2[2,],pal2_refit[2,])

Estimates=c(psi_01,psi_11,psi_02,psi_12)
methods=c(rep('lasso',100),rep('relaxed lasso',100),rep('pdWOLS',100),
          rep('relaxed pdWOLS',100),rep('PAL',100),rep('relaxed PAL',100))
methods=rep(methods,4)

stage=rep(c('Stage 1','Stage 2'),each=1200)

par=rep(c('psi01','psi11','psi02','psi12'),each=600)

ds=data.frame(Estimates,methods,stage,par)
library(latex2exp)
levels(ds$par)=c(TeX('$\\psi_{01}$'),TeX('$\\psi_{02}$'),
                 TeX('$\\psi_{11}$'),TeX('$\\psi_{12}$'))

true=rep(c(0.8,-2,1,-1.5),each=600)
# true=rep(c(-1.46,1.5,-.9,1),each=600)

ds=cbind(ds,true)

library(ggplot2)
ggplot(ds, aes(x=stage,y=Estimates, fill=methods)) +
  geom_boxplot()+
  facet_grid(vars(par),scales='free',labeller = label_parsed)+
  geom_hline(aes(yintercept=true),linetype="dashed")+
  xlab('Stage')+
  ylab('Estimates')

library(kableExtra)
library(dplyr)
load("~/Desktop/weights/dwols.RData")
load("~/Desktop/weights/multiple.RData")

sel_rate=sel_rate[,-c(8,10,12)]

library(kableExtra)
kable(sel_rate*100,  'latex', booktabs = T,
      caption = "Variable Selection Rate (%)")%>%
  add_header_above(c('','Stage 1'=6, 'Stage 2'=3))%>%
  kable_styling(latex_options = "HOLD_position",position = 'center')

err_rate=cbind(err_rate[1:6],err_rate[7:12])
rownames(err_rate)=c('lasso', 'relaxed lasso',
'sail','relaxed sail', 'PAL',
'relaxed PAL')

colnames(err_rate)=c('Error Rate (Stage 1)','Error Rate (Stage 2)')

res=cbind(err_rate,value_function)

kable(res,  'latex', booktabs = T,
      caption = "Value Ratio and Error Rate (%)")%>%
  kable_styling(latex_options = "HOLD_position",position = 'center')




