library(splines)
library(glmnet)
library(xtable)
library(dplyr)
devtools::load_all()
## pdwols
Theta = seq(-2, 2, by=0.2)

knots = seq(-1.8, 1.8, by=0.2)

method = "adaptive lasso" # or "lasso"
criteria = "ERIC"  # or "BIC"
refit = TRUE;n=200;p=10
cn=log(log(n))
kn=cn*n^(1/3)*(log(p))^(2/3)
nsim = 500
nu=0.05
value = matrix(0, nsim, p)
num_knot = matrix(0, nsim, p)
useN = TRUE
pfac=c(0,rep(1,20))
vic1=bic1=vic2=bic2=confounder=list()
rank1=rank2=list()
load("~/Desktop/weights/book_chapter/sim2.RData")
load("pdwols2.RData")
# save.image("pdwols2.RData")
## pdwols sim
for (i in 401:500) {
  Y=ds[[i]][[2]];X=ds[[i]][[1]];A=ds[[i]][[3]]
  m=lm(Y~cbind(X,A*X))
  beta=coef(m)[2:11]
  w=abs(1/beta)^3
  m=glmnet(X,A,'binomial',penalty.factor = w)
  l=length(m$lambda)
  criterion=c()
  for (j in 1:l) {
    pi=exp(cbind(1,X)%*%coef(m)[,j])/(1+exp(cbind(1,X)%*%coef(m)[,j]))
    tau=A/(pi)+(1-A)/(1-pi)
    tauX=sweep(X,1,tau,'*')
    criterion[j]= crossprod(abs(beta),
                            abs(colSums(tauX*A)/sum(tau*A)-colSums(tauX*(1-A))/sum(tau*(1-A))))
  }

  alpha=coef(m)[,which.min(criterion)]
  confounder[[i]]=alpha[-1]
  ind=which(alpha[-1]!=0)
  m=glm(A~X[,ind],binomial)
  ps=fitted(m);w=abs(A-ps)

  ## pdwols
  m=sail(y=Y,e=A,x=X,penalty.factor=pfac,weights = w,
         basis=function(o) o)
  ## rank
  D=as.matrix(coef(m)[13:22,])
  rank1[[i]]=rowMeans(D!=0)

  ## vic
  l=length(m$lambda);v=c()
  for (j in 1:l) {
    theta=coef(m)[,j]
    psi=theta[12:22]
    num=sum(psi[-1]!=0)
    aopt=cbind(A,X)%*%psi>0
    v[j]=sum((A*aopt+(1-A)*(1-aopt))*Y/((A*ps)+(1-A)*(1-ps)))-kn*num
  }
  psi=coef(m)[,which.max(v)]
  ind=which(psi[-1]!=0)
  X=cbind(X,A,A*X)
  linear=lm(Y~X[,ind],weights = w)
  psi[ind+1]=coef(linear)[-1]
  vic1[[i]]=psi[12:22]
  ## bic
  df=m$dfbeta+m$dfalpha+1
  bic=n* log(m$rss/n)+df*log(n)
  psi=coef(m)[,which.min(bic)]
  ind=which(psi[-1]!=0)
  X=cbind(X,A,A*X)
  linear=lm(Y~X[,ind],weights = w)
  psi[ind+1]=coef(linear)[-1]
  bic1[[i]]=psi[12:22]
}
### without confounder selection
for (i in 401:500) {
  Y=ds[[i]][[2]];X=ds[[i]][[1]];A=ds[[i]][[3]]
  m=glm(A~X,binomial)
  ps=fitted(m);w=abs(A-ps)
  ## pdwols
  m=sail(y=Y,e=A,x=X,penalty.factor=pfac,weights = w,
         basis=function(o) o)
  ## rank
  D=as.matrix(coef(m)[13:22,])
  rank2[[i]]=rowMeans(D!=0)
  ## vic
  l=length(m$lambda);v=c()
  for (j in 1:l) {
    theta=coef(m)[,j]
    psi=theta[12:22]
    num=sum(psi[-1]!=0)
    aopt=cbind(A,X)%*%psi>0
    v[j]=sum((A*aopt+(1-A)*(1-aopt))*Y/((A*ps)+(1-A)*(1-ps)))-kn*num
  }
  psi=coef(m)[,which.max(v)]
  ind=which(psi[-1]!=0)
  X=cbind(X,A,A*X)
  linear=lm(Y~X[,ind],weights = w)
  psi[ind+1]=coef(linear)[-1]
  vic2[[i]]=psi[12:22]
  ## bic
  df=m$dfbeta+m$dfalpha+1
  bic=n* log(m$rss/n)+df*log(n)
  psi=coef(m)[,which.min(bic)]
  ind=which(psi[-1]!=0)
  X=cbind(X,A,A*X)
  linear=lm(Y~X[,ind],weights = w)
  psi[ind+1]=coef(linear)[-1]
  bic2[[i]]=psi[12:22]
}
## ranking method

for (i in 16:nsim) {
  Y=ds[[i]][[2]];X=ds[[i]][[1]];A=ds[[i]][[3]]
  dat=data.frame(id=1:n, y=Y, A, X)
  res1 = rank_binary(X[,1],'opposite',i)

  res2 = estimate.regime(dat, tailoring.name="X2", Theta, knots,
                         method = method, criteria = criteria, refit=refit, nu=nu, useN=useN)
  res3 = estimate.regime(dat, tailoring.name="X3", Theta, knots,
                         method = method, criteria = criteria, refit=refit, nu=nu, useN=useN)
  res4 = estimate.regime(dat, tailoring.name="X4", Theta, knots,
                         method = method, criteria = criteria, refit=refit, nu=nu, useN=useN)
  res5 = estimate.regime(dat, tailoring.name="X5", Theta, knots,
                         method = method, criteria = criteria, refit=refit, nu=nu, useN=useN)
  res6 = estimate.regime(dat, tailoring.name="X6", Theta, knots,
                         method = method, criteria = criteria, refit=refit, nu=nu, useN=useN)
  res7 = estimate.regime(dat, tailoring.name="X7", Theta, knots,
                         method = method, criteria = criteria, refit=refit, nu=nu, useN=useN)
  res8 = estimate.regime(dat, tailoring.name="X8", Theta, knots,
                         method = method, criteria = criteria, refit=refit, nu=nu, useN=useN)
  res9 = estimate.regime(dat, tailoring.name="X9", Theta, knots,
                         method = method, criteria = criteria, refit=refit, nu=nu, useN=useN)
  res10 = rank_binary(X[,10],'opposite',i)

  #record the maximum value function of each tailoring variabels, number of knots selected
  value[i,] = c(res1[[2]], max(res2$v[,2]), max(res3$v[,2]),
                max(res4$v[,2]), max(res5$v[,2]), max(res6$v[,2]),max(res7$v[,2]), max(res8$v[,2]),max(res9$v[,2]), res10[[2]])
  num_knot[i,] = c(length(res1$beta), length(res2$beta)-1, length(res3$beta)-1,
                   length(res4$beta)-1, length(res5$beta)-1, length(res6$beta)-1,
                   length(res7$beta)-1, length(res8$beta)-1,length(res9$beta)-1, length(res10$beta))
  #record the estimated regression coefficients
  #so that we can see which coefficients are non-zero
}

avg_value = apply(value, 2, mean)
avg_knot = apply(num_knot, 2, mean)
tailoring_rank = apply(-value, 1, rank)
prop = matrix(0, 10, 10)
for (k in 1:10) {
  for (j in 1:10) {
    prop[k, j] = mean(tailoring_rank[k,] == j)
  }
}
rownames(prop) = c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10")
colnames(prop) = 1:10
prop = prop*100

res=cbind(prop,round(avg_value,2),round(avg_knot,2))


library(kableExtra)
kable(res,  'latex', booktabs = T,align = "c",
      caption = "Rows are tailoring variables and columns are the ranks.
We show the proportion of times that tailoring variables are ranked over 500 replications, together with the
average number of knots and average value functions")%>%
  collapse_rows(columns = 1, latex_hline = "major", valign = "middle")%>%
  kable_styling(latex_options = "HOLD_position",position = 'center')



### pdwols
load("~/Desktop/weights/pdwols.RData")
blip1=do.call(cbind,vic1)
blip2=do.call(cbind,bic1)
blip3=do.call(cbind,vic2)
blip4=do.call(cbind,bic2)
# con_select=do.call(cbind,confounder)
# rowMeans(con_select)
library(latex2exp)
par(mfrow=c(2,2))
boxplot(blip1[1,],blip2[1,],blip3[1,],blip4[1,],
        names =c('VIC+OAL', 'BIC+OAL','VIC','BIC'),
        ylab=TeX('$\\psi_{0}$'),xlab='Methods')
abline(h = .5,lty = 2)

boxplot(blip1[2,],blip2[2,],blip3[2,],blip4[2,],
        names =c('VIC+OAL', 'BIC+OAL','VIC','BIC'),
        ylab=TeX('$\\psi_{1}$'),xlab='Methods')
abline(h = -0.8,lty = 2)

boxplot(blip1[3,],blip2[3,],blip3[3,],blip4[3,],
        names =c('VIC+OAL', 'BIC+OAL','VIC','BIC'),
        ylab=TeX('$\\psi_{2}$'),xlab='Methods')
abline(h = -.5,lty = 2)

boxplot(blip1[4,],blip2[4,],blip3[4,],blip4[4,],
        names =c('VIC+OAL', 'BIC+OAL','VIC','BIC'),
        ylab=TeX('$\\psi_{3}$'),xlab='Methods')
abline(h = -.5,lty = 2)

true_sel=c(T,T,T,T, rep(F,7))
sel_rate1=apply(blip1!=0, 1, mean)
sel_rate2=apply(blip2!=0, 1, mean)
sel_rate3=apply(blip3!=0, 1, mean)
sel_rate4=apply(blip4!=0, 1, mean)

true_sel=c(T,T,T, rep(F,2),T,rep(F,4))
con_sel_rate=apply(con_select!=0, 1, mean)

rate=cbind(sel_rate1[-1],sel_rate2[-1],sel_rate3[-1],sel_rate4[-1])
colnames(rate)=c('VIC+OAL', 'BIC+OAL','VIC','BIC')

library(dplyr)
library(kableExtra)

kable(round(rate,2)*100,  'latex', booktabs = T,
      caption = "Variable Selection Rate of the Blip Parameters (Including the Noise Variables)")%>%
  kable_styling(latex_options = "HOLD_position",position = 'center')

devtools::load_all()
set.seed(999);ds_test=gen_inst(10000,10);Xtest=ds_test[[1]];rm(ds_test);true_psi=c(.5,-.8,-.5,-.5)
opt=as.numeric(cbind(1,Xtest[,1:3])%*%true_psi>=0)

opt1=apply(blip1, 2, function(i) cbind(1,Xtest)%*%i)>0
error1=1-mean(colMeans(apply(opt1, 2, function(i) i==opt)))
opt2=apply(blip2, 2, function(i) cbind(1,Xtest)%*%i)>0
error2=1-mean(colMeans(apply(opt2, 2, function(i) i==opt)))
opt3=apply(blip3, 2, function(i) cbind(1,Xtest)%*%i)>0
error3=1-mean(colMeans(apply(opt3, 2, function(i) i==opt)))
opt4=apply(blip4, 2, function(i) cbind(1,Xtest)%*%i)>0
error4=1-mean(colMeans(apply(opt4, 2, function(i) i==opt)))


tfree=.5-2*Xtest[,1]-0.6*exp(Xtest[,1])-2*Xtest[,2]+Xtest[,3]+2*Xtest[,6]
v_true=mean(tfree+opt*(cbind(1,Xtest[,1:3])%*%true_psi))
v_p1=mean(tfree+opt1*drop((cbind(1,Xtest[,1:3])%*%true_psi)))
v_p2=mean(tfree+opt2*drop((cbind(1,Xtest[,1:3])%*%true_psi)))
v_p3=mean(tfree+opt3*drop((cbind(1,Xtest[,1:3])%*%true_psi)))
v_p4=mean(tfree+opt4*drop((cbind(1,Xtest[,1:3])%*%true_psi)))

v_0=mean(tfree+0*(cbind(1,Xtest[,1:3])%*%true_psi))
v_1=mean(tfree+(cbind(1,Xtest[,1:3])%*%true_psi))

v=rbind(cbind(error1,error2,error3,error4),
        cbind(v_p1,v_p2,v_p3,v_p4))
colnames(v)=c('VIC+OAL', 'BIC+OAL','VIC','BIC')
rownames(v)=c('Error rate','Value')
library(kableExtra)
kable(round(v,2),  'latex', booktabs = T,
      caption = "Error rate and value using pdWOLS with sample size 500 (500 simulations).")%>%
  kable_styling(latex_options = "HOLD_position",position = 'center')

## pdwols rank
for (i in 1:500) {
  Y=ds[[i]][[2]];X=ds[[i]][[1]];A=ds[[i]][[3]]
  m=lm(Y~cbind(X,A*X))
  beta=coef(m)[2:11]
  w=abs(1/beta)^3

  m=glmnet(X,A,'binomial',penalty.factor = w)
  l=length(m$lambda)
  criterion=c()
  for (j in 1:l) {
    pi=exp(cbind(1,X)%*%coef(m)[,j])/(1+exp(cbind(1,X)%*%coef(m)[,j]))
    tau=A/(pi)+(1-A)/(1-pi)
    tauX=sweep(X,1,tau,'*')
    criterion[j]= crossprod(abs(beta),
                            abs(colSums(tauX*A)/sum(tau*A)-colSums(tauX*(1-A))/sum(tau*(1-A))))
  }

  alpha=coef(m)[,which.min(criterion)]
  confounder[[i]]=alpha[-1]
  ind=which(alpha[-1]!=0)
  m=glm(A~X[,ind],binomial)
  ps=fitted(m);w=abs(A-ps)

  ## pdwols
  m=sail(y=Y,e=A,x=X,penalty.factor=pfac,weights = w,
         basis=function(o) o)

  ## rank
  D=as.matrix(coef(m)[13:22,])
  rank1[[i]]=rowMeans(D!=0)

  ## vic
  l=length(m$lambda);v=c()
  for (j in 1:l) {
    theta=coef(m)[,j]
    psi=theta[12:22]
    num=sum(psi[-1]!=0)
    aopt=cbind(A,X)%*%psi>0
    v[j]=sum((A*aopt+(1-A)*(1-aopt))*Y/((A*ps)+(1-A)*(1-ps)))-kn*num
  }
  psi=coef(m)[,which.max(v)]
  ind=which(psi[-1]!=0)
  X=cbind(X,A,A*X)
  linear=lm(Y~X[,ind],weights = w)
  psi[ind+1]=coef(linear)[-1]
  vic1[[i]]=psi[12:22]
  ## bic
  df=m$dfbeta+m$dfalpha+1
  bic=n* log(m$rss/n)+df*log(n)
  psi=coef(m)[,which.min(bic)]
  ind=which(psi[-1]!=0)
  X=cbind(X,A,A*X)
  linear=lm(Y~X[,ind],weights = w)
  psi[ind+1]=coef(linear)[-1]
  bic1[[i]]=psi[12:22]
}

r1=do.call(cbind,rank1)
r2=do.call(cbind,rank2)

r1=rowMeans(r1)
r2=rowMeans(r2)
r=cbind(r1,r2)
kable(round(r,2)*100,  'latex', booktabs = T,
      caption = "Variable selection rate among all the solution paths")%>%
  kable_styling(latex_options = "HOLD_position",position = 'center')


single=do.call(rbind,rank1)

tailoring_rank2 = apply(-single, 1, rank,ties.method = 'first')
prop = matrix(0, 10, 10)
for (k in 1:10) {
  for (j in 1:10) {
    prop[k, j] = mean(tailoring_rank2[k,] == j)
  }
}
rownames(prop) = c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10")
colnames(prop) = 1:10
prop = prop*100

kable(prop,  'latex', booktabs = T,align = "c",
      caption = "Rows are tailoring variables and columns are the ranks.
We show the proportion of times that tailoring variables are ranked over 500 replications, together with the
average number of knots and average value functions")%>%
  collapse_rows(columns = 1, latex_hline = "major", valign = "middle")%>%
  kable_styling(latex_options = "HOLD_position",position = 'center')

D=c()
for (i in 1:500) {
  D[i]=cor(x=tailoring_rank[,i],y=tailoring_rank2[,i])
}


boxplot(apply(tailoring_rank[1:3,], 2, mean),apply(tailoring_rank2[1:3,], 2, mean),
        apply(tailoring_rank[4:10,], 2, mean),apply(tailoring_rank2[4:10,], 2, mean),
        ylab='Ranking',xlab='Methods',names =c(TeX('VR($X_1-X_3$)'),TeX('pdWOLS($X_1-X_3$)'),
                                               TeX('VR($X_4-X_{10}$)'),TeX('pdWOLS($X_4-X_{10}$)') ))




