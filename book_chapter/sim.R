library(splines)
library(glmnet)
library(xtable)
library(dplyr)
devtools::load_all()
# outcome adaptive

## ols outcome

## dwols outcome
# fit_treat=glm(A~X,family = binomial)
# ps=fitted(fit_treat)
# w=abs(A-ps)
#
# m=lm(Y~cbind(X,A*X),weights = w)
# beta=coef(m)[2:11]
# w=abs(1/beta)^3
#
# m=glmnet(X,A,'binomial',penalty.factor = w)
# l=length(m$lambda)
# criterion=c()
# for (i in 1:l) {
#   pi=exp(cbind(1,X)%*%coef(m)[,i])/(1+exp(cbind(1,X)%*%coef(m)[,i]))
#   tau=A/(pi)+(1-A)/(1-pi)
#   tauX=sweep(X,1,tau,'*')
#   criterion[i]= crossprod(abs(beta),
#                           abs(colSums(tauX*A)/sum(tau*A)-colSums(tauX*(1-A))/sum(tau*(1-A))))
# }
# coef(m)[,which.min(criterion)]

## pdwols

Theta = seq(-2, 2, by=0.2)

knots = seq(-1.8, 1.8, by=0.2)

method = "adaptive lasso" # or "lasso"
criteria = "ERIC"  # or "BIC"
refit = TRUE;n=500;p=10
cn=log(log(n))
kn=cn*n^(1/3)*(log(p))^(2/3)
nsim = 500
nu=0.05
value = matrix(0, nsim, p)
num_knot = matrix(0, nsim, p)
useN = TRUE
pfac=c(0,rep(1,20))
vic1=bic1=vic2=bic2=confounder=list()
load("~/Desktop/weights/book_chapter/sim.RData")

## pdwols sim
for (i in 1:50) {
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

for (i in 1:50) {
  Y=ds[[i]][[2]];X=ds[[i]][[1]];A=ds[[i]][[3]]
  m=glm(A~X[,ind],binomial)
  ps=fitted(m);w=abs(A-ps)
  ## pdwols
  m=sail(y=Y,e=A,x=X,penalty.factor=pfac,weights = w,
         basis=function(o) o)
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
coef1 = matrix(0, nsim, length(knots))
coef2 = matrix(0, nsim, length(knots))
coef3 = matrix(0, nsim, length(knots))
coef4 = matrix(0, nsim, length(knots))
coef5 = matrix(0, nsim, length(knots))
coef6 = matrix(0, nsim, length(knots))
coef7 = matrix(0, nsim, length(knots))
coef8 = matrix(0, nsim, length(knots))
coef9 = matrix(0, nsim, length(knots))
coef10 = matrix(0, nsim, length(knots))

for (i in ) {
  Y=ds[[i]][[2]];X=ds[[i]][[1]];A=ds[[i]][[3]]
  dat=data.frame(id=1:n, y=Y, A, X)
  res1 = estimate.regime(dat, tailoring.name="X1", Theta, knots,
                         method = method, criteria = criteria, refit=refit, nu=nu, useN=useN)

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
  res10 = estimate.regime(dat, tailoring.name="X10", Theta, knots,
                          method = method, criteria = criteria, refit=refit, nu=nu, useN=useN)

  #record the maximum value function of each tailoring variabels, number of knots selected
  value[i,] = c(max(res1$v[,2]), max(res2$v[,2]), max(res3$v[,2]),
                max(res4$v[,2]), max(res5$v[,2]), max(res6$v[,2]),max(res7$v[,2]), max(res8$v[,2]),max(res9$v[,2]), max(res10$v[,2]))
  num_knot[i,] = c(length(res1$beta)-1, length(res2$beta)-1, length(res3$beta)-1,
                   length(res4$beta)-1, length(res5$beta)-1, length(res6$beta)-1,
                   length(res7$beta)-1, length(res8$beta)-1,length(res9$beta)-1, length(res10$beta)-1)

  #record the estimated regression coefficients
  #so that we can see which coefficients are non-zero
  knot_location = as.numeric(substr(names(res1$beta[-1]), 3, 10))
  for (k in 1:length(knot_location)) {
    index = which(knots==knot_location[k])
    coef1[i, index] = res1$beta[-1][k]
  }

  knot_location = as.numeric(substr(names(res2$beta[-1]), 3, 10))
  for (k in 1:length(knot_location)) {
    index = which(knots==knot_location[k])
    coef2[i, index] = res2$beta[-1][k]
  }


  knot_location = as.numeric(substr(names(res3$beta[-1]), 3, 10))
  for (k in 1:length(knot_location)) {
    index = which(knots==knot_location[k])
    coef3[i, index] = res3$beta[-1][k]
  }

  knot_location = as.numeric(substr(names(res2$beta[-1]), 3, 10))
  for (k in 1:length(knot_location)) {
    index = which(knots==knot_location[k])
    coef2[i, index] = res2$beta[-1][k]
  }

  knot_location = as.numeric(substr(names(res4$beta[-1]), 3, 10))
  for (k in 1:length(knot_location)) {
    index = which(knots==knot_location[k])
    coef4[i, index] = res4$beta[-1][k]
  }

  knot_location = as.numeric(substr(names(res5$beta[-1]), 3, 10))
  for (k in 1:length(knot_location)) {
    index = which(knots==knot_location[k])
    coef5[i, index] = res5$beta[-1][k]
  }

  knot_location = as.numeric(substr(names(res6$beta[-1]), 3, 10))
  for (k in 1:length(knot_location)) {
    index = which(knots==knot_location[k])
    coef6[i, index] = res6$beta[-1][k]
  }

  knot_location = as.numeric(substr(names(res7$beta[-1]), 3, 10))
  for (k in 1:length(knot_location)) {
    index = which(knots==knot_location[k])
    coef7[i, index] = res7$beta[-1][k]
  }

  knot_location = as.numeric(substr(names(res8$beta[-1]), 3, 10))
  for (k in 1:length(knot_location)) {
    index = which(knots==knot_location[k])
    coef8[i, index] = res8$beta[-1][k]
  }

  knot_location = as.numeric(substr(names(res9$beta[-1]), 3, 10))
  for (k in 1:length(knot_location)) {
    index = which(knots==knot_location[k])
    coef9[i, index] = res9$beta[-1][k]
  }

  knot_location = as.numeric(substr(names(res10$beta[-1]), 3, 10))
  for (k in 1:length(knot_location)) {
    index = which(knots==knot_location[k])
    coef10[i, index] = res10$beta[-1][k]
  }
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

res=cbind(prop,round(avg_value,2))


library(kableExtra)
kable(res,  'latex', booktabs = T,align = "c",
      caption = "Rows are tailoring variables and columns are the ranks.
We show the proportion of times that tailoring variables are ranked over 1000 replications, together with the
average number of knots and average value functions")%>%
  collapse_rows(columns = 1, latex_hline = "major", valign = "middle")%>%
  kable_styling(latex_options = "HOLD_position",position = 'center')

save.image("~/Desktop/weights/book.RData")

blip=do.call(cbind,vic1)
rowMeans(blip)

load("~/Desktop/weights/book.RData")

con_select=do.call(cbind,confounder)
rowMeans(con_select)

library(latex2exp)
boxplot(blip[1,],blip[2,],blip[3,],blip[4,],
        names =c(TeX('$\\psi_{0}$'),TeX('$\\psi_{1}$'),
                 TeX('$\\psi_{2}$'),TeX('$\\psi_{3}$')),
        ylab='Estimates',xlab='Blip parameters')

abline(h = .5,lty = 2)
abline(h = -.5,lty = 2,col='blue')
abline(h = -.8,lty = 2,col='pink')

true_sel=c(T,T,T,T, rep(F,7))
sel_rate=apply(blip!=0, 1, mean)

true_sel=c(T,T,T, rep(F,7))
con_sel_rate=apply(con_select!=0, 1, mean)

rate=cbind(con_sel_rate,sel_rate[-1])

library(dplyr)
library(kableExtra)

kable(round(rate,2)*100,  'latex', booktabs = T,
      caption = "Variable Selection Rate of the Blip Parameters (Including the Noise Variables)")%>%
  kable_styling(latex_options = "HOLD_position",position = 'center')

devtools::load_all()
set.seed(999);ds_test=gen_inst(10000,10);Xtest=ds_test[[1]];rm(ds_test);true_psi=c(.5,-.8,-.5,-.5)
opt=as.numeric(cbind(1,Xtest[,1:3])%*%true_psi>=0)
table(opt)

pdwols_opt=apply(blip, 2, function(i) cbind(1,Xtest)%*%i)>0
error=1-mean(colMeans(apply(pdwols_opt, 2, function(i) i==opt)))

tfree=.5-2*Xtest[,1]-0.6*exp(Xtest[,1])-2*Xtest[,2]+Xtest[,3]
v_true=mean(tfree+opt*(cbind(1,Xtest[,1:3])%*%true_psi))
v_p=mean(tfree+pdwols_opt*drop((cbind(1,Xtest[,1:3])%*%true_psi)))
v_0=mean(tfree+0*(cbind(1,Xtest[,1:3])%*%true_psi))
v_1=mean(tfree+(cbind(1,Xtest[,1:3])%*%true_psi))

v=cbind(error,v_p/v_true)





