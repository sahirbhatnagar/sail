ds=read.csv("~/Desktop//weights/smart.csv")
A=ds$FINAL_GROUP
Y=ds$DASS_score____6_week
index=c(6,8,37,38,41,9, 34,35,30, 10, 18, 23,15,11,21,25)
X=ds[,index]
levels(X$marital)=c(levels(X$marital),'couple','other')
X$marital[X$marital=='Married'|X$marital=="Common law / Living together"]='couple'
X$marital[X$marital!='couple']='other'

levels(X$education)=c(levels(X$education),'other')
X$education[X$education!='University degree']='other'

levels(X$employment)=c(levels(X$employment),'other')
X$employment[X$employment!='Full time employment']='other'
for (i in c(2:5,9:16)) {
  X[,i]=as.integer(factor(X[,i]))
}

for (i in c(2:5,9:16)) {
  X[,i]=X[,i]-1
}

X=as.matrix(X)
Y=-Y
A=as.integer(factor(A))-1
## table one dude
# library(tableone)
# A=as.integer(A)-1
# cindex=c(1,6:8)
# vars=colnames(X)[cindex]
# fac=colnames(X)[-cindex]
# tableOne <- CreateTableOne(vars = colnames(X), strata = "A",
#                            data = as.data.frame(cbind(A,X)),test=F)
# library(kableExtra)
# p <- print(tableOne, printToggle = FALSE, noSpaces = TRUE)
# kableone(p, booktabs = TRUE, format = "latex",
#          caption = "Characteristics of the study population stratified by treatment")%>%
#   kable_styling(latex_options = "HOLD_position",position = 'center')


# devtools::load_all()
# load("~/Desktop/weights/smart.RData")
#
m=glm(A~ds$PROV+ds$DASS_S,binomial)
ps=fitted(m)
w=abs(A-ps)
n=nrow(X);p=ncol(X)
pfac=c(0,rep(1,2*ncol(X)))
cn=log(log(n))
kn=cn*n^(1/3)*(log(p))^(2/3)
m=sail(y=Y,e=A,x=X,penalty.factor=pfac,weights = w,
       basis=function(o) o)
## rank
D=as.matrix(coef(m)[18:34,])
rowMeans(D!=0)
sort(rowMeans(D!=0),decreasing = T)
## model selection
l=length(m$lambda)
v=c()
for (j in 1:l) {
  theta=coef(m)[,j]
  psi=theta[18:34]
  num=sum(psi[-1]!=0)
  aopt=cbind(A,X)%*%psi>0
  v[j]=sum((A*aopt+(1-A)*(1-aopt))*Y/((A*ps)+(1-A)*(1-ps)))-kn*num
}
psi=coef(m)[,which.max(v)]
ind=which(psi[-1]!=0)
x=cbind(X,A,A*X)
linear=lm(Y~x[,ind],weights = w)
psi[ind+1]=coef(linear)[-1]

## bootstrap
value_boot=c()
for (i in 1:4000) {
  ind=sample(1:50,50,replace = T)
  aopt=A[ind]*(cbind(1,X[ind,])%*%psi[18:34])>0
  value_boot[i]=mean(Y[ind]+aopt*(cbind(1,X[ind,]) %*%psi[18:34]))
}

mean(value_boot)
mean(value_boot)+1.96*c(-sd(value_boot),sd(value_boot))

#### OAL
x=cbind(X,as.integer(factor(ds$PROV))-1)
m=lm(Y~cbind(x,A*x))
beta=coef(m)[2:18]
w=abs(1/beta)^3

m=glmnet(x,A,'binomial',penalty.factor = w)
l=length(m$lambda)
criterion=c()
for (j in 1:l) {
  pi=exp(cbind(1,x)%*%coef(m)[,j])/(1+exp(cbind(1,x)%*%coef(m)[,j]))
  tau=A/(pi)+(1-A)/(1-pi)
  tauX=sweep(x,1,tau,'*')
  criterion[j]= crossprod(abs(beta),
                          abs(colSums(tauX*A)/sum(tau*A)-colSums(tauX*(1-A))/sum(tau*(1-A))))
}

alpha=coef(m)[,which.min(criterion)]
ind=which(alpha[-1]!=0)
m=glm(A~x[,ind],binomial)
ps=fitted(m);w=abs(A-ps)

## pdwols
m=sail(y=Y,e=A,x=X,penalty.factor=pfac,weights = w,
       basis=function(o) o)

## vic
l=length(m$lambda);v=c()
for (j in 1:l) {
  theta=coef(m)[,j]
  psi=theta[18:34]
  num=sum(psi[-1]!=0)
  aopt=cbind(A,X)%*%psi>0
  v[j]=sum((A*aopt+(1-A)*(1-aopt))*Y/((A*ps)+(1-A)*(1-ps)))-kn*num
}
psi=coef(m)[,which.max(v)]
ind=which(psi[-1]!=0)
X=cbind(X,A,A*X)
linear=lm(Y~X[,ind],weights = w)
psi[ind+1]=coef(linear)[-1]

####  ranking methods
n=50
dat=data.frame(id=1:n, y=Y, A, X)

method = "adaptive lasso" # or "lasso"
criteria = "ERIC"  # or "BIC"
refit = TRUE;nu=0.05
useN = TRUE

Theta1 = seq(48, 74, by=1)
knots1 = seq(50, 72, by=1)
res1 = estimate.regime(dat, tailoring.name="age", Theta1, knots1,
                       method = method, criteria = criteria, refit=refit, nu=nu, useN=useN)
res2 = rank_binary(X[,2],'opposite')
res3 = rank_binary(X[,3],'opposite')
res4 = rank_binary(X[,4],'opposite')
res5 = rank_binary(X[,5],'opposite')

Theta2 = seq(12, 30, by=.5)
knots2 = seq(13, 29, by=.5)
res6 = estimate.regime(dat, tailoring.name="DASS_S", Theta2, knots2,
                       method = method, criteria = criteria, refit=refit, nu=nu, useN=useN)

Theta3 = seq(26, 58, by=1)
knots3 = seq(27, 57, by=1)
res7 = estimate.regime(dat, tailoring.name="PCS12_b", Theta3, knots3,
                       method = method, criteria = criteria, refit=refit, nu=nu, useN=useN)

Theta4 = seq(29, 54, by=1)
knots4 = seq(30, 53, by=1)
res8 = estimate.regime(dat, tailoring.name="MCS12_b", Theta4, knots4,
                       method = method, criteria = criteria, refit=refit, nu=nu, useN=useN)
res9 = rank_binary(X[,9],'opposite')
res10 = rank_binary(X[,10],'opposite')
res11 = rank_binary(X[,11],'opposite')
res12 = rank_binary(X[,12],'opposite')
res13 = rank_binary(X[,13],'opposite')
res14 = rank_binary(X[,14],'opposite')
res15 = rank_binary(X[,15],'opposite')
res16 = rank_binary(X[,16],'opposite')

res1. = estimate.regime(dat, tailoring.name="age", Theta1, knots1,
                       method = method, criteria = criteria, refit=refit, nu=nu, useN=useN)
res2. = rank_binary(X[,2],'o')
res3. = rank_binary(X[,3],'o')
res4. = rank_binary(X[,4],'o')
res5. = rank_binary(X[,5],'o')

res6. = estimate.regime(dat, tailoring.name="DASS_S", Theta2, knots2,
                       method = method, criteria = criteria, refit=refit, nu=nu, useN=useN)

res7. = estimate.regime(dat, tailoring.name="PCS12_b", Theta3, knots3,
                       method = method, criteria = criteria, refit=refit, nu=nu, useN=useN)

res8. = estimate.regime(dat, tailoring.name="MCS12_b", Theta4, knots4,
                       method = method, criteria = criteria, refit=refit, nu=nu, useN=useN)
res9. = rank_binary(X[,9],'o')
res10. = rank_binary(X[,10],'o')
res11. = rank_binary(X[,11],'o')
res12. = rank_binary(X[,12],'o')
res13. = rank_binary(X[,13],'o')
res14. = rank_binary(X[,14],'o')
res15. = rank_binary(X[,15],'o')
res16. = rank_binary(X[,16],'o')



#record the maximum value function of each tailoring variabels, number of knots selected
res_value = c(max(res1$v[,2]), res2[[2]],res3[[2]],res4[[2]],res5[[2]], max(res6$v[,2]),
              max(res7$v[,2]), max(res8$v[,2]),res9[[2]], res10[[2]],res11[[2]],res12[[2]],
              res13[[2]],res14[[2]],res15[[2]],res16[[2]])

res_value2=c(max(res1.$v[,2]), res2.[[2]],res3.[[2]],res4.[[2]],res5.[[2]], max(res6.$v[,2]),
             max(res7.$v[,2]), max(res8.$v[,2]),res9.[[2]], res10.[[2]],res11.[[2]],res12.[[2]],
             res13.[[2]],res14.[[2]],res15.[[2]],res16.[[2]])

names(res_value)=colnames(X)
names(res_value2)=colnames(X)

rres=cbind(res_value,res_value2)
rres=apply(rres, 1, max)
names(rres)=colnames(X)
names(sort(rres,decreasing = T))



num_knot[i,] = c(length(res1$beta), length(res2$beta)-1, length(res3$beta)-1,
                 length(res4$beta)-1, length(res5$beta)-1, length(res6$beta)-1,
                 length(res7$beta)-1, length(res8$beta)-1,length(res9$beta)-1, length(res10$beta))
#record the estimated regression coefficients



summary(dat$age)
summary(dat$DASS_S)
summary(dat$PCS12_b)
summary(dat$MCS12_b)

value_boot=list()
for (i in 1:4000) {
  ind=sample(1:50,50,replace = T)
  dat=data.frame(id=1:n, y=Y[ind], A[ind], X[ind,])
  res7 = estimate.regime(dat, tailoring.name="PCS12_b", Theta3, knots3,
                         method = method, criteria = criteria, refit=refit, nu=nu, useN=useN)
  res11 = rank_binary(X[ind,11],'opposite')
  res3 = rank_binary(X[ind,3],'opposite')
  res6 = estimate.regime(dat, tailoring.name="DASS_S", Theta2, knots2,
                         method = method, criteria = criteria, refit=refit, nu=nu, useN=useN)

  value_boot[[i]]=c(max(res7$v[,2]),res11[[2]],res3[[2]],max(res6$v[,2]))
}

jude=do.call(cbind,value_boot)

rowMeans(jude)
rowMeans(jude)+1.96*c(sd(jude[1,]),sd(jude[2,]),sd(jude[3,]),sd(jude[4,]))
rowMeans(jude)-1.96*c(sd(jude[1,]),sd(jude[2,]),sd(jude[3,]),sd(jude[4,]))


