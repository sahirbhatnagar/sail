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
save.image("~/Desktop/weights/smart.RData")
for (i in 1:ncol(X)) {
  X[,i]=as.integer(X[,i])
}

X$employment[X$employment==10]=0

for (i in c(2,9:16)) {
  X[,i]=X[,i]-1
}

X$marital=X$marital-7
X$education=X$education-4
X=as.matrix(X)

## table one dude
library(tableone)
A=as.integer(A)-1
cindex=c(1,6:8)
vars=colnames(X)[cindex]
fac=colnames(X)[-cindex]
tableOne <- CreateTableOne(vars = colnames(X), strata = "A",
                           data = as.data.frame(cbind(A,X)),test=F)
library(kableExtra)
p <- print(tableOne, printToggle = FALSE, noSpaces = TRUE)
kableone(p, booktabs = TRUE, format = "latex",
         caption = "Characteristics of the study population stratified by treatment")%>%
  kable_styling(latex_options = "HOLD_position",position = 'center')


devtools::load_all()
load("~/Desktop/weights/smart.RData")
Y=-Y
m=glm(A~ds$PROV+ds$DASS_S,binomial)
ps=fitted(m)
w=abs(A-ps)
n=nrow(X);p=ncol(X)
pfac=c(0,rep(1,2*ncol(X)))
cn=log(log(n))
kn=cn*n^(1/3)*(log(p))^(2/3)
m=sail(y=Y,e=A,x=X,penalty.factor=pfac,weights = w,
       basis=function(o) o)
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







