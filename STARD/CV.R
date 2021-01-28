p=ncol(X2)
pfac2=c(rep(1,p),0,rep(1,p))
pfac2=pfac2/sum(pfac2) * (2*p+1)
pfac2=c(pfac2[p+1],pfac2[-p-1])

## pdwols

m=cv.sail(y=Y2,e=A2,x=X2,nfolds = 4, weights = w2,penalty.factor=pfac2,
          parallel = T,basis=function(i) i)

sail2=coef(m,s='lambda.min')[(p+2):(2*p+2)]

aopt=as.numeric(cbind(1,X2)%*%sail2>0)
Y1[ind]=Y2+(aopt-A2)* cbind(1,X2)%*%sail2

## stage 1
p=ncol(X1)
pfac1=c(rep(1,p),0,rep(1,p))
pfac1=pfac1/sum(pfac1) * (2*p+1)
pfac1=c(pfac1[p+1],pfac1[-p-1])
m=cv.sail(y=Y1,e=A1,x=X1,nfolds = 4,weights = w1,penalty.factor=pfac1,
          parallel = T,basis=function(i) i)

sail1=coef(m,s='lambda.min')



