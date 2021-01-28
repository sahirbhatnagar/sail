g2=function(n,p){
  X1=matrix(rbinom(n*p,size = 1,prob = 0.5),ncol = p)
  X1[X1==0]=-1
  A1=rbinom(n,size = 1,prob = 0.5)
  A2=rbinom(n,size = 1,prob = 0.5)

  expit=function(x,a){
    exp(x+2*a-1)/(1+exp(x+2*a-1))
  }

  prob=expit(X1[,1],A1)
  x2=rbinom(n,1,prob)
  X2=cbind(x2,matrix(rbinom(n*(p-1),size = 1,prob = 0.5),ncol = p-1))
  X2[X2==0]=-1
  mu=-1+2*X1[,1]+X2[,1]-A1+1.5*X1[,1]*A1+A2-1.5*X2[,1]*A2
  Y=rnorm(n,mu,1)
  return(list(Y,X1,X2,A1,A2))
}


