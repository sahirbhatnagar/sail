gen_inst=function(n,p){
  p=p-2
  sig <- matrix(NA, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      sig[i,j]=0.15^(abs(i-j))
    }
  }
  X <- MASS::mvrnorm(n, mu=rep(0,p),Sigma = sig)
  X=cbind(rbinom(n,1,0.5),X,rbinom(n,1,0.5))

  x1=X[,1];x2=X[,2];x3=X[,3];x4=X[,4];x5=X[,5];x6=X[,6]

  expit=function(x) {
    return(1/(1+(exp(-(1+x1+x2+x3+x4+x5)))))
  }
  prob=apply(X[,1:5], 1, expit)
  A=rbinom(n,1,prob)
  ## Generate Y
  tfree=.5-2*x1-0.6*exp(x1)-2*x2+x3+2*x6
  psi=c(.5,-.8,-.5,-.5)
  ymean=tfree+A*(cbind(1,x1,x2,x3)%*%psi)
  y=rnorm(n,mean=ymean,sd = 1)
  return(list(X,y,A))
}


