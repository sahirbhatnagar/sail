gen_high=function(n,p,rho=.15){
  sig <- matrix(NA, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      sig[i,j]=0.15^(abs(i-j))
    }
  }
  X <- MASS::mvrnorm(n, mu=rep(0,p),Sigma = sig)
  x1=X[,1];x2=X[,2]
  A=rbinom(n,1,0.5)
  ## Generate Y
  tfree=1-2*x1+2*sin(pi*x1)-exp(x2)-2*x2
  # psi=c(1,-1.5)
  psi=c(0,0)
  ymean=tfree+A*(cbind(1,X[,1])%*%psi)
  y=rnorm(n,mean=ymean,sd = 1)
  return(list(X,y,A))
}






