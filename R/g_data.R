g_data=function(n,p){
  sig <- matrix(NA, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      sig[i,j]=0.25^(abs(i-j))
    }
  }
  X <- MASS::mvrnorm(n, mu=rep(0,p),Sigma = sig)
  x1=X[,1];x2=X[,2]

  expit=function(x) {
    return(1/(1+(exp(-(1+x[1]+x[2])))))
  }
  prob=apply(X[,1:2], 1, expit)
  A=rbinom(n,1,prob)
  ## Generate Y
  tfree=1-2*(exp(x1)+log(abs(x1)))+2*x2
  psi=c(1,-1.5)
  ymean=tfree+A*(cbind(1,X[,1])%*%psi)
  y=rnorm(n,mean=ymean,sd = 1)
  fit_treat=glm(A~X,family = binomial)
  ps=fitted(fit_treat)
  w=abs(A-ps)
  return(list(X,y,w,A))
}






