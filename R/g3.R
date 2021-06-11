g3=function(n,p){

  X1 <- matrix(rnorm(n*p),nrow = n)

  x1=X1[,1];x2=X1[,2]
  expit=function(x) {
    return(exp(x[1]-x[2])/(1+exp(x[1]-x[2])))
  }

  prob1=apply(X1[,1:2], 1, expit)

  A1=rbinom(n,1,prob1)

  X2=rnorm(n,mean=0.5*A1+0.8*x1)

  for (j in 2:p) {
    X2=cbind(X2,rnorm(n,mean = 0.8*X1[,j]))
  }

  prob2=apply(X2[,1:2], 1, expit)
  A2=rbinom(n,1,prob2)
  opt1=0.8-2*x1>0
  opt2=1-1.5*X2[,1]>0

  ## Generate Y
  yopt=.5+2*x1+2*x2
  ymean=yopt-(opt1-A1)*(.8-2*x1)-(opt2-A2)*(1-1.5*X2[,1])
  y=rnorm(n,ymean,1)

  fit_treat=glm(A1~X1,family = binomial)
  ps=fitted(fit_treat)
  w1=abs(A1-ps)
  fit_treat=glm(A2~X2,family = binomial)
  ps=fitted(fit_treat)
  w2=abs(A2-ps)
  return(list(X1,X2,y,w1,w2,A1,A2))

}






