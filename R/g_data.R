
g_data=function(){
  X=matrix(rnorm(10000),1000)
  x1=X[,1];x2=X[,2];x3=X[,3];x4=X[,4]

  ## treatment model
  # b0=0.4
  # b1=b2=b3=b4=1
  # b=c(b1,b2,b3,b4)

  ##  A~X1-X4
  expit=function(x) {
    return(1/(1+(exp(-(0.4+x[1]+x[2]+x[3]+x[4])))))
  }

  prob=apply(X[,1:4], 1, expit)

  A=rbinom(1000,1,prob)


  # library(tableone)
  # ds=data.frame(cbind(A,X))
  # colnames(ds)[-1]=colnames(ds)[-1]=c('X1',"X2"  ,"X3" , "X4" , "X5" , "X6" , "X7" , "X8" , "X9" , "X10")
  # print(CreateTableOne(names(ds)[-1],strata = 'A',data = ds,test = F),smd = T)
  ### Propensity score distribution
  # treat_model=glm(A~X,family = binomial)
  # ps=fitted(treat_model)
  # par(mfrow=c(2,1))
  # par(mar=c(0,5,3,3))
  # hist(ps[A==1] , main="" , xlim=c(0,1),
  #      ylab="Treatment", xlab="",
  #      xaxt="n", las=1 , col="slateblue1")
  # par(mar=c(5,5,0,3))
  # hist(ps[A==0], main="" , xlim=c(0,1),ylim = c(100,0), ylab="Non-Treatment",
  #      xlab="Propensity Score",  las=1 , col="tomato3")

  ## Generate Y
  tfree=1+exp(x1)+exp(x2)-2*x3+4*x4
  psi=c(5,5,5,-1,-1)
  ymean=tfree+A*(cbind(1,X[,1:4])%*%psi)
  y=rnorm(1000,ymean,1)

  ## Both Correct

  fit_treat=glm(A~X,family = binomial)
  ps=fitted(fit_treat)
  w=abs(A-ps)
  w=w*length(w)/sum(w)
  return(list(X,y,w,A))
}






