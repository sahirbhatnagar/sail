cv.dtr=function(x,y,e,w){
  train=sample(length(y),round(length(y)*.6))
  y=y-weighted.mean(y,w)
  A=A-weighted.mean(A,w)
  ff=function(s){
    return(s-weighted.mean(s,w))
  }
  x=apply(x, 2, ff)
  xe=apply(x*A, 2, ff)

  x_train=x[train,];x_test=x[-train,]
  y_train=y[train];y_test=y[-train]
  A_train=A[train];A_test=A[-train]
  w_train=w[train]

  xe_test=xe[-train,]
  model=sail(x = x_train, alpha=0.5 , weights = w_train,
             maxit=150, y = y_train, e = A_train, basis=function(i) i)

  par=list()
  test_error=c()

  for (i in 1:length(model$lambda)) {
    par[[i]]=coef(model)[,i]
    pred=x_test%*%model$beta[,i]+model$bE[i]*A_test+xe_test%*%model$alpha[,i]
    test_error[i]=sum((y_test-pred)^2)
  }

  return(par[[which.min(test_error)]])

}
