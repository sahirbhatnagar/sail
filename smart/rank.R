create.pseudo = function(dat, tailoring.name, Theta_vec,direction="smaller") {
  n = nrow(dat)
  dat.pseudo = cbind(dat, 0)
  names(dat.pseudo) = c(names(dat),"regime")
  dat.pseudo = dat.pseudo[NULL,]
  A0 = dat$A
  P0 = dat[,tailoring.name]
  for(ii in 1:(length(Theta_vec))){
    if (direction == "smaller") {
      temp.index = (A0 == 1 & P0 < Theta_vec[ii]) | (A0 == 0 & P0 >= Theta_vec[ii])
    } else if (direction == "greater") {
      temp.index = (A0 == 1 & P0 > Theta_vec[ii]) | (A0 == 0 & P0 <= Theta_vec[ii])
    }
    temp.index[is.na(temp.index)] = FALSE
    if (sum(temp.index) == 0) next
    dat.p = cbind( dat[temp.index,], rep(Theta_vec[ii], sum(temp.index)) )
    names(dat.p) = c(names(dat), "regime")
    dat.pseudo = rbind(dat.pseudo, dat.p)
  }
  dat.pseudo = dat.pseudo[order(dat.pseudo$id, dat.pseudo$regime),]
  return(dat.pseudo)
}

create.analysis = function(dat, tailoring.name, Theta_vec,i) {
  dat$weights.trt = w[[i]]
  #create pseudo data set
  dat.pseudo = create.pseudo(dat, tailoring.name, Theta_vec)
  return(dat.pseudo)
}

linear.spline = function(x, knots) {
  n = length(x)
  p = length(knots)
  newx = matrix(NA, n, p+1)
  newx[,1] = x
  for (j in 1:p) newx[,j+1] = (x - knots[j])*(x > knots[j])
  return(newx)
}

estimate.regime = function(dat, tailoring.name, Theta_vec, knots,
                           method="adaptive lasso", criteria="ERIC",
                           refit=TRUE, nu=0.5, useN=FALSE)
{
  if (method == "lasso" & criteria == "ERIC") stop("ERIC cannot be used with LASSO")
  dat.pseudo = create.analysis(dat, tailoring.name, Theta_vec,i)
  n = nrow(dat)
  N = nrow(dat.pseudo)
  X = linear.spline(dat.pseudo$regime, knots)
  #X = cbind(X, dat.pseudo[,tailoring.name]) #add baseline variable
  for (j in 1:ncol(X)) X[,j] = X[,j]/max(abs(X[,j]))
  colnames(X) = c("linear", knots)

  ###### method #####
  if (method == "adaptive lasso") {
    fit0 = lm(dat.pseudo$y ~ X, weights=dat.pseudo$weights.trt)
    w = 1/abs(fit0$coefficients[c(-1,-2)]^2)
    w[w>10^2] = 10^2
    penalty.factor = c(0, w)
    penalty.factor = penalty.factor/sum(penalty.factor)*ncol(X)
  }
  else if (method == "lasso") {
    penalty.factor = c(0, rep(1, ncol(X)-1))
  }
  else if (method == "nonparametric") {
    #### nonparametric estimates ####
    v = dat.pseudo %>% group_by(regime) %>%
      summarise(value=mean(weights.trt*Y)/mean(weights.trt),
                n=n(), n0=sum(A0==0), n1=sum(A0==1))
    return(v)
  }

  fit = glmnet(x=X, y=dat.pseudo$y, weights=dat.pseudo$weights.trt,
               lambda=NULL, penalty.factor=penalty.factor, standardize = FALSE)
  pred = predict(fit, newx=X)

  ##### model selection criteria #####
  if(useN==TRUE)  num = N else num = n
  aic = numeric(ncol(pred))
  bic = numeric(ncol(pred))
  eric = numeric(ncol(pred))
  for (j in 1:ncol(pred)) {
    wmse = mean(dat.pseudo$weights.trt*(dat.pseudo$y-pred[,j])^2)
    p = fit$df[j] - 1
    aic[j] = log(wmse)*num + 2*p
    bic[j] = log(wmse)*num + log(num)*p
    eric[j] = log(wmse)*num + 2*nu*log(num*wmse/fit$lambda[j])*p
  }
  if (criteria == "AIC")
    min.lambda = which.min(aic)
  else if (criteria == "BIC")
    min.lambda = which.min(bic)
  else if (criteria == "ERIC")
    min.lambda = which.min(eric)

  ##### whether OLS refit is performed #####
  if (refit == FALSE) {
    beta = fit$beta[,min.lambda]
    fitted = predict(fit, newx=X, s=fit$lambda[min.lambda])
    dat.pseudo$fitted = fitted
    v = dat.pseudo %>% group_by(regime) %>% summarise(value=mean(weights.trt*fitted)/mean(weights.trt))
  }
  else if (refit == TRUE) {
    var = abs(fit$beta[,min.lambda])>1e-10
    if (sum(var) == 0) {
      res = lm(dat.pseudo$y ~ 1, weights=dat.pseudo$weights.trt)
    } else {
      sX = X[,var]
      res = lm(dat.pseudo$y ~ sX, weights=dat.pseudo$weights.trt)
    }
    beta = res$coef[-1]
    dat.pseudo$fitted = res$fitted
    v = dat.pseudo %>% group_by(regime) %>% summarise(value=mean(weights.trt*fitted)/mean(weights.trt))
  }
  return(list(beta=beta, v=v))
}

rank_binary=function(x,choice,i){
  aopt=rep(0,length(x))
  if(choice =='opposite'){
    aopt[x==0]=1
  } else{aopt[x==1]=1}
  va=sum((Y* w[[i]])[A==aopt])/sum(w[[i]][A==aopt])
  beta=0
  list(beta,va)
}




