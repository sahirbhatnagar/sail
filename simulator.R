library(simulator)

devtools::load_all()
library(doParallel)
registerDoParallel(cores = 4)
get_data=function(n, p, rho=0.15) {
  sig <- matrix(rho, p, p)
  diag(sig) <- 1
  X <- MASS::mvrnorm(n, mu=rep(0,p),Sigma = sig)
  x1=X[,1];x2=X[,2];x3=X[,3];x4=X[,4];x5=X[,5]

  expit=function(x) {
    return(1/(1+(exp(-(1+x[1]+x[2]+x[3]+x[4]+x[5])))))
  }

  prob=apply(X[,1:5], 1, expit)

  A=rbinom(n,1,prob)

  ## Generate Y
  tfree=3+3*exp(x1)+3*x2+3*x3+3*x4+3*x5
  psi=c(1,-1.5,-1.5,1)
  ymean=tfree+A*(cbind(1,X[,1:3])%*%psi)

  fit_treat=glm(A~X,family = binomial)
  ps=fitted(fit_treat)
  w=abs(A-ps)

  new_model(name = "model", label = sprintf("n = %s, p = %s", n, p),
            params = list(x=X, A=A,w=w,n = n,p = p,ymean=ymean),
            simulate = function(ymean, nsim) {
              y <- as.numeric(ymean)+ matrix(rnorm(nsim * n), n, nsim)
              return(split(y, col(y))) # make each col its own list element
            })
}

## The Methods

lasso=new_method("lasso", "Lasso",
                 method = function(model, draw, lambda = NULL) {


                   fit <- glmnet::cv.glmnet(x = cbind(model$x,model$A,model$x*model$A),
                                    y = draw)
                   list(psi = coef(fit,s='lambda.min')[12:22])})

sail=new_method("sail", "Sail",
                method = function(model, draw, lambda = NULL) {
                  fit <-sail:: cv.sail(x = model$x,e=model$A,y = draw,
                                 basis=function(i) i,parallel = T)
                list(psi = coef(fit,s='lambda.min')[12:22])})



lasso_weights=new_method("lasso_weights", "Lasso with Weights",
                 method = function(model, draw, lambda = NULL) {


                   fit <- glmnet::cv.glmnet(x = cbind(model$x,model$A,model$x*model$A),
                                            weights = model$w,y = draw)
                 list(psi = coef(fit,s='lambda.min')[12:22])})

sail_weights=new_method("sail_weights", "Sail with Weights",
                method = function(model, draw, lambda = NULL) {
                  fit <-cv.sail(x = model$x,e=model$A,y = draw,weights = model$w,
                                       basis=function(i) i,parallel = T)
                list(psi = coef(fit,s='lambda.min')[12:22])})


### The metrics

# select_rate=new_metric("sl_rate", "slection rate",
#                        metric = function(model, out) {
#                          colSums(as.matrix(out$psi!=0))
#                        })


# value_function=new_metric("v_function", "value function",
#                        metric = function(model, out) {
#                          set.seed(24)
#                          sig <- matrix(0.15, model$p, model$p)
#                          diag(sig) <- 1
#                          X_test <- MASS::mvrnorm(10000, mu=rep(0,model$p),Sigma = sig)
#                          opt=as.numeric(cbind(1,X_test)%*%out$psi>=0)
#
#                          tfree=3+3*exp(X_test[,1])+3*X_test[,2]+3*X_test[,3]+
#                            3*X_test[,4]+3*X_test[,5]
#                          psi=c(3,1,1,1)
#                          value=tfree+opt*(cbind(1,X_test[,1:3])%*%psi)
#                          mean(value)
#                        })


estimation=new_metric("estimation", "Estimation",
                       metric = function(model, out) {
                         out$psi
                       })

name_of_simulation <- "dtr"

sim <- new_simulation("DTR", "DTR") %>%
  generate_model(get_data, n = 800, p = 10) %>%
  simulate_from_model(nsim = 200, index = 1) %>%
  run_method(list(lasso,sail,lasso_weights,sail_weights)) %>%
  evaluate(list(estimation))

library(simulator)

# sim <- sim %>%
#   simulate_from_model(nsim = 200, index = 1) %>%
#   run_method(list(lasso,sail,lasso_weights,sail_weights)) %>%
#   evaluate(list(estimation))

## nsim=200

sim=load_simulation('dtr')
simulator::save_simulation(sim)

par=data.frame(evals(sim))

lasso_sel=subset(par,Method=='lasso')$estimation%>%
  matrix(nrow = 11)
sail_sel=subset(par,Method=='sail')$estimation%>%
  matrix(nrow = 11)
lasso_w_sel=subset(par,Method=='lasso_weights')$estimation%>%
  matrix(nrow = 11)
sail_w_sel=subset(par,Method=='sail_weights')$estimation%>%
  matrix(nrow = 11)

### Selection Rate
lasso_rate=paste(as.character(rowMeans(lasso_sel!=0)*100),'%')

sail_rate=paste(as.character(rowMeans(sail_sel!=0)*100),'%')

lasso_w_rate=paste(as.character(rowMeans(lasso_w_sel!=0)*100),'%')

sail_w_rate=paste(as.character(rowMeans(sail_w_sel!=0)*100),'%')

select_rate=rbind(lasso_rate,sail_rate,lasso_w_rate,sail_w_rate)

## Bias (SD)

mean_value=cbind(rowMeans(lasso_sel),rowMeans(sail_sel),
                 rowMeans(lasso_w_sel),rowMeans(sail_w_sel))

bias=mean_value[1:4,]-c(1,-1.5,-1.5,1)
names(bias)=c('lasso','sail','lasso','sail')
sd_value=cbind(apply(lasso_sel, 1, sd),apply(sail_sel, 1, sd),
         apply(lasso_w_sel, 1, sd),apply(sail_w_sel, 1, sd))

## Error Rate

set.seed(24)
sig <- matrix(0.15, 10, 10)
true_psi=c(1,-1.5,-1.5,1)
diag(sig) <- 1
X_test <- MASS::mvrnorm(10000, mu=rep(0,10),Sigma = sig)
opt=as.numeric(cbind(1,X_test[,1:3])%*%true_psi>=0)

error_rate=function(psi){
  est_opt=as.numeric(cbind(1,X_test)%*%psi>=0)
  mean(est_opt==opt)
}

error_lasso=mean(apply(lasso_sel, 2, error_rate))
error_w_lasso=mean(apply(lasso_w_sel, 2, error_rate))
error_sail=mean(apply(sail_sel, 2, error_rate))
error_w_sail=mean(apply(sail_w_sel, 2, error_rate))

mean_error=1-c(error_lasso,error_sail,error_w_lasso,error_w_sail)


## Value Function
tfree=3+3*exp(X_test[,1])+3*X_test[,2]+3*X_test[,3]+
  3*X_test[,4]+3*X_test[,5]
value=function(psi){
  est_opt=as.numeric(cbind(1,X_test)%*%psi>=0)
  value=tfree+est_opt*(cbind(1,X_test[,1:3])%*%true_psi)
  mean(value)
  }
v_true=mean(tfree+opt*(cbind(1,X_test[,1:3])%*%true_psi))

value_lasso=mean(apply(lasso_sel, 2, value))
value_w_lasso=mean(apply(lasso_w_sel, 2, value))
value_sail=mean(apply(sail_sel, 2, value))
value_w_sail=mean(apply(sail_w_sel, 2, value))

mean_value=c(value_lasso,value_sail,value_w_lasso,value_w_sail)

## Correct sparsity
true=c(1,1,1,1,rep(0,7))

lasso_sp=lasso_sel!=0
cp_lasso=true==lasso_sp
cp_lasso=apply(cp_lasso, 2, mean)
cp_lasso=mean(cp_lasso)

sail_sp=sail_sel!=0
cp_sail=true==sail_sp
cp_sail=apply(cp_sail, 2, mean)
cp_sail=mean(cp_sail)

lasso_w_sp=lasso_w_sel!=0
cp_w_lasso=true==lasso_w_sp
cp_w_lasso=apply(cp_w_lasso, 2, mean)
cp_w_lasso=mean(cp_w_lasso)

sail_w_sp=sail_w_sel!=0
cp_w_sail=true==sail_w_sp
cp_w_sail=apply(cp_w_sail, 2, mean)
cp_w_sail=mean(cp_w_sail)

cp=c(cp_lasso,cp_sail,cp_w_lasso,cp_w_sail)

## TP and FP

tp_lasso=mean(apply(lasso_sp[1:4,]==c(1,1,1,1),2,mean))
fp_lasso=1-mean(apply(lasso_sp[-(1:4),]==rep(0,7),2,mean))

tp_sail=mean(apply(sail_sp[1:4,]==c(1,1,1,1),2,mean))
fp_sail=1-mean(apply(sail_sp[-(1:4),]==rep(0,7),2,mean))

tp_w_lasso=mean(apply(lasso_w_sp[1:4,]==c(1,1,1,1),2,mean))
fp_w_lasso=1-mean(apply(lasso_w_sp[-(1:4),]==rep(0,7),2,mean))

tp_w_sail=mean(apply(sail_w_sp[1:4,]==c(1,1,1,1),2,mean))
fp_w_sail=1-mean(apply(sail_w_sp[-(1:4),]==rep(0,7),2,mean))

tp=c(tp_lasso,tp_sail,tp_w_lasso,tp_w_sail)

fp=c(fp_lasso,fp_sail,fp_w_lasso,fp_w_sail)




