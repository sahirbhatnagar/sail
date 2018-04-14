rm(list=ls())
devtools::load_all("/home/sahir/git_repositories/sail/")
pacman::p_load(magrittr)
pacman::p_load(dplyr)
pacman::p_load(splines)
pacman::p_load(foreach)
pacman::p_load(doMC)

f.basis <- function(i) splines::bs(i, degree = 3)
data("oasis")
system.time(fit <- sail(x = oasis$x, y = oasis$y, e = oasis$e,
                        basis = f.basis,
                        # dfmax = 10,
                        # thresh = 1e-3,
                        # expand = FALSE, group = attr(Xm,"assign"),
                        # fdev = 1e-10,
                        nlambda = 100, verbose = 2, alpha = 0.2))

plot(fit)
fit

rm(list=ls())
devtools::load_all("/home/sahir/git_repositories/sail/")
pacman::p_load(magrittr)
pacman::p_load(dplyr)
pacman::p_load(splines)
pacman::p_load(foreach)
pacman::p_load(doMC)
f.basis <- function(i) splines::bs(i, degree = 5)
data("oasis")
registerDoMC(cores = 8)
system.time(cvfit <- cv.sail(x = oasis$x, y = oasis$y, e = oasis$e,
                             basis = f.basis,
                             dfmax = 10,
                             fdev = 1e-10,
                             # expand = FALSE, group = attr(Xm,"assign"),
                             parallel = TRUE, nfolds = 20,
                             alpha = 0.2,
                             nlambda = 100, verbose = 1))

data("oasis")
f.basis <- function(i) splines::bs(i, degree = 3)
Xm <- model.matrix(~0+bs(Age, 3)+EDUC+bs(MMSE,3)+bs(eTIV,3)+bs(nWBV,3)+bs(ASF,3)+
                     bs(noise1, 3) + bs(noise2, 3)+ bs(noise3, 3)+ bs(noise4, 3)+ bs(noise5, 3)+
                     bs(noise6, 3) + bs(noise7, 3)+ bs(noise8, 3)+ bs(noise9, 3)+ bs(noise10, 3)+
                     bs(noise11, 3) + bs(noise12, 3)+ bs(noise13, 3)+ bs(noise14, 3),
                   data = as.data.frame(oasis$x))
head(Xm)
system.time(fit <- sail(x = oasis$x, y = oasis$y, e = oasis$e,
                        basis = f.basis,
                        # expand = FALSE, group = attr(Xm,"assign"),
                        # fdev = 1e-10,
                        nlambda = 100, verbose = 2, alpha = 0.1))
plot(fit)
fit


rm(list=ls())
devtools::load_all("/home/sahir/git_repositories/sail/")
pacman::p_load(magrittr)
pacman::p_load(dplyr)
pacman::p_load(splines)
pacman::p_load(foreach)
pacman::p_load(doMC)
data("oasis")
head(oasis$x)
Xm <- model.matrix(~0+bs(Age,3)+EDUC+bs(MMSE,3)+bs(eTIV,3)+bs(nWBV,3)+bs(ASF,3)+
                     bs(noise1, 3) + bs(noise2, 3)+ bs(noise3, 3)+ bs(noise4, 3)+ bs(noise5, 3)+
                     bs(noise6, 3) + bs(noise7, 3)+ bs(noise8, 3)+ bs(noise9, 3)+ bs(noise10, 3)+
                     bs(noise11, 3) + bs(noise12, 3)+ bs(noise13, 3)+ bs(noise14, 3),
                   data = as.data.frame(oasis$x))
registerDoMC(cores = 8)
system.time(cvfit <- cv.sail(x = Xm, y = oasis$y, e = oasis$e,
                              basis = f.basis,
                              # fdev = 1e-10,
                              # expand = FALSE, group = attr(Xm,"assign"),
                              parallel = TRUE, nfolds = 10,
                              alpha = 0.1,
                              nlambda = 100, verbose = 1))

fit
plot(fit)
coef(fit)
coef(cvfit, s="lambda.min")
coef(cvfit, s="lambda.1se")
plot(cvfit)


system.time(fit <- sail(x = X, y = Y, e = E,
                        nlambda = 100, verbose = 2, alpha = 0.5))
system.time(cvfit <- cv.sail(x = Xm, y = Y, e = E, expand = FALSE, group = attr(Xm, "assign"),
                             fdev = 1e-10,
                             parallel = TRUE, nfolds = 10, alpha = 0.1,
                             nlambda = 100, verbose = 1))

plot(cvfit)
plot(fit)
str(fit)
coef(fit)
fit
registerDoMC(cores = 8)
# f.basis <- function(i) splines::bs(i, degree = 3)
f.basis <- function(i) splines::bs(i, degree = 3)
foldid <- sample(seq_along(Y2),size=length(Y2),replace=TRUE)
system.time(cvfit <- cv.sail(x = X2, y = Y2, e = E2, basis = f.basis,
                             # fdev = 1e-10,
                             parallel = TRUE, nfolds = 5, alpha = 0.5,
                             nlambda = 100, verbose = 1))
plot(fit)
coef(cvfit, s="lambda.min")
coef(cvfit, s="lambda.1se")
coef(fit)
plot(cvfit)



cor(predict(cvfit, s = "lambda.min", newx = X2, newe = E2),
    Y2)^2
plot(predict(cvfit, s = "lambda.min", newx = X2, newe = E2), Y2)
abline(a=0,b=1)

cor(predict(cvfit, s = "lambda.1se", newx = X2, newe = E2),
    Y2)^2
plot(predict(cvfit, s = "lambda.1se", newx = X2, newe = E2), Y2)
abline(a=0,b=1)

head(X)
X.m <- model.matrix(~ 0+(Age+EDUC+MMSE+eTIV+nWBV+ASF+noise1+noise2+noise3+
                           noise4+noise5+noise6+noise7+noise8+noise9+noise10+
                         noise11+noise12+noise13+noise14)*E,
                    data = as.data.frame(cbind(X,E)))
head(X.m)
glmnetfit <- cv.glmnet(x = X.m, y = Y)
plot(glmnetfit)
coef(glmnetfit)

X.mt <- model.matrix(~ 0+(Age+EDUC+MMSE+eTIV+nWBV+ASF+noise1+noise2+noise3+
                           noise4+noise5+noise6+noise7+noise8+noise9+noise10+
                           noise11+noise12+noise13+noise14)*E2,
                    data = as.data.frame(cbind(X2,E2)))

cor(predict(glmnetfit, newx = X.mt, s = "lambda.min"),
    Y2)^2





SAheart %>% str
X %>% str
newData <- lapply(as.data.frame(X),
                  function(x) {
                    r <- range(as.numeric(x))
                    y <- seq(r[1], r[2], length.out = 400)
                    if(is.factor(x))
                      y <- levels(x)[as.integer(y)]
                    return(y)
                  }
)

newData
plot_apoe_inter <- function(cv_obj, X, original_name = "X60", sail_name = "X19", lambda_type = "lambda.min",
                            xlab = "supramarginal gyrus right", ylab = "Mini-Mental State Examination",
                            apoe = FALSE, main = "APOE e4 = 0",
                            color = RSkittleBrewer::RSkittleBrewer(flavor = "trop"),
                            legend = TRUE, ylim =  c(15,30)) {

  cv_obj = cvfit; X =X;
  lambda_type = "lambda.min";
  sail_name = "MMSE";
  color = RColorBrewer::brewer.pal(9,"Set1");
  xlab =  "supramarginal gyrus right";
  ylab = "Mini-Mental State Examination";
  main = "APOE e4 = 0"
  legend = TRUE; ylim =  c(15,30)
  # ==================
  dfs <- cv_obj$sail.fit$ncols
  lin_pred <- coef(cv_obj, s = lambda_type)["(Intercept)",,drop=T] +
    # cv_obj$sail.fit$design[,paste("X97",seq_len(dfs), sep = "_")] %*%
    # coef(cv_obj, s = lambda_type)[paste("X97",seq_len(dfs), sep = "_"),,drop=F] +
    cv_obj$sail.fit$design[,paste(sail_name,seq_len(dfs), sep = "_")] %*%
    coef(cv_obj, s = lambda_type)[paste(sail_name,seq_len(dfs), sep = "_"),,drop=F] +
    cv_obj$sail.fit$design[,"E"] %*%
    coef(cv_obj, s = lambda_type)["E",,drop=F] +
    cv_obj$sail.fit$design[,paste0(paste(sail_name,seq_len(dfs), sep = "_"),":E")] %*%
    coef(cv_obj, s = lambda_type)[paste0(paste(sail_name,seq_len(dfs), sep = "_"),":E"),,drop=F] #+
    # cv_obj$sail.fit$design[,paste("X98",seq_len(dfs), sep = "_")] %*%
    # coef(cv_obj, s = lambda_type)[paste("X98",seq_len(dfs), sep = "_"),,drop=F] +
    # cv_obj$sail.fit$design[,paste0(paste("X98",seq_len(dfs), sep = "_"),":X_E")] %*%
    # coef(cv_obj, s = lambda_type)[paste0(paste("X98",seq_len(dfs), sep = "_"),":X_E"),,drop=F]

  # all(rownames(coef(cv_obj, s = lambda_type))[-1] ==  colnames(cv_obj$sail.fit$design))
  # lin_pred <- cbind(1,cv_obj$sail.fit$design) %*% coef(cv_obj, s = lambda_type)

  ees <- unique(cv_obj$sail.fit$design[,"E"])
  unexposed_index <- which(cv_obj$sail.fit$design[,"E"]==ees[1])
  exposed_index <- which(cv_obj$sail.fit$design[,"E"]==ees[2])
  # browser()
  # 1=control, 2=MCI (Mild Cognitive Impairment) and 3=Alzeimer Disease
  # c(cont_index0, cont_index1, mci_index0, mci_index1, az_index0, az_index1) %>% length()

  e0 <- lin_pred[unexposed_index,]
  e1 <- lin_pred[exposed_index,]

  cont_pred0 <- lin_pred[cont_index0,]
  cont_pred1 <- lin_pred[cont_index1,]
  mci_pred0 <- lin_pred[mci_index0,]
  mci_pred1 <- lin_pred[mci_index1,]
  az_pred0 <- lin_pred[az_index0,]
  az_pred1 <- lin_pred[az_index1,]

  min.length.top <- range(lin_pred)[1] ; max.length.top <- range(lin_pred)[2]
  par(mai=c(1,1,1,0.2))
  plot(X[,sail_name], lin_pred,
       pch = 19,
       ylab = ylab,
       xlab = xlab,
       col = color[1],
       bty="n",
       xaxt="n",
       type = "n",
       cex.lab = 2,
       cex.axis = 2,
       cex = 2,
       main = main,
       cex.main = 2.5,
       ylim = c(min.length.top-3, max.length.top+3))#,
       # ylim = ylim)
  axis(1, labels = T, cex.axis = 2)


  x1 <- X[exposed_index,sail_name]
  x0 <- X[unexposed_index,sail_name]

  points(x1[order(x1)], e1[order(x1)], pch = 19, col = color[2], cex = 1.5, type = "l")
  points(x0[order(x0)], e0[order(x0)], pch = 19, col = color[1], cex = 1.5, type = "l")



  if (apoe) {
    points(X[cont_index1,original_name], cont_pred1, pch = 19, col = color[1], cex = 1.5)
    points(X[mci_index1,original_name], mci_pred1, pch = 19, col = color[2], cex = 1.5)
    points(X[az_index1,original_name], az_pred1, pch = 19, col = color[3], cex = 1.5)
    rug(X[exposed_index,original_name], side = 1)

  } else {
    points(X[cont_index0,original_name], cont_pred0, pch = 19, col = color[1], cex = 1.5)
    points(X[mci_index0,original_name], mci_pred0, pch = 19, col = color[2], cex = 1.5)
    points(X[az_index0,original_name], az_pred0, pch = 19, col = color[3], cex = 1.5)
    rug(X[unexposed_index,original_name], side = 1)
  }
  # if (legend) legend("bottomright", c("APOE = 1","APOE = 0"), col = color[2:1], pch = 19, cex = 2, bty = "n")
  # text(3, main, cex = 2)
  if (legend) legend("bottom", c("Control", "Mild Cognitive Impairment","Alzeimer Disease"),
                     col = color[1:3], pch = 19, cex = 2, bty = "n")

  rug(X[,sail_name], side = 1)
}










system.time(
  fit <- sail(x = X, y = Y, e = E, df = 3, degree = 3, thresh = 1e-4,
              maxit = 1000,
              alpha = .3,
              # dfmax = 15,
              verbose = TRUE, nlambda = 100)
)


plot(fit)
fit

registerDoMC(cores = 8)
foldid=sample(1:5,size=length(Y),replace=TRUE)
foldid %>% table
sample(rep(seq(10), length = length(Y)))

# cv1=cv.glmnet(x,y,foldid=foldid,alpha=1)
# cv.5=cv.glmnet(x,y,foldid=foldid,alpha=.5)
# cv0=cv.glmnet(x,y,foldid=foldid,alpha=0)

system.time(
  cvfit.9 <- cv.sail(x = X, y = Y, e = E, df = 3, degree = 3, thresh = 1e-3,
                   maxit = 1000,
                   # foldid = foldid,
                   alpha = .1,
                   parallel = TRUE,
                   nfolds = 3,
                   # dfmax = 15,
                   verbose = TRUE, nlambda = 100)
)

plot(cvfit.9)
coef(cvfit, s = "lambda.min")[nonzero(coef(cvfit, s = "lambda.min")),,drop=F]
coef(cvfit, s = "lambda.1se")[nonzero(coef(cvfit, s = "lambda.1se")),,drop=F]
xv <- "MMSE"
ind <- cvfit$sail.fit$group == which(cvfit$sail.fit$vnames == xv)
design.mat <- cvfit$sail.fit$design[,cvfit$sail.fit$main.effect.names[ind],drop = FALSE]
# f.truth <- design.mat %*% DT$b4
# f.truth <- DT$f2
plotMain(object = cvfit$sail.fit, xvar = xv, s = cvfit$lambda.min, legend.position = "topleft")


