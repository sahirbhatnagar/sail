pacman::p_load(splines)
pacman::p_load(magrittr)
pacman::p_load(foreach)
pacman::p_load(methods)
pacman::p_load(doMC)
pacman::p_load(glmnet)
pacman::p_load(LassoBacktracking)
pacman::p_load(glinternet)

# rm(list=ls())
# dev.off()
devtools::load_all("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/")

parameterIndex <- as.numeric(as.character(commandArgs(trailingOnly = T)[1]))
# parameterIndex = 1

if (parameterIndex == 1) { # 1a
  hierarchy = "strong" ; nonlinear = TRUE ; interactions = TRUE
} else if (parameterIndex == 2) { # 1b
  hierarchy = "weak" ; nonlinear = TRUE ; interactions = TRUE
} else if (parameterIndex == 3) { # 1c
  hierarchy = "none" ; nonlinear = TRUE ; interactions = TRUE
} else if (parameterIndex == 4) { # 2
  hierarchy = "strong"; nonlinear = FALSE; interactions = TRUE
} else if (parameterIndex == 5) { # 3
  hierarchy = "strong" ; nonlinear = TRUE ; interactions = FALSE
}

lambda.type <- "lambda.min"
# hierarchy = "strong", nonlinear = TRUE, interactions = TRUE, # scenario 1a
# hierarchy = "weak", nonlinear = TRUE, interactions = TRUE, # scenario 1b
# hierarchy = "none", nonlinear = TRUE, interactions = TRUE, # scenario 1c
# hierarchy = "strong", nonlinear = FALSE, interactions = TRUE, # scenario 2
# hierarchy = "strong", nonlinear = TRUE, interactions = FALSE, # scenario 3

# Simulate Data -----------------------------------------------------------

n = 200
p = 1000

DT <- gendataPaper(n = n, p = p, SNR = 2, betaE = 1,
                   hierarchy = hierarchy, nonlinear = nonlinear, interactions = interactions,
                   corr = 0,
                   E = truncnorm::rtruncnorm(n, a = -1, b = 1))

# used for glmnet and lasso backtracking
X <- design_sail(x = DT$x, e = DT$e, nvars = p,
                 vnames = paste0("X",1:p), degree = 1,
                 center.x = FALSE, basis.intercept = FALSE)$design

EX <- cbind(E = DT$e, DT$x)

Y.scaled <- scale(DT$y, center = TRUE, scale = FALSE)

DT_test <- gendataPaper(n = n, p = p, SNR = 2, betaE = 1,
                   hierarchy = hierarchy, nonlinear = nonlinear, interactions = interactions,
                   corr = 0,
                   E = truncnorm::rtruncnorm(n, a = -1, b = 1))

# used for glmnet and lasso backtracking
X_test <- design_sail(x = DT_test$x, e = DT_test$e, nvars = p,
                 vnames = paste0("X",1:p), degree = 1,
                 center.x = FALSE, basis.intercept = FALSE)$design

Y_test.scaled <- scale(DT_test$y, center = TRUE, scale = FALSE)

EX_test <- cbind(E = DT_test$e, DT_test$x)

# Lasso -------------------------------------------------------------------

cvfitlasso <- cv.glmnet(x = X, y = DT$y, alpha = 1, nfolds = 10)
plot(cvfitlasso)

# on the training
coef(cvfitlasso, s = cvfitlasso[["lambda.min"]])[nonzeroCoef(coef(cvfitlasso,s = cvfit[["lambda.min"]])),,drop=F]
coef(cvfitlasso, s = cvfitlasso[["lambda.1se"]])[nonzeroCoef(coef(cvfitlasso,s = cvfit[["lambda.1se"]])),,drop=F]
cvfitlasso$cvm[which(cvfitlasso[["lambda.min"]]==cvfitlasso$lambda)]
cor(predict(cvfitlasso, newx = X, s = cvfitlasso$lambda.min), DT$y)^2
crossprod(predict(cvfitlasso, newx = X, s = cvfitlasso$lambda.min) - DT$y)/n

# on the test
cor(predict(cvfitlasso, newx = X_test, s = cvfitlasso$lambda.min), DT_test$y)^2
crossprod(predict(cvfitlasso, newx = X_test, s = cvfitlasso$lambda.min) - DT_test$y) / n


# Glinternet --------------------------------------------------------------

cvglinternet <- glinternet.cv(X = EX, Y = DT$y,
                              numLevels = rep(1, ncol(EX)),
                              nLambda = 100, interactionCandidates = c(1),
                              verbose = T)

cvglinternet$activeSet
tc <- coef(cvglinternet)

colnames(EX)[tc$mainEffects$cont]


paste0(colnames(EX)[tc$interactions$contcont[,2]],":E")

# on the training
crossprod(predict(cvglinternet$glinternetFit, X = EX, lambda = cvglinternet$lambdaHat) - DT$y) / n
crossprod(predict(cvglinternet, X = EX) - DT$y) / n

plot(cvglinternet)

cvglinternet$cvErr[which(cvglinternet$lambdaHat==cvglinternet$lambda)]
cor(predict(cvglinternet, X = EX), DT$y)^2
registerDoMC(cores = 10)


numLevels = sample(1:5, 10, replace=TRUE)
X = sapply(numLevels, function(x) if (x==1) rnorm(200) else sample(0:(x-1), 200, replace=TRUE))
Y = X[,1] + X[,2] + X[,3]*X[,4] + rnorm(200)

fit = glinternet.cv(X, Y, numLevels, nFolds=5, nLambda = 50)
fit$activeSet
fit$cvErr[which(fit$lambdaHat==fit$lambda)]
crossprod(predict(fit, X = X) - Y) / n

# Lasso BT ----------------------------------------------------------------


out <- cvLassoBT(EX, DT$y, iter_max=1, nperms=1, nfolds = 10, mc.cores = 10, nlambda = 100)
out <- cvLassoBT(model@params$EX, draws@draws$r1.2, iter_max=10, nperms=1, nfolds = 10, mc.cores = 10, nlambda = 100)
coef(out$BT_fit)
out$cv_opt_err
out$cv_opt

coefBT <- as.matrix(predict(out$BT_fit, type = "coef", s = out$lambda[out$cv_opt[1]], iter = out$cv_opt[2]))

nz <- coefBT[nonzero(coefBT),,drop=F]
# as.numeric(as.character(rownames(nz)[-1]))

coefBT[,1] <- 0

nonzero(coefBT) %>% length()

yhat <- tryCatch({
  # tf <- as.matrix(predict(out$BT_fit, s = out$lambda[out$cv_opt[1]])[,out$cv_opt[1], out$cv_opt[2], drop=TRUE])
  preds <- predict(out$BT_fit, s = out$lambda[out$cv_opt[1]])
  tf <- as.matrix(preds[,out$cv_opt[1], min(dim(preds)[3],out$cv_opt[2]), drop=TRUE])
  a
},
error = function(err) {
  return(matrix(0, nrow = nrow(EX), ncol = 1))
} # return NULL on error
)

rm(yhat)

ints <- grep(":",rownames(nz)[-1], value=T)
if (length(ints) != 0) {
inters <- sapply(stringr::str_split(ints, ":"), function(i) paste(colnames(EX)[as.numeric(as.character(i))], collapse = ":"))
mains <- colnames(EX)[as.numeric(as.character(setdiff(rownames(nz)[-1], ints)))]
} else {
  mains <- colnames(EX)[as.numeric(as.character(rownames(nz)[-1]))]
}

mains

coef(out$BT_fit) %>% dim
predict(out$BT_fit, s = out$lambda[out$cv_opt[1]]) %>% dim
as.matrix(predict(out$BT_fit, s = out$lambda[out$cv_opt[1]])[,out$cv_opt[1], out$cv_opt[2], drop=F]) %>% length()
plot(predict(out$BT_fit, s = out$lambda[out$cv_opt[1]])[,out$cv_opt[1], out$cv_opt[2]], DT$y)
abline(a=0,b=1)
cor(predict(out$BT_fit)[,out$cv_opt[1], out$cv_opt[2]], DT$y)^2


# foldid <- sample(1:10,size=length(DT$y),replace=TRUE)

registerDoMC(cores = 10)

cvfit <- cv.sail(x = DT$x, y = DT$y, e = DT$e, df = 5, degree = 3, basis.intercept = FALSE,
                 thresh = 1e-4,
                 maxit = 1000,
                 # group.penalty = "SCAD",
                 alpha = .2,
                 parallel = TRUE,
                 # foldid = foldid,
                 nfolds = 10, verbose = T, nlambda = 100)
plot(cvfit)
cvfit$cvm[which(cvfit$lambda.min==cvfit$lambda)]
cor(predict(cvfit, s = cvfit$lambda.min), DT$y)^2
coef(cvfit, s = "lambda.min")[nonzero(coef(cvfit, s = "lambda.min")),,drop=F]
coef(cvfit, s = "lambda.1se")[nonzero(coef(cvfit, s = "lambda.1se")),,drop=F]

# par(mfrow=c(2,2))
# for (i in 1:4){
#   xv <- paste0("X",i)
#   ind <- cvfit$sail.fit$group == which(cvfit$sail.fit$vnames == xv)
#   # design.mat <- cvfit$sail.fit$design[,cvfit$sail.fit$main.effect.names[ind],drop = FALSE]
#   # f.truth <- design.mat %*% DT$b1
#   f.truth <- DT[[paste0("f",i)]]
#   plotMain(object = cvfit$sail.fit, xvar = xv, s = cvfit$lambda.min, f.truth = f.truth, legend.position = "topleft")
# }
#
#
# f3.persp = function(X, E) {
#   # E * as.vector(DT$betaE) + DT$f3.f(X) +
#   E * DT$f3.f(X)
# }
#
# f4.persp = function(X, E) {
#   # E * as.vector(DT$betaE) + DT$f4.f(X) +
#   E * DT$f4.f(X)
# }
#
# i = 3
# xv <- paste0("X",i)
# plotInter(object = cvfit$sail.fit, xvar = xv, s = cvfit$lambda.min,
#           f.truth = f3.persp,
#           simulation = T,
#           npoints = 40)
#
# i = 4
# xv <- paste0("X",4)
# plotInter(object = cvfit$sail.fit, xvar = xv, s = cvfit$lambda.min,
#           f.truth = f4.persp,
#           simulation = T,
#           npoints = 40)



