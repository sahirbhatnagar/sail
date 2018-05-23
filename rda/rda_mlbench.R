rm(list=ls())
pacman::p_load(mlbench)
# devtools::document()
devtools::load_all()
# source("R/ksfilter.R")
pacman::p_load(data.table)
pacman::p_load(tidyr)
pacman::p_load(dplyr)
pacman::p_load(magrittr)
pacman::p_load(tidyverse)
pacman::p_load(doParallel)
pacman::p_load(splines)
pacman::p_load(foreach)
pacman::p_load(doMC)
pacman::p_load(glmnet)
pacman::p_load_current_gh('sahirbhatnagar/glinternet')
pacman::p_load(DataExplorer)
pacman::p_load(LassoBacktracking)


pacman::p_data("mlbench")
data("BostonHousing")
data("BostonHousing2")
?BostonHousing2

Y <- BostonHousing2$cmedv
hist(Y)
Xorig <- BostonHousing2 %>% select(-town, -medv, -cmedv, -chas) %>% as.matrix
X <- sail:::standardize(Xorig, center = TRUE, normalize = TRUE)$x
colnames(X)
Xc <- X[, -which(colnames(X)=="ptratio")]
E <- X[, which(colnames(X)=="ptratio")]
# E <- BostonHousing2 %>% pull(zn) %>% as.numeric
# E <- BostonHousing2 %>% pull(chas) %>% as.numeric
# E <- BostonHousing2 %>% pull(age) %>% as.numeric # dont use age as E
DataExplorer::plot_histogram(X)
BostonHousing2 %>% str

f.basis <- function(i) splines::bs(i, degree = 3)
fit <- sail(x = Xc, y = Y, e = E, basis = f.basis, verbose = 2,
            thresh = 5e-3, maxit = 250,
            strong = FALSE)
plot(fit)


registerDoMC(cores = 10)
set.seed(123456)
cvfit <- cv.sail(x = Xc, y = Y, e = E, basis = f.basis, verbose = 2,
            thresh = 5e-3,
            # alpha = 0.3,
            maxit = 250,
            strong = FALSE, parallel = TRUE, nfolds = 10,
            center.x = T, center.e = T)
plot(cvfit)
cvfit03=cvfit
# dev.off()
predict(cvfit, type="non", s = "lambda.min")
predict(cvfit, type="non", s = "lambda.1se")
sprintf("%.2f (%.2f, %.2f)",
        cvfit03$cvm[which(cvfit03$lambda.min==cvfit03$lambda)],
        cvfit03$cvlo[which(cvfit03$lambda.min==cvfit03$lambda)],
        cvfit03$cvup[which(cvfit03$lambda.min==cvfit03$lambda)])

cor(predict(cvfit03, s="lambda.min"), Y)^2


set.seed(123456)
which(colnames(X)=="ptratio")
cvglinternet <- glinternet.cv(X = Xorig, Y = Y, nFolds = 10,
                              numLevels = rep(1, ncol(X)),
                              family = "gaussian",
                              interactionCandidates = which(colnames(Xorig)=="ptratio"),
                              verbose = T)
tc <- coef(cvglinternet, lambdaType = "lambdaHat")


head(X)
(mains <- colnames(Xorig)[tc$mainEffects$cont])
(inters <- paste(colnames(Xorig)[tc$interactions$contcont[,1]],colnames(Xorig)[tc$interactions$contcont[,2]], sep=":"))
(cvmse <- cvglinternet$cvErr[which(cvglinternet$lambdaHat==cvglinternet$lambda)])
plot(cvglinternet)


set.seed(123456)
cvfitlasso <- cv.glmnet(Xorig, Y, alpha = 1, nfolds = 10)
plot(cvfitlasso)
# coef(cvfitlasso, s = "lambda.min")
sprintf("%.2f (%.2f, %.2f)",
        cvfitlasso$cvm[which(cvfitlasso$lambda.min==cvfitlasso$lambda)],
        cvfitlasso$cvlo[which(cvfitlasso$lambda.min==cvfitlasso$lambda)],
        cvfitlasso$cvup[which(cvfitlasso$lambda.min==cvfitlasso$lambda)])
