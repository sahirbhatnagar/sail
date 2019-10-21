#' This script runs the 200 bootstraps for the support analysis
#' Results are saved in a .RData file and then subsequently used
#' in support_plots.R. The code in support_plots.R generates
#' the figures which compare the methods in terms of variable selection and
#' AUC


library(sail)
library(splines)
library(glmnet)
library(LassoBacktracking)
library(HierBasis)
library(glinternet)
library(pROC)

dat <- read.csv("support2.csv", header = T, row.names = 1)

library(dplyr)

dat %>%
  select(age, sex, dzclass, num.co, diabetes,
         dementia, meanbp, wblc, hrt, resp, temp, crea, sod, adlsc) -> usedatcomplete

usedatcomplete$survive6M <- ifelse(dat$death == 0, 1,
                                   ifelse(dat$d.time >= 180, 1, -1))

usedatcomplete <- na.omit(usedatcomplete)

# ARF/MOSF: 1
# Cancer: 2
# Coma: 3
# Coma COPD/CHF/Cirrhosis: 4
usedatcomplete$dzclass <- as.numeric(usedatcomplete$dzclass)
usedatcomplete$dzclass <- ifelse(usedatcomplete$dzclass==1,1,0)
usedatcomplete$sex <- ifelse(usedatcomplete$sex=="female",0,1)

cat("done pre-processing", collapse="\n")

# training-validation-test split
set.seed(8451)
sailROC <- c()
lassoROC <- c()
lassoBTROC <- c()
HBROC <- c()
GLROC <- c()
sailOR <- c()
lassoOR <- c()
lassoBTOR <- c()
HBOR <- c()
GLOR <- c()
sailCoef <- c()
lassoCoef <- c()
lassoBTCoef <- c()
HBCoef <- c()
GLCoef <- c()
COEFlist <- list()
for (i in 1:200) {

  cat("iteration: ",i, collapse="\n")
  # split dataset
  trainID <- sample(1:nrow(usedatcomplete),nrow(usedatcomplete)*0.34)
  valID <- sample((1:nrow(usedatcomplete))[! 1:nrow(usedatcomplete) %in% trainID],nrow(usedatcomplete)*0.33)
  testID <- setdiff(setdiff((1:nrow(usedatcomplete)),trainID),valID)
  traincomplete <- usedatcomplete[trainID,]
  valcomplete <- usedatcomplete[valID,]
  testcomplete <- usedatcomplete[testID,]

  # sail
  cat("iteration: ",i," - sail", collapse="\n")
  designmattrain <- model.matrix(~ 0
                                 + bs(age, degree = 3)
                                 + sex
                                 + bs(num.co, degree = 3)
                                 + diabetes
                                 + dementia
                                 + bs(meanbp, degree = 3)
                                 + bs(wblc, degree = 3)
                                 + bs(hrt, degree = 3)
                                 + bs(resp, degree = 3)
                                 + bs(temp, degree = 3)
                                 + bs(crea, degree = 3)
                                 + bs(sod, degree = 3)
                                 + bs(adlsc, degree = 3), data = traincomplete)
  sailfittrain <- sail(x = designmattrain,
                       y = traincomplete$survive6M,
                       e = traincomplete$dzclass,
                       expand = FALSE,
                       group = attr(designmattrain, "assign"),
                       verbose = 0,
                       alpha = 0.1,
                       strong = FALSE)
  designmatval <- model.matrix(~ 0
                               + bs(age, degree = 3)
                               + sex
                               + bs(num.co, degree = 3)
                               + diabetes
                               + dementia
                               + bs(meanbp, degree = 3)
                               + bs(wblc, degree = 3)
                               + bs(hrt, degree = 3)
                               + bs(resp, degree = 3)
                               + bs(temp, degree = 3)
                               + bs(crea, degree = 3)
                               + bs(sod, degree = 3)
                               + bs(adlsc, degree = 3), data = valcomplete)
  valscore <- predict(sailfittrain, newx = designmatval, newe = valcomplete$dzclass)
  apply(valscore[,-1],2,function(x) exp(coef(glm(valcomplete$survive6M~scale(x)))[2])) -> valOR
#  sum(as.vector(coef(sailfittrain, s = sailfittrain$lambda[which.max(valOR)+1])) != 0) - 1 -> Ncoef
  as.matrix(coef(sailfittrain, s = sailfittrain$lambda[which.max(valOR)+1])) -> COEF
  COEFlist[[i]] <- COEF
  names(COEF[COEF[,1]!=0,]) -> Ncoef
  length(Ncoef) - length(grep("bs",Ncoef)) - 1 + length(grep("bs",Ncoef))/3 -> Ncoef
  sailCoef <- c(sailCoef, Ncoef)
  designmattest <- model.matrix(~ 0
                                + bs(age, degree = 3)
                                + sex
                                + bs(num.co, degree = 3)
                                + diabetes
                                + dementia
                                + bs(meanbp, degree = 3)
                                + bs(wblc, degree = 3)
                                + bs(hrt, degree = 3)
                                + bs(resp, degree = 3)
                                + bs(temp, degree = 3)
                                + bs(crea, degree = 3)
                                + bs(sod, degree = 3)
                                + bs(adlsc, degree = 3), data = testcomplete)
  pROC::auc(testcomplete$survive6M, as.vector(predict(sailfittrain, newx = designmattest, newe = testcomplete$dzclass, s = sailfittrain$lambda[which.max(valOR)+1])))[1] -> testROC
  sailROC <- c(sailROC, testROC)
  exp(coef(glm(testcomplete$survive6M~scale(as.vector(predict(sailfittrain, newx = designmattest, newe = testcomplete$dzclass, s = sailfittrain$lambda[which.max(valOR)+1])))))[2]) -> testOR
  sailOR <- c(sailOR, testOR)

  # Lasso
  cat("iteration: ",i," - lasso", collapse="\n")
  lassofittrain <- glmnet(x = as.matrix(traincomplete[,1:14]), y = traincomplete$survive6M, family = "binomial", alpha = 1)
  lassovalscore <- predict(lassofittrain, newx = as.matrix(valcomplete[,1:14]))
  apply(lassovalscore[,-1],2,function(x) exp(coef(glm(valcomplete$survive6M~scale(x)))[2])) -> valOR
  sum(as.vector(coef(lassofittrain, s = lassofittrain$lambda[which.max(valOR)+1])) != 0) - 1 -> Ncoef
  lassoCoef <- c(lassoCoef, Ncoef)
  pROC::auc(testcomplete$survive6M, as.vector(predict(lassofittrain, newx = as.matrix(testcomplete[,1:14]), s = lassofittrain$lambda[which.max(valOR)+1])))[1] -> testROC
  lassoROC <- c(lassoROC, testROC)
  exp(coef(glm(testcomplete$survive6M~scale(as.vector(predict(lassofittrain, newx = as.matrix(testcomplete[,1:14]), s = lassofittrain$lambda[which.max(valOR)+1])))))[2]) -> testOR
  lassoOR <- c(lassoOR, testOR)

  # LassoBT
  cat("iteration: ",i," - lassoBT", collapse="\n")
  lassoBTfittrain <- LassoBT(x = as.matrix(traincomplete[,1:14]), y = traincomplete$survive6M, iter_max = 20)
  lassoBTvalscore <- predict(lassoBTfittrain, newx = as.matrix(valcomplete[,1:14]))
  apply(lassoBTvalscore,3,function(x) apply(x[,-1],2,function(xcol) exp(coef(glm(valcomplete$survive6M ~ scale(xcol)))[2]))) -> lassoBTORrun
  which.max(apply(lassoBTORrun,2,max)) -> iterID
  lassoBTvalscore <- predict(lassoBTfittrain, newx = as.matrix(valcomplete[,1:14]))[,,iterID]
  apply(lassoBTvalscore[,-1],2,function(x) exp(coef(glm(valcomplete$survive6M~scale(x)))[2])) -> valOR
  sum(as.vector(predict(lassoBTfittrain,type = "coefficients",s = lassoBTfittrain$lambda[which.max(valOR)+1], iter = iterID)) != 0) -1 -> Ncoef
  lassoBTCoef <- c(lassoBTCoef, Ncoef)
  lassoBTtestscore <- predict(lassoBTfittrain, newx = as.matrix(testcomplete[,1:14]))[,,iterID]
  pROC::auc(testcomplete$survive6M, as.vector(lassoBTtestscore[,which.max(valOR)+1]))[1] -> testROC
  lassoBTROC <- c(lassoBTROC, testROC)
  exp(coef(glm(testcomplete$survive6M~scale(as.vector(lassoBTtestscore[,which.max(valOR)+1]))))[2]) -> testOR
  lassoBTOR <- c(lassoBTOR, testOR)

  # HierBasis
  cat("iteration: ",i," - HierBasis", collapse="\n")
  HBfittrain <- AdditiveHierBasis(x = as.matrix(traincomplete[,1:14]), y = ifelse(traincomplete$survive6M==1,1,0), type = "binomial", nlam = 100, max.iter = 100)
  HBvalscore <- predict(HBfittrain, new.x = as.matrix(valcomplete[,1:14]))
  apply(HBvalscore[,apply(HBvalscore,2,function(x) length(unique(x)))!=1],2,function(x) exp(coef(glm(valcomplete$survive6M~scale(x)))[2])) -> valOR
  sum(HBfittrain$beta[,apply(HBvalscore,2,function(x) length(unique(x)))!=1][,which.max(valOR)] != 0) -> Ncoef
  HBCoef <- c(HBCoef, Ncoef)
  HBtestscore <- predict(HBfittrain, new.x = as.matrix(testcomplete[,1:14]))
  pROC::auc(testcomplete$survive6M, HBtestscore[,apply(HBvalscore,2,function(x) length(unique(x)))!=1][,which.max(valOR)])[1] -> testROC
  HBROC <- c(HBROC, testROC)
  exp(coef(glm(testcomplete$survive6M~scale(HBtestscore[,apply(HBvalscore,2,function(x) length(unique(x)))!=1][,which.max(valOR)])))[2]) -> testOR
  HBOR <- c(HBOR, testOR)

  # GLinternet
  cat("iteration: ",i," - GLinternet", collapse="\n")
  glinternetfittrain <- glinternet(X = as.matrix(traincomplete[,1:14]), Y = ifelse(traincomplete$survive6M==1,1,0), numLevels = c(1,2,2,1,2,2,1,1,1,1,1,1,1,1), nLambda = 100, interactionCandidates = c(3), family = "binomial")
  glinternetvalscore <- predict(glinternetfittrain, X = as.matrix(valcomplete[,1:14]))
  apply(glinternetvalscore[,-1],2,function(x) exp(coef(glm(valcomplete$survive6M~scale(x)))[2])) -> valOR
  sum(glinternetfittrain$betahat[[which.max(valOR) + 1]] != 0) -> Ncoef
  GLCoef <- c(GLCoef, Ncoef)
  pROC::auc(testcomplete$survive6M, as.vector(predict(glinternetfittrain, X = as.matrix(testcomplete[,1:14]), lambda = glinternetfittrain$lambda[which.max(valOR)+1])))[1] -> testROC
  GLROC <- c(GLROC, testROC)
  exp(coef(glm(testcomplete$survive6M~scale(as.vector(predict(glinternetfittrain, X = as.matrix(testcomplete[,1:14]), lambda = glinternetfittrain$lambda[which.max(valOR)+1])))))[2]) -> testOR
  GLOR <- c(GLOR, testOR)
}

# split dataset
trainID <- sample(1:nrow(usedatcomplete),nrow(usedatcomplete)*0.34)
valID <- sample((1:nrow(usedatcomplete))[! 1:nrow(usedatcomplete) %in% trainID],nrow(usedatcomplete)*0.33)
testID <- setdiff(setdiff((1:nrow(usedatcomplete)),trainID),valID)
traincomplete <- usedatcomplete[trainID,]
valcomplete <- usedatcomplete[valID,]
testcomplete <- usedatcomplete[testID,]

# sail
designmattrain <- model.matrix(~ 0
                               + bs(age, degree = 3)
                               + sex
                               + bs(num.co, degree = 3)
                               + diabetes
                               + dementia
                               + bs(meanbp, degree = 3)
                               + bs(wblc, degree = 3)
                               + bs(hrt, degree = 3)
                               + bs(resp, degree = 3)
                               + bs(temp, degree = 3)
                               + bs(crea, degree = 3)
                               + bs(sod, degree = 3)
                               + bs(adlsc, degree = 3), data = traincomplete)
sailfittrain <- sail(x = designmattrain,
                     y = traincomplete$survive6M,
                     e = traincomplete$dzclass,
                     expand = FALSE,
                     group = attr(designmattrain, "assign"),
                     verbose = 0,
                     alpha = 0.1,
                     strong = FALSE)
designmatval <- model.matrix(~ 0
                             + bs(age, degree = 3)
                             + sex
                             + bs(num.co, degree = 3)
                             + diabetes
                             + dementia
                             + bs(meanbp, degree = 3)
                             + bs(wblc, degree = 3)
                             + bs(hrt, degree = 3)
                             + bs(resp, degree = 3)
                             + bs(temp, degree = 3)
                             + bs(crea, degree = 3)
                             + bs(sod, degree = 3)
                             + bs(adlsc, degree = 3), data = valcomplete)
valscore <- predict(sailfittrain, newx = designmatval, newe = valcomplete$dzclass)
apply(valscore[,-1],2,function(x) exp(coef(glm(valcomplete$survive6M~scale(x)))[2])) -> valOR
#  sum(as.vector(coef(sailfittrain, s = sailfittrain$lambda[which.max(valOR)+1])) != 0) - 1 -> Ncoef
sailfittrain$lambda[which.max(valOR)+1] -> optLambda
designmattrain -> designmat2

save(designmat2,usedatcomplete,sailfittrain,optLambda,COEFlist,file="support_sail200Bootstrap.RData")

permRes <- data.frame(OR=c(sailOR,lassoOR,lassoBTOR,HBOR,GLOR),
                      AUC=c(sailROC,lassoROC,lassoBTROC,HBROC,GLROC),
                      Ncoef=c(sailCoef,lassoCoef,lassoBTCoef,HBCoef,GLCoef),
                      Method=rep(c("sail","lasso","lassoBT","HierBasis","GLinternet"),each=200))
describeBy(permRes[,c("AUC","Ncoef")], permRes$Method, mat=TRUE) -> errorSummary

save(errorSummary, file="support_methods_comparison_error_summary.RData")


