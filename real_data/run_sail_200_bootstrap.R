setwd("~/Desktop/sail/real_data/")
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
set.seed(20000101)
seedpool <- sample(1:10000,200)
sailROC <- c()
sailOR <- c()
sailCoef <- c()
COEFlist <- list()
for (i in 1:200) {
  set.seed(seedpool[i])
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
}  

save(sailCoef,sailROC,sailOR,COEFlist,file="sailRes200Bootstrap.RData")
