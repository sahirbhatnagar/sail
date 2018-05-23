rm(list=ls())
library(data.table)
library(magrittr)
devtools::load_all()

DT <- fread("/home/sahir/git_repositories/sail/data-nogit/sebastien/dataSahir.txt")

setkey(DT, "IID_new")
key(DT)
DT[, table(CarrierStatus)]

DT %>% colnames() %>% grep("^sum", ., value = TRUE)

DT[, `:=`(sum_Size_Ok.COMB = sum(sum_Size_Ok.DEL, sum_Size_Ok.DUP, na.rm = TRUE),
          sum_nb_genes_prot.COMB = sum(sum_nb_genes_prot.DEL + sum_nb_genes_prot.DUP, na.rm = TRUE),
          sum_compte_PSD.COMB = sum(sum_compte_PSD.DEL + sum_compte_PSD.DUP, na.rm = TRUE),
          sum_pLI.COMB = sum(sum_pLI.DEL + sum_pLI.DUP, na.rm = TRUE),
          sum_del.score_3.COMB = sum(sum_del.score_3.DEL + sum_dup.score_3.DUP, na.rm = TRUE),
          sum_ExAC_005.popn.COMB = sum(sum_ExAC_005.popn.DEL + sum_ExAC_005.popn.DUP, na.rm = TRUE),
          sum_PPI.COMB = sum(sum_PPI.DEL + sum_PPI.DUP, na.rm = TRUE),
          sum_eQTL.COMB = sum(sum_eQTL.DEL + sum_eQTL.DUP, na.rm = TRUE),
          sum_Pearson_Min.COMB = sum(sum_Pearson_Min.DEL + sum_Pearson_Min.DUP, na.rm = TRUE),
          sum_compte_FMRP.COMB = sum(sum_compte_FMRP.DEL + sum_compte_FMRP.DUP, na.rm = TRUE)), by = IID_new]
DT[, Enum := ifelse(CarrierStatus=="no", 0, ifelse(CarrierStatus=="del", 1, ifelse(CarrierStatus=="dup", 2, 3)))]
DT[, table(Enum, CarrierStatus)]

DT %>% colnames() %>% grep("COMB", ., value = TRUE)

DT[, DT %>% colnames() %>% grep("COMB", ., value = TRUE) , with=F] %>% complete.cases() %>% sum

DT_cc_PIQ <- DT[!is.na(Combine_PIQ)]
DT_cc_PIQ %>% colnames() %>% grep("COMB", ., value = TRUE) %>% dput
x_vars <- c(DT_cc_PIQ %>% colnames() %>% grep("COMB", ., value = TRUE),
            paste0("C",1:6),"Age")
x_vars <- c(#"sum_Size_Ok.COMB",
            "sum_nb_genes_prot.COMB",
            "sum_compte_PSD.COMB",
            "sum_pLI.COMB",
            "sum_del.score_3.COMB",
            # "sum_ExAC_005.popn.COMB",
            # "sum_PPI.COMB",
            "sum_eQTL.COMB",
            "sum_Pearson_Min.COMB",
            "sum_compte_FMRP.COMB",paste0("C",1:6),"Age")
DT$sum_ExAC_005.popn.COMB %>% table

X <- DT_cc_PIQ[,x_vars, with=F] %>% as.matrix()
X <- standardize(X, center = TRUE, normalize = TRUE)$x
E <- DT_cc_PIQ$Enum
E <- standardize(E, center = TRUE, normalize = TRUE)$x
Y <- DT_cc_PIQ$Combine_PIQ
any(is.na(Y))
hist(Y)
DataExplorer::plot_histogram(X)


f.basis <- function(x) splines::bs(x)
fit <- sail(x = X, y = Y, e = E, basis = f.basis, center.e = FALSE,
            verbose = 2, strong = FALSE, center.x = FALSE, thresh = 5e-03,
            fdev = 1e-10)
coef(fit)
plot(fit)
fit
DataExplorer::plot_histogram(X)

library(doMC)
registerDoMC(cores = 8)
fit <- cv.sail(x = X, y = Y, e = E, basis = f.basis, center.e = FALSE,
            verbose = 2, strong = FALSE, center.x = FALSE, thresh = 5e-03,
            fdev = 1e-10, nfolds = 10, parallel = TRUE)
plot(fit)
predict(fit, type="non", s = "lambda.min")
DT <- fread("~/Downloads/abalone/Dataset.data")
setnames(DT, c("Sex","Length","Diameter", "Height", "Whole_weight", "Shucked_weight",
               "Viscera_weight","Shell_weight","Rings"))
hist(DT$Rings)

DT[, Gender := ifelse(Sex=="M", 0, ifelse(Sex=="F", 1, 2))]
DT[, table(Gender, Sex)]

train <- caret::createDataPartition(DT$Rings)[[1]]
# train <- seq(nrow(DT))
test <- setdiff(seq(nrow(DT)), train)

X <- DT[train, c("Length","Diameter", "Height", "Whole_weight", "Shucked_weight",
           "Viscera_weight","Shell_weight"), with = F] %>% as.matrix()
X_test <- DT[test, c("Length","Diameter", "Height", "Whole_weight", "Shucked_weight",
                      "Viscera_weight","Shell_weight"), with = F] %>% as.matrix()
E <- DT[train]$Gender
E_test <- DT[test]$Gender
Y <- DT[train]$Rings
Y_test <- DT[test]$Rings

f.basis <- function(x) splines::bs(x, degree = 3)
fit <- sail(x = X, y = Y, e = E, basis = f.basis,
            # center.e = FALSE,
            verbose = 2, strong = FALSE)
plot(fit)

colMeans((Y_test - predict(fit, newx = X_test, newe = E_test))^2)


library(doMC)
registerDoMC(cores = 8)
plot(fit)
coef(fit)
cvfit <- cv.sail(x = X, y = Y, e = E, basis = f.basis,
                 # center.e = FALSE,
                 nfolds = 10, parallel = TRUE,
                 verbose = 2, strong = TRUE)
plot(cvfit)
dev.off()
