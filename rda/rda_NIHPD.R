rm(list=ls())

pacman::p_load(data.table)
pacman::p_load(dplyr)
pacman::p_load(zoo)
`%ni%` = Negate(`%in%`)
source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_git_v2/sail/rda/rda_NIHPD_data_cleaning.R")
#
DT <- nihdata(nprobes = 100, phenoVariable = "WASI_Full_Scale_IQ", exposure = "E", filter = "pvalue", data = "AAL")

# topn <- 5000 # how many probes (filtered based on highest sd)
# # depends on if youre using 80k cortical thickness.. or 80
# phenotypeVariable <- "WASI_Full_Scale_IQ"
# # phenotypeVariable <- "WASI_Verbal_IQ"
# # phenotypeVariable <- "WASI_Performance_IQ"
# # exposureVariable <- "E2"
# exposureVariable <- "E" # binary age
# # exposureVariable <- "gender_num"
# # exposureVariable <- "age"
# # exposureVariable <- "income_binary2"
# filter_var <- "sd"
# dataset <- "notAAL"

# test methods ------------------------------------------------------------

devtools::load_all("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_git_v2/sail/")
pacman::p_load(doParallel)
pacman::p_load(splines)
pacman::p_load(foreach)
pacman::p_load(doMC)
pacman::p_load(glmnet)
pacman::p_load(LassoBacktracking)
pacman::p_load(glinternet)
pacman::p_load_gh('asadharis/HierBasis')
pacman::p_load(SAM)
pacman::p_load(gamsel)



# simulator style ---------------------------------------------------------

pacman::p_load(simulator) # this file was created under simulator version 0.2.1
source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_git_v2/sail/my_sims/model_functions.R")
source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_git_v2/sail/my_sims/method_functions_rda.R")
source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_git_v2/sail/my_sims/eval_functions.R")


sim <- new_simulation(name = "rda_may_11_2018_p2k",
                      label = "rda_may_11_2018_p2k",
                      dir = ".") %>%
  generate_model(make_nihpd_data_split, seed = 12345,
                 nprobes = 2000,
                 # phenoVariable = list("WASI_Full_Scale_IQ","WASI_Verbal_IQ","WASI_Performance_IQ"),
                 phenoVariable = "WASI_Full_Scale_IQ",
                 exposure = "E",
                 filter = "pvalue",
                 data = "notAAL"
                 # vary_along = "phenoVariable",
  ) %>%
  simulate_from_model(nsim = 5, index = 1:20) %>%
  run_method(list(lassosplit,lassosplitadaptive, sailsplitweak, #lassoBTsplit, sailsplitweak, sailsplitadaptiveweak
                  sailsplit, sailsplitadaptive),
             parallel = list(socket_names = 20,
                             libraries = c("LassoBacktracking", "glinternet","glmnet","splines",
                                           "magrittr","sail","gamsel","SAM","HierBasis","simulator", "parallel")))
simulator::save_simulation(sim)
sim

sim <- load_simulation("rda_may_11_2018_p2k")
# rm(sim)
# sim <- load_simulation("rda_may_11_2018_AAL")
sim <- sim %>% evaluate(list(msevalid, nactive, r2))
sim %>% plot_eval(metric_name = "mse")
sim %>% plot_eval(metric_name = "nactive")

sim <- sim %>% run_method(list(lassoBTsplit, Hiersplit, SPAMsplit, gamselsplit), #GLinternetsplit, ),
                   parallel = list(socket_names = 20,
                                   libraries = c("LassoBacktracking", "glinternet","glmnet","splines",
                                                 "magrittr","sail","gamsel","SAM","HierBasis","simulator", "parallel")))

sim <- new_simulation(name = "rda_may_11_2018_p1k",
                      label = "rda_may_11_2018_p1k",
                      dir = ".") %>%
  generate_model(make_nihpd_data_split, seed = 12345,
                 nprobes = 1000,
                 # phenoVariable = list("WASI_Full_Scale_IQ","WASI_Verbal_IQ","WASI_Performance_IQ"),
                 phenoVariable = "WASI_Full_Scale_IQ",
                 exposure = "E",
                 filter = "pvalue",
                 data = "notAAL"
                 # vary_along = "phenoVariable",
  ) %>%
  simulate_from_model(nsim = 5, index = 1:20) %>%
  run_method(list(lassosplit,lassosplitadaptive, lassoBTsplit,
                  sailsplit, sailsplitadaptive, sailsplitweak, sailsplitadaptiveweak),
             parallel = list(socket_names = 20,
                             libraries = c("LassoBacktracking", "glinternet","glmnet","splines",
                                           "magrittr","sail","gamsel","SAM","HierBasis","simulator", "parallel")))
simulator::save_simulation(sim)
sim

sim <- load_simulation("rda_may_11_2018_p1k")
sim <- sim %>% evaluate(list(msevalid, nactive, r2))
sim %>% plot_eval(metric_name = "mse")
sim %>% plot_eval(metric_name = "nactive")


sim <- new_simulation(name = "rda_may_11_2018_AAL",
                      label = "rda_may_11_2018_AAL",
                      dir = ".") %>%
  generate_model(make_nihpd_data_split, seed = 12345,
                 nprobes = 1000,
                 # phenoVariable = list("WASI_Full_Scale_IQ","WASI_Verbal_IQ","WASI_Performance_IQ"),
                 phenoVariable = "WASI_Full_Scale_IQ",
                 exposure = "E",
                 filter = "pvalue",
                 data = "AAL"
                 # vary_along = "phenoVariable",
  ) %>%
  simulate_from_model(nsim = 5, index = 1:20) %>%
  run_method(list(lassosplit,lassosplitadaptive, lassoBTsplit,
                  sailsplit, sailsplitadaptive, sailsplitweak, sailsplitadaptiveweak),
             parallel = list(socket_names = 20,
                             libraries = c("LassoBacktracking", "glinternet","glmnet","splines",
                                           "magrittr","sail","gamsel","SAM","HierBasis","simulator", "parallel")))
simulator::save_simulation(sim)
sim


sim <- new_simulation(name = "rda_may_11_2018_p200",
                      label = "rda_may_11_2018_p200",
                      dir = ".") %>%
  generate_model(make_nihpd_data_split, seed = 12345,
                 nprobes = 200,
                 # phenoVariable = list("WASI_Full_Scale_IQ","WASI_Verbal_IQ","WASI_Performance_IQ"),
                 phenoVariable = "WASI_Full_Scale_IQ",
                 exposure = "E",
                 filter = "pvalue",
                 data = "notAAL"
                 # vary_along = "phenoVariable",
  ) %>%
  simulate_from_model(nsim = 1, index = 1:2) %>%
  run_method(list(lassosplit,lassosplitadaptive, sailsplit, sailsplitadaptive, sailsplitweak))#,
# parallel = list(socket_names = 10,
# libraries = c("LassoBacktracking", "glinternet","glmnet","splines",
# "magrittr","sail","gamsel","SAM","HierBasis","simulator", "parallel")))



sim <- load_simulation("testrda_may_11_2018")

sim <- sim %>% run_method(list(sailsplitadaptive))
# ,
# parallel = list(socket_names = 20,
#                 libraries = c("LassoBacktracking", "glinternet","glmnet","splines",
#                               "magrittr","sail","gamsel","SAM","HierBasis","simulator", "parallel")))

sim <- sim %>% evaluate(list(msevalid, nactive))
sim %>% plot_eval(metric_name = "mse")


# testing glinternet

sim <- new_simulation(name = "rda_may_13_2018_AAL_test_glinternet",
                      label = "rda_may_13_2018_AAL_test_glinternet",
                      dir = ".") %>%
  generate_model(make_nihpd_data_split, seed = 12345,
                 nprobes = 200,
                 # phenoVariable = list("WASI_Full_Scale_IQ","WASI_Verbal_IQ","WASI_Performance_IQ"),
                 phenoVariable = "WASI_Full_Scale_IQ",
                 exposure = "E",
                 filter = "pvalue",
                 data = "notAAL"
                 # vary_along = "phenoVariable",
  ) %>%
  simulate_from_model(nsim = 1, index = 1:2) %>%
  run_method(list(lassosplit,GLinternetsplit, sailsplit))


sim <- sim %>% evaluate(list(msevalid, nactive))
sim %>% plot_eval(metric_name = "mse")












# lasso -------------------------------------------------------------------

# fitlasso <- glmnet(x = xtrain_lasso, y = ytrain, alpha = 1)
# # plot(fitlasso)
# ytest_hat <- predict(fitlasso, newx = xtest_lasso)
# msetest <- colMeans((ytest - ytest_hat)^2)
# lambda.min.index <- as.numeric(which.min(msetest))
# lambda.min <- fitlasso$lambda[which.min(msetest)]
#
# yvalid_hat <- predict(fitlasso, newx = xvalid_lasso, s = lambda.min)
# (msevalid <- mean((yvalid - drop(yvalid_hat))^2))
#
# (nzcoef <- coef(fitlasso, s = lambda.min)[nonzeroCoef(coef(fitlasso, s = lambda.min)),,drop=F])
# nrow(nzcoef)

train_test_ind <- caret::createDataPartition(DT$ytrain, p = 200/338)[[1]]
validate_ind <- seq(length(DT$ytrain))[-train_test_ind]
train_ind <- sample(train_test_ind, floor(length(train_test_ind)/2))
test_ind <- setdiff(train_test_ind, train_ind)

c(train_ind, test_ind, validate_ind) %>% duplicated() %>% any


fitlasso <- glmnet(x = DT$xtrain_lasso[train_ind,], y = DT$ytrain[train_ind], alpha = 1)
ytest_hat <- predict(fitlasso, newx = DT$xtrain_lasso[test_ind,])
msetest <- colMeans((DT$ytrain[test_ind] - ytest_hat)^2)
lambda.min.index <- as.numeric(which.min(msetest))
lambda.min <- fitlasso$lambda[which.min(msetest)]

yvalid_hat <- predict(fitlasso, newx = DT$xtrain_lasso[validate_ind,], s = lambda.min)
(msevalid <- mean((DT$ytrain[validate_ind] - drop(yvalid_hat))^2))
cor(predict(fitlasso, s=lambda.min, newx = DT$xtrain_lasso[train_ind,]), DT$ytrain[train_ind])^2

(nzcoef <- coef(fitlasso, s = lambda.min)[nonzeroCoef(coef(fitlasso, s = lambda.min)),,drop=F])
nrow(nzcoef)















set.seed(1234)
cvfitlasso <- cv.glmnet(xtrain_lasso, ytrain, alpha = 1, nfolds = 10)
plot(cvfitlasso)
(nzcoef <- coef(cvfitlasso, s = "lambda.min")[nonzeroCoef(coef(cvfitlasso, s = "lambda.min")),,drop=F])
nrow(nzcoef)
# coef(cvfitlasso, s = "lambda.min")
sprintf("%.2f (%.2f, %.2f)",
        cvfitlasso$cvm[which(cvfitlasso$lambda.min==cvfitlasso$lambda)],
        cvfitlasso$cvlo[which(cvfitlasso$lambda.min==cvfitlasso$lambda)],
        cvfitlasso$cvup[which(cvfitlasso$lambda.min==cvfitlasso$lambda)])
cor(predict(cvfitlasso, s="lambda.min", newx = xtrain_lasso), ytrain)^2



# glinternet --------------------------------------------------------------
DT %>% names
DT$xtrain_lasso %>% dim
cvglinternet <- glinternet.cv(X = DT[["xtrain_lasso"]], Y = DT[["ytrain"]], nFolds = 10,
                              numLevels = c(2,rep(1, ncol(DT[["xtrain_lasso"]])-1)),
                              family = "gaussian",
                              nLambda = 100, interactionCandidates = c(1),
                              verbose = T)

fitGL <- glinternet(X = DT$xtrain_lasso[train_ind,], Y = DT$ytrain[train_ind],
                    numLevels = c(2,rep(1, ncol(DT$xtrain_lasso[train_ind,])-1)),
                    nLambda = 100, interactionCandidates = c(1),
                    verbose = T)


# sail --------------------------------------------------------------------

# data split entire data
f.basis <- function(i) splines::bs(i, degree = 5)
# f.basis <- function(i) i
fit <- sail(x = xtrain, y = ytrain, e = etrain,
            basis = f.basis, alpha = 0.5, fdev = 1e-07,
            maxit = 250, strong = FALSE, verbose = 2, dfmax = 30,
            nlambda = 100, center.e = TRUE)

plot(fit)

ytest_hat <- predict(fit, newx = xtest, newe = etest)
msetest <- colMeans((ytest - ytest_hat)^2)
plot(log(fit$lambda), msetest)
lambda.min.index <- as.numeric(which.min(msetest))
lambda.min <- fit$lambda[which.min(msetest)]

yvalid_hat <- predict(fit, newx = xvalid, newe = evalid, s = lambda.min)
(msevalid <- mean((yvalid - drop(yvalid_hat))^2))
cor(yvalid_hat, yvalid)^2
(nzcoef <- predict(fit, s = lambda.min, type = "nonzero"))
fit$active[lambda.min.index]

# data split training only

# set.seed(123456)

train_test_ind <- caret::createDataPartition(DT$ytrain, p = 200/338)[[1]]
validate_ind <- seq(length(DT$ytrain))[-train_test_ind]
train_ind <- sample(train_test_ind, floor(length(train_test_ind)/2))
test_ind <- setdiff(train_test_ind, train_ind)

f.basis <- function(i) splines::bs(i, degree = 3)
fit <- sail(x = DT$xtrain[train_ind,], y = DT$ytrain[train_ind], e = DT$etrain[train_ind],
            basis = f.basis, alpha = 0.5, fdev = 1e-07,
            maxit = 250, strong = FALSE, verbose = 2, dfmax = 40,
            nlambda = 100, center.e = TRUE)

ytest_hat <- predict(fit, newx = DT$xtrain[test_ind,], newe = DT$etrain[test_ind])
msetest <- colMeans((DT$ytrain[test_ind] - ytest_hat)^2)
# plot(log(fit$lambda)[5:30], msetest[5:30])
lambda.min.index <- as.numeric(which.min(msetest))
lambda.min <- fit$lambda[which.min(msetest)]

yvalid_hat <- predict(fit, newx = DT$xtrain[validate_ind,], newe = DT$etrain[validate_ind], s = lambda.min)
(msevalid <- mean((DT$ytrain[validate_ind] - drop(yvalid_hat))^2))
cor(yvalid_hat, DT$ytrain[validate_ind])^2
(nzcoef <- predict(fit, s = lambda.min, type = "nonzero"))
fit$active[lambda.min.index]






# cv

f.basis <- function(i) splines::bs(i, degree = 5)
registerDoMC(cores = 10)
set.seed(123)
cvfit05 <- cv.sail(x = xtrain, y = ytrain, e = etrain, basis = f.basis, alpha = 0.5,
                   maxit = 500, thresh = 1e-04,
                   fdev = 1e-8,
                   grouped = FALSE,
                   # dfmax = 30,
                   strong = FALSE, verbose = 2, nfolds = 10, parallel = TRUE)
# saveRDS(cvfit03, file = "rda/scott/cvfit_alpha01_groupedFalse_strongFalse_nfolds5_bs3_Estd2_pop_100m.rds")
plot(cvfit03)
plot(cvfit03$sail.fit)
predict(cvfit05, type="non", s = "lambda.min")
cor(predict(cvfit05, s="lambda.min"), ytrain)^2


sprintf("%.2f (%.2f, %.2f)",
        cvfit05$cvm[which(cvfit05$lambda.min==cvfit05$lambda)],
        cvfit05$cvlo[which(cvfit05$lambda.min==cvfit05$lambda)],
        cvfit05$cvup[which(cvfit05$lambda.min==cvfit05$lambda)])
cvfit05$sail.fit$active[which(cvfit05$lambda.min==cvfit05$lambda)]

