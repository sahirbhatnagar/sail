######################################
# R Source code file for real data analysis of ADNI data
# Notes:
# I used the simulator framework to run this real data analysis
# the method functions were slightly modified for this
# because the real data input are slightly different
# see my_sims/method_functions_rda.R for details
# we didnt fit spam or gamsel because it wouldnt let expansion of binary predictors
# see my_sims/plot_functions_rda_ADNI.R for functions used for plotting
# Author: Sahir Bhatnagar
# Created: 2016
# Updated: May 17, 2018
#####################################


# load packages ---------------------------------------------------------

# rm(list=ls())
pacman::p_load(simulator) # this file was created under simulator version 0.2.1
devtools::load_all()
pacman::p_load(data.table)
pacman::p_load(magrittr)
pacman::p_load(genefilter)
pacman::p_load(tidyverse)
pacman::p_load(doParallel)
pacman::p_load(splines)
pacman::p_load(foreach)
pacman::p_load(doMC)
pacman::p_load(glmnet)
pacman::p_load_current_gh('sahirbhatnagar/glinternet')
pacman::p_load(LassoBacktracking)
pacman::p_load(SAM)
pacman::p_load(gamsel)
pacman::p_load_gh('asadharis/HierBasis')
pacman::p_load(psych) # for error.crosses plot
pacman::p_load(ggplot2)



# load source code helper files -------------------------------------------

source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_git_v2/sail/my_sims/model_functions.R")
source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_git_v2/sail/my_sims/method_functions_rda.R")
source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_git_v2/sail/my_sims/eval_functions.R")
source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_git_v2/sail/my_sims/plot_functions_rda_ADNI.R")


# load ADNI data and merge covariates with phenotypes ---------------------

amy_mat <- read.csv("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_git_v2/sail/rda/csf_amyloid_final.csv", stringsAsFactors = FALSE)
covr <- read.csv("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_git_v2/sail/rda/covariates.csv", stringsAsFactors = FALSE, sep = ";")

# these are used for plotting. we use the entire data set for the plots --------------------
DT <- dplyr::inner_join(amy_mat, covr, by = c("PTID" = "IID")) %>%
  select(-AV45_path_bl)
X <- DT %>% select(starts_with("X"), diag_3bl.x, APOE_bin) %>%
  mutate(diag_3bl.x = diag_3bl.x - 1) %>%
  as.matrix()
Xnorm <- sail:::standardize(X, center = TRUE, normalize = TRUE)$x
E <- Xnorm[, "diag_3bl.x"]


# Run real data analysis --------------------------------------------------

sim <- new_simulation(name = "rda_ADNI_may_17_2018v2",
                      label = "rda_ADNI_may_17_2018v2",
                      dir = ".") %>%
  generate_model(make_ADNI_data_split, seed = 12345,
                 phenoVariable = "MMSCORE_bl",
                 exposure = "diag_3bl.x",n_train_test = 250) %>%
  simulate_from_model(nsim = 6, index = 1:35)  %>%
  run_method(list(lassosplit, sailsplit, GLinternetsplit,lassoBTsplit, Hiersplit),
             parallel = list(socket_names = 35,
                             libraries = c("LassoBacktracking", "glinternet","glmnet","splines",
                                           "magrittr","sail","gamsel","SAM","HierBasis","simulator", "parallel")))

simulator::save_simulation(sim)


# load simulator object for plots -----------------------------------------

sim <- simulator::load_simulation("rda_ADNI_may_17_2018v2")
sim <- sim %>% evaluate(list(msevalid, nactive, r2))
df <- as.data.frame(evals(sim))
saveRDS(df, file = "my_sims/rda_results/rda_ADNI_may_17_2018v2.rds")



# plot results ------------------------------------------------------------

sim %>% plot_eval(metric_name = "nactive")
sim %>% plot_eval(metric_name = "mse")
sim %>% plot_eval(metric_name = "r2")


## ---- load-results-adni ----
df <- readRDS("my_sims/rda_results/rda_ADNI_may_17_2018v2.rds")
DT_res <- df %>% as.data.table()


## ---- get-best-sail-model ----

draw_ind <- DT_res[Method=="sail"][which.min(mse)]$Draw
# dat <- draws(sim)@draws[[draw_ind]]
# saveRDS(dat, "my_sims/rda_results/rda_ADNI_may_17_2018v2_data.rds")
dat <- readRDS("my_sims/rda_results/rda_ADNI_may_17_2018v2_data.rds")
fit <- sail(x = dat$xtrain, y = dat$ytrain, e = dat$etrain,
            expand = FALSE,
            center.x = F,
            center.e = T,
            group = dat$group,
            alpha = 0.1,
            maxit = 250,
            strong = TRUE,
            verbose = 1)
ytest_hat <- predict(fit, newx = dat$xtest, newe = dat$etest)
msetest <- colMeans((dat$ytest - ytest_hat)^2)
lambda.min.index <- as.numeric(which.min(msetest))
lambda.min <- fit$lambda[which.min(msetest)]
yvalid_hat <- predict(fit, newx = dat$xvalid, newe = dat$evalid, s = lambda.min)
msevalid <- mean((dat$yvalid - drop(yvalid_hat))^2)
nzcoef <- predict(fit, s = lambda.min, type = "nonzero")

design <- do.call(rbind, list(dat$xtrain, dat$xtest, dat$xvalid))



## ---- plot-main-adni ----

# names taken from email from lai (search for ADNI in gmail)
plotMainADNI(fit, x = DT[as.numeric(as.character(rownames(design))), "X175"],
             xvar = paste0("bs(X175, 3)",1:3),design = design, s = lambda.min,
             ylab = "f(X175)", xlab = "Cuneus right")
plotMainADNI(fit, x = DT[as.numeric(as.character(rownames(design))), "X196"],
             xvar = paste0("bs(X196, 3)",1:3),design = design, s = lambda.min,
             ylab = "f(X196)", xlab = "Lateral occipitotemporal gyrus left")
plotMainADNI(fit, x = DT[as.numeric(as.character(rownames(design))), "X154"],
             xvar = paste0("bs(X154, 3)",1:3),design = design, s = lambda.min,
             ylab = "f(X154)", xlab =  "Middle occipital gyrus left")
plotMainADNI(fit, x = DT[as.numeric(as.character(rownames(design))), "X24"],
             xvar = paste0("bs(X24, 3)",1:3),design = design, s = lambda.min,
             ylab = "f(X24)", xlab =  "Middle occipital gyrus left")

## ---- plot-inter-adni-175 ----

par(mfrow=c(1,2), tcl=-0.5, family="serif",
    omi=c(0.2,0.2,0,0))
plotInterADNI(fit, x = DT[as.numeric(as.character(rownames(design))), "X175"],
              xvar = paste0("bs(X175, 3)",1:3), design = design, s = lambda.min,
              e = E, apoe = FALSE, legend = T, legend.position = "bottomleft",
              ylab = "f(X175)", xlab = "Cuneus right", main = "APOE e4 = 0")
plotInterADNI(fit, x = DT[as.numeric(as.character(rownames(design))), "X175"],
              xvar = paste0("bs(X175, 3)",1:3), design = design, s = lambda.min,
              e = E, apoe = TRUE, legend = F,
              ylab = "f(X175)", xlab = "Cuneus right", main = "APOE e4 = 1")


## ---- plot-inter-adni-154 ----

par(mfrow=c(1,2), tcl=-0.5, family="serif",
    omi=c(0.2,0.2,0,0))
plotInterADNI(fit, x = DT[as.numeric(as.character(rownames(design))), "X154"],
              xvar = paste0("bs(X154, 3)",1:3), design = design, s = lambda.min,
              e = E, apoe = FALSE, legend = F,
              ylab = "f(X154)", xlab = "Middle occipital gyrus left", main = "APOE e4 = 0")
plotInterADNI(fit, x = DT[as.numeric(as.character(rownames(design))), "X154"],
              xvar = paste0("bs(X154, 3)",1:3), design = design, s = lambda.min,
              e = E, apoe = TRUE, legend = T,
              ylab = "f(X154)", xlab = "Middle occipital gyrus left", main = "APOE e4 = 1")


## ----error-crosses ----

affect.mat2 <- describeBy(DT_res[, c("mse","nactive")], DT_res$Method, mat = TRUE)

par(family="serif")
error.crosses(affect.mat2[c(6:10),],
              affect.mat2[c(1:5),],
              labels=unique(affect.mat2$group1),
              xlab="Number of Active Variables",
              main = "ADNI Data: Means (+/- 1 SD) from 200 Train/Validate/Test Splits",
              sd = TRUE,
              cex.lab = 1.4,
              cex.axis = 1.4,
              cex.main = 1.5,
              xlim = c(0, 34),
              ylab="Test Set MSE",
              colors = sail:::cbbPalette[c(1,3,4,7,2)],
              pch=16,cex=2)

# error.crosses(affect.mat2[c(1:5),],
#               affect.mat2[c(6:10),],
#               labels=unique(affect.mat2$group1),
#               ylab="Number of active variables",
#               sd = TRUE,
#               # xlim = c(0, 35),
#               xlab="Test set MSE",
#               colors = sail:::cbbPalette[c(1,3,4,7,2)],
#               pch=16,cex=2)


## ---- not used under this line ----



# sim <- new_simulation(name = "rda_ADNI_may_17_2018",
#                       label = "rda_ADNI_may_17_2018",
#                       dir = ".") %>%
#   generate_model(make_ADNI_data_split, seed = 12345,
#                  phenoVariable = "MMSCORE_bl",
#                  exposure = "diag_3bl.x",n_train_test = 200) %>%
#   simulate_from_model(nsim = 1, index = 1:2)  %>%
#   run_method(list(lassosplit,lassosplitadaptive, sailsplitweak, #lassoBTsplit, sailsplitweak, sailsplitadaptiveweak
#                   sailsplit, sailsplitadaptive),
#              parallel = list(socket_names = 20,
#                              libraries = c("LassoBacktracking", "glinternet","glmnet","splines",
#                                            "magrittr","sail","gamsel","SAM","HierBasis","simulator", "parallel")))

#
# draws(sim)@draws$r2.1$xtrain_lasso[, -which(colnames(draws(sim)@draws$r2.1$xtrain_lasso) %in% c("APOE_bin", "E"))]
# draws(sim)@draws$r2.1$group


rm(list=ls())
# devtools::document()
devtools::load_all()
# source("R/ksfilter.R")
pacman::p_load(data.table)
pacman::p_load(magrittr)
pacman::p_load(genefilter)
pacman::p_load(tidyverse)
pacman::p_load(doParallel)
pacman::p_load(splines)
pacman::p_load(foreach)
pacman::p_load(doMC)
pacman::p_load(glmnet)
pacman::p_load_current_gh('sahirbhatnagar/glinternet')
# pacman::p_load(mice)


# amy_pheno <- xlsx::read.xlsx("~/Downloads/DrCelia_data.xlsx", sheetIndex = 1)
# amy_mat <- read.csv("~/git_repositories/sail/data-nogit/adni_new/csf_amyloid_final.csv", stringsAsFactors = FALSE)
amy_mat <- read.csv("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_git_v2/sail/rda/csf_amyloid_final.csv", stringsAsFactors = FALSE)
covr <- read.csv("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_git_v2/sail/rda/covariates.csv", stringsAsFactors = FALSE, sep = ";")
# surv <- read.csv("~/git_repositories/sail/data/adni_new/fdg_info.csv", stringsAsFactors = FALSE, sep = ",")

# sum(as.character(amy_pheno$PTID) %in% amy_mat$PTID)
# sum(as.character(covr$IID) %in% amy_mat$PTID)
# as.character(covr$IID) %>% unique() %>% length()

DT <- dplyr::inner_join(amy_mat, covr, by = c("PTID" = "IID")) %>%
  select(-AV45_path_bl)
colnames(DT)
DT$diag_3bl.x %>% table
# DT <- DT %>% mutate(diag_3bl.x = ifelse(diag_3bl.x==3, 2, diag_3bl.x))
DT$diag_3bl.x %>% table
# DT <- DT %>% dplyr::filter(diag_3bl.x == 3)

brain_regions <- grep("X", colnames(DT), value=T)
# fmla <- reformulate(c(sapply(brain_regions, function(i) sprintf("bs(%s,3)",i)),
#                       "diag_3bl.x"), intercept = FALSE)
fmla <- reformulate(c(sapply(brain_regions, function(i) sprintf("bs(%s,3)",i)),
                      "APOE_bin"), intercept = FALSE)
# fmla <- reformulate(c(sapply(brain_regions, function(i) sprintf("bs(%s, degree = 5)",i)),
#                       "diag_3bl.y"), intercept = FALSE)
# fmla <- reformulate(c(sapply(brain_regions, function(i) sprintf("bs(%s)",i)),
#                       "bs(Age_bl)"), intercept = FALSE)
# fmla <- reformulate(c(sapply(brain_regions, function(i) sprintf("bs(%s)",i)),
#                       "bs(Age_bl)"), intercept = FALSE)
# fmla <- reformulate(c(sapply(brain_regions, function(i) sprintf("bs(%s)",i))),
#                     intercept = FALSE)

# train <- caret::createDataPartition(DT$MMSCORE_bl)[[1]]
# train <- seq(nrow(DT))
# test <- setdiff(seq(nrow(DT)), train)


# DT$diag_3bl.x[train] %>% table
# DT$diag_3bl.x[test] %>% table

X <- DT %>% select(starts_with("X"), diag_3bl.x, APOE_bin) %>%
  mutate(diag_3bl.x = diag_3bl.x - 1) %>%
  as.matrix()
# X <- DT %>% select(starts_with("X"), Age_bl) %>%
#   as.matrix()

Xnorm <- sail:::standardize(X, center = TRUE, normalize = TRUE)$x

model_mat <- model.matrix(fmla, data = as.data.frame(Xnorm))

# head(model_mat)

# X <- DT %>% select(starts_with("X"), Age_bl, diag_3bl.x) %>% as.matrix()

# X <- DT %>% select(starts_with("X"), Age_bl, EDUCAT) %>% as.matrix()
# X <- DT %>% select(starts_with("X")) %>% as.matrix()
# X <- DT %>% select(starts_with("X"), diag_3bl.x) %>% as.matrix()
# dimnames(X)[[1]] <- DT$PTID
# X <- X[,80:97]
# colnames(X)
# ind <- which(DT$diag_3bl.x==2)
# X <- X[ind,,drop = F]


# E <- DT[train,] %>% pull(APOE_bin) %>% as.numeric
E <- DT[train,] %>% pull(APOE_bin) %>% as.numeric
# E <- DT %>% pull(APOE_bin) %>% as.numeric

# E <- E[ind]
# E <- DT %>% pull(EDUCAT) %>% as.numeric
# E <- Xnorm[, "EDUCAT"]
# E <- DT[train,] %>% pull(diag_3bl.x) %>% as.numeric
# E <- DT[train,] %>% pull(Age_bl) %>% as.numeric

# the minus 1 for diag is for glinternet which needs 0,1,2,
# X <- DT %>% select(starts_with("X"), Age_bl, diag_3bl.x) %>%
#   mutate(diag_3bl.x = diag_3bl.x - 1) %>%
#   as.matrix()
X[,"diag_3bl.x"] %>% table
# E <- DT %>% pull(APOE_bin) %>% as.numeric
E <- Xnorm[, "diag_3bl.x"]
Y <- DT %>% pull(MMSCORE_bl) %>% as.numeric

set.seed(12345)
dat_lasso <- sail:::partition_data(x = X[,-which(colnames(X) %in% c("diag_3bl.x"))],
                                   y = Y, e = X[,"diag_3bl.x"], p = 200/length(Y),
                                   partition_on = Xnorm[, "diag_3bl.x"], type = "train_test_val")
dat <- sail:::partition_data(x = model_mat, y = Y, e = E, p = 200/length(Y),
                             partition_on = Xnorm[,"diag_3bl.x"], type = "train_test_val")

# dat_lasso <- sail:::partition_data(x = X, y = Y, e = E, p = 1,
#                                    partition_on = Y, type = "train_test")
# dat <- sail:::partition_data(x = model_mat, y = Y, e = E, p = 1,
#                              partition_on = Y, type = "train_test")

f.basis <- function(i) splines::bs(i, df = 5)
# fit <- sail(x = dat$xtrain, y = dat$ytrain, e = dat$etrain,
#             basis = f.basis,
#             alpha = 0.1,
#             thresh = 1e-02,
#             strong = FALSE,
#             verbose = 2)

#strong=TRUE works better for ADNI
fit <- sail(x = dat$xtrain, y = dat$ytrain, e = dat$etrain,
            # basis = f.basis,
            expand = FALSE,
            center.x = F,
            center.e = T,
            group = attr(model_mat, "assign"),
            alpha = 0.1,
            maxit = 250,
            # thresh = 1e-02,
            strong = TRUE,
            verbose = 2)
plot(fit)
ytest_hat <- predict(fit, newx = dat$xtest, newe = dat$etest)
msetest <- colMeans((dat$ytest - ytest_hat)^2)
plot(log(fit$lambda), msetest)
lambda.min.index <- as.numeric(which.min(msetest))
lambda.min <- fit$lambda[which.min(msetest)]

yvalid_hat <- predict(fit, newx = dat$xvalid, newe = dat$evalid, s = lambda.min)
(msevalid <- mean((dat$yvalid - drop(yvalid_hat))^2))

(nzcoef <- predict(fit, s = lambda.min, type = "nonzero"))


plotMainADNI(fit, x = X[dat$train_ind,"X175"], xvar = paste0("bs(X175, 3)",1:3), s = lambda.min,
         ylab = "f(Distance to Airport (YYZ) (m))", xlab = "Distance to Airport (YYZ) (m)")
plotMainADNI(fit, x = X[dat$train_ind,"X243"], xvar = paste0("bs(X243, 3)",1:3), s = lambda.min,
             ylab = "f(Distance to Airport (YYZ) (m))", xlab = "Distance to Airport (YYZ) (m)")

dev.off()

plotInterADNI(fit, x = X[,"X175"], xvar = paste0("bs(X175, 3)",1:3), s = lambda.min,
              model_matrix = model_mat, e = E,
              ylab = "f(X175)", xlab = "X175", main = "")
plotInterADNI(fit, x = X[,"X243"], xvar = paste0("bs(X243, 3)",1:3), s = lambda.min,
              model_matrix = model_mat, e = E,
              ylab = "f(X243)", xlab = "X243", main = "")

View(plotInter)





registerDoMC(cores = 10)
dat$xtrain %>% dim
cvfit <- cv.sail(x = dat$xtrain, y = dat$ytrain, e = dat$etrain,
                 # basis = f.basis,
               expand = FALSE,
               center.x = F,
               center.e = T,
               group = attr(model_mat, "assign"),
               alpha = 0.1,
               maxit = 250,
               parallel = TRUE,
               nfolds = 10,
               # thresh = 1e-02,
               strong = TRUE,
               verbose = 2)
plot(cvfit)
predict(cvfit, type="non", s = "lambda.min")

yvalid_hat <- predict(fit, newx = dat$xtest, newe = dat$etest, s = fit$lambda.min)
(msevalid <- mean((dat$ytest - drop(yvalid_hat))^2))

(nzcoef <- predict(fit, s = fit$lambda.min, type = "nonzero"))


fitGL <- glinternet(X = dat_lasso$xtrain_lasso, Y = dat_lasso$ytrain,
                    numLevels = c(3, rep(1, ncol(dat_lasso$xtrain_lasso)-2), 2),
                    # numLevels = c(2, rep(1, ncol(dat_lasso$xtrain_lasso)-1)),
                    nLambda = 100, interactionCandidates = c(1),
                    verbose = T)

ytest_hat <- predict(fitGL, X = dat_lasso$xtest_lasso)
msetest <- colMeans((dat_lasso$ytest - ytest_hat)^2)
lambda.min.index <- as.numeric(which.min(msetest))
lambda.min <- fitGL$lambda[which.min(msetest)]

yvalid_hat <- predict(fitGL, X = dat_lasso$xvalid_lasso, lambda = lambda.min)
(msevalid <- mean((dat_lasso$yvalid - drop(yvalid_hat))^2))
tc <- coef(fitGL, lambdaIndex = lambda.min.index)

(mains <- c(colnames(dat_lasso$xtrain_lasso)[tc[[1]]$mainEffects$cont],
  colnames(dat_lasso$xtrain_lasso)[tc[[1]]$mainEffects$cat]))

(inters <- c(paste(colnames(dat_lasso$xtrain_lasso)[tc[[1]]$interactions$catcont[,1]],
                 colnames(dat_lasso$xtrain_lasso)[tc[[1]]$interactions$catcont[,2]], sep=":"),
             paste(colnames(dat_lasso$xtrain_lasso)[tc[[1]]$interactions$contcont[,1]],
                   colnames(dat_lasso$xtrain_lasso)[tc[[1]]$interactions$contcont[,2]], sep=":"),
             paste(colnames(dat_lasso$xtrain_lasso)[tc[[1]]$interactions$catcat[,1]],
                   colnames(dat_lasso$xtrain_lasso)[tc[[1]]$interactions$catcat[,2]], sep=":")))


fitGL <- glinternet.cv(X = dat_lasso$xtrain_lasso, Y = dat_lasso$ytrain,
                       family = "gaussian",
                       nFolds = 10,
                       # numLevels = c(2, rep(1, ncol(dat_lasso$xtrain_lasso)-2), 3),
                       numLevels = c(2, rep(1, ncol(dat_lasso$xtrain_lasso)-1)),
                       nLambda = 100, interactionCandidates = c(1),
                       verbose = T)
plot(fitGL)
yvalid_hat <- predict(fitGL, X = dat_lasso$xtest_lasso, lambda = "lambdaHat")
(msevalid <- mean((dat_lasso$ytest - drop(yvalid_hat))^2))


tc <- coef(fitGL, lambdaType = "lambdaHat")
(mains <- colnames(dat_lasso$xtrain_lasso)[tc[[1]]$mainEffects$cont])
(inters <- paste(colnames(dat_lasso$xtrain_lasso)[tc[[1]]$interactions$catcont[,1]],
                 colnames(dat_lasso$xtrain_lasso)[tc[[1]]$interactions$catcont[,2]], sep=":"))



fit <- glmnet(x = dat_lasso[["xtrain_lasso"]], y = dat_lasso[["ytrain"]],
              alpha = 1)

ytest_hat <- predict(fit, newx = dat_lasso[["xtest_lasso"]])
msetest <- colMeans((dat_lasso[["ytest"]] - ytest_hat)^2)
lambda.min.index <- as.numeric(which.min(msetest))
lambda.min <- fit$lambda[which.min(msetest)]

yvalid_hat <- predict(fit, newx = dat_lasso[["xvalid_lasso"]], s = lambda.min)
(msevalid <- mean((dat_lasso[["yvalid"]] - drop(yvalid_hat))^2))

nzcoef <- coef(fit, s = lambda.min)[nonzeroCoef(coef(fit, s = lambda.min)),,drop=F]


fit <- cv.glmnet(x = dat_lasso[["xtrain_lasso"]], y = dat_lasso[["ytrain"]],
              alpha = 1, nfolds = 10)
plot(fit)
yvalid_hat <- predict(fit, newx = dat_lasso[["xtest_lasso"]], s = fit$lambda.min)
(msevalid <- mean((dat_lasso[["ytest"]] - drop(yvalid_hat))^2))

nzcoef <- coef(fit, s = lambda.min)[nonzeroCoef(coef(fit, s = lambda.min)),,drop=F]









# not-used ----------------------------------------------------------------


vars <- c("bs(X6, 3)1","bs(X6, 3)2","bs(X6, 3)3","bs(X154, 3)1","bs(X154, 3)2",
  "bs(X154, 3)3","bs(X175, 3)1","bs(X175, 3)2","bs(X175, 3)3","APOEbin",
  "bs(X6, 3)1:E","bs(X6, 3)2:E","bs(X6, 3)3:E","APOEbin:E","E")

(ints <- grep(":", vars, value=T))
(mains <- setdiff(vars, ints))
apoe <- if(any(mains=="APOEbin")) "APOEbin" else NULL
environ <- if(any(mains=="E")) "E" else NULL
mains_wo_apoe <- setdiff(mains, c("APOEbin","E"))

orig_mains <- unique(stringr::str_extract(mains_wo_apoe, "X\\d*"))
c(orig_mains, apoe, environ)



f.basis <- function(i) i
fit <- sail(x = X, y = Y, e = E, basis = f.basis, alpha = 0.5, thresh = 1e-02, strong = FALSE,
            verbose = 2)
plot(fit)
=======
system.time(
  fit <- sail(x = X, y = Y, e = E, basis = f.basis, alpha = 0.1, strong = FALSE)
)

fit <- sail(
  x = model_mat,
  # x = as.matrix(DT[,brain_regions[1:10]]),
  y = Y, e = E,
  expand = FALSE,
  thresh = 5e-3,
  strong = FALSE,
  # expand = FALSE,
  # center.e = FALSE,
  alpha = 0.2,
  # fdev = 1e-8,
  group = attr(model_mat, "assign"),
  verbose = 2)

plot(fit)
fit
as.matrix(coef(fit)[nonzero(coef(fit)),,])

as.matrix(fit$alpha[nonzero(fit$alpha),,])

help(sail)

fit
any(!fit$converged)
plot(fit)
fit$active
coef(fit)

tt <- KKT(b0 = fit$a0, betaE = fit$bE, beta = fit$beta, gamma = fit$gamma,
          alpha = fit$alpha, y = Y, phij = fit$Phi_j, xe_phij = fit$XE_Phi_j,
          e = E, df = fit$df,
          lambda = fit$lambda, lambda2 = fit$lambda2, group = fit$group,
          we = fit$we, wj = fit$wj, wje = fit$wje, thr = 1e-1, loss = "ls")

registerDoMC(cores = 8)
system.time(
  cvfit <- cv.sail(x = X, y = Y, e = E, df = 3, degree = 3, maxit = 1000,
                   alpha = 0.1,
                   nlambda = 100,
                   nfolds = 10,
                   verbose = TRUE,
                   thresh = 1e-4, parallel = TRUE)
)

plot(cvfit)
coef(cvfit, s = "lambda.min")[nonzero(coef(cvfit, s = "lambda.min")),,drop=F]
coef(cvfit, s = "lambda.1se")[nonzero(coef(cvfit, s = "lambda.1se")),,drop=F]
