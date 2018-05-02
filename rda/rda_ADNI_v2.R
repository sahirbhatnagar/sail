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
# pacman::p_load(mice)


# amy_pheno <- xlsx::read.xlsx("~/Downloads/DrCelia_data.xlsx", sheetIndex = 1)
amy_mat <- read.csv("~/git_repositories/sail/data-nogit/adni_new/csf_amyloid_final.csv", stringsAsFactors = FALSE)
covr <- read.csv("~/git_repositories/sail/data-nogit/adni_new/covariates.csv", stringsAsFactors = FALSE, sep = ";")
# surv <- read.csv("~/git_repositories/sail/data/adni_new/fdg_info.csv", stringsAsFactors = FALSE, sep = ",")

# sum(as.character(amy_pheno$PTID) %in% amy_mat$PTID)
# sum(as.character(covr$IID) %in% amy_mat$PTID)
# as.character(covr$IID) %>% unique() %>% length()

DT <- dplyr::inner_join(amy_mat, covr, by = c("PTID" = "IID")) %>%
  select(-AV45_path_bl)
colnames(DT)
DT$diag_3bl.x %>% table

DT <- DT %>% dplyr::filter(diag_3bl.x == 2)

brain_regions <- grep("X", colnames(DT), value=T)
fmla <- reformulate(c(sapply(brain_regions, function(i) sprintf("bs(%s)",i)),
                      "Age_bl", "diag_3bl.y"), intercept = FALSE)
fmla <- reformulate(c(sapply(brain_regions, function(i) sprintf("bs(%s)",i)),
                      "Age_bl"), intercept = FALSE)

train <- caret::createDataPartition(DT$diag_3bl.x)[[1]]
train <- seq(nrow(DT))
test <- setdiff(seq(nrow(DT)), train)
DT$diag_3bl.x[train] %>% table
DT$diag_3bl.x[test] %>% table

model_mat <- model.matrix(fmla, data = DT[train,])
head(model_mat)

# X <- DT %>% select(starts_with("X"), Age_bl, diag_3bl.x) %>% as.matrix()
# X <- DT %>% select(starts_with("X"), Age_bl, EDUCAT) %>% as.matrix()
# X <- DT %>% select(starts_with("X")) %>% as.matrix()
# X <- DT %>% select(starts_with("X"), diag_3bl.x) %>% as.matrix()
# dimnames(X)[[1]] <- DT$PTID
# X <- X[,80:97]

# ind <- which(DT$diag_3bl.x==2)
# X <- X[ind,,drop = F]

E <- DT[train,] %>% pull(APOE_bin) %>% as.numeric
# E <- DT %>% pull(APOE_bin) %>% as.numeric
# E <- E[ind]
# E <- DT[train,] %>% pull(EDUCAT) %>% as.numeric
# E <- DT[train,] %>% pull(diag_3bl.x) %>% as.numeric
# E <- DT[train,] %>% pull(Age_bl) %>% as.numeric

Y <- DT[train,] %>% pull(MMSCORE_bl) %>% as.numeric
# Y <- Y[ind]
dev.off()
hist(DT$MMSCORE_bl)
hist(Y)
hist(E)
table(Y)
hist(scale(Y, scale = F))
# summary(lm(Y ~ X*E))
f.basis <- function(i) splines::bs(i, df = 5)
system.time(
  fit <- sail(x = X, y = Y, e = E, basis = f.basis, alpha = 0.1, strong = FALSE)
)

fit <- sail(
  x = model_mat,
  # x = as.matrix(DT[,brain_regions[1:10]]),
  y = Y, e = E,
  strong = FALSE,
  # expand = FALSE,
  # center.e = FALSE,
  # alpha = 0.2
  # fdev = 1e-8,
  # group = attr(model_mat, "assign"),
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
