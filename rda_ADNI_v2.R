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
amy_mat <- read.csv("~/git_repositories/sail/data/adni_new/csf_amyloid_final.csv", stringsAsFactors = FALSE)
covr <- read.csv("~/git_repositories/sail/data/adni_new/covariates.csv", stringsAsFactors = FALSE, sep = ";")
# surv <- read.csv("~/git_repositories/sail/data/adni_new/fdg_info.csv", stringsAsFactors = FALSE, sep = ",")

# sum(as.character(amy_pheno$PTID) %in% amy_mat$PTID)
# sum(as.character(covr$IID) %in% amy_mat$PTID)
# as.character(covr$IID) %>% unique() %>% length()

DT <- dplyr::inner_join(amy_mat, covr, by = c("PTID" = "IID")) %>%
  select(-AV45_path_bl)

# X <- DT %>% select(starts_with("X"), Age_bl, diag_3bl.x) %>% as.matrix()
X <- DT %>% select(starts_with("X"), Age_bl, EDUCAT) %>% as.matrix()
# X <- DT %>% select(starts_with("X")) %>% as.matrix()
# X <- DT %>% select(starts_with("X"), diag_3bl.x) %>% as.matrix()
dimnames(X)[[1]] <- DT$PTID
# X <- X[,80:97]

# ind <- which(DT$diag_3bl.x==2)
# X <- X[ind,,drop = F]

# E <- DT %>% pull(APOE_bin) %>% as.numeric
# E <- DT %>% pull(APOE_bin) %>% as.numeric
# E <- E[ind]
# E <- DT %>% pull(EDUCAT) %>% as.numeric
E <- DT %>% pull(diag_3bl.x) %>% as.numeric

Y <- DT %>% pull(MMSCORE_bl) %>% as.numeric
# Y <- Y[ind]

hist(Y)
hist(E)
table(Y)
hist(scale(Y, scale = F))
# summary(lm(Y ~ X*E))
system.time(
fit <- sail(x = X, y = Y, e = E, df = 3, degree = 3, maxit = 1000,
            thresh = 1e-4,
            nlambda = 100,
            verbose = TRUE, alpha = 0.05)
)

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
