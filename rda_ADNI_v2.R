rm(list=ls())
devtools::document()
devtools::load_all()
# source("R/ksfilter.R")
pacman::p_load(data.table)
pacman::p_load(magrittr)
pacman::p_load(genefilter)
pacman::p_load(tidyverse)
pacman::p_load(doParallel)
pacman::p_load(mice)


amy_pheno <- xlsx::read.xlsx("~/Downloads/DrCelia_data.xlsx", sheetIndex = 1)
amy_mat <- read.csv("~/git_repositories/sail/data/adni_new/csf_amyloid_final.csv", stringsAsFactors = FALSE)
covr <- read.csv("~/git_repositories/sail/data/adni_new/covariates.csv", stringsAsFactors = FALSE, sep = ";")
surv <- read.csv("~/git_repositories/sail/data/adni_new/fdg_info.csv", stringsAsFactors = FALSE, sep = ",")

# sum(as.character(amy_pheno$PTID) %in% amy_mat$PTID)
# sum(as.character(covr$IID) %in% amy_mat$PTID)
# as.character(covr$IID) %>% unique() %>% length()

DT <- dplyr::inner_join(amy_mat, covr, by = c("PTID" = "IID")) %>%
  select(-AV45_path_bl)

X <- DT %>% select(starts_with("X"), Age_bl) %>% as.matrix()
dimnames(X)[[1]] <- DT$PTID
X <- X[,80:97]

E <- DT %>% pull(APOE_bin) %>% as.numeric

Y <- DT %>% pull(MMSCORE_bl) %>% as.numeric

summary(lm(Y ~ X*E))
fit <- sail(x = X, y = Y, e = E, df = 3, maxit = 500,
            nlambda.gamma = 5, nlambda.beta = 5,
            nlambda = 25,
            # group.penalty = "SCAD",
            lambda.factor = 0.0001,
            thresh = 1e-4, center=T, normalize=F, verbose = T)
coef(fit)
cvfit <- cv.sail(x = X, y = Y, e = E, df = 3, maxit = 500,
                 nlambda.gamma = 5, nlambda.beta = 5,
                 nlambda = 25,
                 group.penalty = "SCAD",
                 # lambda.beta = exp(seq(log(0.01 * 1000), log(1000), length.out = 25)),
                 # lambda.gamma =rep(1000,25),
                 lambda.factor = 0.0001,
                 nfolds = 5,
                 thresh = 1e-5, center=T, normalize=F, verbose = T)

