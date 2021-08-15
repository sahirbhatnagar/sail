## @knitr models

# devtools::load_all("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/")

make_easy_sail_model <- function(n, p, df, SNR, betaE) {

  # n = 400;p=10;df=5;SNR=4;betaE=2

  #==============

  # group membership
  group <- rep(seq_len(p), each = df)

  # covariates
  X <- replicate(n = p, runif(n))
  E <- rnorm(n = n, sd = 0.5)

  # Expand X's
  Phi_j <- do.call(cbind, lapply(seq_len(p), function(j) splines::bs(X[,j], df = df)))
  XE_Phi_j <- E * Phi_j

  main_effect_names <- paste(paste0("X", group), rep(seq_len(df), times = p), sep = "_")
  interaction_names <- paste(main_effect_names, "X_E", sep = ":")

  dimnames(Phi_j)[[2]] <- main_effect_names
  dimnames(XE_Phi_j)[[2]] <- interaction_names

  design <- cbind(Phi_j, "X_E" = E, XE_Phi_j)

  true_beta <- matrix(rep(0, ncol(design), ncol = 1))
  dimnames(true_beta)[[1]] <- colnames(design)
  # the first 5 main effects and the first 2 interactions are active
  true_beta[c(1:(5*df),(p * df + 2):(p * df + 1 + 2 * df) ),1] <- rnorm(n = 7 * df)
  true_beta["X_E",] <- betaE

  causal <- rownames(true_beta[which(true_beta[,1]!=0),,drop=F])
  not_causal <- setdiff(rownames(true_beta), causal)


  new_model(name = "sail_easy",
            label = sprintf("n = %s, p = %s, df = %s, SNR = %s, betaE = %s", n, p, df, SNR, betaE),
            params = list(n = n, p = p, df = df, SNR = SNR, betaE = betaE,
                          true_beta = true_beta, causal = causal, design = design, X = X, E = E,
                          not_causal = not_causal),
            simulate = function(n, design, true_beta, SNR, nsim) {
              # error <- stats::rnorm(n)
              error <- MASS::mvrnorm(nsim, mu = rep(0, n), Sigma = diag(n))
              Y.star <- as.numeric(design %*% true_beta)
              k <- sqrt(stats::var(Y.star) / (SNR * apply(error, 1, var)))
              error2 <- sweep(t(error), 2, k, FUN = "*")
              y <- Y.star + error2
              return(split(y, col(y))) # make each col its own list element
            })
}


make_gendata_Paper <- function(n, p, corr, betaE, SNR, lambda.type, parameterIndex) {

  main <- paste0("X", seq_len(p))
  vnames <- c(main, "E", paste0(main,":E"))
  # used for glmnet and lasso backtracking
  # f.identity <- function(i) i

  new_model(name = "gendata_thesis_test",
            label = sprintf("n = %s, p = %s, corr = %s, betaE = %s, SNR = %s, index = %s, lambda = %s",
                            n, p, corr, betaE, SNR, parameterIndex, lambda.type),
            params = list(n = n, p = p, corr = corr, betaE = betaE, SNR = SNR,
                          lambda.type = lambda.type, parameterIndex = parameterIndex),
            simulate = function(n, p, corr, betaE, SNR, parameterIndex, lambda.type, nsim) {
              models <- list()
              for(i in seq(nsim)) {
                # test and training set
                DT <- sail::gendata(n = n, p = p, corr = corr, betaE = betaE,
                                    SNR = SNR, parameterIndex = parameterIndex)

                trainIndex <- drop(caret::createDataPartition(DT$y, p = 0.5, list = FALSE, times = 1))
                testIndex <- setdiff(seq_along(DT$y), trainIndex)

                xtrain <- DT$x[trainIndex,,drop=FALSE]
                xtest <- DT$x[testIndex,,drop=FALSE]

                etrain <- DT$e[trainIndex]
                etest <- DT$e[testIndex]

                ytrain <- DT$y[trainIndex]
                ytest <- DT$y[testIndex]

                main <- colnames(DT$x)
                vnames <- c(main, "E", paste0(main,":E"))
                vnames_lasso <- c("E", main) # needs to be in this order for glinternet

                # X_linear_design <- sail::design_sail(x = DT$x, e = DT$e, expand = TRUE, basis = function(i) i,
                #                                      nvars = ncol(DT$x),
                #                                      vnames = colnames(DT$x),
                #                                      center.x = FALSE, center.e = FALSE)$design

                # as is done in pliable lasso, only feed (X,E) to glmnet
                X_main <- cbind(E = etrain, xtrain)

                # test set
                # DT_test <- sail::gendata(n = n, p = p, corr = corr, betaE = betaE,
                #                          SNR = SNR, parameterIndex = parameterIndex)

                # X_main_test <- sail::design_sail(x = DT_test$x, e = DT_test$e, expand = TRUE, basis = function(i) i,
                #                                           nvars = ncol(DT_test$x),
                #                                           vnames = colnames(DT_test$x),
                #                                           center.x = FALSE, center.e = FALSE)$design
                X_main_test <- cbind(E = etest, xtest)

                # models[[i]] <- list(xtrain = DT$x, etrain = DT$e, ytrain = DT$y, xtrain_lasso = X_main,
                #                     xtest = DT_test$x, etest = DT_test$e, ytest = DT_test$y, xtest_lasso = X_main_test,
                #                     causal = DT$causal, not_causal = DT$not_causal, vnames = vnames, vnames_lasso = vnames_lasso)
                models[[i]] <- list(xtrain = xtrain, etrain = etrain, ytrain = ytrain, xtrain_lasso = X_main,
                                    xtest = xtest, etest = etest, ytest = ytest, xtest_lasso = X_main_test,
                                    causal = DT$causal, not_causal = DT$not_causal, vnames = vnames, vnames_lasso = vnames_lasso)
              }
              return(models)
            })

}



# n should be the total of train and test. 2*n is chosen as the validation
# this was use in CSDA paper
make_gendata_Paper_data_split <- function(n, p, corr, betaE, SNR, lambda.type, parameterIndex) {

  main <- paste0("X", seq_len(p))
  vnames <- c(main, "E", paste0(main,":E"))
  # used for glmnet and lasso backtracking
  # f.identity <- function(i) i

  new_model(name = "gendata_thesis_split_v4",
            label = sprintf("n = %s, p = %s, corr = %s, betaE = %s, SNR = %s, index = %s, lambda = %s",
                            n, p, corr, betaE, SNR, parameterIndex, lambda.type),
            params = list(n = n, p = p, corr = corr, betaE = betaE, SNR = SNR,
                          lambda.type = lambda.type, parameterIndex = parameterIndex),
            simulate = function(n, p, corr, betaE, SNR, parameterIndex, lambda.type, nsim) {
              models <- list()
              for(i in seq(nsim)) {
                # test and training set
                DT <- sail::gendata(n = n, p = p, corr = corr, betaE = betaE,
                                    SNR = SNR, parameterIndex = parameterIndex)

                #validation set
                DT_val <- sail::gendata(n = 2*n, p = p, corr = corr, betaE = betaE,
                                        SNR = SNR, parameterIndex = parameterIndex)

                trainIndex <- drop(caret::createDataPartition(DT$y, p = 0.5, list = FALSE, times = 1))
                testIndex <- setdiff(seq_along(DT$y), trainIndex)

                xtrain <- DT$x[trainIndex,,drop=FALSE]
                xtest <- DT$x[testIndex,,drop=FALSE]
                xvalid <- DT_val$x

                etrain <- DT$e[trainIndex]
                etest <- DT$e[testIndex]
                evalid <- DT_val$e

                ytrain <- DT$y[trainIndex]
                ytest <- DT$y[testIndex]
                yvalid <- DT_val$y

                main <- colnames(DT$x)
                vnames <- c(main, "E", paste0(main,":E"))
                vnames_lasso <- c("E", main) # needs to be in this order for glinternet

                # as is done in pliable lasso, only feed (E,X) to glmnet
                X_main <- cbind(E = etrain, xtrain)

                # test set
                X_main_test <- cbind(E = etest, xtest)

                # validation set
                X_main_valid <- cbind(E = evalid, xvalid)

                models[[i]] <- list(xtrain = xtrain, etrain = etrain, ytrain = ytrain, xtrain_lasso = X_main,
                                    xtest = xtest, etest = etest, ytest = ytest, xtest_lasso = X_main_test,
                                    xvalid = xvalid, evalid = evalid, yvalid = yvalid, xvalid_lasso = X_main_valid,
                                    causal = DT$causal, not_causal = DT$not_causal,
                                    vnames = vnames, vnames_lasso = vnames_lasso)
              }
              return(models)
            })

}





make_gendata_Paper_data_split_not_simulator <- function(n, p, corr, betaE, SNR, lambda.type, parameterIndex) {

  main <- paste0("X", seq_len(p))
  vnames <- c(main, "E", paste0(main,":E"))
  # used for glmnet and lasso backtracking
  # f.identity <- function(i) i

  DT <- sail::gendata(n = n, p = p, corr = corr, betaE = betaE,
                      SNR = SNR, parameterIndex = parameterIndex)

  #validation set
  DT_val <- sail::gendata(n = 2*n, p = p, corr = corr, betaE = betaE,
                          SNR = SNR, parameterIndex = parameterIndex)

  trainIndex <- drop(caret::createDataPartition(DT$y, p = 0.5, list = FALSE, times = 1))
  testIndex <- setdiff(seq_along(DT$y), trainIndex)

  xtrain <- DT$x[trainIndex,,drop=FALSE]
  xtest <- DT$x[testIndex,,drop=FALSE]
  xvalid <- DT_val$x

  etrain <- DT$e[trainIndex]
  etest <- DT$e[testIndex]
  evalid <- DT_val$e

  ytrain <- DT$y[trainIndex]
  ytest <- DT$y[testIndex]
  yvalid <- DT_val$y

  main <- colnames(DT$x)
  vnames <- c(main, "E", paste0(main,":E"))
  vnames_lasso <- c("E", main) # needs to be in this order for glinternet

  # as is done in pliable lasso, only feed (E,X) to glmnet
  X_main <- cbind(E = etrain, xtrain)

  # test set
  X_main_test <- cbind(E = etest, xtest)

  # validation set
  X_main_valid <- cbind(E = evalid, xvalid)

  models <- list(xtrain = xtrain, etrain = etrain, ytrain = ytrain, xtrain_lasso = X_main,
                      xtest = xtest, etest = etest, ytest = ytest, xtest_lasso = X_main_test,
                      xvalid = xvalid, evalid = evalid, yvalid = yvalid, xvalid_lasso = X_main_valid,
                      causal = DT$causal, not_causal = DT$not_causal,
                      vnames = vnames, vnames_lasso = vnames_lasso)
  return(models)
}






make_gendata_Paper_not_simulator <- function(n, p, corr, betaE, SNR, lambda.type, parameterIndex) {

  main <- paste0("X", seq_len(p))
  vnames <- c(main, "E", paste0(main,":E"))

  if (parameterIndex == 1) { # 1a
    hierarchy = "strong" ; nonlinear = TRUE ; interactions = TRUE
    causal <- c("X1","X2","X3","X4","E","X3:E","X4:E")
  } else if (parameterIndex == 2) { # 1b
    hierarchy = "weak" ; nonlinear = TRUE ; interactions = TRUE
    causal <- c("X1","X2","E","X3:E","X4:E")
  } else if (parameterIndex == 3) { # 1c
    hierarchy = "none" ; nonlinear = TRUE ; interactions = TRUE
    causal <- c("X3:E","X4:E")
  } else if (parameterIndex == 4) { # 2
    hierarchy = "strong"; nonlinear = FALSE; interactions = TRUE
    causal <- c("X1","X2","X3","X4","E","X3:E","X4:E")
  } else if (parameterIndex == 5) { # 3
    hierarchy = "strong" ; nonlinear = TRUE ; interactions = FALSE
    causal <- c("X1","X2","X3","X4","E")
  }

  not_causal <- setdiff(vnames, causal)

  DT <- gendataPaper(n = n, p = p, SNR = SNR, betaE = betaE,
                     hierarchy = hierarchy, nonlinear = nonlinear, interactions = interactions,
                     corr = corr)
                     # , E = truncnorm::rtruncnorm(n, a = -1, b = 1))

  return(DT)
  # # used for glmnet and lasso backtracking
  # X_linear_design <- design_sail(x = DT$x, e = DT$e, nvars = p,
  #                                vnames = paste0("X",1:p), degree = 1,
  #                                center.x = FALSE, basis.intercept = FALSE)$design
  #
  # new_model(name = "gendata_Paper",
  #           label = sprintf("n = %s, p = %s, corr = %s, betaE = %s, SNR = %s, hierarchy = %s,
  #                           nonlinear = %s, interactions = %s, scenario = %s",
  #                           n, p, corr, betaE, SNR, hierarchy,
  #                           nonlinear, interactions, parameterIndex),
  #           params = list(n = n, p = p, corr = corr, betaE = betaE, SNR = SNR, lambda.type = lambda.type,
  #                         hierarchy = hierarchy, nonlinear = nonlinear, vnames = vnames,
  #                         interactions = interactions, causal = causal, X_linear_design = X_linear_design,
  #                         not_causal = not_causal, x = DT$x, e = DT$e, Y.star = DT$Y.star, EX = cbind(E=DT$e, DT$x)),
  #           simulate = function(n, Y.star, nsim) {
  #             error <- MASS::mvrnorm(nsim, mu = rep(0, n), Sigma = diag(n))
  #             # Y.star <- as.numeric(design %*% true_beta)
  #             k <- sqrt(stats::var(Y.star) / (SNR * apply(error, 1, var)))
  #             error2 <- sweep(t(error), 2, k, FUN = "*")
  #             y <- Y.star + error2
  #             return(split(y, col(y))) # make each col its own list element
  #           })

}




# source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_git_v2/sail/rda/rda_NIHPD_data_cleaning.R")

make_nihpd_data_split <- function(nprobes, phenoVariable, exposure, filter, data) {

  new_model(name = "nihpd_split",
            label = sprintf("p = %s, pheno = %s, exposure = %s, filter = %s, data = %s",
                            nprobes, phenoVariable, exposure, filter, data),
            params = list(nprobes = nprobes,
                          phenoVariable = phenoVariable,
                          exposure = exposure,
                          filter = filter,
                          data = data),
            simulate = function(nprobes, phenoVariable, exposure, filter, data, nsim) {

              DT <- nihdata(nprobes = nprobes,
                            phenoVariable = phenoVariable,
                            exposure = exposure,
                            filter = filter,
                            data = data)

              models <- list()

              for(i in seq(nsim)) {

                train_test_ind <- caret::createDataPartition(DT$ytrain, p = 200/338)[[1]]

                validate_ind <- seq(length(DT$ytrain))[-train_test_ind]
                train_ind <- sample(train_test_ind, floor(length(train_test_ind)/2))
                test_ind <- setdiff(train_test_ind, train_ind)

                xtrain <- DT$xtrain[train_ind, , drop=FALSE]
                xtest <- DT$xtrain[test_ind, , drop=FALSE]
                xvalid <- DT$xtrain[validate_ind, , drop=FALSE]

                xtrain_lasso <- DT$xtrain_lasso[train_ind, , drop=FALSE]
                xtest_lasso <- DT$xtrain_lasso[test_ind, , drop=FALSE]
                xvalid_lasso <- DT$xtrain_lasso[validate_ind, , drop=FALSE]

                etrain <- DT$etrain[train_ind]
                etest <- DT$etrain[test_ind]
                evalid <- DT$etrain[validate_ind]

                ytrain <- DT$ytrain[train_ind]
                ytest <- DT$ytrain[test_ind]
                yvalid <- DT$ytrain[validate_ind]

                # main <- colnames(DT$x)
                # vnames <- c(main, "E", paste0(main,":E"))
                # vnames_lasso <- c("E", main) # needs to be in this order for glinternet

                models[[i]] <- list(xtrain = xtrain, etrain = etrain, ytrain = ytrain, xtrain_lasso = xtrain_lasso,
                                    xtest = xtest, etest = etest, ytest = ytest, xtest_lasso = xtest_lasso,
                                    xvalid = xvalid, evalid = evalid, yvalid = yvalid, xvalid_lasso = xvalid_lasso)
              }
              return(models)
            })

}



make_ADNI_data_split <- function(phenoVariable = "MMSCORE_bl", exposure = "diag_3bl.x", n_train_test = 200) {

  amy_mat <- read.csv("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_git_v2/sail/rda/csf_amyloid_final.csv", stringsAsFactors = FALSE)
  covr <- read.csv("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_git_v2/sail/rda/covariates.csv", stringsAsFactors = FALSE, sep = ";")
  DT <- dplyr::inner_join(amy_mat, covr, by = c("PTID" = "IID")) %>%
    select(-AV45_path_bl)

  brain_regions <- grep("X", colnames(DT), value=T)

  fmla <- reformulate(c(sapply(brain_regions, function(i) sprintf("bs(%s,3)",i)),
                        "APOE_bin"), intercept = FALSE)

  X <- DT %>% select(starts_with("X"), diag_3bl.x, APOE_bin) %>%
    mutate(diag_3bl.x = diag_3bl.x - 1) %>%
    as.matrix()

  Xnorm <- sail:::standardize(X, center = TRUE, normalize = TRUE)$x

  model_mat <- model.matrix(fmla, data = as.data.frame(Xnorm))
  group = attr(model_mat, "assign")

  E <- Xnorm[, "diag_3bl.x"]
  Y <- DT %>% pull(MMSCORE_bl) %>% as.numeric


  new_model(name = "ADNI_split",
            label = sprintf("traintest = %s, pheno = %s, exposure = %s",
                            n_train_test, phenoVariable, exposure),
            params = list(Xnorm = Xnorm,
                          phenoVariable = phenoVariable,
                          exposure = exposure,
                          E = E,
                          Y = Y,
                          X = X,
                          model_mat = model_mat,
                          group = group,
                          n_train_test = n_train_test),
            simulate = function(Xnorm, phenoVariable, exposure, E, Y, X, model_mat, group, n_train_test, nsim) {

              models <- list()

              for(i in seq(nsim)) {

                #need to do it seperately, because for sail we need to normalize externally
                dat_lasso <- partition_data(x = X[,-which(colnames(X) %in% c("diag_3bl.x"))],
                                                   y = Y, e = X[,"diag_3bl.x"], p = n_train_test/length(Y),
                                                   partition_on = Xnorm[, "diag_3bl.x"], type = "train_test_val")
                dat <- partition_data(x = model_mat, y = Y, e = E, p = n_train_test/length(Y),
                                             partition_on = Xnorm[,"diag_3bl.x"], type = "train_test_val")

                xtrain <- dat[["xtrain"]]
                xtest <- dat[["xtest"]]
                xvalid <- dat[["xvalid"]]

                xtrain_lasso <- dat_lasso[["xtrain_lasso"]]
                xtest_lasso <- dat_lasso[["xtest_lasso"]]
                xvalid_lasso <- dat_lasso[["xvalid_lasso"]]

                etrain <- dat[["etrain"]]
                etest <- dat[["etest"]]
                evalid <- dat[["evalid"]]

                etrain_lasso <- dat_lasso[["etrain"]]
                etest_lasso <- dat_lasso[["etest"]]
                evalid_lasso <- dat_lasso[["evalid"]]

                ytrain <- dat[["ytrain"]]
                ytest <- dat[["ytest"]]
                yvalid <- dat[["yvalid"]]

                ytrain_lasso <- dat_lasso[["ytrain"]]
                ytest_lasso <- dat_lasso[["ytest"]]
                yvalid_lasso <- dat_lasso[["yvalid"]]

                # main <- colnames(DT$x)
                # vnames <- c(main, "E", paste0(main,":E"))
                # vnames_lasso <- c("E", main) # needs to be in this order for glinternet

                models[[i]] <- list(xtrain = xtrain, xtrain_lasso = xtrain_lasso,
                                    etrain = etrain, etrain_lasso = etrain_lasso,
                                    ytrain = ytrain, ytrain_lasso = ytrain_lasso,
                                    xtest = xtest, xtest_lasso = xtest_lasso,
                                    etest = etest, etest_lasso = etest_lasso,
                                    ytest = ytest, ytest_lasso = ytest_lasso,
                                    xvalid = xvalid, xvalid_lasso = xvalid_lasso,
                                    evalid = evalid, evalid_lasso = evalid_lasso,
                                    yvalid = yvalid, yvalid_lasso = yvalid_lasso, group = group)
              }
              return(models)
            })

}



partition_data <- function(x, y, e, p, partition_on, type = c("train_test_val", "train_test")) {

  type <- match.arg(type)

  if (type == "train_test_val") {
    ex <- cbind(E = e, x)
    train_test_ind <- caret::createDataPartition(partition_on, p = p)[[1]]

    validate_ind <- seq(length(y))[-train_test_ind]
    train_ind <- sample(train_test_ind, floor(length(train_test_ind)/2))
    test_ind <- setdiff(train_test_ind, train_ind)

    xtrain <- x[train_ind, , drop=FALSE]
    xtest <- x[test_ind, , drop=FALSE]
    xvalid <- x[validate_ind, , drop=FALSE]

    xtrain_lasso <- ex[train_ind, , drop=FALSE]
    xtest_lasso <- ex[test_ind, , drop=FALSE]
    xvalid_lasso <- ex[validate_ind, , drop=FALSE]

    etrain <- e[train_ind]
    etest <- e[test_ind]
    evalid <- e[validate_ind]

    ytrain <- y[train_ind]
    ytest <- y[test_ind]
    yvalid <- y[validate_ind]

    return(list(xtrain = xtrain, etrain = etrain, ytrain = ytrain, xtrain_lasso = xtrain_lasso,
                xtest = xtest, etest = etest, ytest = ytest, xtest_lasso = xtest_lasso,
                xvalid = xvalid, evalid = evalid, yvalid = yvalid, xvalid_lasso = xvalid_lasso,
                train_ind = train_ind, test_ind = test_ind, validate_ind = validate_ind))
  } else if (type == "train_test"){
    ex <- cbind(E = e, x)
    train_ind <- caret::createDataPartition(partition_on, p = p)[[1]]

    xtrain <- x[train_ind, , drop=FALSE]
    xtest <- x[-train_ind, , drop=FALSE]

    xtrain_lasso <- ex[train_ind, , drop=FALSE]
    xtest_lasso <- ex[-train_ind, , drop=FALSE]

    etrain <- e[train_ind]
    etest <- e[-train_ind]

    ytrain <- y[train_ind]
    ytest <- y[-train_ind]

    return(list(xtrain = xtrain, etrain = etrain, ytrain = ytrain, xtrain_lasso = xtrain_lasso,
                xtest = xtest, etest = etest, ytest = ytest, xtest_lasso = xtest_lasso,
                train_ind = train_ind))
  }

}


# nsim = 10;n=100;SNR=3
# error <- MASS::mvrnorm(nsim, mu = rep(0, n), Sigma = diag(n))
# Y.star <- rnorm(n)
# k <- sqrt(stats::var(Y.star) / (SNR * apply(error, 1, var)))
# error2 <- sweep(t(error), 2, k, FUN = "*")
# y <- Y.star + error2
# split(y, col(y))
#
# all.equal(tr[,1], t(error)[,1]*k[1])
