## @knitr models

# devtools::load_all("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/")


######With raw data

make_data_split <- function(phenoVariable, DT, X, E, Y, exposure, n_train_test ) {


  PRS <- grep("PRS", colnames(DT), value=T)


  fmla <- reformulate(c(sapply(PRS, function(i) sprintf("bs(%s,3)",i))
  ), intercept = FALSE)

  #Code if we want to normalize the data for values that are not on a similar scale
  #Xnorm <- sail:::standardize(X, center = TRUE, normalize = TRUE)$x

  model_mat <- model.matrix(fmla, data = as.data.frame(X))

  #We group together the columns that have been separated into spline components corresponding to the same PRS variable
  group = attr(model_mat, "assign")


  new_model(name = "PRS_IQ4y_data_split_Cond2_split2_4",
            label = sprintf("traintest = %s, pheno = %s, exposure = %s",
                            n_train_test, phenoVariable, exposure),
            params = list(
                          phenoVariable = phenoVariable,
                          exposure = exposure,
                          E = E,
                          Y = Y,
                          X = X,
                          model_mat = model_mat,
                          group = group,
                          n_train_test = n_train_test),
            simulate = function( phenoVariable, exposure, E, Y, X, model_mat, group, n_train_test, nsim) {

              models <- list()

              for(i in seq(nsim)) {

                #need to do it seperately, because for sail we need to normalize externally
                dat_lasso <- partition_data(x = X[,-which(colnames(X) %in% c("Tx_group_bin","IQ_4yrs"))],
                                            y = X[,"IQ_4yrs"], e = X[,"Tx_group_bin"], p = n_train_test/length(Y),
                                            partition_on = X[, "Tx_group_bin"], type = "train_test_val")
                dat <- partition_data(x = X[,-which(colnames(X) %in% c("Tx_group_bin","IQ_4yrs"))], y = X[,"IQ_4yrs"], e = X[,"Tx_group_bin"], p = n_train_test/length(Y),
                                      partition_on = X[,"Tx_group_bin"], type = "train_test_val")

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
    train_ind <- sample(train_test_ind, ceiling((0.73333)*length(train_test_ind)))
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







####################################################
#EXTRA - NOT USED

######With normalized data

make_data_splitN <- function(phenoVariable , exposure , n_train_test ) {



  PRS <- grep("PRS", colnames(DT), value=T)


  fmla <- reformulate(c(sapply(PRS, function(i) sprintf("bs(%s,3)",i))
  ), intercept = FALSE)

  #We will normalize the datra since IQ is not on the same scale as the PRS scores #Amanda
  Xnorm <- sail:::standardize(X, center = TRUE, normalize = TRUE)$x

  model_mat <- model.matrix(fmla, data = as.data.frame(Xnorm))

  #We group together the columns that have been separated into spline components corresponding to the same PRS variable
  group = attr(model_mat, "assign")

  E <- Xnorm[, "Tx_group_bin"]



  new_model(name = "testm_E",
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
                dat_lasso <- partition_data(x = Xnorm[,-which(colnames(X) %in% c("Tx_group_bin","IQ_4yrs"))],
                                            y = Xnorm[,"IQ_4yrs"], e = Xnorm[,"Tx_group_bin"], p = n_train_test/length(Y),
                                            partition_on = Xnorm[, "Tx_group_bin"], type = "train_test_val")
                dat <- partition_data(x = model_mat, y = Xnorm[,"IQ_4yrs"], e = Xnorm[,"Tx_group_bin"], p = n_train_test/length(Y),
                                      partition_on = Xnorm[,"Tx_group_bin"], type = "train_test_val")

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


