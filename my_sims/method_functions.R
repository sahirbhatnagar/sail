## @knitr methods
# library(doMC)
# registerDoMC(cores = 10)
# error_return <- list(beta = NA,
#                      # cvfit = NULL,
#                      vnames = NA,
#                      nonzero_coef = NA,
#                      active = c(""),
#                      not_active = NA,
#                      yhat_test = NA,
#                      cvmse = NA,
#                      causal = NA,
#                      not_causal = NA,
#                      ytest = NA)

error_return <- list(beta = NA,
                     vnames = NA,
                     nonzero_coef = NA,
                     active = c(""),
                     not_active = NA,
                     yvalid_hat = NA,
                     msevalid = NA,
                     causal = NA,
                     not_causal = NA,
                     yvalid = NA)

sail <- new_method("sail", "Sail",
                   method = function(model, draw) {
                     tryCatch({
                     cvfit <- cv.sail(x = draw[["xtrain"]], y = draw[["ytrain"]], e = draw[["etrain"]],
                                      basis = function(i) splines::bs(i, degree = 5),
                                      # verbose = 1,
                                      # basis = function(i) i,
                                      # dfmax = 20,
                                      parallel = FALSE, nfolds = 10, nlambda = 100)

                     nzcoef <- coef(cvfit, s = model$lambda.type)[sail:::nonzero(coef(cvfit, s = model$lambda.type)),,drop=F]

                     list(beta = coef(cvfit, s = model$lambda.type)[-1,,drop=F],
                          # cvfit = cvfit,
                          vnames = draw[["vnames"]],
                          nonzero_coef = nzcoef,
                          active = cvfit$sail.fit$active[[which(cvfit[[model$lambda.type]]==cvfit$lambda)]],
                          not_active = setdiff(draw[["vnames"]], cvfit$sail.fit$active[[which(cvfit[[model$lambda.type]]==cvfit$lambda)]]),
                          yhat_test = predict(cvfit, s = cvfit[[model$lambda.type]], newx = draw[["xtest"]], newe = draw[["etest"]]),
                          cvmse = cvfit$cvm[which(cvfit[[model$lambda.type]]==cvfit$lambda)],
                          causal = draw[["causal"]],
                          not_causal = draw[["not_causal"]],
                          ytest = draw[["ytest"]])
                     },
                     error = function(err) {
                       return(error_return)
                     }
                     )
                   })

sailsplit <- new_method("sail", "Sail",
                   method = function(model, draw) {
                     tryCatch({
                       fit <- sail(x = draw[["xtrain"]], y = draw[["ytrain"]], e = draw[["etrain"]],
                                   basis = function(i) splines::bs(i, degree = 5))

                       ytest_hat <- predict(fit, newx = draw[["xtest"]], newe = draw[["etest"]])
                       msetest <- colMeans((draw[["ytest"]] - ytest_hat)^2)
                       lambda.min.index <- as.numeric(which.min(msetest))
                       lambda.min <- fit$lambda[which.min(msetest)]

                       yvalid_hat <- predict(fit, newx = draw[["xvalid"]], newe = draw[["evalid"]], s = lambda.min)
                       msevalid <- mean((draw[["yvalid"]] - drop(yvalid_hat))^2)

                       nzcoef <- predict(fit, s = lambda.min, type = "nonzero")

                       return(list(beta = coef(fit, s = lambda.min)[-1,,drop=F],
                            # fit = fit,
                            vnames = draw[["vnames"]],
                            nonzero_coef = nzcoef,
                            active = fit$active[[lambda.min.index]],
                            not_active = setdiff(draw[["vnames"]], fit$active[[lambda.min.index]]),
                            yvalid_hat = yvalid_hat,
                            msevalid = msevalid,
                            causal = draw[["causal"]],
                            not_causal = draw[["not_causal"]],
                            yvalid = draw[["yvalid"]]))
                     },
                     error = function(err) {
                     return(error_return)
                     }
                     )
                   })

sailsplitlinear <- new_method("linearsail", "Linear Sail",
                        method = function(model, draw) {
                          tryCatch({
                            fit <- sail(x = draw[["xtrain"]], y = draw[["ytrain"]], e = draw[["etrain"]],
                                        basis = function(i) i)

                            ytest_hat <- predict(fit, newx = draw[["xtest"]], newe = draw[["etest"]])
                            msetest <- colMeans((draw[["ytest"]] - ytest_hat)^2)
                            lambda.min.index <- as.numeric(which.min(msetest))
                            lambda.min <- fit$lambda[which.min(msetest)]

                            yvalid_hat <- predict(fit, newx = draw[["xvalid"]], newe = draw[["evalid"]], s = lambda.min)
                            msevalid <- mean((draw[["yvalid"]] - drop(yvalid_hat))^2)

                            nzcoef <- predict(fit, s = lambda.min, type = "nonzero")

                            return(list(beta = coef(fit, s = lambda.min)[-1,,drop=F],
                                        # fit = fit,
                                        vnames = draw[["vnames"]],
                                        nonzero_coef = nzcoef,
                                        active = fit$active[[lambda.min.index]],
                                        not_active = setdiff(draw[["vnames"]], fit$active[[lambda.min.index]]),
                                        yvalid_hat = yvalid_hat,
                                        msevalid = msevalid,
                                        causal = draw[["causal"]],
                                        not_causal = draw[["not_causal"]],
                                        yvalid = draw[["yvalid"]]))
                          },
                          error = function(err) {
                            return(error_return)
                          }
                          )
                        })


sailsplitadaptive <- new_method("Adaptivesail", "Adaptive Sail",
                        method = function(model, draw) {
                          tryCatch({
                            fit <- sail(x = draw[["xtrain"]], y = draw[["ytrain"]], e = draw[["etrain"]],
                                        basis = function(i) splines::bs(i, degree = 5))

                            ytest_hat <- predict(fit, newx = draw[["xtest"]], newe = draw[["etest"]])
                            msetest <- colMeans((draw[["ytest"]] - ytest_hat)^2)
                            lambda.min.index <- as.numeric(which.min(msetest))
                            lambda.min <- fit$lambda[which.min(msetest)]

                            pfs <- coef(fit, s = lambda.min)[-1,]
                            pfe <- 1/abs(pfs["E"])

                            pfmain <- pfs[fit$main.effect.names]
                            pfmain <- 1/sapply(split(pfmain, fit$group), sail:::l2norm)
                            pfmain[which(pfmain==Inf)] <- 50

                            pfinter <- pfs[fit$interaction.names]
                            pfinter <- 1/sapply(split(pfinter, fit$group), sail:::l2norm)
                            pfinter[which(pfinter==Inf)] <- 50

                            fit <- sail(x = draw[["xtrain"]], y = draw[["ytrain"]], e = draw[["etrain"]],
                                        basis = function(i) splines::bs(i, degree = 5),
                                        penalty.factor = c(pfe, pfmain, pfinter))

                            ytest_hat <- predict(fit, newx = draw[["xtest"]], newe = draw[["etest"]])
                            msetest <- colMeans((draw[["ytest"]] - ytest_hat)^2)
                            lambda.min.index <- as.numeric(which.min(msetest))
                            lambda.min <- fit$lambda[which.min(msetest)]

                            yvalid_hat <- predict(fit, newx = draw[["xvalid"]], newe = draw[["evalid"]], s = lambda.min)
                            msevalid <- mean((draw[["yvalid"]] - drop(yvalid_hat))^2)

                            nzcoef <- predict(fit, s = lambda.min, type = "nonzero")

                            return(list(beta = coef(fit, s = lambda.min)[-1,,drop=F],
                                        # fit = fit,
                                        vnames = draw[["vnames"]],
                                        nonzero_coef = nzcoef,
                                        active = fit$active[[lambda.min.index]],
                                        not_active = setdiff(draw[["vnames"]], fit$active[[lambda.min.index]]),
                                        yvalid_hat = yvalid_hat,
                                        msevalid = msevalid,
                                        causal = draw[["causal"]],
                                        not_causal = draw[["not_causal"]],
                                        yvalid = draw[["yvalid"]]))
                          },
                          error = function(err) {
                            return(error_return)
                          }
                          )
                        })


sailsplitweak <- new_method("sailweak", "Sail Weak",
                        method = function(model, draw) {
                          tryCatch({
                            fit <- sail(x = draw[["xtrain"]], y = draw[["ytrain"]], e = draw[["etrain"]],
                                        basis = function(i) splines::bs(i, degree = 5), strong = FALSE)

                            ytest_hat <- predict(fit, newx = draw[["xtest"]], newe = draw[["etest"]])
                            msetest <- colMeans((draw[["ytest"]] - ytest_hat)^2)
                            lambda.min.index <- as.numeric(which.min(msetest))
                            lambda.min <- fit$lambda[which.min(msetest)]

                            yvalid_hat <- predict(fit, newx = draw[["xvalid"]], newe = draw[["evalid"]], s = lambda.min)
                            msevalid <- mean((draw[["yvalid"]] - drop(yvalid_hat))^2)

                            nzcoef <- predict(fit, s = lambda.min, type = "nonzero")

                            return(list(beta = coef(fit, s = lambda.min)[-1,,drop=F],
                                        # fit = fit,
                                        vnames = draw[["vnames"]],
                                        nonzero_coef = nzcoef,
                                        active = fit$active[[lambda.min.index]],
                                        not_active = setdiff(draw[["vnames"]], fit$active[[lambda.min.index]]),
                                        yvalid_hat = yvalid_hat,
                                        msevalid = msevalid,
                                        causal = draw[["causal"]],
                                        not_causal = draw[["not_causal"]],
                                        yvalid = draw[["yvalid"]]))
                          },
                          error = function(err) {
                            return(error_return)
                          }
                          )
                        })

sailsplitadaptiveweak <- new_method("Adaptivesailweak", "Adaptive Sail Weak",
                                method = function(model, draw) {
                                  tryCatch({
                                    fit <- sail(x = draw[["xtrain"]], y = draw[["ytrain"]], e = draw[["etrain"]],
                                                strong = FALSE,
                                                basis = function(i) splines::bs(i, degree = 5))

                                    ytest_hat <- predict(fit, newx = draw[["xtest"]], newe = draw[["etest"]])
                                    msetest <- colMeans((draw[["ytest"]] - ytest_hat)^2)
                                    lambda.min.index <- as.numeric(which.min(msetest))
                                    lambda.min <- fit$lambda[which.min(msetest)]

                                    pfs <- coef(fit, s = lambda.min)[-1,]
                                    pfe <- 1/abs(pfs["E"])

                                    pfmain <- pfs[fit$main.effect.names]
                                    pfmain <- 1/sapply(split(pfmain, fit$group), sail:::l2norm)
                                    pfmain[which(pfmain==Inf)] <- 50

                                    pfinter <- pfs[fit$interaction.names]
                                    pfinter <- 1/sapply(split(pfinter, fit$group), sail:::l2norm)
                                    pfinter[which(pfinter==Inf)] <- 50

                                    fit <- sail(x = draw[["xtrain"]], y = draw[["ytrain"]], e = draw[["etrain"]],
                                                strong = FALSE,
                                                basis = function(i) splines::bs(i, degree = 5),
                                                penalty.factor = c(pfe, pfmain, pfinter))

                                    ytest_hat <- predict(fit, newx = draw[["xtest"]], newe = draw[["etest"]])
                                    msetest <- colMeans((draw[["ytest"]] - ytest_hat)^2)
                                    lambda.min.index <- as.numeric(which.min(msetest))
                                    lambda.min <- fit$lambda[which.min(msetest)]

                                    yvalid_hat <- predict(fit, newx = draw[["xvalid"]], newe = draw[["evalid"]], s = lambda.min)
                                    msevalid <- mean((draw[["yvalid"]] - drop(yvalid_hat))^2)

                                    nzcoef <- predict(fit, s = lambda.min, type = "nonzero")

                                    return(list(beta = coef(fit, s = lambda.min)[-1,,drop=F],
                                                # fit = fit,
                                                vnames = draw[["vnames"]],
                                                nonzero_coef = nzcoef,
                                                active = fit$active[[lambda.min.index]],
                                                not_active = setdiff(draw[["vnames"]], fit$active[[lambda.min.index]]),
                                                yvalid_hat = yvalid_hat,
                                                msevalid = msevalid,
                                                causal = draw[["causal"]],
                                                not_causal = draw[["not_causal"]],
                                                yvalid = draw[["yvalid"]]))
                                  },
                                  error = function(err) {
                                    return(error_return)
                                  }
                                  )
                                })



gbm <- new_method("gbm", "GBM",
                   method = function(model, draw) {

                     tryCatch({

                     gbmfit <- gbm::gbm.fit(x = draw[["xtrain_lasso"]], y = draw[["ytrain"]],
                                            distribution = "gaussian", verbose = FALSE,
                                            n.trees = 100, #default
                                            interaction.depth = 1)#default)

                     yhat_test <- gbm::predict.gbm(gbmfit,
                                                   newdata = draw[["xtest"]],
                                                   n.trees = 100)

                     topgbm <- summary(gbmfit, n.trees = 100)
                     active <- rownames(topgbm[which(topgbm$rel.inf>0),,drop=FALSE])


                     return(list(beta = NA,
                          # cvfit = cvfit,
                          vnames = draw[["vnames_lasso"]],
                          nonzero_coef = NA,
                          active = active,
                          not_active = setdiff(draw[["vnames_lasso"]], active),
                          yhat_test = yhat_test,
                          cvmse = NA,
                          causal = draw[["causal"]],
                          not_causal = draw[["not_causal"]],
                          ytest = draw[["ytest"]]))
                     },
                     error = function(err) {
                       return(error_return)
                     }
                     )
                   })


lasso <- new_method("lasso", "Lasso",
                    method = function(model, draw) {

                      tryCatch({

                      fitglmnet <- cv.glmnet(x = draw[["xtrain_lasso"]], y = draw[["ytrain"]],
                                             alpha = 1, nfolds = 10)

                      nzcoef <- coef(fitglmnet, s = model$lambda.type)[nonzeroCoef(coef(fitglmnet, s = model$lambda.type)),,drop=F]

                      return(list(beta = coef(fitglmnet, s = model$lambda.type)[-1,,drop=F],
                           # cvfit = fitglmnet,
                           vnames = draw[["vnames_lasso"]],
                           nonzero_coef = nzcoef,
                           active = setdiff(rownames(nzcoef), c("(Intercept)")),
                           not_active = setdiff(colnames(draw[["xtrain_lasso"]]), setdiff(rownames(nzcoef), c("(Intercept)"))),
                           yhat_test = predict(fitglmnet, newx = draw[["xtest_lasso"]], s = model$lambda.type),
                           cvmse = fitglmnet$cvm[which(fitglmnet[[model$lambda.type]]==fitglmnet$lambda)],
                           causal = draw[["causal"]],
                           not_causal = draw[["not_causal"]],
                           ytest = draw[["ytest"]]))
                      },
                      error = function(err) {
                        return(error_return)
                      }
                      )
                    })


lassosplit <- new_method("lasso", "Lasso",
                    method = function(model, draw) {

                      tryCatch({

                        fit <- glmnet(x = draw[["xtrain_lasso"]], y = draw[["ytrain"]],
                                      alpha = 1)

                        ytest_hat <- predict(fit, newx = draw[["xtest_lasso"]])
                        msetest <- colMeans((draw[["ytest"]] - ytest_hat)^2)
                        lambda.min.index <- as.numeric(which.min(msetest))
                        lambda.min <- fit$lambda[which.min(msetest)]

                        yvalid_hat <- predict(fit, newx = draw[["xvalid_lasso"]], s = lambda.min)
                        msevalid <- mean((draw[["yvalid"]] - drop(yvalid_hat))^2)

                        nzcoef <- coef(fit, s = lambda.min)[nonzeroCoef(coef(fit, s = lambda.min)),,drop=F]

                        return(list(beta = coef(fit, s = lambda.min)[-1,,drop=F],
                                    vnames = draw[["vnames_lasso"]],
                                    nonzero_coef = nzcoef,
                                    active = setdiff(rownames(nzcoef), c("(Intercept)")),
                                    not_active = setdiff(colnames(draw[["xtrain_lasso"]]), setdiff(rownames(nzcoef), c("(Intercept)"))),
                                    yvalid_hat = yvalid_hat,
                                    msevalid = msevalid,
                                    causal = draw[["causal"]],
                                    not_causal = draw[["not_causal"]],
                                    yvalid = draw[["yvalid"]]))
                      },
                      error = function(err) {
                        return(error_return)
                      }
                      )
                    })


lassosplitadaptive <- new_method("Adaptivelasso", "Adaptive Lasso",
                         method = function(model, draw) {

                           tryCatch({

                             fit <- glmnet(x = draw[["xtrain_lasso"]], y = draw[["ytrain"]],
                                           alpha = 1)

                             ytest_hat <- predict(fit, newx = draw[["xtest_lasso"]])
                             msetest <- colMeans((draw[["ytest"]] - ytest_hat)^2)
                             lambda.min.index <- as.numeric(which.min(msetest))
                             lambda.min <- fit$lambda[which.min(msetest)]

                             pfs <- coef(fit, s = lambda.min)[-1,]

                             fit <- glmnet(x = draw[["xtrain_lasso"]], y = draw[["ytrain"]],
                                           alpha = 1, penalty.factor = 1/abs(pfs))

                             ytest_hat <- predict(fit, newx = draw[["xtest_lasso"]])
                             msetest <- colMeans((draw[["ytest"]] - ytest_hat)^2)
                             lambda.min.index <- as.numeric(which.min(msetest))
                             lambda.min <- fit$lambda[which.min(msetest)]

                             yvalid_hat <- predict(fit, newx = draw[["xvalid_lasso"]], s = lambda.min)
                             msevalid <- mean((draw[["yvalid"]] - drop(yvalid_hat))^2)

                             nzcoef <- coef(fit, s = lambda.min)[nonzeroCoef(coef(fit, s = lambda.min)),,drop=F]

                             return(list(beta = coef(fit, s = lambda.min)[-1,,drop=F],
                                         vnames = draw[["vnames_lasso"]],
                                         nonzero_coef = nzcoef,
                                         active = setdiff(rownames(nzcoef), c("(Intercept)")),
                                         not_active = setdiff(colnames(draw[["xtrain_lasso"]]), setdiff(rownames(nzcoef), c("(Intercept)"))),
                                         yvalid_hat = yvalid_hat,
                                         msevalid = msevalid,
                                         causal = draw[["causal"]],
                                         not_causal = draw[["not_causal"]],
                                         yvalid = draw[["yvalid"]]))
                           },
                           error = function(err) {
                             return(error_return)
                           }
                           )
                         })


# lassoBT only gives the minimum CV error.. doesnt have lambda.1se
lassoBT <- new_method("lassoBT", "LassoBT",
                      method = function(model, draw) {

                        tryCatch({

                        # fitglmnet <- cv.glmnet(x = model$X_linear_design, y = draw, alpha = 1, nfolds = 10)
                        fitBT <- cvLassoBT(x = draw[["xtrain_lasso"]],
                                           y = draw[["ytrain"]], iter_max=10, nperms=1, nfolds = 10, nlambda = 100)

                        coefBT <- as.matrix(predict(fitBT$BT_fit, type = "coef",
                                                    s = fitBT$lambda[fitBT$cv_opt[1]], iter = fitBT$cv_opt[2]))
                        nzcoef <- coefBT[sail:::nonzero(coefBT),,drop=F]

                        ints <- grep(":",rownames(nzcoef)[-1], value=T)
                        if (length(ints) != 0) {
                          inters <- sapply(stringr::str_split(ints, ":"), function(i) paste(colnames(draw[["xtrain_lasso"]])[as.numeric(as.character(i))], collapse = ":"))
                          mains <- colnames(draw[["xtrain_lasso"]])[as.numeric(as.character(setdiff(rownames(nzcoef)[-1], ints)))]
                        } else {
                          inters <- character(0)
                          mains <- colnames(draw[["xtrain_lasso"]])[as.numeric(as.character(rownames(nzcoef)[-1]))]
                        }

                        active <- if (length(c(mains, inters)) == 0) " " else c(mains, inters)
                        yhat <- tryCatch({
                          preds <- predict(fitBT$BT_fit, newx = draw[["xtest_lasso"]], s = fitBT$lambda[fitBT$cv_opt[1]], type = "response")
                          tf <- as.matrix(preds[,fitBT$cv_opt[1], min(dim(preds)[3],fitBT$cv_opt[2]),drop=TRUE])
                        },
                        error = function(err) {
                          return(matrix(0, nrow = nrow(draw[["xtrain_lasso"]]), ncol = 1))
                        } # return NULL on error
                        )

                        return(list(beta = coefBT,
                             # cvfit = fitBT,
                             vnames = draw[["vnames"]],
                             nonzero_coef = nzcoef,
                             active = active,
                             not_active = setdiff(model$vnames, active),
                             yhat_test = yhat,
                             cvmse = fitBT$cv_opt_err,
                             causal = draw[["causal"]],
                             not_causal = draw[["not_causal"]],
                             ytest = draw[["ytest"]]))
                        },
                        error = function(err) {
                          return(error_return)
                        }
                        )
                      })



lassoBTsplit <- new_method("lassoBT", "LassoBT",
                           method = function(model, draw) {

                             tryCatch({

                               # fitBT <- cvLassoBT(x = draw[["xtrain_lasso"]],
                               #                    y = draw[["ytrain"]], iter_max=10, nperms=1, nfolds = 10, nlambda = 100)
                               fitBT <- LassoBT(x = draw[["xtrain_lasso"]],
                                              y = draw[["ytrain"]], iter_max=10, nlambda = 100)

                               ytest_hat <- predict(fitBT, newx = draw[["xtest_lasso"]])
                               msetest <- colMeans((draw[["ytest"]] - ytest_hat)^2)
                               lambda.min.index <- which(msetest==min(msetest), arr.ind = TRUE)
                               lambda.min <- fitBT$lambda[lambda.min.index[1,1]]
                               iter.min <- lambda.min.index[1,2]

                               coefBT <- as.matrix(predict(fitBT, type = "coef",
                                                           s = lambda.min, iter = iter.min))
                               nzcoef <- coefBT[sail:::nonzero(coefBT),,drop=F]

                               ints <- grep(":",rownames(nzcoef)[-1], value=T)
                               if (length(ints) != 0) {
                                 inters <- sapply(stringr::str_split(ints, ":"), function(i) paste(colnames(draw[["xtrain_lasso"]])[as.numeric(as.character(i))], collapse = ":"))
                                 mains <- colnames(draw[["xtrain_lasso"]])[as.numeric(as.character(setdiff(rownames(nzcoef)[-1], ints)))]
                               } else {
                                 inters <- character(0)
                                 mains <- colnames(draw[["xtrain_lasso"]])[as.numeric(as.character(rownames(nzcoef)[-1]))]
                               }

                               active <- if (length(c(mains, inters)) == 0) " " else c(mains, inters)
                               yvalid_hat <- tryCatch({
                                 as.matrix(predict(fitBT, newx = draw[["xvalid_lasso"]], s = lambda.min, iter = iter.min, type = "response"))
                               },
                               error = function(err) {
                                 return(matrix(0, nrow = nrow(draw[["xvalid_lasso"]]), ncol = 1))
                               } # return NULL on error
                               )

                               msevalid <- mean((draw[["yvalid"]] - drop(yvalid_hat))^2)

                               return(list(beta = coefBT,
                                           # cvfit = fitBT,
                                           vnames = draw[["vnames_lasso"]],
                                           nonzero_coef = nzcoef,
                                           active = active,
                                           not_active = setdiff(draw[["vnames"]], active),
                                           yvalid_hat = yvalid_hat,
                                           msevalid = msevalid,
                                           causal = draw[["causal"]],
                                           not_causal = draw[["not_causal"]],
                                           yvalid = draw[["yvalid"]]))
                             },
                             error = function(err) {
                               return(error_return)
                             }
                             )
                           })

GLinternet <- new_method("GLinternet", "GLinternet",
                         method = function(model, draw) {

                           tryCatch({
                             cvglinternet <- glinternet.cv(X = draw[["xtrain_lasso"]], Y = draw[["ytrain"]],
                                                           numLevels = rep(1, ncol(draw[["xtrain_lasso"]])),
                                                           nLambda = 100, interactionCandidates = c(1),
                                                           verbose = F)
                             tc <- coef(cvglinternet, lambdaType = ifelse(model$lambda.type=="lambda.min", "lambdaHat","lambdaHat1Std"))

                             mains <- colnames(draw[["xtrain_lasso"]])[tc$mainEffects$cont]
                             inters <- paste0(colnames(draw[["xtrain_lasso"]])[tc$interactions$contcont[,2]],":E")

                             # on the training
                             # crossprod(predict(cvglinternet$glinternetFit, X = EX, lambda = cvglinternet$lambdaHat) - DT$y) / n
                             # crossprod(predict(cvglinternet, X = EX) - DT$y) / n

                             if (model$lambda.type=="lambda.min") {
                               cvmse <- cvglinternet$cvErr[which(cvglinternet$lambdaHat==cvglinternet$lambda)]
                             } else {
                               cvmse <- cvglinternet$cvErr[which(cvglinternet$lambdaHat1Std==cvglinternet$lambda)]
                             }

                             return(list(beta = NULL,
                                         # cvfit = cvglinternet,
                                         nonzero_coef = NULL,
                                         active = c(mains, inters),
                                         not_active = setdiff(draw[["vnames"]], c(mains, inters)),
                                         yhat_test = predict(cvglinternet, X = draw[["xtest_lasso"]], type = "response",
                                                        lambdaType = ifelse(model$lambda.type=="lambda.min", "lambdaHat","lambdaHat1Std")),
                                         cvmse = cvmse,
                                         causal = draw[["causal"]],
                                         not_causal = draw[["not_causal"]],
                                         ytest = draw[["ytest"]]))
                           },
                           error = function(err) {
                             return(error_return)
                           }
                           )
                         })


GLinternetsplit <- new_method("GLinternet", "GLinternet",
                         method = function(model, draw) {

                           tryCatch({
                             fitGL <- glinternet(X = draw[["xtrain_lasso"]], Y = draw[["ytrain"]],
                                                 numLevels = rep(1, ncol(draw[["xtrain_lasso"]])),
                                                 nLambda = 100, interactionCandidates = c(1),
                                                 verbose = F)

                             ytest_hat <- predict(fitGL, X = draw[["xtest_lasso"]])
                             msetest <- colMeans((draw[["ytest"]] - ytest_hat)^2)
                             lambda.min.index <- as.numeric(which.min(msetest))
                             lambda.min <- fitGL$lambda[which.min(msetest)]

                             yvalid_hat <- predict(fitGL, X = draw[["xvalid_lasso"]], lambda = lambda.min)
                             msevalid <- mean((draw[["yvalid"]] - drop(yvalid_hat))^2)

                             # nzcoef <- coef(fit, s = lambda.min)[nonzeroCoef(coef(fit, s = lambda.min)),,drop=F]

                             tc <- coef(fitGL, lambdaIndex = lambda.min.index)
                             mains <- colnames(draw[["xtrain_lasso"]])[tc[[1]]$mainEffects$cont]
                             inters <- paste0(colnames(draw[["xtrain_lasso"]])[tc[[1]]$interactions$contcont[,2]],":E")

                             return(list(beta = NULL,
                                         vnames = draw[["vnames_lasso"]],
                                         nonzero_coef = NULL,
                                         active = c(mains, inters),
                                         not_active = setdiff(draw[["vnames"]], c(mains, inters)),
                                         yvalid_hat = yvalid_hat,
                                         msevalid = msevalid,
                                         causal = draw[["causal"]],
                                         not_causal = draw[["not_causal"]],
                                         yvalid = draw[["yvalid"]]))
                           },
                           error = function(err) {
                             return(error_return)
                           }
                           )
                         })





Hiersplit <- new_method("HierBasis", "HierBasis",
                         method = function(model, draw) {

                           tryCatch({

                             fitHier <- AdditiveHierBasis(x = draw[["xtrain_lasso"]],
                                                          y = draw[["ytrain"]],
                                                          type = "gaussian",
                                                          nlam = 100)

                             ytest_hat <- predict(fitHier, new.x = draw[["xtest_lasso"]])
                             msetest <- colMeans((draw[["ytest"]] - ytest_hat)^2)
                             lambda.min.index <- as.numeric(which.min(msetest))
                             lambda.min <- fitHier$lam[which.min(msetest),]

                             yvalid_hat <- predict(fitHier, new.x = draw[["xvalid_lasso"]])[,lambda.min.index]
                             msevalid <- mean((draw[["yvalid"]] - drop(yvalid_hat))^2)

                             components <- view.model.addHierBasis(fitHier, lam.index = lambda.min.index)
                             nonzeroInd <- sapply(components, function(i) if (i[1]=="zero function") FALSE else TRUE)
                             active <- colnames(draw[["xtrain_lasso"]])[nonzeroInd]

                             return(list(beta = NULL,
                                         vnames = draw[["vnames_lasso"]],
                                         nonzero_coef = NULL,
                                         active = active,
                                         not_active = setdiff(colnames(draw[["xtrain_lasso"]]), active),
                                         yvalid_hat = yvalid_hat,
                                         msevalid = msevalid,
                                         causal = draw[["causal"]],
                                         not_causal = draw[["not_causal"]],
                                         yvalid = draw[["yvalid"]]))
                           },
                           error = function(err) {
                             return(error_return)
                           }
                           )
                         })



SPAMsplit <- new_method("SPAM", "SPAM",
                        method = function(model, draw) {
                          tryCatch({
                            fitspam <- samQL(X = draw[["xtrain_lasso"]],
                                             y = draw[["ytrain"]], p = 5, nlambda = 100)

                            ytest_hat <- predict(fitspam, newdata = draw[["xtest_lasso"]])$values
                            msetest <- colMeans((draw[["ytest"]] - ytest_hat)^2)
                            lambda.min.index <- as.numeric(which.min(msetest))
                            lambda.min <- fitspam$lambda[which.min(msetest)]

                            yvalid_hat <- predict(fitspam, newdata = draw[["xvalid_lasso"]])$values[,lambda.min.index]
                            msevalid <- mean((draw[["yvalid"]] - drop(yvalid_hat))^2)

                            coefs <- fitspam$w[,lambda.min.index,drop=FALSE]
                            dimnames(coefs)[[1]] <- paste(rep(draw[["vnames_lasso"]], each = fitspam$p),
                                                          rep(seq_len(fitspam$p), times = length(draw[["vnames_lasso"]])),
                                                          sep = "_")
                            active <- unique(gsub("\\_\\d*", "", names(which(abs(coefs[, 1]) > 0))))


                            return(list(beta = coefs,
                                        # fit = fit,
                                        vnames = draw[["vnames"]],
                                        nonzero_coef = NULL,
                                        active = active,
                                        not_active = setdiff(colnames(draw[["xtrain_lasso"]]), active),
                                        yvalid_hat = yvalid_hat,
                                        msevalid = msevalid,
                                        causal = draw[["causal"]],
                                        not_causal = draw[["not_causal"]],
                                        yvalid = draw[["yvalid"]]))
                          },
                          error = function(err) {
                            return(error_return)
                          }
                          )
                        })


gamselsplit <- new_method("gamsel", "gamsel",
                        method = function(model, draw) {
                          tryCatch({
                            bases <- pseudo.bases(draw[["xtrain_lasso"]])
                            fitgamsel <- gamsel(draw[["xtrain_lasso"]],
                                                draw[["ytrain"]],
                                                bases = bases, num_lambda = 100, family = "gaussian")

                            ytest_hat <- predict(fitgamsel, newdata = draw[["xtest_lasso"]], type = "response")
                            msetest <- colMeans((draw[["ytest"]] - ytest_hat)^2)
                            lambda.min.index <- as.numeric(which.min(msetest))
                            lambda.min <- fitgamsel$lambda[which.min(msetest)]

                            yvalid_hat <- predict(fitgamsel, newdata = draw[["xvalid_lasso"]], index = lambda.min.index,
                                                  type = "response")
                            msevalid <- mean((draw[["yvalid"]] - drop(yvalid_hat))^2)

                            active <- colnames(draw[["xtrain_lasso"]])[predict(fitgamsel,
                                                                               index = lambda.min.index, type = "nonzero")[[1]]]

                            return(list(beta = NA,
                                        # fit = fit,
                                        vnames = draw[["vnames"]],
                                        nonzero_coef = NULL,
                                        active = active,
                                        not_active = setdiff(colnames(draw[["xtrain_lasso"]]), active),
                                        yvalid_hat = yvalid_hat,
                                        msevalid = msevalid,
                                        causal = draw[["causal"]],
                                        not_causal = draw[["not_causal"]],
                                        yvalid = draw[["yvalid"]]))
                          },
                          error = function(err) {
                            return(error_return)
                          }
                          )
                        })

# pacman::p_load(mgcv)
# their_gam <- new_method("gam", "Gam",
#                            method = function(model, draw) {
#                              dat <- data.frame(Y = draw, model$X, E = model$E)
#                              gamfit <- gam(Y ~ E + s(X1) + s(X2) + s(X3)+ s(X4)+ s(X5)+ s(X6)+ s(X7) +
#                                              s(X8) + s(X9) + s(X10) +
#                                              ti(X1, E) + ti(X2, E) + ti(X3, E) + ti(X4, E) + ti(X5, E) + ti(X6, E) +
#                                              ti(X7, E) + ti(X8, E) + ti(X9, E) + ti(X10, E),
#                                            data = dat, select=TRUE, method="REML")
#                              list(gamfit = gamfit,
#                                   yhat = predict(gamfit),
#                                   beta = coef(gamfit),
#                                   y = draw)
#                            })




