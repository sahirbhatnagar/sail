## @knitr methods
# library(doMC)
# registerDoMC(cores = 8)

sail <- new_method("sail", "Sail",
                   method = function(model, draw) {

                     cvfit <- cv.sail(x = draw[["xtrain"]], y = draw[["ytrain"]], e = draw[["etrain"]],
                                      basis = function(i) splines::bs(i, degree = 5),
                                      # basis = function(i) i,
                                      dfmax = 20,
                                      parallel = TRUE, nfolds = 10, nlambda = 100)

                     nzcoef <- coef(cvfit, s = model$lambda.type)[nonzero(coef(cvfit, s = model$lambda.type)),,drop=F]

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
                   })


gbm <- new_method("gbm", "GBM",
                   method = function(model, draw) {

                     gbmfit <- gbm::gbm.fit(x = draw[["xtrain_lasso"]], y = draw[["ytrain"]],
                                            distribution = "gaussian", verbose = FALSE,
                                            n.trees = 100, #default
                                            interaction.depth = 1)#default)

                     yhat_test <- gbm::predict.gbm(gbmfit,
                                                   newdata = draw[["xtest"]],
                                                   n.trees = 100)

                     topgbm <- summary(gbmfit, n.trees = 100)
                     active <- rownames(topgbm[which(topgbm$rel.inf>0),,drop=FALSE])


                     list(beta = NA,
                          # cvfit = cvfit,
                          vnames = draw[["vnames_lasso"]],
                          nonzero_coef = NA,
                          active = active,
                          not_active = setdiff(draw[["vnames_lasso"]], active),
                          yhat_test = yhat_test,
                          cvmse = NA,
                          causal = draw[["causal"]],
                          not_causal = draw[["not_causal"]],
                          ytest = draw[["ytest"]])
                   })


lasso <- new_method("lasso", "Lasso",
                    method = function(model, draw) {
                      fitglmnet <- cv.glmnet(x = draw[["xtrain_lasso"]], y = draw[["ytrain"]],
                                             alpha = 1, nfolds = 10)

                      nzcoef <- coef(fitglmnet, s = model$lambda.type)[nonzeroCoef(coef(fitglmnet, s = model$lambda.type)),,drop=F]

                      list(beta = coef(fitglmnet, s = model$lambda.type)[-1,,drop=F],
                           # cvfit = fitglmnet,
                           vnames = draw[["vnames_lasso"]],
                           nonzero_coef = nzcoef,
                           active = setdiff(rownames(nzcoef), c("(Intercept)")),
                           not_active = setdiff(colnames(draw[["xtrain_lasso"]]), setdiff(rownames(nzcoef), c("(Intercept)"))),
                           yhat_test = predict(fitglmnet, newx = draw[["xtest_lasso"]], s = model$lambda.type),
                           cvmse = fitglmnet$cvm[which(fitglmnet[[model$lambda.type]]==fitglmnet$lambda)],
                           causal = draw[["causal"]],
                           not_causal = draw[["not_causal"]],
                           ytest = draw[["ytest"]])
                    })

# lassoBT only gives the minimum CV error.. doesnt have lambda.1se
lassoBT <- new_method("lassoBT", "LassoBT",
                      method = function(model, draw) {

                        # fitglmnet <- cv.glmnet(x = model$X_linear_design, y = draw, alpha = 1, nfolds = 10)
                        fitBT <- cvLassoBT(x = draw[["xtrain_lasso"]],
                                           y = draw[["ytrain"]], iter_max=10, nperms=1, nfolds = 10, nlambda = 100)

                        coefBT <- as.matrix(predict(fitBT$BT_fit, type = "coef",
                                                    s = fitBT$lambda[fitBT$cv_opt[1]], iter = fitBT$cv_opt[2]))
                        nzcoef <- coefBT[nonzero(coefBT),,drop=F]

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

                        list(beta = coefBT,
                             # cvfit = fitBT,
                             vnames = draw[["vnames"]],
                             nonzero_coef = nzcoef,
                             active = active,
                             not_active = setdiff(model$vnames, active),
                             yhat_test = yhat,
                             cvmse = fitBT$cv_opt_err,
                             causal = draw[["causal"]],
                             not_causal = draw[["not_causal"]],
                             ytest = draw[["ytest"]])
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
                             return(list(beta = NULL,
                                         # cvfit = NULL,
                                         nonzero_coef = NULL,
                                         active = c(""),
                                         not_active = draw[["vnames"]],
                                         yhat_test = NA,
                                         cvmse = NA,
                                         causal = draw[["causal"]],
                                         not_causal = draw[["not_causal"]],
                                         ytest = draw[["ytest"]]))

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




