## @knitr methods

sail <- new_method("sail", "Sail",
                   method = function(model, draw) {

                     draw <- scale(draw, center = TRUE, scale = FALSE)

                     cvfit <- cv.sail(x = model$x, y = draw, e = model$e,
                                      df = 5, degree = 3, basis.intercept = FALSE,
                                      thresh = 1e-3, maxit = 1000, alpha = .2,
                                      dfmax = 75,
                                      parallel = TRUE, nfolds = 10, verbose = T, nlambda = 100)

                     nzcoef <- coef(cvfit, s = model$lambda.type)[nonzero(coef(cvfit, s = model$lambda.type)),,drop=F]

                     list(beta = coef(cvfit, s = model$lambda.type)[-1,,drop=F],
                          cvfit = cvfit,
                          nonzero_coef = nzcoef,
                          active = cvfit$sail.fit$active[[which(cvfit[[model$lambda.type]]==cvfit$lambda)]],
                          not_active = setdiff(model$vnames, cvfit$sail.fit$active[[which(cvfit[[model$lambda.type]]==cvfit$lambda)]]),
                          yhat = predict(cvfit, s = cvfit[[model$lambda.type]]),
                          cvmse = cvfit$cvm[which(cvfit[[model$lambda.type]]==cvfit$lambda)],
                          causal = model$causal,
                          not_causal = model$not_causal,
                          y = draw)
                   })

lasso <- new_method("lasso", "Lasso",
                    method = function(model, draw) {
                      fitglmnet <- cv.glmnet(x = model$X_linear_design, y = draw, alpha = 1, nfolds = 10)

                      nzcoef <- coef(fitglmnet, s = model$lambda.type)[nonzeroCoef(coef(fitglmnet, s = model$lambda.type)),,drop=F]

                      list(beta = coef(fitglmnet, s = model$lambda.type)[-1,,drop=F],
                           cvfit = fitglmnet,
                           nonzero_coef = nzcoef,
                           active = setdiff(rownames(nzcoef), c("(Intercept)")),
                           not_active = setdiff(model$vnames, setdiff(rownames(nzcoef), c("(Intercept)"))),
                           yhat = predict(fitglmnet, newx = model$X_linear_design, s = model$lambda.type),
                           cvmse = fitglmnet$cvm[which(fitglmnet[[model$lambda.type]]==fitglmnet$lambda)],
                           causal = model$causal,
                           not_causal = model$not_causal,
                           y = draw)
                    })

lassoBT <- new_method("lassoBT", "LassoBT",
                      method = function(model, draw) {

                        # fitglmnet <- cv.glmnet(x = model$X_linear_design, y = draw, alpha = 1, nfolds = 10)
                        fitBT <- cvLassoBT(model$EX, draw, iter_max=10, nperms=1, nfolds = 10, nlambda = 100)

                        coefBT <- as.matrix(predict(fitBT$BT_fit, type = "coef", s = fitBT$lambda[fitBT$cv_opt[1]], iter = fitBT$cv_opt[2]))
                        nzcoef <- coefBT[nonzero(coefBT),,drop=F]

                        ints <- grep(":",rownames(nzcoef)[-1], value=T)
                        if (length(ints) != 0) {
                          inters <- sapply(stringr::str_split(ints, ":"), function(i) paste(colnames(model$EX)[as.numeric(as.character(i))], collapse = ":"))
                          mains <- colnames(model$EX)[as.numeric(as.character(setdiff(rownames(nzcoef)[-1], ints)))]
                        } else {
                          inters <- character(0)
                          mains <- colnames(model$EX)[as.numeric(as.character(rownames(nzcoef)[-1]))]
                        }

                        active <- if (length(c(mains, inters)) == 0) " " else c(mains, inters)
                        yhat <- tryCatch({
                          preds <- predict(fitBT$BT_fit, s = fitBT$lambda[fitBT$cv_opt[1]])
                          tf <- as.matrix(preds[,fitBT$cv_opt[1], min(dim(preds)[3],fitBT$cv_opt[2]),drop=TRUE])
                        },
                        error = function(err) {
                          return(matrix(0, nrow = nrow(model$EX), ncol = 1))
                        } # return NULL on error
                        )

                        list(beta = coefBT,
                             cvfit = fitBT,
                             nonzero_coef = nzcoef,
                             active = active,
                             not_active = setdiff(model$vnames, active),
                             yhat = yhat,
                             cvmse = fitBT$cv_opt_err,
                             causal = model$causal,
                             not_causal = model$not_causal,
                             y = draw)
                      })


GLinternet <- new_method("GLinternet", "GLinternet",
                         method = function(model, draw) {

                           tryCatch({
                             cvglinternet <- glinternet.cv(X = model$EX, Y = draw,
                                                           numLevels = rep(1, ncol(model$EX)),
                                                           nLambda = 100, interactionCandidates = c(1),
                                                           verbose = F)
                             tc <- coef(cvglinternet)

                             mains <- colnames(model$EX)[tc$mainEffects$cont]
                             inters <- paste0(colnames(model$EX)[tc$interactions$contcont[,2]],":E")

                             # on the training
                             # crossprod(predict(cvglinternet$glinternetFit, X = EX, lambda = cvglinternet$lambdaHat) - DT$y) / n
                             # crossprod(predict(cvglinternet, X = EX) - DT$y) / n

                             return(list(beta = NULL,
                                         cvfit = cvglinternet,
                                         nonzero_coef = NULL,
                                         active = c(mains, inters),
                                         not_active = setdiff(model$vnames, c(mains, inters)),
                                         yhat = predict(cvglinternet, X = model$EX),
                                         cvmse = cvglinternet$cvErr[which(cvglinternet$lambdaHat==cvglinternet$lambda)],
                                         causal = model$causal,
                                         not_causal = model$not_causal,
                                         y = draw))
                           },
                           error = function(err) {
                             return(list(beta = NULL,
                                         cvfit = NULL,
                                         nonzero_coef = NULL,
                                         active = c(""),
                                         not_active = model$vnames,
                                         yhat = NA,
                                         cvmse = NA,
                                         causal = model$causal,
                                         not_causal = model$not_causal,
                                         y = draw))

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




