## @knitr methods

devtools::load_all()
pacman::p_load(glmnet)
pacman::p_load(doParallel)
doParallel::registerDoParallel(cores = 5)
getDoParWorkers()

my_sail <- new_method("sail", "Sail",
                        method = function(model, draw) {
                          cvfit <- cv.sail(x = model$X, y = draw, e = model$E,
                                           df = 5, maxit = 500,
                                           nlambda.gamma = 10, nlambda.beta = 10,
                                           nlambda = 100,
                                           # group.penalty = "SCAD",
                                           # lambda.beta = exp(seq(log(0.01 * 1000), log(1000), length.out = 25)),
                                           # lambda.gamma =rep(1000,25),
                                           lambda.factor = 0.0001,
                                           nfolds = 5,
                                           thresh = 1e-3, center=T, normalize=F, verbose = T)
                          list(beta = coef(cvfit, s = "lambda.min")[-1,,drop=F],
                               cvfit = cvfit,
                               nonzero = coef(cvfit, s = "lambda.min")[nonzeroCoef(coef(cvfit, s = "lambda.min")),,drop=F],
                               nonzero_names = setdiff(rownames(coef(cvfit, s = "lambda.min")[nonzeroCoef(coef(cvfit, s = "lambda.min")),,drop=F]),c("(Intercept)")),
                               yhat = cbind(1,model$design) %*% coef(cvfit, s = "lambda.min"),
                               y = draw)
                        })


pacman::p_load(mgcv)
their_gam <- new_method("gam", "Gam",
                           method = function(model, draw) {
                             dat <- data.frame(Y = draw, model$X, E = model$E)
                             gamfit <- gam(Y ~ E + s(X1) + s(X2) + s(X3)+ s(X4)+ s(X5)+ s(X6)+ s(X7) +
                                             s(X8) + s(X9) + s(X10) +
                                             ti(X1, E) + ti(X2, E) + ti(X3, E) + ti(X4, E) + ti(X5, E) + ti(X6, E) +
                                             ti(X7, E) + ti(X8, E) + ti(X9, E) + ti(X10, E),
                                           data = dat, select=TRUE, method="REML")
                             list(gamfit = gamfit,
                                  yhat = predict(gamfit),
                                  beta = coef(gamfit),
                                  y = draw)
                           })


lasso <- new_method("lasso", "Lasso",
                    method = function(model, draw) {
                      fitglmnet <- cv.glmnet(x = model$design, y = draw, alpha = 1, standardize = F, nfolds = 5)
                      list(beta = coef(fitglmnet, s = "lambda.min")[-1,,drop=F],
                           yhat = predict(fitglmnet, newx = model$design, s = "lambda.min"),
                           nonzero = coef(fitglmnet)[nonzeroCoef(coef(fitglmnet)),,drop=F],
                           nonzero_names = setdiff(rownames(coef(fitglmnet)[nonzeroCoef(coef(fitglmnet)),,drop=F]),c("(Intercept)")),
                           y = draw)
                    })

