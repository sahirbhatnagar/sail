## @knitr metrics

sqrerr <- new_metric("sqrerr", "squared error",
                     metric = function(model, out) {
                       colMeans(as.matrix(out$beta - model$true_beta)^2)
                     })

mse <- new_metric("mse", "MSE",
                   metric = function(model, out) {
                     # as.numeric(sqrt(crossprod(out$y - out$yhat)))
                     as.numeric(crossprod(out$y - out$yhat) / model$n)
                   })


cvmse <- new_metric("cvmse", "10-Fold CV MSE",
                   metric = function(model, out) {
                     as.numeric(out$cvmse)
                   })

tpr <- new_metric("tpr", "True Positive Rate",
                  metric = function(model, out) {
                    length(intersect(out$active, model$causal))/length(model$causal)
                  })

r2 <- new_metric("r2", "R squared",
                 metric = function(model, out) {
                   cor(out$y,as.vector(out$yhat))^2
                 })

"%ni%" <- Negate("%in%")

fpr <- new_metric("fpr", "False Positive Rate",
                  metric = function(model, out){
                    active <- out$active
                    FPR <- sum(active %ni% model$causal) / length(model$not_causal)
                    FPR
                  })

correct_sparsity <- new_metric("cs", "Correct Sparsity",
                               metric = function(model, out){
                                 correct_nonzeros <- sum(out$active %in% model$causal)
                                 correct_zeros <- length(setdiff(model$not_causal,out$active))
                                 #correct sparsity
                                 (1 / length(model$vnames)) * (correct_nonzeros + correct_zeros)

                               })

nactive <- new_metric("nactive", "Number of Active Variables",
                    metric = function(model, out) {
                      length(out$active)
                    })
