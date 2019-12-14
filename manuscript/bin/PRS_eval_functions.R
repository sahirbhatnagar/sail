## @knitr metrics

sqrerr <- new_metric("sqrerr", "squared error",
                     metric = function(model, out) {
                       colMeans(as.matrix(out$beta - model$true_beta)^2)
                     })

# mse <- new_metric("mse", "Test Set MSE",
#                   metric = function(model, out) {
#                     # as.numeric(sqrt(crossprod(out$y - out$yhat)))
#                     as.numeric(crossprod(out$ytest - out$yhat_test) / (length(out$ytest)))
#                   })

msevalid <- new_metric("mse", "Validation Set MSE",
                       metric = function(model, out) {
                         # as.numeric(sqrt(crossprod(out$y - out$yhat)))
                         as.numeric(out$msevalid)
                         # as.numeric(crossprod(out$yvalid - out$yvalid_hat) / (length(out$yvalid)))
                       })


cvmse <- new_metric("cvmse", "10-Fold CV MSE",
                    metric = function(model, out) {
                      as.numeric(out$cvmse)
                    })

tpr <- new_metric("tpr", "True Positive Rate",
                  metric = function(model, out) {
                    length(intersect(out$active, out$causal))/length(out$causal)
                  })

r2 <- new_metric("r2", "R squared",
                 metric = function(model, out) {
                   cor(out$yvalid,as.vector(out$yvalid_hat))^2
                 })

"%ni%" <- Negate("%in%")

fpr <- new_metric("fpr", "False Positive Rate",
                  metric = function(model, out){
                    active <- out$active
                    FPR <- sum(active %ni% out$causal) / length(out$not_causal)
                    FPR
                  })

correct_sparsity <- new_metric("cs", "Correct Sparsity",
                               metric = function(model, out){
                                 correct_nonzeros <- sum(out$active %in% out$causal)
                                 correct_zeros <- length(setdiff(out$not_causal,out$active))
                                 #correct sparsity
                                 (1 / length(model$vnames)) * (correct_nonzeros + correct_zeros)

                               })

nactive <- new_metric("nactive", "Number of Active Variables",
                      metric = function(model, out) {
                        length(out$active)
                      })
