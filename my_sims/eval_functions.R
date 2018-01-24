## @knitr metrics

sqrerr <- new_metric("sqrerr", "squared error",
                     metric = function(model, out) {
                       colMeans(as.matrix(out$beta - model$true_beta)^2)
                     })

rmse <- new_metric("rmse", "Root mean squared error",
                   metric = function(model, out) {
                     as.numeric(sqrt(crossprod(out$y - out$yhat)))
                   })


tpr <- new_metric("tpr", "True Positive Rate",
                  metric = function(model, out) {
                    length(intersect(out$nonzero_names, model$causal))/length(model$causal)
                  })

r2 <- new_metric("r2", "R squared",
                 metric = function(model, out) {
                   cor(out$y,as.vector(out$yhat))^2
                 })

"%ni%" <- Negate("%in%")

fpr <- new_metric("fpr", "False Positive Rate",
                  metric = function(model, out){
                    # these are the terms which the model identified as zero
                    modelIdentifyZero <- setdiff(colnames(model$design), out$nonzero_names)
                    FPR <- sum(out$nonzero_names %ni% model$causal) / (sum(out$nonzero_names %ni% model$causal) + sum(modelIdentifyZero %in% model$not_causal))
                    FPR
                  })
