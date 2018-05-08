# This is the main simulator file
# I ran this on tmux by copy pasting into an rsession on hydra, and chaning the index argument
# in simulate from model

# use this code to save the results. turns out you cant run simulator in parallel "by hand"
# you need to load the simulator object prior to adding new simulations to the object.
# so for now I just saved the evals object in separate .rds file in mcgillsims/sim_results
# source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/my_sims/eval_functions.R")
# sim <- sim %>% evaluate(list(mse, cvmse, r2, tpr, fpr, correct_sparsity,nactive))
# res <- as.data.frame(evals(sim))
# saveRDS(res, file = "sail/sail_lambda_branch/mcgillsims/sim_results/res7.rds")


# rm(list=ls())
pacman::p_load(simulator) # this file was created under simulator version 0.2.1
pacman::p_load(splines)
pacman::p_load(magrittr)
# pacman::p_load(foreach)
pacman::p_load(methods)
# pacman::p_load(doMC)
pacman::p_load(glmnet)
pacman::p_load(LassoBacktracking)
pacman::p_load(glinternet)
pacman::p_load(gbm)
# pacman::p_load_gh('sahirbhatnagar/sail', dependencies = FALSE)
library(sail)
pacman::p_load_gh('asadharis/HierBasis')
pacman::p_load(SAM)
pacman::p_load(gamsel)
pacman::p_load(cowplot)
pacman::p_load(dplyr)
pacman::p_load(tidyr)
pacman::p_load(data.table)
pacman::p_load(ggplot2)
pacman::p_load(latex2exp)
pacman::p_load(lemon)
# registerDoMC(cores = 10)
# Index <- as.numeric(as.character(commandArgs(trailingOnly = T)[1]))

# source("/home/sahir/git_repositories/sail/my_sims/model_functions.R")
# source("/home/sahir/git_repositories/sail/my_sims/method_functions.R")
# source("/home/sahir/git_repositories/sail/my_sims/eval_functions.R")

source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_git_v2/sail/my_sims/model_functions.R")
source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_git_v2/sail/my_sims/method_functions.R")
source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_git_v2/sail/my_sims/eval_functions.R")


# run simulation in parallel on tmux --------------------------------------

sim <- new_simulation(name = "apr_25_2018",
                      label = "apr_25_2018",
                      dir = ".") %>%
  generate_model(make_gendata_Paper_data_split, seed = 1234,
                 n = 400, p = 1000, corr = 0, betaE = 2, SNR = 2, lambda.type = "lambda.min",
                 parameterIndex = list(1,2,3,4,5),
                 vary_along = "parameterIndex") %>%
  simulate_from_model(nsim = 6, index = 1:35) %>%
  run_method(list(sailsplit, lassosplit, lassoBTsplit, GLinternetsplit, Hiersplit, SPAMsplit, gamselsplit),
             parallel = list(socket_names = 35,
                             libraries = c("LassoBacktracking", "glinternet","glmnet","splines",
                                           "magrittr","sail","gamsel","SAM","HierBasis","simulator", "parallel")))
simulator::save_simulation(sim)

s2 <- new_simulation(name = "may_1_2018",
                      label = "may_1_2018",
                      dir = ".") %>%
  generate_model(make_gendata_Paper_data_split, seed = 1234,
                 n = 400, p = 50, corr = 0, betaE = 2, SNR = 2, lambda.type = "lambda.min",
                 parameterIndex = list(1,2),
                 vary_along = "parameterIndex") %>%
  simulate_from_model(nsim = 2, index = 1:2) %>%
  run_method(list(sailsplitadaptive, sailsplit, sailsplitweak, sailsplitadaptiveweak))


s2 <- s2 %>% evaluate(list(msevalid, tpr, fpr, nactive, r2))
s2 %>% plot_eval_by(metric_name = "fpr", varying = "parameterIndex")
# load simulation ---------------------------------------------------------

sim <- load_simulation("apr_25_2018")
sim <- sim %>% evaluate(list(msevalid, tpr, fpr, nactive, r2))
# simulator::save_simulation(sim)
# sim <- sim %>% run_method(list(sailsplitweak, sailsplitadaptiveweak),
#                           parallel = list(socket_names = 35,
#                                           libraries = c("LassoBacktracking", "glinternet","glmnet","splines",
#                                                         "magrittr","sail","gamsel","SAM","HierBasis","simulator", "parallel")))
# simulator::save_simulation(sim)
# ls()

# analyze results ---------------------------------------------------------

# df <- as.data.frame(evals(sim))
# saveRDS(df, file = "my_sims/simulation_results/apr_25_2018_results.rds")
df <- readRDS("my_sims/simulation_results/apr_25_2018_results.rds")
df <- df %>% separate(Model, into = c("simnames","betaE","corr","lambda.type","n","p","parameterIndex","SNR_2"),
                      sep = "/")
DT <- as.data.table(df)

trop <- RSkittleBrewer::RSkittleBrewer("trop")

DT[parameterIndex=="parameterIndex_1", table(Method)]


cbbPalette <- c("#8720B6","#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
trop <- RSkittleBrewer::RSkittleBrewer("trop")
gg_sy <- theme(legend.position = "bottom", axis.text = element_text(size = 20),
               axis.title = element_text(size = 20), legend.text = element_text(size = 20),
               legend.title = element_text(size = 20),plot.title = element_text(size = 20) )

appender <- function(string) TeX(paste(string))

DT[, scenario:= as.numeric(as.character(stringr::str_extract_all(parameterIndex, "\\d", simplify = T)))]
DT[, scen:=ifelse(scenario==1,"Strong Hierarchy",ifelse(scenario==2, "Weak Hierarchy", ifelse(scenario==3,"Interactions Only",ifelse(scenario==4, "Strong Hierarchy (Linear)", "Main Effects Only"))))]
DT[, scen:=factor(scen, levels = c("Strong Hierarchy", "Weak Hierarchy","Interactions Only","Strong Hierarchy (Linear)", "Main Effects Only"))]
DT$scen %>% table
#Truth obeys strong hierarchy (parameterIndex = 1)
#Truth obeys weak hierarchy (parameterIndex = 2)
#Truth only has interactions (parameterIndex = 3)
#Truth is linear (parameterIndex = 4)
#Truth only has main effects (parameterIndex = 5)


(p1_mse <- ggplot(DT, aes(Method, mse, fill = Method)) +
    ggplot2::geom_boxplot() +
    facet_rep_wrap(~scen, scales = "free", ncol = 2,repeat.tick.labels = 'left',
               labeller = as_labeller(appender,
                                      default = label_parsed))+
    scale_fill_manual(values=c(cbbPalette, "red","pink"), guide=guide_legend(ncol=2)) +
    ggplot2::labs(y = "Test Set MSE", title = "") + xlab("") + panel_border()+background_grid()+
    theme(legend.position = "right", legend.text=element_text(size=18)) )

reposition_legend(p1_mse, 'center', panel='panel-2-3')



cowplot::background_grid()
save_plot("mcgillsims/figures/p1_cvmse.pdf", p1_cvmse,
          base_height = 7, base_width = 9)






# Not used underneath this line -------------------------------------------








sim %>% plot_eval_by(metric_name = "mse", "parameterIndex")



evals(sim)@evals



subset_simulation(sim, parameterIndex == 1) %>% plot_evals("fpr", "tpr")
simulator::load_outputs()
o <- output(sim)
o[[1]][[1]]@out$r1.1$beta

# sim %>% plot_eval(metric_name = "cvmse")
# sim %>% plot_eval(metric_name = "nactive")
# sim %>% plot_eval(metric_name = "tpr")

  # run_method(list(sail, lasso, GLinternet, lassoBT, gbm),
  #            parallel = list(socket_names = 5,
  #                            libraries = c("LassoBacktracking", "glinternet","glmnet","splines","magrittr","sail","gbm")))

simulator::save_simulation(sim)


# Data split simulation ---------------------------------------------------
rm(sim)

sim <- new_simulation(name = "sail_lasso_glinternet_lassoBT_gbm_split_testing",
                      label = "sail v7_split",
                      dir = ".") %>%
  generate_model(make_gendata_Paper_data_split, seed = 1234,
                 n = 400, p = 20, corr = 0, betaE = 2, SNR = 2, lambda.type = "lambda.min",
                 parameterIndex = list(1),
                 vary_along = "parameterIndex") %>%
  simulate_from_model(nsim = 2, index = 1:10) %>%
  # run_method(list(SPAMsplit)) %>%
  # run_method(list(sailsplit, lassosplit, lassoBTsplit, GLinternetsplit, Hiersplit, SPAMsplit, gamselsplit)) %>%
  run_method(list(sailsplit, lassosplit, lassoBTsplit, GLinternetsplit, Hiersplit, SPAMsplit, gamselsplit),
             parallel = list(socket_names = 10,
                             libraries = c("LassoBacktracking", "glinternet","glmnet","splines","magrittr","sail","gamsel","SAM","HierBasis"))) %>%
  evaluate(list(msevalid)) #%>%

sim <- sim %>%
  run_method(list(lassosplit, GLinternetsplit, lassoBTsplit)) %>%
  evaluate(list(msevalid, tpr, fpr, nactive,r2))

sim <- sim %>%
  run_method(list(gamselsplit)) %>%
  evaluate(list(msevalid, tpr, fpr, nactive,r2))

sim %>%
  plot_eval(metric_name = "r2")

as.data.frame(evals(sim))
  # plot_eval(metric_name = "msevalid")
  #            parallel = list(socket_names = 5,
  #                            libraries = c("LassoBacktracking", "glinternet","glmnet","splines","magrittr","sail","gbm")))

# simulator::save_simulation(sim)


xtrain <- draws(sim)@draws$r2.2[["xtrain"]]
xtrain_lasso <- draws(sim)@draws$r2.2[["xtrain_lasso"]]
etrain <- draws(sim)@draws$r2.2[["etrain"]]
ytrain <- draws(sim)@draws$r2.2[["ytrain"]]

xtest <- draws(sim)@draws$r2.2[["xtest"]]
xtest_lasso <- draws(sim)@draws$r2.2[["xtest_lasso"]]
etest <- draws(sim)@draws$r2.2[["etest"]]
ytest <- draws(sim)@draws$r2.2[["ytest"]]

xvalid <- draws(sim)@draws$r2.2[["xvalid"]]
xvalid_lasso <- draws(sim)@draws$r2.2[["xvalid_lasso"]]
evalid <- draws(sim)@draws$r2.2[["evalid"]]
yvalid <- draws(sim)@draws$r2.2[["yvalid"]]
vnames <- draws(sim)@draws$r2.2[["vnames_lasso"]]
# no built-in tuning
pacman::p_load_gh('asadharis/HierBasis')
fitHier <- AdditiveHierBasis(x = xtrain_lasso,
                         y = ytrain,
                         type = "gaussian",
                         nlam = 100)
ytest_hat <- predict(fitHier, new.x = xtest_lasso)
msetest <- colMeans((ytest - ytest_hat)^2)
plot(log(fitHier$lam), msetest)

lambda.min.index <- as.numeric(which.min(msetest))
lambda.min <- fitHier$lam[which.min(msetest),]

yvalid_hat <- predict(fitHier, new.x = xvalid_lasso)[,lambda.min.index]
msevalid <- mean((yvalid - drop(yvalid_hat))^2)

components <- view.model.addHierBasis(fitHier, lam.index = lambda.min.index)
nonzeroInd <- sapply(components, function(i) if (i[1]=="zero function") FALSE else TRUE)
active <- colnames(xtrain_lasso)[nonzeroInd]

# Spam - also no built-in tuning
pacman::p_load(SAM)
fitspam <- samQL(X = xtrain_lasso,
                 y = ytrain, p = 5, nlambda = 100)
plot(fitspam)
print(fitspam)

class(fitspam)
ytest_hat <- predict(fitspam, newdata = xtest_lasso)$values
msetest <- colMeans((ytest - ytest_hat)^2)
plot(log(fitspam$lambda), msetest)

lambda.min.index <- as.numeric(which.min(msetest))
lambda.min <- fitspam$lambda[which.min(msetest)]

# SAM::predict.samQL()
yvalid_hat <- predict(fitspam, newdata = xvalid_lasso)$values[,lambda.min.index]
msevalid <- mean((yvalid - drop(yvalid_hat))^2)

vnames <- draws(sim)@draws$r2.2[["vnames_lasso"]]
coefs <- fitspam$w[,lambda.min.index,drop=FALSE]
dimnames(coefs)[[1]] <- paste(rep(vnames, each = 5), rep(seq_len(5), times = length(vnames)), sep = "_")
active <- unique(gsub("\\_\\d*", "", names(which(abs(coefs[, 1]) > 0))))






# gamsel ------------------------------------------------------------------

pacman::p_load(gamsel)
bases=pseudo.bases(xtrain_lasso)
# Gaussian gam
fitgamsel <- gamsel(xtrain_lasso, ytrain, bases = bases, num_lambda = 100, family = "gaussian")
gamsel::predict.gamsel()
ytest_hat <- predict(fitgamsel, newdata = xtest_lasso, type = "response")
msetest <- colMeans((ytest - ytest_hat)^2)
plot(log(fitgamsel$lambdas), msetest)

lambda.min.index <- as.numeric(which.min(msetest))
lambda.min <- fitgamsel$lambdas[which.min(msetest)]

yvalid_hat <- predict(fitgamsel, newdata = xvalid_lasso, index = lambda.min.index, type = "response")
msevalid <- mean((yvalid - drop(yvalid_hat))^2)

active <- colnames(xtrain_lasso)[predict(fitgamsel, index = lambda.min.index, type = "nonzero")[[1]]]


par(mfrow=c(1,2),mar=c(5,4,3,1))
summary(gamsel.out)
gamsel::getActive(gamsel.out)

par(mfrow=c(2,2))
plot(gamsel.out, newx=X[,getActive(gamsel.out, index = 24, type = "nonlinear")$l24],index=24)
preds<-predict(gamsel.out, X, index = 24, type = "terms")









SAM::predict.samQL(fitspam, newdata = xtest_lasso)
fitspam
predict(fitspam, draws(sim)@draws[["r1.1"]][["xtest_lasso"]])



fitBT <- LassoBT(x = xtrain_lasso,
                 y = ytrain, iter_max=10, nlambda = 100)

ytest_hat <- predict(fitBT, newx = xtest_lasso)

fitGL <- glinternet(X = xtrain_lasso, Y = ytrain,
                           numLevels = rep(1, ncol(xtrain_lasso)),
                           nLambda = 100, interactionCandidates = c(1),
                           verbose = F)

ytest_hat <- predict(fitGL, X = xtest_lasso)
msetest <- colMeans((ytest - ytest_hat)^2)
lambda.min.index <- as.numeric(which.min(msetest))
lambda.min <- fitGL$lambda[which.min(msetest)]
plot(log(fitGL$lambda), msetest)
yvalid_hat <- predict(fitGL, X = xvalid_lasso, lambda = lambda.min)
msevalid <- mean((yvalid - drop(yvalid_hat))^2)
cor(yvalid, as.vector(yvalid_hat))^2
coef(fitGL, lambdaIndex = lambda.min.index)

tc <- coef(fitGL, lambdaIndex = lambda.min.index)

mains <- colnames(xtrain_lasso)[tc[[1]]$mainEffects$cont]
inters <- paste0(colnames(xtrain_lasso)[tc[[1]]$interactions$contcont[,2]],":E")


fit <- sail(x = xtrain, y = ytrain, e = etrain, basis = function(i) splines::bs(i, degree = 5),
            verbose = 1)
plot(fit)

fit <- glmnet(x = xtrain_lasso, y = ytrain, alpha = 1)

DT <- sail::gendata(n = 200, p = 40, corr = 0, SNR = 2, betaE = 2, parameterIndex = 1)
DT_test <- sail::gendata(n = 200, p = 40, corr = 0, SNR = 2, betaE = 2, parameterIndex = 1)
x_lasso <- cbind(E = DT$e, DT$x)
x_lasso_test <- cbind(E = DT_test$e, DT_test$x)

fit <- LassoBT(x = x_lasso,
               y = DT$y, iter_max=10, nlambda = 100)
predict(fit)

ytest_hat <- predict(fit, newx = x_lasso, newe = etest)
ytest_hat <- predict(fit, newx = xtest_lasso)


ytest_hat <- predict(fit, newx = x_lasso+rnorm(x_lasso, sd = 2))
ytest_hat <- predict(fit, newx = x_lasso)
dy <- dim(ytest_hat)
res <- matrix(nrow = dy[3], ncol = dy[2])
for (j in seq(dy[3])) {
  res[j, ] <- colMeans((DT$y - ytest_hat[,,j])^2)
}
all.equal(t(res),msetest)


msetest <- colMeans((DT$y - ytest_hat)^2)
which(msetest == min(msetest), arr.ind = TRUE)
lambda.min.index <- which(msetest==min(msetest), arr.ind = TRUE)
lambda.min <- fitBT$lambda[lambda.min.index[1,1]]
iter.min <- lambda.min.index[1,2]

coefBT <- as.matrix(predict(fitBT, type = "response",
                            s = lambda.min, iter = iter.min))
nzcoef <- coefBT[sail:::nonzero(coefBT),,drop=F]

yvalid_hat <- tryCatch({
  as.matrix(predict(fitBT, newx = x_lasso, s = lambda.min, iter = iter.min, type = "response"))
},
error = function(err) {
  return(matrix(0, nrow = 200, ncol = 1))
} # return NULL on error
)

# m2 <- colMeans(sweep(ytest_hat, 1, ytest, FUN = "-")^2)
# all.equal(msetest, m2)
plot(log(fit$lambda), msetest)
msetest[which.min(msetest)]
lambda.min.index <- as.numeric(which.min(msetest))
lambda.min <- fit$lambda[which.min(msetest)]
abline(v=log(lambda.min))
coef(fit, s = lambda.min)

yvalid_hat <- predict(fit, newx = xvalid_lasso, s = lambda.min)
msevalid <- mean((yvalid - drop(yvalid_hat))^2)

nzcoef <- predict(fit, s = lambda.min, type = "nonzero")
nzcoef <- coef(fit, s = lambda.min)[nonzeroCoef(coef(fit, s = lambda.min)),,drop=F]






# sim <- load_simulation("sail_lasso_glinternet_lassoBT_gbm")

# sim <- sim %>% evaluate(list(mse, cvmse, tpr, fpr, nactive))
# sim %>% plot_eval(metric_name = "mse")
# sim %>% plot_eval(metric_name = "cvmse")
# sim %>% plot_eval(metric_name = "nactive")
# sim %>% plot_eval(metric_name = "tpr")


# sim <- new_simulation(name = "sail_lassoBT_glinternet",
#                        label = "Sail McGill Other", dir = "/mnt/GREENWOOD_SCRATCH/sahir.bhatnagar/sail_simulations/") %>%
#   generate_model(make_gendata_Paper, seed = 123,
#                  n = 200, p = 1000, corr = 0, betaE = 1, SNR = 2, lambda.type = "lambda.min",
#                  parameterIndex = list(1,4,5),
#                  vary_along = "parameterIndex") %>%
#   simulate_from_model(nsim = 2, index = 21) %>%
#   run_method(list(lassoBT, GLinternet),
#              parallel = list(socket_names = 5,
#                         libraries = c("LassoBacktracking", "glinternet","glmnet","splines","magrittr")))


# sim <- sim %>% evaluate(list(mse, cvmse, r2, tpr, fpr, correct_sparsity,nactive))
# res <- as.data.frame(evals(sim))
# # saveRDS(res, file = "sail/sail_lambda_branch/mcgillsims/sim_results/res7.rds")
#
# saveRDS(object = res,
#         file = tempfile(pattern = sprintf("res_%s_",2),
#                         tmpdir = "/mnt/GREENWOOD_SCRATCH/sahir.bhatnagar/sail_simulations/glinternet_lassoBT/",
#                         fileext = ".rds")
# )

# ls()


## @knitr init

# name_of_simulation <- "sail_mcgill_talk"

## @knitr main

# sim <- new_simulation(name = name_of_simulation,
#                       label = "Sail easy p25") %>%
#   generate_model(make_easy_sail_model, seed = 123,
#                  n = 200, p = 25, df = 5, SNR = 2, betaE = 1) %>%
#   simulate_from_model(nsim = 25) %>%
#   run_method(list(lasso, my_sail))

# getOption("simulator.files")

# sim <- new_simulation(name = name_of_simulation,
#                       label = "Sail McGill", dir = "/mnt/GREENWOOD_SCRATCH/sahir.bhatnagar/sail_simulations/") %>%
#   generate_model(make_gendata_Paper, seed = 123,
#                  n = 200, p = 1000, corr = 0, betaE = 1, SNR = 2, lambda.type = "lambda.min",
#                  # parameterIndex = 1
#                  parameterIndex = list(1,2,3,4,5),
#                  vary_along = "parameterIndex"
#                  ) %>%
#   simulate_from_model(nsim = 10, index = 20) %>%
#   run_method(list(lasso, sail))
#
#
# sim <- load_simulation(name_of_simulation, dir = "/mnt/GREENWOOD_SCRATCH/sahir.bhatnagar/sail_simulations/")


# Simulation for LassoBT and Glinternet -----------------------------------

# sim2 <- new_simulation(name = "sail_lassoBT_glinternet",
#                        label = "Sail McGill Other", dir = "/mnt/GREENWOOD_SCRATCH/sahir.bhatnagar/sail_simulations/") %>%
#   generate_model(make_gendata_Paper, seed = 123,
#                  n = 200, p = 1000, corr = 0, betaE = 1, SNR = 2, lambda.type = "lambda.min",
#                  parameterIndex = list(1,2,3,4,5),
#                  vary_along = "parameterIndex") %>%
#   simulate_from_model(nsim = 2, index = 1) %>%
#   run_method(list(lassoBT, GLinternet))
#
#
#
# sim2 <- sim2 %>% evaluate(list(mse,cvmse,tpr, fpr,nactive))
# sim2 %>% plot_eval_by("nactive","parameterIndex")
#
# as.data.frame(evals(sim2))
# sim2 %>% evaluate(list(mse)) %>% evals() %>% as.data.frame()
#
# sim2 %>% plot_eval(metric_name = "cvmse")
# Plot Results ------------------------------------------------------------
# pacman::p_load(cowplot)
# pacman::p_load(tidyverse)
#
# df <- do.call(rbind, lapply(list.files(path = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/mcgillsims/sim_results/",
#                                        pattern = "*.rds", full.names = TRUE), readRDS))
#
# df <- df %>% separate(Model, into = c("simnames","betaE","corr","lambda.type","n","p","parameterIndex","SNR_2"),
#                 sep = "/")
#
#
# trop = RSkittleBrewer::RSkittleBrewer("trop")
#
# df[which(df$parameterIndex=="parameterIndex_1"),]$Method  %>% table
#
# (p1_mse <- ggplot(df[which(df$parameterIndex=="parameterIndex_5"),], aes(Method, fpr, fill = Method)) +
#     ggplot2::geom_boxplot() +
#     # geom_jitter()+
#     scale_fill_manual(values=trop[1:2]) +
#     theme(legend.position = "none",
#           axis.text=element_text(size=16),
#           # axis.text.x=element_text(size=1),
#           axis.title=element_text(size=16)) +
#     ggplot2::labs(y = "10-Fold CV MSE", title = "Simulation Scenario 1") + xlab(""))
#
#
# (p1_mse <- ggplot(df, aes(Method, mse, fill = Method)) +
#     ggplot2::geom_boxplot() +
#     scale_fill_manual(values=trop[1:2]) +
#     theme(legend.position = "none",
#           axis.text=element_text(size=16),
#           # axis.text.x=element_text(size=1),
#           axis.title=element_text(size=16)) +
#     ggplot2::labs(y = "10-Fold CV MSE", title = "Simulation Scenario 1") + xlab(""))

# sim <- sim %>% evaluate(list(cvmse))
#
# sim %>% plot_eval(metric_name = "cvmse")
# sim <- sim %>% evaluate(list(rmse))
#
# sim %>% subset_simulation(methods = c("lasso","sail","gam")) %>% evaluate(list(rmse, tpr,fpr)) %>% plot_eval(metric_name = "rmse")
#
#
# ## @knitr plots
#
# plot_eval_by(sim, "hisloss", varying = "prob")
#
# ## @knitr tables
#
# tabulate_eval(sim, "herloss", output_type = "markdown",
#               format_args = list(digits = 1))
