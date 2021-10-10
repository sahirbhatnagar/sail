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
# pacman::p_load_gh('sahirbhatnagar/sail')
devtools::load_all()
pacman::p_load(splines)
pacman::p_load(magrittr)
pacman::p_load(foreach)
pacman::p_load(methods)
pacman::p_load(doMC)
pacman::p_load(glmnet)
pacman::p_load(LassoBacktracking)
pacman::p_load(glinternet)
pacman::p_load(gbm)

registerDoMC(cores = 8)
# Index <- as.numeric(as.character(commandArgs(trailingOnly = T)[1]))

source("/home/sahir/git_repositories/sail/my_sims/model_functions.R")
source("/home/sahir/git_repositories/sail/my_sims/method_functions.R")
source("/home/sahir/git_repositories/sail/my_sims/eval_functions.R")
# source("/home/sahir/git_repositories/sail/R/functions.R")
# devtools::load_all("/home/sahir/git_repositories/sail/")

source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/my_sims/model_functions.R")
source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/my_sims/method_functions.R")
source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/my_sims/eval_functions.R")
# source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/R/functions.R")

# devtools::load_all("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/")

sim <- new_simulation(name = "sail_lasso",
                      label = "sail v2",
                      dir = ".") %>%
  generate_model(make_gendata_Paper, seed = 1234,
                 n = 200, p = 20, corr = 0, betaE = 2, SNR = 2, lambda.type = "lambda.min",
                 parameterIndex = list(3),
                 vary_along = "parameterIndex") %>%
  simulate_from_model(nsim = 2, index = 1:2) %>%
  # run_method(list(GLinternet))
  run_method(list(sail, lasso, GLinternet, lassoBT, gbm))#,
             # parallel = list(socket_names = 30,
                             # libraries = c("LassoBacktracking", "glinternet","glmnet","splines","magrittr")))

sim <- sim %>% evaluate(list(mse, cvmse, tpr, fpr, nactive))
sim %>% plot_eval(metric_name = "mse")
sim %>% plot_eval(metric_name = "cvmse")
sim %>% plot_eval(metric_name = "nactive")
sim %>% plot_eval(metric_name = "tpr")


# no built-in tuning
pacman::p_load_gh('asadharis/HierBasis')
fit <- AdditiveHierBasis(x = draws(sim)@draws[["r1.1"]][["xtrain_lasso"]],
                 y = draws(sim)@draws[["r1.1"]][["ytrain"]],
                 type = "gaussian")


# Spam - also no built-in tuning
pacman::p_load(SAM)
fitspam <- samQL(X = draws(sim)@draws[["r1.1"]][["xtrain_lasso"]],
                 y = draws(sim)@draws[["r1.1"]][["ytrain"]])
plot(fitspam)
fitspam
predict(fitspam, draws(sim)@draws[["r1.1"]][["xtest_lasso"]])


cvglinternet <- glinternet.cv(X = draws(sim)@draws[["r1.1"]][["xtrain_lasso"]],
                              Y = draws(sim)@draws[["r1.1"]][["ytrain"]],
                              numLevels = rep(1, ncol(draws(sim)@draws[["r1.1"]][["xtrain_lasso"]])),
                              nLambda = 100, interactionCandidates = c(1),
                              verbose = T)
tc <- coef(cvglinternet)

simulator::model(sim)@params$lambda.type

mains <- colnames(draws(sim)@draws[["r1.1"]][["xtrain_lasso"]])[tc$mainEffects$cont]
inters <- paste0(colnames(draws(sim)@draws[["r1.1"]][["xtrain_lasso"]])[tc$interactions$contcont[,2]],":E")
predict(cvglinternet, X = draws(sim)@draws[["r1.1"]][["xtest_lasso"]],
        type = "response",
        # lambdaType = ifelse(simulator::model(sim)@params$lambda.type=="lambda.min", "lambdaHat","lambdaHat1Std"),
        lambdaType = "lambdaHat")

gbm1 <- gbm::gbm.fit(x = draws(sim)@draws[["r1.1"]][["xtrain_lasso"]],
                     y = draws(sim)@draws[["r1.1"]][["ytrain"]],
                     distribution = "gaussian",
                     n.trees = 1000,
                     verbose = TRUE)

predict(gbm1, )

gbm1
summary(gbm1)
best.iter <- gbm.perf(gbm1,method="cv")
print(best.iter)

gbm::predict.gbm(gbmfit, newdata = draws(sim)@draws[["r1.1"]][["xtrain"]],
                 n.trees = 100)
topgbm <- summary(gbm1, n.trees = 100)
rownames(topgbm[which(topgbm$rel.inf>0),,drop=FALSE])

draws(sim)@draws[["r1.1"]][["xtrain_lasso"]]
simulator::model(sim)@params




sim <- sim %>% evaluate(list(mse, cvmse, tpr, fpr, correct_sparsity,nactive))
e <- evals(sim)
res <- as(e, "data.frame")
# saveRDS(res, file = "sail/sail_lambda_branch/mcgillsims/sim_results/res7.rds")

saveRDS(object = res,
        file = tempfile(pattern = sprintf("res_lassoBT_glinternet"),
                        tmpdir = "/mnt/GREENWOOD_SCRATCH/sahir.bhatnagar/sail_simulations_other",
                        fileext = ".rds")
)



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
