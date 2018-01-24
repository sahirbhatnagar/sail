# This is the main simulator file
rm(list=ls())
pacman::p_load(simulator) # this file was created under simulator version 0.2.1
source("/home/sahir/git_repositories/sail/my_sims/model_functions.R")
source("/home/sahir/git_repositories/sail/my_sims/method_functions.R")
source("/home/sahir/git_repositories/sail/my_sims/eval_functions.R")

## @knitr init

name_of_simulation <- "sail_udem_talk"

## @knitr main

sim <- new_simulation(name = name_of_simulation,
                      label = "Sail easy p25") %>%
  generate_model(make_easy_sail_model, seed = 123,
                 n = 200, p = 25, df = 5, SNR = 2, betaE = 1) %>%
  simulate_from_model(nsim = 25) %>%
  run_method(list(lasso, my_sail))


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
