######################################
# R Source code file for simulation of sail paper
# Notes:
# This is the main simulator file
# I ran this on tmux by copy pasting into an rsession on hydra.
# parallel doesnt work in rstudio
# use this code to save the results. turns out you cant run simulator in parallel "by hand"
# you need to load the simulator object prior to adding new simulations to the object.
#
# Since the simulations take a long time to run, I save the results to a .rds file
# and then create figures based on that
# Author: Sahir Bhatnagar
# Created: 2016
# Updated: September 7, 2021
#####################################


# load packages -----------------------------------------------------------

pacman::p_load(simulator) # this file was created under simulator version 0.2.1
pacman::p_load(splines)
pacman::p_load(magrittr)
pacman::p_load(methods)
pacman::p_load(glmnet)
pacman::p_load(LassoBacktracking)
pacman::p_load_gh('sahirbhatnagar/glinternet')
pacman::p_load(gbm)
# remotes::install_github('asadharis/HierBasis')
library("HierBasis")
# devtools::load_all()
pacman::p_load(SAM)
pacman::p_load(gamsel)
pacman::p_load(cowplot)
pacman::p_load(dplyr)
pacman::p_load(tidyr)
pacman::p_load(data.table)
pacman::p_load(ggplot2)
pacman::p_load(latex2exp)
pacman::p_load(lemon)
pacman::p_load(here)
# pacman::p_load(sail)
pacman::p_load(truncnorm)
# remotes::install_local()
library(sail)
# remotes::install_github("sahirbhatnagar/sail@csda-review")

# source helper functions -------------------------------------------------

source(here::here("my_sims/model_functions.R"))
source(here::here("my_sims/method_functions.R"))
source(here::here("my_sims/eval_functions.R"))

# run simulation in parallel on tmux on NEWTON --------------------------------------
# this is the simulation study being used in the CSDA paper now (Aug14_2021 results)
sim <- new_simulation(name = "aug_14_2021",
                      label = "aug_14_2021",
                      dir = ".") %>%
  generate_model(make_gendata_Paper_data_split, seed = 1234,
                 n = 400, p = 1000, corr = 0, betaE = 2, SNR = 2, lambda.type = "lambda.min",
                 parameterIndex = list(1,2,3,4,5),
                 vary_along = "parameterIndex") %>%
  simulate_from_model(nsim = 6, index = 1:35)

simulator::save_simulation(sim)
sim

sim <- sim %>% run_method(list(sailsplit, sailsplitlinear,
                               sailsplitweak,
                               lassosplitadaptive, sailsplitadaptive,
                               lassosplit, lassoBTsplit, GLinternetsplit, Hiersplit, SPAMsplit, gamselsplit),
                          parallel = list(socket_names = 35,
                                          libraries = c("LassoBacktracking", "glinternet","glmnet","splines",
                                                        "magrittr","sail","gamsel","SAM","HierBasis","simulator", "parallel")))

sim <- sim %>% evaluate(list(msevalid, tpr, fpr, nactive, r2))
simulator::save_simulation(sim)
sim


# save evals to data.frame --------------------------------------------------

# aug14_2021 results are the revised simulation results after fixing intercept issue
sim <- load_simulation("aug_14_2021")
df <- as.data.frame(evals(sim))
saveRDS(df, file = "my_sims/simulation_results/aug_14_2021_results.rds")


# plot MSE ----------------------------------------------------------------
# This is just to quickly check results. All plots for paper are in the manuscript folder
df <- readRDS("my_sims/simulation_results/aug_14_2021_results.rds")
df <- df %>% separate(Model, into = c("simnames","betaE","corr","lambda.type","n","p","parameterIndex","SNR_2"),
                      sep = "/")
DT <- as.data.table(df)

DT[parameterIndex=="parameterIndex_1", table(Method)]
DT[, table(parameterIndex)]

cbbPalette <- c("#8720B6","#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
gg_sy <- theme(legend.position = "bottom", axis.text = element_text(size = 20),
               axis.title = element_text(size = 20), legend.text = element_text(size = 20),
               legend.title = element_text(size = 20),plot.title = element_text(size = 20) )

appender <- function(string) TeX(paste(string))

DT[, scenario:= as.numeric(as.character(stringr::str_extract_all(parameterIndex, "\\d", simplify = T)))]
DT$scenario %>% table
DT[, scen:=ifelse(scenario==1,"Strong Hierarchy",ifelse(scenario==2, "Weak Hierarchy", ifelse(scenario==3,"Interactions Only",ifelse(scenario==4, "Strong Hierarchy (Linear)", ifelse(scenario==5, "Main Effects Only", "Linear v2")))))]
DT$scen %>% table
DT[, scen:=factor(scen, levels = c("Strong Hierarchy", "Weak Hierarchy","Interactions Only","Strong Hierarchy (Linear)","Main Effects Only"))]
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
    scale_fill_manual(values=c(cbbPalette, "red","pink","darkblue","black"), guide=guide_legend(ncol=2)) +
    ggplot2::labs(y = "Test Set MSE", title = "") + xlab("") + panel_border()+background_grid()+
    theme(legend.position = "right", legend.text=element_text(size=18)) )

reposition_legend(p1_mse, 'center', panel='panel-2-3')

