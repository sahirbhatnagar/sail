#' This code separates the IQ 4yrs data and related PRS scores for the different
#' individuals into a train, test and validate set as a tool to obtain the most
#' appropriate SAIL model.This will be used as a "real" data example for Sahir's
#' manuscript to illustrate the uses of the SAIL algorithm.

# load packages ---------------------------------------------------------

pacman::p_load(simulator) # this file was created under simulator version 0.2.1
pacman::p_load(data.table)
pacman::p_load(magrittr)
pacman::p_load(tidyverse)
pacman::p_load(doParallel)
pacman::p_load(splines)
pacman::p_load(foreach)
pacman::p_load(doMC)
pacman::p_load(glmnet)
pacman::p_load_current_gh('sahirbhatnagar/glinternet')
pacman::p_load_gh('sahirbhatnagar/sail')
pacman::p_load(LassoBacktracking)
pacman::p_load(SAM)
pacman::p_load(gamsel)
pacman::p_load_gh('asadharis/HierBasis')
pacman::p_load(psych) # for error.crosses plot
pacman::p_load(ggplot2)
pacman::p_load(mice)


# load source code helper files -------------------------------------------

# Uneven split of train and validate datasets   (train=110,valid=40,test=20)
source("rda/PRS/05_model_functions_2.R")
source("rda/PRS/05_method_functions.R")
source("rda/PRS/05_eval_functions.R")
source("rda/PRS/05_plot_functions_rda.R")


# Load data and format appropriately ----------------------------------------

gen3pc <- read.table("rda/PRS/Gen_3PC_scores.txt", header = TRUE)
gen3pc <- cbind(gen3pc,rownames(gen3pc))
colnames(gen3pc)[4] <- "SentrixID"
iq_md <- read.table("rda/PRS/IQ_and_mental_development_variables_for_Sahir_with_study_ID.txt",
                    header=TRUE)
snp_prs_na <- read.table("rda/PRS/NFP_170614_INFO08_nodup_hard09_noambi_GWAS_EduYears_Pooled_beta_withaf_5000pruned_noambi_16Jan2018.score",
                         header = TRUE, sep = ",")

# Merge the iq_md (SentrixID), snp_prs_na (Label.1) and gen3pc (SentrixID) datasets together
m1 <- merge(iq_md,snp_prs_na, by.x = "SentrixID",by.y = "Label.1")
IQdat <- merge(m1,gen3pc, by.x = "SentrixID", by.y = "SentrixID")

# Keep relavant columns
IQ4y <- subset(IQdat,
               select=c("IQ_4yrs", #Y variable
                        "Tx_group_bin", #E variable
                        "PRS_0.0001","PRS_0.001","PRS_0.01",
                        "PRS_0.05","PRS_0.1","PRS_0.2","PRS_0.3",
                        "PRS_0.4","PRS_0.5"))   #X variables

pt <- mice::mice(IQ4y, m = 1, seed = 1234)
IQ4y <- mice::complete(pt)
IQ4y <- as.matrix(IQ4y)
rownames(IQ4y) <- IQdat[,"SentrixID"]
IQ4y <- IQ4y[!is.na(IQ4y[,"IQ_4yrs"]),]

# cleanup
rm(list = c("IQdat","m1","snp_prs_na","iq_md","pt","gen3pc"))


# Comparison of methods over 200 random splits of the data ----------------
# Glinternet did not converge
# Create the 210 simulations (different splits of train/validate/test)
sim <- new_simulation(name = "PRS_IQ_4yrs_cond2_split2_4",
                      label = "PRS_IQ_4yrs_cond2_split2_4",
                      dir = "my_sims/") %>%
  generate_model(make_data_split,
                 DT = IQ4y,
                 X = IQ4y,
                 E = IQ4y[,"Tx_group_bin"],
                 Y = IQ4y[,"IQ_4yrs"],
                 seed = 1234,
                 phenoVariable = "IQ_4yrs",
                 exposure = "Tx_group_bin",
                 n_train_test = 150) %>%
  simulate_from_model(nsim = 6, index = 1:35)

# Now, we run the 4 methods  save the results
sim <- run_method(sim,
                  list(lassosplit,
                       sailsplit,
                       sailsplitadaptiveweak,
                       sailsplitweak,
                       lassoBTsplit,
                       Hiersplit),
                  parallel = list(socket_names = 40,# refers to the number of cores on which to split the simulation
                                  libraries = c("LassoBacktracking",
                                                "glmnet","splines",
                                                "magrittr","sail","gamsel",
                                                "SAM","HierBasis","simulator",
                                                "parallel"))) %>%
  evaluate(list(msevalid, nactive, r2))

save_simulation(sim)
simulation_results_PRS <- as.data.frame(evals(sim))
simulation_results_PRS <- as.data.table(simulation_results_PRS)

# this isnt used anymore. use sail_PRS_results.RData instead
# saveRDS(simulation_results_PRS, file = "my_sims/rda_results/rda_IQ_4yrs_cond2_s2_4.rds")
# simulation_results_PRS <- readRDS(file = "my_sims/rda_results/rda_IQ_4yrs_cond2_s2_4.rds")

index <- paste0("r",rep(1:35,each=6),".",1:6)

# the 4 is the index for sail weak
sail_coef_list_PRS <- lapply(index, function(i) {
  simulator::output(sim)[[4]]@out[[i]][["beta"]]
})

save(simulation_results_PRS, sail_coef_list_PRS, file = "rda/PRS/sail_PRS_results.RData")
# ls()

# sim <- load_simulation(name = "PRS_IQ_4yrs_cond2_split2_4", dir = "my_sims/")
# simulator::output(sim)[[4]]@out[[3]][["beta"]]
# simulator::output(sim)[[2]]@out[[]][["msevalid"]]
# sim %>% subset_simulation(methods = c("Adaptivesailweak", "lasso", "lassoBT", "sail",
#                                      "sailweak")) %>% plot_eval("mse")
# sim %>% subset_simulation(methods = c("Adaptivesailweak", "lasso", "lassoBT", "sail",
#                                       "sailweak")) %>% plot_eval("nactive")
