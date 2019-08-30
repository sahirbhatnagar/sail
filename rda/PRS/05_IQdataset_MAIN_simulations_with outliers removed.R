#This code separates the IQ 4yrs data and related PRS scores for the different individuals into a train, test and validate set
#as a tool to obtain the most appropriate SAIL model.This will be used as a "real" data example for Sahir's manuscript to illustrate
#the uses of the SAIL algorithm. We also remove any outliers in order to improve legibility of the graphs.
#This is the code used to generate the figures for the manuscript




# load packages ---------------------------------------------------------

# rm(list=ls())
pacman::p_load(simulator) # this file was created under simulator version 0.2.1
# devtools::load_all()
pacman::p_load(data.table)
pacman::p_load(magrittr)
# pacman::p_load(genefilter)
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
# pacman::p_load(sortinghat)
# pacman::p_load(sail)
# pacman::p_load(ggplot2)







# load source code helper files -------------------------------------------

#Splits data evenly between train and validate (train=60,valid=60,test=50)
  #source("/mnt/GREENWOOD_JBOD1/GREENWOOD_BACKUP/home/amanda.lovato/DNA methylation/Sahir/Script/05_model_functions.R")
#Uneven split of train and validate datasets   (train=110,valid=40,test=20)
  # source("/mnt/GREENWOOD_JBOD1/GREENWOOD_BACKUP/home/amanda.lovato/DNA methylation/Sahir/Script/05_model_functions_2.R")
  # source("/mnt/GREENWOOD_JBOD1/GREENWOOD_BACKUP/home/amanda.lovato/DNA methylation/Sahir/Script/05_method_functions.R")
  # source("/mnt/GREENWOOD_JBOD1/GREENWOOD_BACKUP/home/amanda.lovato/DNA methylation/Sahir/Script/05_eval_functions.R")
  # source("/mnt/GREENWOOD_JBOD1/GREENWOOD_BACKUP/home/amanda.lovato/DNA methylation/Sahir/Script/05_plot_functions_rda.R")

source("~/git_repositories/sail/rda/PRS/05_model_functions_2.R")
source("~/git_repositories/sail/rda/PRS/05_method_functions.R")
source("~/git_repositories/sail/rda/PRS/05_eval_functions.R")
source("~/git_repositories/sail/rda/PRS/05_plot_functions_rda.R")


#######################################################################################################################################################
######Load data and format appropriately


gen3pc=read.table("~/git_repositories/sail/rda/PRS/Gen_3PC_scores.txt", header=TRUE)
#Missing one column of info:  ### #eigvals:	8.654	1.217	1.198     ###   (might be important)
gen3pc=cbind(gen3pc,rownames(gen3pc))
colnames(gen3pc)[4]="SentrixID"

iq_md=read.table("~/git_repositories/sail/rda/PRS/IQ_and_mental_development_variables_for_Sahir_with_study_ID.txt", header=TRUE)

snp_prs_na=read.table("~/git_repositories/sail/rda/PRS/NFP_170614_INFO08_nodup_hard09_noambi_GWAS_EduYears_Pooled_beta_withaf_5000pruned_noambi_16Jan2018.score", header=TRUE,sep=",")



######Merge the iq_md (SentrixID), snp_prs_na (Label.1) and gen3pc (SentrixID) datasets together

m1=merge(iq_md,snp_prs_na,by.x="SentrixID",by.y="Label.1")
IQdat=merge(m1,gen3pc,by.x="SentrixID",by.y="SentrixID")

###Keep relavant columns

IQ4y=
  subset(IQdat, select=c("IQ_4yrs",                                                                                           #Y variable
                         "Tx_group_bin",                                                                                      #E variable
                         "PRS_0.0001","PRS_0.001","PRS_0.01","PRS_0.05","PRS_0.1","PRS_0.2","PRS_0.3","PRS_0.4","PRS_0.5"))   #X variables

IQ4y=as.matrix(IQ4y)
rownames(IQ4y)=IQdat[,"SentrixID"]
IQ4y=IQ4y[!is.na(IQ4y[,"IQ_4yrs"]),]  #remove rows with missing IQ values




# these are used for plotting. we use the entire data set for the plots --------------------
DT <- IQ4y
X <- IQ4y %>% as.matrix()
Xnorm <- sail:::standardize(X, center = TRUE, normalize = TRUE)$x
e <- IQ4y[, "Tx_group_bin"]
Y <- DT[,"IQ_4yrs"] %>% as.numeric




#######################################################################################################################################################
######Generate the 200 simulations and test them against sail, lasso, lassobt, hierbasis and glinternetsplit



#-----------------
#4 methods were run here (sail, lasso, lassobt, hierbasis)
#Glinternetsplit did not converge (you can add glinternetsplit as one of the factors in run_method to see for yourself)
#When plots of mse and number of avtive variables were fun for these 4 methods, we noticed that there were some outliers.
#These outliers made the plots more difficult to interpret. As such, the next step was to identify which of the 210 simulation
#resulted in these results and remove them from the analysis

#sim <- new_simulation(name = "IQ_4yrs_cond2_split2_4",
#                      label = "IQ_4yrs_cond2_split2_4",
#                      dir = ".") %>%
#  generate_model(make_data_split, seed = 12254,
#                 phenoVariable = "IQ_4yrs",
#                 exposure = "Tx_group_bin",n_train_test = 150) %>%
#  simulate_from_model(nsim = 6, index = 1:35)  %>%
#  run_method(list(lassosplit, sailsplit, lassoBTsplit, Hiersplit),
#             parallel = list(socket_names = 35,#refers to the number of cores on which to split the simulation
#                             libraries = c("LassoBacktracking", "glinternet","glmnet","splines",
#                                           "magrittr","sail","gamsel","SAM","HierBasis","simulator", "parallel")))
#
#sim <- sim %>% evaluate(list(msevalid, nactive, r2))
#
#save(sim,file="/mnt/GREENWOOD_JBOD1/GREENWOOD_BACKUP/home/amanda.lovato/files/IQ_4yrs_cond2_s2_4.Rdata")



##########Here we run the 4 simulations (but wait before running the 4 methods)

#Create the 210 simulations (different splits of train/validate/test)
sim <- new_simulation(name = "IQ_4yrs_cond2_split2_4",
                      label = "IQ_4yrs_cond2_split2_4",
                      dir = "my_sims/") %>%
  generate_model(make_data_split, seed = 12254,
                 phenoVariable = "IQ_4yrs",
                 exposure = "Tx_group_bin",n_train_test = 150) %>%
  simulate_from_model(nsim = 6, index = 1:35)





#To see which simulations were leading to outlier values for MSE and number of active variables,
#we created the following code that explicitly tests that:

#"/mnt/GREENWOOD_JBOD1/GREENWOOD_BACKUP/home/amanda.lovato/DNA methylation/Sahir/Script/Finding the outliers.R"

#This code found r8.5 to be outlier for mse HierBasis and r25.6 and r30.2 to be an outlier in terms of number of variables
#We will remove these draws and re-run sail, lasso,lassobt and hierbasis on the subsetted list. Unfortunately, we were unable to
#access the results for each draw separately. Thus, we removed all of r8,r25 and r30 (each of which has 6 draws).
#This removes a total of 18 draws from the 210 leaving 192. #This is done right below:


# load("/mnt/GREENWOOD_JBOD1/GREENWOOD_BACKUP/home/amanda.lovato/files/sim-IQ_4yrs_cond2_split2_4.Rdata")
# sim@draws_refs[[1]][[30]]=NULL
# sim@draws_refs[[1]][[25]]=NULL
# sim@draws_refs[[1]][[8]]=NULL
# save(sim,file="/mnt/GREENWOOD_JBOD1/GREENWOOD_BACKUP/home/amanda.lovato/files/sim-IQ_4yrs_cond2_split2_4.Rdata")
# load("/mnt/GREENWOOD_JBOD1/GREENWOOD_BACKUP/home/amanda.lovato/files/sim-IQ_4yrs_cond2_split2_4.Rdata")


#Now, we run the 4 methods on the 192 remaining simulations and save the results:
sim <- run_method(sim,list(lassosplit, sailsplit,sailsplitadaptive, sailsplitadaptiveweak,sailsplitweak, lassoBTsplit, Hiersplit),
                  parallel = list(socket_names = 40,#refers to the number of cores on which to split the simulation
                                  libraries = c("LassoBacktracking", "glmnet","splines",
                                                "magrittr","sail","gamsel","SAM","HierBasis","simulator", "parallel"))) %>%
  evaluate(list(msevalid, nactive, r2))

save_simulation(sim)
as.data.frame(evals(sim))
ls()

sim <- load_simulation(name = "IQ_4yrs_cond2_split2_4",
                       dir = "my_sims/")

plot_eval(sim, "mse")

sim %>% subset_simulation(methods = c("lasso", "sail","sailweak",
                                      "Adaptivesail", "Adaptivesailweak",
                                      "lassoBT")) %>% plot_eval("nactive")

sim %>% subset_simulation(methods = c("lasso", "sail","sailweak",
                                      "Adaptivesail", "Adaptivesailweak",
                                      "lassoBT")) %>% plot_eval("mse")
# save(sim,file="/mnt/GREENWOOD_JBOD1/GREENWOOD_BACKUP/home/amanda.lovato/files/IQ_4yrs_cond2_s2_4.Rdata")




#######################################################################################################################################################
######Plot the results


# load simulator object for plots -----------------------------------------

# load("/mnt/GREENWOOD_JBOD1/GREENWOOD_BACKUP/home/amanda.lovato/files/IQ_4yrs_cond2_s2_4.Rdata")

df <- as.data.frame(evals(sim))
# saveRDS(df, file = "~/git_repositories/sail/my_sims/rda_results/rda_IQ_4yrs_cond2_s2_4.rds")




# plot results ------------------------------------------------------------

# sim %>% plot_eval(metric_name = "nactive")
# sim %>% plot_eval(metric_name = "mse")
# sim %>% plot_eval(metric_name = "r2")


## ---- load-results
#df <- readRDS("/mnt/GREENWOOD_JBOD1/GREENWOOD_BACKUP/home/amanda.lovato/files/IQ_4yrs_cond2_s2_4.Rdata")
df <- readRDS(file = "~/git_repositories/sail/my_sims/rda_results/rda_IQ_4yrs_cond2_s2_4.rds")
DT_res <- df %>% as.data.table()


## ---- get-best-sail-model ---- (the one that results in the lowest MSE)

draw_ind <- DT_res[Method=="sail"][which.min(mse)]$Draw
draw_ind <- sample(DT_res[Method=="sail"]$Draw, size = 1)
dat <- draws(sim)@draws[[draw_ind]]
saveRDS(dat, file = "~/git_repositories/sail/my_sims/rda_results/rda_IQ_4yrs_cond2_s2_4_data.rds")
dat <- readRDS("~/git_repositories/sail/my_sims/rda_results/rda_IQ_4yrs_cond2_s2_4_data.rds")
# devtools::load_all()
fit <- sail(x = dat$xtrain, y = dat$ytrain, e = dat$etrain,
            basis = function(i) splines::bs(i, degree = 3),
            expand = TRUE,
            center.x = T,
            center.e = T,
            #group = dat$group,
            alpha = 0.1,
            #maxit = 250,
            strong = FALSE,
            verbose = 2
)
plot(fit)
fit
ytest_hat <- predict(fit, newx = dat$xtest, newe = dat$etest)
msetest <- colMeans((dat$ytest - ytest_hat)^2)
lambda.min.index <- as.numeric(which.min(msetest))
lambda.min <- fit$lambda[which.min(msetest)]
yvalid_hat <- predict(fit, newx = dat$xvalid, newe = dat$evalid, s = lambda.min)
msevalid <- mean((dat$yvalid - drop(yvalid_hat))^2)
nzcoef <- predict(fit, s = lambda.min, type = "nonzero")

design <- do.call(rbind, list(dat$xtrain, dat$xtest, dat$xvalid))
#design2=sail:::design_sail(design,e,expand=TRUE,basis = function(i) splines::bs(i, degree = 3),center.e=TRUE,center.x=TRUE,nvars=ncol(design),vnames=colnames(design))$design




#Get the fit for the entire dataset (no special selection of the best model- just simple application of sail to the dataset)
fit <- sail(x = X[,-c(1,2)], y = Y, e = IQ4y[, "Tx_group_bin"],
              basis = function(i) splines::bs(i, degree = 3),
              expand = TRUE,
              center.x = T,
              center.e = T,
              #group = dat$group,
              alpha = 0.1,
              #maxit = 250,
              strong = FALSE,
              verbose = 1
)

library(doMC)
registerDoMC(cores = 10)
fit <- cv.sail(x = X[,-c(1,2)], y = Y, e = IQ4y[, "Tx_group_bin"],
            basis = function(i) splines::bs(i, degree = 3),
            expand = TRUE,
            center.x = T,
            center.e = T,
            #group = dat$group,
            alpha = 0.1,
            #maxit = 250,
            strong = FALSE,
            verbose = 1,
            parallel = TRUE
)


plot(fit)



## ----error-crosses plots----

affect.mat2 <- describeBy(DT_res[, c("mse","nactive")], DT_res$Method, mat = TRUE)
affect.mat2 <- affect.mat2[which(affect.mat2$group1 %in% c("sail","lasso","lassoBT","HierBasis")),]
par(family="serif")
error.crosses(affect.mat2[grep("nactive", rownames(affect.mat2)),],
              affect.mat2[grep("mse", rownames(affect.mat2)),],
              labels=unique(affect.mat2$group1),
              xlab="Number of Active Variables",
              main = "IQ Data: Means (+/- 1 SE) from 200 Train/Validate/Test Splits",
              sd = FALSE,
              cex.lab = 1.4,
              cex.axis = 1.4,
              cex.main = 1.5,
              xlim = c(0, 10),
              # ylim=c(100,250),
              #ylim=c(-1000,1000),
              ylab="Test Set MSE",
              colors = sail:::cbbPalette[c(3,4,7,2)],
              pch=16, cex=2)





#########################

#Plot interactions separately for the dichotomous Environment variable






#To extrapolate the entire dataset (not just the training set on which we fitted)
#Uncomment the following lines

#object=fit_t   #fit on the entire dataset
#x=X[,-c(1,2)]

#To extrapolate the  just the training set on which we fitted
#Uncomment the following lines
object=fit      #fit on the training dataset
x=dat$xtrain


#Define important variables
i=1 #(since we are interested in PRS_0.0001 which is the first of the covariates, i=1)
cv=c("IQ4y")
xv=colnames(get(cv[1]))[-c(1,2)]  #Remove the first 2 columns that contain the response and environment variable
xvar = xv[i]
s = lambda.min
xlab=xv[i]
ylab=paste0("f(",xv[i],")")
color=c("orange","blue")

#Extract relevant coefficients
ind <- object$group == which(object$vnames == xvar)
allCoefs <- coef(object, s = s)
a0 <- allCoefs[1, ]
betas <- as.matrix(allCoefs[object$main.effect.names[ind], , drop = FALSE])
alphas <- as.matrix(allCoefs[object$interaction.names[ind], , drop = FALSE])
betaE <- as.matrix(allCoefs["E", , drop = FALSE])



#####Extract the design matrix from the fit object (from the training dataset - 110 individuals, from the entire dataset-170 individuals)
design.mat.main <- object$design[, object$main.effect.names[ind], drop = FALSE]
design.mat.int <- object$design[, object$interaction.names[ind], drop = FALSE]



ind2 <- which(object$vnames %in% xvar)

originalE <- object$design[, "E", drop = FALSE]
originalX <- x[,ind2,drop=FALSE]


#Obtain the associated IQ value for the corresponding PRS score and environment status
f.hat= drop(originalE * as.vector(betaE) +
              design.mat.main %*% betas + design.mat.int %*% alphas)

ylims <- range(f.hat)


#Separate data according to the Environement variable(control or intervention)

intervention = max(drop(unique(originalE)))
control = min(drop(unique(originalE)))


cont_index <- which(originalE==control)
int_index <- which(originalE==intervention)


cont_pred <- f.hat[cont_index]
int_pred <- f.hat[int_index]


#Plot the effects

min.length.top <- range(f.hat)[1] ; max.length.top <- range(f.hat)[2]
# par(mai=c(1,1,1,0.2))
plot(originalX, f.hat,
     pch = 19,
     ylab = ylab,
     xlab = xlab,
     bty="n",
     xaxt="n",
     type = "n",
     # cex.lab = 2,
     # cex.axis = 2,
     # cex = 2,
     main = "Effect of Intervention and PRS
     on IQ at 4 years of age"
     # cex.main = 2.5,
     # ylim = c(min.length.top-3, max.length.top+3),
     # ylim = ylims
)


# axis(1, labels = T, cex.axis = 2)
axis(1, labels = T)

lines(originalX[cont_index][order(originalX[cont_index])], cont_pred[order(originalX[cont_index])], col = color[1], lwd = 3)
lines(originalX[int_index][order(originalX[int_index])], int_pred[order(originalX[int_index])], col = color[2], lwd = 3)
legend("bottomright",c("Intervention","Control"),col=c(color[2],color[1]),lwd=c(3,3))











#####################$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#####################$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#####################$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#########Extra... looking at main and interaction effects for all variables

############
####Plot Main effects (fit=fit with training dataset giving lowest MSE   ; fit_t=fit with the entire dataset)

######At lambda.min
par(mfrow=c(2,2))
cv=c("IQ4y")
xv=colnames(get(cv[1]))[-c(1,2)]

for(j in 1:length(cv)){
  pdf(paste0("/mnt/GREENWOOD_JBOD1/GREENWOOD_BACKUP/home/amanda.lovato/DNA methylation/Sahir/DATA/IQ and mental development/",cv[j],"_main_c2_2.pdf"))
  for(i in 1:length(xv)){


    plotMain(fit, x = dat$xtrain, e = dat$etrain_lasso,
             xvar = xv[i], s = lambda.min,ylab = paste0("f(",xv[i],")"), xlab = xv[i],legend.position = "bottomleft", main="Train plot")


    plotMain(fit, x = dat$xtest, e = dat$etest_lasso,
             xvar = xv[i], s = lambda.min,ylab = paste0("f(",xv[i],")"), xlab = xv[i],legend.position = "bottomleft", main="Test plot")


    plotMain(fit, x = dat$xvalid, e = dat$evalid,
             xvar = xv[i], s = lambda.min,ylab = paste0("f(",xv[i],")"), xlab = xv[i],legend.position = "bottomleft", main="Validate plot")


    plotMain(fit_t, x = X[,-c(1,2)], e = e,
             xvar = xv[i], s = lambda.min,ylab = paste0("f(",xv[i],")"), xlab = xv[i],legend.position = "bottomleft", main="Full dataset plot")

  }
  dev.off()
}





############
####Plot Inter effects (fit=fit with training dataset giving lowest MSE   ; fit_t=fit with the entire dataset)

######At lambda.min (based on the best sail model)

cv=c("IQ4y")
xv=colnames(get(cv[1]))[-c(1,2)]

for(j in 1:length(cv)){
  pdf(paste0("/mnt/GREENWOOD_JBOD1/GREENWOOD_BACKUP/home/amanda.lovato/DNA methylation/Sahir/DATA/IQ and mental development/",cv[j],"_inter_c2_2.pdf"))
  for(i in 1:length(xv)){


    try( plotInter(object = fit, x = dat$xtrain, e = dat$etrain_lasso,
                   xvar = xv[i], s = lambda.min,ylab = paste0("f(",xv[i],")"), xlab = xv[i],legend.position = "bottomleft", main="Train plot") )


    try( plotInter(object = fit, x = dat$xtest, e = dat$etest_lasso,
                   xvar = xv[i], s = lambda.min,ylab = paste0("f(",xv[i],")"), xlab = xv[i],legend.position = "bottomleft", main="Test plot") )


    try( plotInter(object = fit,  x = dat$xvalid, e = dat$evalid,
                   xvar = xv[i], s = lambda.min,ylab = paste0("f(",xv[i],")"), xlab = xv[i],legend.position = "bottomleft", main="Validate plot") )


    try( plotInter(fit_t, x = X[,-c(1,2)], xvar=xv[i], e = e,
                   s = lambda.min, ylab = paste0("f(",xv[i],")"), xlab = xv[i], legend.position = "bottomleft", main="Full dataset plot") )



  }
  dev.off()
}












