#################################
# Script to analyze ADNI data
# I used the data given by Lai. The plots generated in this script were
# used for the UdeM Job Talk. I used  MMSCORE ~ AMY_beta * APOE, and also tried
# MMSCORE ~ AMY_beta * EDUCAT
#################################3

#' 1) Visit_BL is for the visit. Subjects had several visits: baseline, 48
#' months, 36 months, 60 months…However, I believe that you have the combined
#' information (A-beta and genotypes) only for the baseline. 2) Diag_3bl is the
#' diagnosis: 1=control, 2=MCI (Mild Cognitive Impairment) and 3=Alzeimer
#' Disease 3) APOE_bin is an indicator if the individual has the « bad » allele
#' for the APOE gene…APOE gene is THE gene for Alzheimer and usually one wants
#' to control for that in the analysis. 4) ABETA142 is an overall measure of
#' A-beta for the whole brain.

message("loading packages")

rm(list=ls())
# devtools::document()
devtools::load_all()
# source("R/ksfilter.R")
pacman::p_load(data.table)
pacman::p_load(magrittr)
pacman::p_load(genefilter)
pacman::p_load(tidyverse)
pacman::p_load(doParallel)
# pacman::p_load(mice)


message("loading data")

# amy_pheno <- xlsx::read.xlsx("~/Downloads/DrCelia_data.xlsx", sheetIndex = 1)
amy_mat <- read.csv("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_git/sail/data/csf_amyloid_final.csv", stringsAsFactors = FALSE)
covr <- read.csv("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_git/sail/data/covariates.csv", stringsAsFactors = FALSE, sep = ";")
# surv <- read.csv("~/git_repositories/sail/data/adni_new/fdg_info.csv", stringsAsFactors = FALSE, sep = ",")

# sum(as.character(amy_pheno$PTID) %in% amy_mat$PTID)
# sum(as.character(covr$IID) %in% amy_mat$PTID)
# as.character(covr$IID) %>% unique() %>% length()

DT <- dplyr::inner_join(amy_mat, covr, by = c("PTID" = "IID")) %>%
  select(-AV45_path_bl)

X <- DT %>% select(starts_with("X"), Age_bl, diag_3bl.x) %>% as.matrix()
dimnames(X)[[1]] <- DT$PTID
# X <- X[,80:97]

E <- DT %>% pull(APOE_bin) %>% as.numeric
# E <- DT %>% pull(EDUCAT) %>% as.numeric
# E <- DT %>% pull(diag_3bl.x) %>% as.numeric

Y <- DT %>% pull(MMSCORE_bl) %>% as.numeric


# pacman::p_load(pheatmap)
# pacman::p_load(viridis)
# ab_mat <- DT %>% select(starts_with("X")) %>% as.matrix()
# dimnames(ab_mat)[[1]] <- DT$PTID
# annotation_row <- DT %>% select(APOE_bin) %>%
#   dplyr::rename(APOE = APOE_bin) %>%
#   mutate(APOE = factor(APOE, levels = 0:1, labels = c("APOE = 0", "APOE = 1")))
# rownames(annotation_row) <- DT$PTID
#
# pdf("results/figures/amy_mat_heatmap.pdf", height = 7, width = 10)
# pheatmap::pheatmap(ab_mat,
#                    show_rownames = F, show_colnames = F,
#                    color = viridis(100),
#                    # annotation_col = annotation_col,
#                    # annotation_row = annotation_row,
#                    # annotation_colors = ann_colors,
#                    annotation_names_row = FALSE,
#                    annotation_names_col = TRUE,
#                    drop_levels = FALSE,
#                    annotation_legend = TRUE)
# dev.off()
#
# str(annotation_row)
# summary(lm(Y ~ X*E))
# fit <- sail(x = X, y = Y, e = E, df = 3, maxit = 500,
#             nlambda.gamma = 5, nlambda.beta = 5,
#             nlambda = 25,
#             # group.penalty = "SCAD",
#             lambda.factor = 0.0001,
#             thresh = 1e-4, center=T, normalize=F, verbose = T)
# coef(fit)

message("fitting cv.sail")

pacman::p_load(cowplot)
dp <- data.frame(X,Y, E) %>% mutate(Status = factor(diag_3bl.x, levels = 1:3,
                                                    labels = c("Control", "Mild Cognitive Impairment","Alzeimer Disease")),
                                    APOE = factor(E, levels = 0:1, labels = c("APOE e4 = 0", "APOE e4 = 1")))

ggplot(dp, aes(x = Status, y = Y, fill = Status)) + geom_boxplot() + theme(legend.position = "none") +
  ylab("Mini-Mental State Examination")+scale_fill_manual(values = trop[1:3]) + xlab("")

trop <- RSkittleBrewer::RSkittleBrewer("trop")
ggplot(dp, aes(x = APOE, y = Y, fill = APOE)) + geom_boxplot() +
  theme(legend.position = "none",axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20), axis.title = element_text(size = 20)) +
  ylab("Mini-Mental State Examination") + scale_fill_manual(values = trop[1:2]) + xlab("")

# ran this using tmux (name of seesion is adni on hydra)
doParallel::registerDoParallel(cores = 5)
getDoParWorkers()
cvfit <- cv.sail(x = X, y = Y, e = E, df = 3, maxit = 500,
                 nlambda.gamma = 10, nlambda.beta = 10,
                 nlambda = 100,
                 # group.penalty = "SCAD",
                 # lambda.beta = exp(seq(log(0.01 * 1000), log(1000), length.out = 25)),
                 # lambda.gamma =rep(1000,25),
                 lambda.factor = 0.0001,
                 nfolds = 5,
                 thresh = 1e-3, center=T, normalize=F, verbose = T)

message("cv.sail DONE")


# plot(cvfit)
# plot(Y,cbind(1,cvfit$sail.fit$design) %*% coef(cvfit, s = "lambda.min"))
# coef(cvfit, s = "lambda.min") %>% dim

saveRDS(cvfit,
        file = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_git/sail/results/cvfit_adni_n_340_p_98_control_for_diag_APOE_E_df_3.rds")
# cvfit <- readRDS(file = "results/cvfit_adni_n_340_p_97_E_diag_df_3.rds")
# cvfit <- readRDS(file = "results/cvfit_adni_n_340_p_97_E_apoe_df_3.rds")
cvfit <- readRDS(file = "results/cvfit_adni_n_340_p_98_control_for_diag_APOE_E_df_3.rds")
plot(cvfit)
# cvfit$sail.fit$design[,"X_E"] %>% table
#
# yhat <- cbind(1,cvfit$sail.fit$design) %*% coef(cvfit, s = "lambda.1se")

yhat <- cbind(1,cvfit$sail.fit$design) %*% coef(cvfit$sail.fit, s = "s65")
which.min(cvfit$cvm)
cvfit$cvm[which.min(cvfit$cvm)]

which(cvfit$cvm<2.46)
coef(cvfit$sail.fit, s = "s65")[nonzeroCoef(coef(cvfit$sail.fit, s = "s65")),,drop=F]
cvfit$lambda.min.name <- "s65"

cvfit$cvm[65]
cvfit$nz.main
cvfit$sail.fit
par(mfrow = c(1,2))
hist(yhat[,1])
hist(Y)

(r2_sail <- cor(yhat[,1],Y)^2)
(rmse <- sqrt(crossprod(as.matrix(as.matrix(Y) - yhat))) / nrow(X) )
#
fmla <- as.formula(paste0("~0+",paste(colnames(X), collapse = "+"), "+ E +",paste(colnames(X),":E", collapse = "+") ))
Xmat <- model.matrix(fmla, data = data.frame(X,E))
#
pacman::p_load(glmnet)
# glfit <- cv.glmnet(x = cvfit$sail.fit$design, y = Y, nfolds = 5)
set.seed(123456)
glfit <- cv.glmnet(x = Xmat, y = Y, nfolds = 5)
dev.off()
plot(glfit)
coef(glfit)
coef(glfit, s = "lambda.min")[nonzeroCoef(coef(glfit, s = "lambda.min")),,drop=F]

# yhat_glmnet <- cbind(1, cvfit$sail.fit$design) %*% coef(glfit, s = "lambda.min")
yhat_glmnet <- cbind(1, Xmat) %*% coef(glfit, s = "lambda.min")
(r2_glmnet <- cor(drop(yhat_glmnet),Y)^2)
(rmse <- sqrt(crossprod(as.matrix(Y) - yhat_glmnet)) / nrow(X))

# cross validated mean MSE
(cv_mse <- data.frame(model = c("sail","lasso"),
                      mean = c(cvfit$df[cvfit$lambda.min.name,,drop=F]$mse,
                               glfit$cvm[which(glfit$lambda.min==glfit$lambda)]),
                      lower = c(cvfit$df[cvfit$lambda.min.name,,drop=F]$lower,
                                glfit$cvlo[which(glfit$lambda.min==glfit$lambda)]),
                      upper = c(cvfit$df[cvfit$lambda.min.name,,drop=F]$upper,
                                glfit$cvup[which(glfit$lambda.min==glfit$lambda)]),
                      r2 = c(r2_sail, r2_glmnet)))

pacman::p_load(latex2exp)

pdf(file="results/figures/apoe_interactions_5foldcv_mse_with_lasso_main_int_design.pdf",width=11,height=8)
par(mai=c(1,1,0.1,0.2))
par(tcl=-0.5, family="serif", omi=c(0.2,0.2,0,0))
plot(2:3, seq(2,4.4, length.out = 2),
     pch = 19,
     ylab = "5-Fold CV error",
     xlab = "",
     # col = color[1],
     bty="n",
     xaxt="n",
     type = "n",
     cex.lab = 2,
     cex.axis = 2,
     cex = 2)
# ylim = c(min.length.top-3, max.length.top+3),
# ylim = ylim)
axis(1, labels = c(TeX("sail"), "lasso"), at = c(2.1,2.3), cex.axis = 2)
# text(2.1, 4,
#      labels = sprintf("%.1f [%.1f, %.1f]", cv_mse[1,"mean"], cv_mse[1,"lower"],cv_mse[1,"upper"]), cex = 1.5)
text(2.1, 4.40,
     labels = TeX("$R^2 = 0.61$"), cex = 1.5)
text(2.3, 4.4,
     labels = TeX("$R^2 = 0.44$"), cex = 1.5)
points(c(2.1,2.3), cv_mse[,"mean"], col = "red", pch = 19, cex = 1.5 )
glmnet:::error.bars(c(2.1,2.3), lower = cv_mse[,"lower"], upper = cv_mse[,"upper"])
dev.off()



group <- rep(seq_len(ncol(X)), each = cvfit$sail.fit$df)
orig_names <- paste(rep(colnames(X), each = cvfit$sail.fit$df), rep(seq_len(cvfit$sail.fit$df), times = ncol(X)), sep = "_")

map_names <- data.frame(orig = c(orig_names, "E", paste0(orig_names,":X_E")),
                        sail = rownames(coef(cvfit, s = "lambda.min"))[-1], stringsAsFactors = FALSE)

active_betas <- coef(cvfit, s = "lambda.min")[nonzero(coef(cvfit, s = "lambda.min")),,drop = F][-1,,drop=F] %>% rownames()
map_names[map_names$sail %in% active_betas,,drop=F]
(active_betas_lasso <- coef(glfit, s = "lambda.min")[nonzeroCoef(coef(glfit, s = "lambda.min"))[-1],,drop=F] %>% rownames)

unique(unlist(stringr::str_extract_all(active_betas_lasso, "^[^_]+(?=_)")))
unique(active_betas_lasso) %>% length()
unlist(stringr::str_extract_all(active_betas, "^[^_]+(?=_)"))
unique(unlist(stringr::str_extract_all(active_betas, "^[^_]+(?=_)"))) %>% length()
length(active_betas)

sail_coef <- as.data.frame(as.matrix(coef(cvfit, s = "lambda.min")[nonzero(coef(cvfit, s = "lambda.min")),,drop = F][-1,,drop=F]))
str(sail_coef)
sail_coef$sail <- rownames(sail_coef)
sail_coef <- left_join(sail_coef, map_names)
str(sail_coef)
sail_coef$orig <- replace(sail_coef$orig, which(sail_coef$orig=="E"),"X_E")
sail_coef$region <- unlist(stringr::str_extract_all(sail_coef$orig, "^[^_]+(?=_)"))

region_names <- read.csv("data/brain_regions.csv", header = F, stringsAsFactors = FALSE)
region_names$V1 <- paste0("X",region_names$V1)

sail_coef <- left_join(sail_coef, region_names, by = c("region" = "V1"))

dput(sail_coef$V2)

sail_coef$mod_names <- c(paste("lateral front-orbital gyrus right",1:3),
                         paste("medial front-orbital gyrus right",1:3),
                         paste("supramarginal gyrus right",1:3),
                         paste("occipital pole left",1:3),
                         paste("Age",1:3),
                         paste("Diagnosis",1:3),
                         "APOE e4 = 1",
                         paste(paste("lateral front-orbital gyrus right",1:3),"x APOE"),
                         paste(paste("medial front-orbital gyrus right",1:3),"x APOE"),
                         paste(paste("supramarginal gyrus right",1:3),"x APOE"),
                         paste(paste("occipital pole left",1:3),"x APOE"),
                         paste(paste("Age",1:3),"x APOE"),
                         paste(paste("Diagnosis",1:3),"x APOE"))

sail_coef[order(sail_coef$mod_names),,drop=F]
dput(sail_coef[order(sail_coef$mod_names),"mod_names",drop=F])

order_sail_names <- c("APOE e4 = 1",
                      "Age 1",
                      "Age 2",
                      "Age 3",
                      "Age 1 x APOE",
                      "Age 2 x APOE",
                      "Age 3 x APOE",
                      "Diagnosis 1",
                      "Diagnosis 2",
                      "Diagnosis 3",
                      "Diagnosis 1 x APOE",
                      "Diagnosis 2 x APOE",
                      "Diagnosis 3 x APOE",
                      "lateral front-orbital gyrus right 1",
                      "lateral front-orbital gyrus right 2",
                      "lateral front-orbital gyrus right 3",
                      "lateral front-orbital gyrus right 1 x APOE",
                      "lateral front-orbital gyrus right 2 x APOE",
                      "lateral front-orbital gyrus right 3 x APOE",
                      "medial front-orbital gyrus right 1",
                      "medial front-orbital gyrus right 2",
                      "medial front-orbital gyrus right 3",
                      "medial front-orbital gyrus right 1 x APOE",
                      "medial front-orbital gyrus right 2 x APOE",
                      "medial front-orbital gyrus right 3 x APOE",
                      "occipital pole left 1",
                      "occipital pole left 2",
                      "occipital pole left 3",
                      "occipital pole left 1 x APOE",
                      "occipital pole left 2 x APOE",
                      "occipital pole left 3 x APOE",
                      "supramarginal gyrus right 1",
                      "supramarginal gyrus right 2",
                      "supramarginal gyrus right 3",
                      "supramarginal gyrus right 1 x APOE",
                      "supramarginal gyrus right 2 x APOE",
                      "supramarginal gyrus right 3 x APOE")


trop <- RColorBrewer::brewer.pal(7,"Set1")
par(mar=c(5,13,4,2)) # increase y-axis margin.

pdf(file="results/figures/apoe_interactions_sail_coefs_with_diagnosis.pdf",width=11,height=8)#,units="in",res=150)
par(mai=c(1,3.7,0.3,0.2), family="serif")
barplot(sail_coef[match(order_sail_names,sail_coef$mod_names),"s65"],
        horiz=TRUE,
        las=1,
        names.arg=order_sail_names,
        cex.names = 1.2,
        cex.axis = 1.5,
        cex.lab = 1.5,
        cex.main = 1.5,
        col=c(trop[1], rep(trop[2:7], each = 6)),
        main="sail", ylab = "",xlab = "Coefficient")
dev.off()


# Coef plot for the lasso results -----------------------------------------

# this lasso result is when we used glmnet(x = Xmat) where
# Xmat is the design matrix in that contains all main effects and
# their interaction with E

(lasso_coef <- as.data.frame(as.matrix(
  coef(glfit, s = "lambda.min")[nonzeroCoef(coef(glfit, s = "lambda.min"))[-1],,drop=F]), stringsAsFactors = FALSE))

lasso_coef$lasso <- rownames(lasso_coef)
lasso_coef <- left_join(lasso_coef, map_names, by = c("lasso" = "sail"))
str(lasso_coef)
# lasso_coef$orig <- replace(lasso_coef$orig, which(lasso_coef$orig=="E"),"X_E")
# lasso_coef$region <- unlist(stringr::str_extract_all(lasso_coef$orig, "^[^_]+(?=_)"))
lasso_coef$region <- lasso_coef$lasso # since we used the original design matrix

region_names <- read.csv("data/brain_regions.csv", header = F, stringsAsFactors = FALSE)
region_names$V1 <- paste0("X",region_names$V1)

lasso_coef <- left_join(lasso_coef, region_names, by = c("region" = "V1"))
dput(lasso_coef$V2)
# lasso_coef$mod_names <- c("superior frontal gyrus right 2", "superior frontal gyrus right 3",
#                           "lateral front-orbital gyrus right 1",
#                           "medial front-orbital gyrus right 1",
#                           "cingulate region right 1",
#                           "supramarginal gyrus right 1", "supramarginal gyrus right 2",
#                           "angular gyrus left 2", "angular gyrus left 3",
#                           "superior temporal gyrus left 3",
#                           "lateral occipitotemporal gyrus right 1",
#                           "amygdala right 1", "amygdala right 3",
#                           "occipital pole left 1", "occipital pole left 3",
#                           "middle occipital gyrus left 1",
#                           "caudate nucleus left 3",
#                           "putamen right 2",
#                           "lateral ventricle right 1",
#                           "occipital lobe WM right 2",
#                           "nucleus accumbens left 1",
#                           "anterior limb of internal capsule right 2",
#                           "cerebellum right 1",
#                           "Age 1", "Age 2",
#                           "medial front-orbital gyrus left 3 x APOE",
#                           "caudate nucleus left 1 x APOE",
#                           "subarachnoid cerebro-spinal fluid 3 x APOE",
#                           "background 3 x APOE",
#                           "nucleus accumbens left 1 x APOE", "nucleus accumbens left 2 x APOE")
lasso_coef$mod_names <- c("medial front-orbital gyrus left", "supramarginal gyrus right",
                          "angular gyrus left", "lateral occipitotemporal gyrus right",
                          "middle occipital gyrus left", "occipital lobe WM right", "background",
                          "scalp", "nucleus accumbens left", "anterior limb of internal capsule right", "Age", "Diagnosis", "Diagnosis x APOE")

library(RColorBrewer)
# n <- 24
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# col_for_barplot_lasso <- c(rep(col_vector[1],2),rep(col_vector[2],1),rep(col_vector[3],1),rep(col_vector[4],1),rep(col_vector[5],2),rep(col_vector[6],2),
#                            rep(col_vector[7],1),rep(col_vector[8],1),rep(col_vector[9],2),rep(col_vector[10],2),rep(col_vector[11],1),rep(col_vector[12],1),
#                            rep(col_vector[13],1),rep(col_vector[14],1),rep(col_vector[15],1),rep(col_vector[16],1),rep(col_vector[17],1),rep(col_vector[18],1),
#                            rep(col_vector[19],2),rep(col_vector[20],1),rep(col_vector[21],1),rep(col_vector[22],1),rep(col_vector[23],1),rep(col_vector[24],2))
n <- 13
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_for_barplot_lasso <- c(brewer.pal(12,"Paired"),"black" )

pdf(file="results/figures/apoe_interactions_lasso_coefs_with_diagnosis.pdf",width=11,height=8)#,units="in",res=150)
par(mai=c(1,3.85,0.3,0.2), family="serif")
barplot(lasso_coef$`1`,
        horiz=TRUE,
        las=1,
        names.arg=lasso_coef$mod_names,
        cex.names = 1.2,
        cex.axis = 1.5,
        cex.lab = 1.5,
        cex.main = 1.5,
        col=col_for_barplot_lasso,
        main="lasso", ylab = "",xlab = "Coefficient")
dev.off()

RColorBrewer::display.brewer.all()
#
# Interactions with APOE --------------------------------------------------

# i is the real names of the amy beta in the data, j is the corresponding name in sail
# i = "X60";j = "X19";
plot_apoe_inter <- function(cv_obj, X, original_name = "X60", sail_name = "X19", lambda_type = "lambda.min",
                            xlab = "supramarginal gyrus right", ylab = "Mini-Mental State Examination",
                            apoe = FALSE, main = "APOE e4 = 0",
                            color = RSkittleBrewer::RSkittleBrewer(flavor = "trop"), legend = TRUE, ylim =  c(15,30)) {

  # cv_obj = cvfit; original_name = "X60"; sail_name = "X19"; xlab =  "supramarginal gyrus right";
  # lambda_type = "lambda.min";ylab = "Mini-Mental State Examination";
  # color = RColorBrewer::brewer.pal(9,"Set1"); legend = TRUE; ylim =  c(15,30)
  # ==================
  dfs <- cv_obj$sail.fit$df
  lin_pred <- coef(cv_obj, s = lambda_type)["(Intercept)",,drop=T] +
    # cv_obj$sail.fit$design[,paste("X97",seq_len(dfs), sep = "_")] %*%
    # coef(cv_obj, s = lambda_type)[paste("X97",seq_len(dfs), sep = "_"),,drop=F] +
    cv_obj$sail.fit$design[,paste(sail_name,seq_len(dfs), sep = "_")] %*%
    coef(cv_obj, s = lambda_type)[paste(sail_name,seq_len(dfs), sep = "_"),,drop=F] +
    cv_obj$sail.fit$design[,"X_E"] %*%
    coef(cv_obj, s = lambda_type)["X_E",,drop=F] +
    cv_obj$sail.fit$design[,paste0(paste(sail_name,seq_len(dfs), sep = "_"),":X_E")] %*%
    coef(cv_obj, s = lambda_type)[paste0(paste(sail_name,seq_len(dfs), sep = "_"),":X_E"),,drop=F] +
    cv_obj$sail.fit$design[,paste("X98",seq_len(dfs), sep = "_")] %*%
    coef(cv_obj, s = lambda_type)[paste("X98",seq_len(dfs), sep = "_"),,drop=F] +
    cv_obj$sail.fit$design[,paste0(paste("X98",seq_len(dfs), sep = "_"),":X_E")] %*%
    coef(cv_obj, s = lambda_type)[paste0(paste("X98",seq_len(dfs), sep = "_"),":X_E"),,drop=F]

  # all(rownames(coef(cv_obj, s = lambda_type))[-1] ==  colnames(cv_obj$sail.fit$design))
  # lin_pred <- cbind(1,cv_obj$sail.fit$design) %*% coef(cv_obj, s = lambda_type)

  unexposed_index <- which(cv_obj$sail.fit$design[,"X_E"]==0)
  exposed_index <- which(cv_obj$sail.fit$design[,"X_E"]==1)
  # browser()
  # 1=control, 2=MCI (Mild Cognitive Impairment) and 3=Alzeimer Disease
  cont_index0 <- which(X[,"diag_3bl.x"]==1 & cv_obj$sail.fit$design[,"X_E"]==0)
  cont_index1 <- which(X[,"diag_3bl.x"]==1 & cv_obj$sail.fit$design[,"X_E"]==1)
  mci_index0 <- which(X[,"diag_3bl.x"]==2 & cv_obj$sail.fit$design[,"X_E"]==0)
  mci_index1 <- which(X[,"diag_3bl.x"]==2 & cv_obj$sail.fit$design[,"X_E"]==1)
  az_index0 <- which(X[,"diag_3bl.x"]==3 & cv_obj$sail.fit$design[,"X_E"]==0)
  az_index1 <- which(X[,"diag_3bl.x"]==3 & cv_obj$sail.fit$design[,"X_E"]==1)

  # c(cont_index0, cont_index1, mci_index0, mci_index1, az_index0, az_index1) %>% length()

  e0 <- lin_pred[unexposed_index,]
  e1 <- lin_pred[exposed_index,]

  cont_pred0 <- lin_pred[cont_index0,]
  cont_pred1 <- lin_pred[cont_index1,]
  mci_pred0 <- lin_pred[mci_index0,]
  mci_pred1 <- lin_pred[mci_index1,]
  az_pred0 <- lin_pred[az_index0,]
  az_pred1 <- lin_pred[az_index1,]

  min.length.top <- range(lin_pred)[1] ; max.length.top <- range(lin_pred)[2]
  par(mai=c(1,1,1,0.2))
  plot(X[,original_name], lin_pred,
       pch = 19,
       ylab = ylab,
       xlab = xlab,
       col = color[1],
       bty="n",
       xaxt="n",
       type = "n",
       cex.lab = 2,
       cex.axis = 2,
       cex = 2,
       main = main,
       cex.main = 2.5,
       # ylim = c(min.length.top-3, max.length.top+3),
       ylim = ylim)
  axis(1, labels = T, cex.axis = 2)

  # points(X[exposed_index,original_name], e1, pch = 19, col = color[2], cex = 1.5)
  # points(X[unexposed_index,original_name], e0, pch = 19, col = color[1], cex = 1.5)
  if (apoe) {
    points(X[cont_index1,original_name], cont_pred1, pch = 19, col = color[1], cex = 1.5)
    points(X[mci_index1,original_name], mci_pred1, pch = 19, col = color[2], cex = 1.5)
    points(X[az_index1,original_name], az_pred1, pch = 19, col = color[3], cex = 1.5)
    rug(X[exposed_index,original_name], side = 1)

  } else {
    points(X[cont_index0,original_name], cont_pred0, pch = 19, col = color[1], cex = 1.5)
    points(X[mci_index0,original_name], mci_pred0, pch = 19, col = color[2], cex = 1.5)
    points(X[az_index0,original_name], az_pred0, pch = 19, col = color[3], cex = 1.5)
    rug(X[unexposed_index,original_name], side = 1)
  }
  # if (legend) legend("bottomright", c("APOE = 1","APOE = 0"), col = color[2:1], pch = 19, cex = 2, bty = "n")
  # text(3, main, cex = 2)
  if (legend) legend("bottom", c("Control", "Mild Cognitive Impairment","Alzeimer Disease"),
                     col = color[1:3], pch = 19, cex = 2, bty = "n")

  # rug(X[,original_name], side = 1)
}


dev.off()
# trop <- RSkittleBrewer::RSkittleBrewer("trop")
pdf(file="results/figures/apoe_interactions_mmscore_amy_beta.pdf",width=13,height=11)#,units="in",res=150)
par(mfrow=c(2,2), tcl=-0.5, family="serif", omi=c(0.2,0.2,0,0))
plot_apoe_inter(cv_obj = cvfit, X = X, original_name = "X6", sail_name = "X9", xlab = "lateral front-orbital gyrus right", legend = T)
plot_apoe_inter(cv_obj = cvfit, X = X, original_name = "X1", sail_name = "X11", xlab =  "medial front-orbital gyrus right", legend = F)
plot_apoe_inter(cv_obj = cvfit, X = X, original_name = "X60", sail_name = "X19", xlab =  "supramarginal gyrus right", legend = F)
plot_apoe_inter(cv_obj = cvfit, X = X, original_name = "X251", sail_name = "X42", xlab =  "occipital pole left")
dev.off()


pdf(file="results/figures/apoe_interactions_mmscore_amy_beta_supramarginal_gyrus_right_diag.pdf",width=11,height=8)#,units="in",res=150)
par(mfrow=c(1,2), tcl=-0.5, family="serif",
    omi=c(0.2,0.2,0,0))
# mai = c(1,1,2,0))
# oma = c(1,1,2,1))
plot_apoe_inter(cv_obj = cvfit, X = X, original_name = "X60", sail_name = "X19", xlab =  "supramarginal gyrus right", legend = T,
                apoe = FALSE)
plot_apoe_inter(cv_obj = cvfit, X = X, original_name = "X60", sail_name = "X19", xlab =  "supramarginal gyrus right", legend = F,
                apoe = TRUE, ylab = "", main = "APOE e4 = 1")
dev.off()
#
#
# # Interactions with diag --------------------------------------------------
#
# # i is the real names of the amy beta in the data, j is the corresponding name in sail
# # i = "X60";j = "X19";
# plot_diag_inter <- function(cv_obj, X, original_name = "X60", sail_name = "X19", lambda_type = "lambda.min",
#                             xlab = "supramarginal gyrus right", ylab = "Mini-Mental State Examination",
#                             color = trop[1:3], legend = TRUE, ylim =  c(15,30)) {
#   dfs <- cv_obj$sail.fit$df; lambda_type = "lambda.min"
#   lin_pred <- coef(cv_obj, s = lambda_type)["(Intercept)",,drop=T] +
#     # cv_obj$sail.fit$design[,paste("X97",seq_len(dfs), sep = "_")] %*%
#     # coef(cv_obj, s = lambda_type)[paste("X97",seq_len(dfs), sep = "_"),,drop=F] +
#     cv_obj$sail.fit$design[,paste(sail_name,seq_len(dfs), sep = "_")] %*%
#     coef(cv_obj, s = lambda_type)[paste(sail_name,seq_len(dfs), sep = "_"),,drop=F] +
#     cv_obj$sail.fit$design[,"X_E"] %*%
#     coef(cv_obj, s = lambda_type)["X_E",,drop=F] +
#     cv_obj$sail.fit$design[,paste0(paste(sail_name,seq_len(dfs), sep = "_"),":X_E")] %*%
#     coef(cv_obj, s = lambda_type)[paste0(paste(sail_name,seq_len(dfs), sep = "_"),":X_E"),,drop=F]
#
#   # 1=control, 2=MCI (Mild Cognitive Impairment) and 3=Alzeimer Disease
#
#   control_index <- which(cv_obj$sail.fit$design[,"X_E"]==1)
#   mci_index <- which(cv_obj$sail.fit$design[,"X_E"]==2)
#   ad_index <- which(cv_obj$sail.fit$design[,"X_E"]==3)
#
#   e1 <- lin_pred[control_index,]
#   e2 <- lin_pred[mci_index,]
#   e3 <- lin_pred[ad_index,]
#
#   min.length.top <- range(lin_pred)[1] ; max.length.top <- range(lin_pred)[2]
#   par(mai=c(1,1,0.1,0.2))
#   plot(X[,original_name], lin_pred,
#        pch = 19,
#        ylab = ylab,
#        xlab = xlab,
#        col = color[1],
#        bty="n",
#        xaxt="n",
#        type = "n",
#        cex.lab = 2,
#        cex.axis = 2,
#        cex = 2,
#        # ylim = c(min.length.top-3, max.length.top+3),
#        ylim = ylim)
#   axis(1, labels = T, cex.axis = 2)
#   points(X[control_index,original_name], e1, pch = 19, col = color[1], cex = 1.5)
#   points(X[mci_index,original_name], e2, pch = 19, col = color[2], cex = 1.5)
#   points(X[ad_index,original_name], e3, pch = 19, col = color[3], cex = 1.5)
#   if (legend) legend("bottomright", c("Control","Mild Cognitive Impairment","Alzeimers"), col = color[1:3],
#                      pch = 19, cex = 2, bty = "n")
#   rug(X[,original_name], side = 1)
# }
#
#
# dev.off()
# trop <- RSkittleBrewer::RSkittleBrewer("trop")
# pdf(file="results/figures/apoe_interactions_mmscore_amy_beta.pdf",width=13,height=11)#,units="in",res=150)
# par(mfrow=c(2,2), tcl=-0.5, family="serif", omi=c(0.2,0.2,0,0))
# plot_diag_inter(cv_obj = cvfit, X = X, original_name = "X6", sail_name = "X9", xlab = "lateral front-orbital gyrus left", legend = T)
# plot_diag_inter(cv_obj = cvfit, X = X, original_name = "X60", sail_name = "X19", xlab =  "putamen right", legend = T)
# plot_diag_inter(cv_obj = cvfit, X = X, original_name = "Age_bl", sail_name = "X97", xlab =  "thalamus left", legend = T)
# plot_apoe_inter(cv_obj = cvfit, X = X, original_name = "X60", sail_name = "X19", xlab =  "supramarginal gyrus right")
# dev.off()
#
#
# pdf(file="results/figures/apoe_interactions_mmscore_amy_beta_supramarginal_gyrus_right.pdf",width=11,height=8)#,units="in",res=150)
# plot_apoe_inter(cv_obj = cvfit, X = X, original_name = "X60", sail_name = "X19", xlab =  "supramarginal gyrus right",
#                 ylim = c(20, 30))
# dev.off()
#
#
#
#
#
# #
# #
# #
# #
# #
# # # par(mfrow=c(2,5), tcl=-0.5, family="serif", omi=c(0.2,0.2,0,0))
# # dev.off()
# # par(mai=c(0.55,0.61,0.1,0.2))
# # plot(X[,i], coef(cvfit, s = lambda_type)["(Intercept)",,drop=T] +
# #        # cvfit$sail.fit$design[,paste("X97",seq_len(dfs), sep = "_")] %*%
# #        # coef(cvfit, s = lambda_type)[paste("X97",seq_len(dfs), sep = "_"),,drop=F] +
# #        cvfit$sail.fit$design[,paste(j,seq_len(dfs), sep = "_")] %*%
# #        coef(cvfit, s = lambda_type)[paste(j,seq_len(dfs), sep = "_"),,drop=F],
# #      pch = 19,
# #      ylab = TeX(sprintf("$f(x_{%s})$",as.numeric(gsub("X","",i)))),
# #      xlab = TeX(sprintf("$x_{%s}$",as.numeric(gsub("X","",i)))),
# #      col = trop[1],
# #      bty="n",
# #      xaxt="n",
# #      cex.lab = 2,
# #      cex.axis = 2)#,
# #      # ylim = c(min.length.top, max.length.top))
# # axis(1, labels = F)
# # points(X[,i], lin_pred, pch = 19, col = trop[2])
# # legend(0.0,-1.3, c("Truth","Estimated"), col = trop[2:1], pch = 19, cex = 2, bty = "n")
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# # i = "X60";j = "X19"
# # x <- seq(range(X[,i])[1], range(X[,i])[2], length.out = 30)
# # e <- seq(range(E)[1], range(E)[2], length.out = 30)
# # f.est.0 <- function(x, e) { e * splines::bs(x, cvfit$sail.fit$df) %*%
# #     as.matrix(coef(cvfit, s = "lambda.min")[paste0(j,"_", seq_len(cvfit$sail.fit$df),":X_E"),]) }
# # f.est.1 <- function(x, e) { 1 * splines::bs(x, cvfit$sail.fit$df) %*%
# #     as.matrix(coef(cvfit, s = "lambda.min")[paste0(j,"_", seq_len(cvfit$sail.fit$df),":X_E"),]) }
# # z.est <- outer(x, e, f.est)
# # # z.truth <- outer(x, e, f.truth)
# #
# # op <- par(bg = "white")
# # z_range <- c(min(z.est), max(z.est))
# # # z_range <- c(-1.868998, 2.118003)
# # trop <- RSkittleBrewer::RSkittleBrewer("trop")
# #
# # dev.off()
# # pdf(file="~/Dropbox/jobs/hec/talk/gendata_inter_X1.pdf",width=11,height=8)#,units="in",res=150)
# # # par(mfrow=c(1,2), tcl=-0.5, family="serif", omi=c(0.2,0.2,0,0))
# # par(mai=c(0.,0.,0.3,0.))
# # persp(x, e, z.est,
# #       theta=30, phi=30,
# #       ltheta = 120, expand = 0.5,
# #       r=2, shade=0.3, axes=TRUE,scale=TRUE, box=T,
# #       nticks=5,
# #       zlim = z_range,
# #       # ticktype="detailed",
# #       col=trop[1],
# #       xlab=sprintf("f(x_%s)",as.numeric(gsub("X","",i))),
# #       ylab="X_E",
# #       zlab="Y", main="Estimated X_E*f(X_1)")
# # dev.off()
# #
# #
# #
