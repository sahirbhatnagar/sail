rm(list=ls())
# devtools::document()
devtools::load_all()
# source("R/ksfilter.R")
pacman::p_load(data.table)
pacman::p_load(magrittr)
pacman::p_load(genefilter)
pacman::p_load(tidyverse)
pacman::p_load(doParallel)
pacman::p_load(splines)
pacman::p_load(foreach)
pacman::p_load(doMC)
pacman::p_load(glmnet)
pacman::p_load_current_gh('sahirbhatnagar/glinternet')
# pacman::p_load(mice)


# amy_pheno <- xlsx::read.xlsx("~/Downloads/DrCelia_data.xlsx", sheetIndex = 1)
# amy_mat <- read.csv("~/git_repositories/sail/data-nogit/adni_new/csf_amyloid_final.csv", stringsAsFactors = FALSE)
amy_mat <- read.csv("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_git_v2/sail/rda/csf_amyloid_final.csv", stringsAsFactors = FALSE)
covr <- read.csv("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_git_v2/sail/rda/covariates.csv", stringsAsFactors = FALSE, sep = ";")
# surv <- read.csv("~/git_repositories/sail/data/adni_new/fdg_info.csv", stringsAsFactors = FALSE, sep = ",")

# sum(as.character(amy_pheno$PTID) %in% amy_mat$PTID)
# sum(as.character(covr$IID) %in% amy_mat$PTID)
# as.character(covr$IID) %>% unique() %>% length()

DT <- dplyr::inner_join(amy_mat, covr, by = c("PTID" = "IID")) %>%
  select(-AV45_path_bl)
colnames(DT)
DT$diag_3bl.x %>% table
# DT <- DT %>% mutate(diag_3bl.x = ifelse(diag_3bl.x==3, 2, diag_3bl.x))
DT$diag_3bl.x %>% table
# DT <- DT %>% dplyr::filter(diag_3bl.x == 3)

brain_regions <- grep("X", colnames(DT), value=T)
# fmla <- reformulate(c(sapply(brain_regions, function(i) sprintf("bs(%s,3)",i)),
#                       "diag_3bl.x"), intercept = FALSE)
fmla <- reformulate(c(sapply(brain_regions, function(i) sprintf("bs(%s,3)",i)),
                      "APOE_bin"), intercept = FALSE)
# fmla <- reformulate(c(sapply(brain_regions, function(i) sprintf("bs(%s, degree = 5)",i)),
#                       "diag_3bl.y"), intercept = FALSE)
# fmla <- reformulate(c(sapply(brain_regions, function(i) sprintf("bs(%s)",i)),
#                       "bs(Age_bl)"), intercept = FALSE)
# fmla <- reformulate(c(sapply(brain_regions, function(i) sprintf("bs(%s)",i)),
#                       "bs(Age_bl)"), intercept = FALSE)
# fmla <- reformulate(c(sapply(brain_regions, function(i) sprintf("bs(%s)",i))),
#                     intercept = FALSE)

# train <- caret::createDataPartition(DT$MMSCORE_bl)[[1]]
# train <- seq(nrow(DT))
# test <- setdiff(seq(nrow(DT)), train)


# DT$diag_3bl.x[train] %>% table
# DT$diag_3bl.x[test] %>% table

X <- DT %>% select(starts_with("X"), diag_3bl.x, APOE_bin) %>%
  mutate(diag_3bl.x = diag_3bl.x - 1) %>%
  as.matrix()
# X <- DT %>% select(starts_with("X"), Age_bl) %>%
#   as.matrix()

Xnorm <- sail:::standardize(X, center = TRUE, normalize = TRUE)$x

model_mat <- model.matrix(fmla, data = as.data.frame(Xnorm))

# head(model_mat)

# X <- DT %>% select(starts_with("X"), Age_bl, diag_3bl.x) %>% as.matrix()

# X <- DT %>% select(starts_with("X"), Age_bl, EDUCAT) %>% as.matrix()
# X <- DT %>% select(starts_with("X")) %>% as.matrix()
# X <- DT %>% select(starts_with("X"), diag_3bl.x) %>% as.matrix()
# dimnames(X)[[1]] <- DT$PTID
# X <- X[,80:97]
# colnames(X)
# ind <- which(DT$diag_3bl.x==2)
# X <- X[ind,,drop = F]

# E <- DT[train,] %>% pull(APOE_bin) %>% as.numeric
# E <- E[ind]
# E <- DT %>% pull(EDUCAT) %>% as.numeric
# E <- Xnorm[, "EDUCAT"]
# E <- DT[train,] %>% pull(diag_3bl.x) %>% as.numeric
# E <- DT[train,] %>% pull(Age_bl) %>% as.numeric

# the minus 1 for diag is for glinternet which needs 0,1,2,
# X <- DT %>% select(starts_with("X"), Age_bl, diag_3bl.x) %>%
#   mutate(diag_3bl.x = diag_3bl.x - 1) %>%
#   as.matrix()
X[,"diag_3bl.x"] %>% table
# E <- DT %>% pull(APOE_bin) %>% as.numeric
E <- Xnorm[, "diag_3bl.x"]
Y <- DT %>% pull(MMSCORE_bl) %>% as.numeric

set.seed(12345)
dat_lasso <- sail:::partition_data(x = X[,-which(colnames(X) %in% c("diag_3bl.x"))],
                                   y = Y, e = X[,"diag_3bl.x"], p = 200/length(Y),
                                   partition_on = Xnorm[, "diag_3bl.x"], type = "train_test_val")
dat <- sail:::partition_data(x = model_mat, y = Y, e = E, p = 200/length(Y),
                             partition_on = Xnorm[,"diag_3bl.x"], type = "train_test_val")

# dat_lasso <- sail:::partition_data(x = X, y = Y, e = E, p = 1,
#                                    partition_on = Y, type = "train_test")
# dat <- sail:::partition_data(x = model_mat, y = Y, e = E, p = 1,
#                              partition_on = Y, type = "train_test")

f.basis <- function(i) splines::bs(i, df = 5)
# fit <- sail(x = dat$xtrain, y = dat$ytrain, e = dat$etrain,
#             basis = f.basis,
#             alpha = 0.1,
#             thresh = 1e-02,
#             strong = FALSE,
#             verbose = 2)

#strong=TRUE works better for ADNI
fit <- sail(x = dat$xtrain, y = dat$ytrain, e = dat$etrain,
            # basis = f.basis,
            expand = FALSE,
            center.x = F,
            center.e = T,
            group = attr(model_mat, "assign"),
            alpha = 0.1,
            maxit = 250,
            # thresh = 1e-02,
            strong = TRUE,
            verbose = 2)
plot(fit)
ytest_hat <- predict(fit, newx = dat$xtest, newe = dat$etest)
msetest <- colMeans((dat$ytest - ytest_hat)^2)
plot(log(fit$lambda), msetest)
lambda.min.index <- as.numeric(which.min(msetest))
lambda.min <- fit$lambda[which.min(msetest)]

yvalid_hat <- predict(fit, newx = dat$xvalid, newe = dat$evalid, s = lambda.min)
(msevalid <- mean((dat$yvalid - drop(yvalid_hat))^2))

(nzcoef <- predict(fit, s = lambda.min, type = "nonzero"))


plotMainADNI(fit, x = X[dat$train_ind,"X175"], xvar = paste0("bs(X175, 3)",1:3), s = lambda.min,
         ylab = "f(Distance to Airport (YYZ) (m))", xlab = "Distance to Airport (YYZ) (m)")
plotMainADNI(fit, x = X[dat$train_ind,"X243"], xvar = paste0("bs(X243, 3)",1:3), s = lambda.min,
             ylab = "f(Distance to Airport (YYZ) (m))", xlab = "Distance to Airport (YYZ) (m)")

dev.off()

plotInterADNI(fit, x = X[,"X175"], xvar = paste0("bs(X175, 3)",1:3), s = lambda.min,
              model_matrix = model_mat, e = E,
              ylab = "f(X175)", xlab = "X175", main = "")
plotInterADNI(fit, x = X[,"X243"], xvar = paste0("bs(X243, 3)",1:3), s = lambda.min,
              model_matrix = model_mat, e = E,
              ylab = "f(X243)", xlab = "X243", main = "")

View(plotInter)

# this extrapolates to the entire sample. not just the training set
plotMainADNI <- function(object, x, design, xvar, s, f.truth, col = c("#D55E00", "#009E73"),
                         legend.position = "bottomleft", rug = TRUE, ...) {

  # browser()

  if (length(s) > 1) {
    s <- s[[1]]
    warning("More than 1 s value provided. Only first element will be used for the estimated coefficients.")
  }

  ind <- which(object$vnames %in% xvar)
  allCoefs <- coef(object, s = s)
  a0 <- allCoefs[1, ]
  betas <- as.matrix(allCoefs[object$main.effect.names[ind], , drop = FALSE])
  design.mat <- design[, object$main.effect.names[ind], drop = FALSE]
  originalX <- x

  # f.hat <- drop(a0 + design.mat %*% betas)
  f.hat <- drop(design.mat %*% betas)
  if (!missing(f.truth)) {
    seqs <- seq(range(originalX)[1], range(originalX)[2], length.out = 100)
    f.truth.eval <- f.truth(seqs)
    ylims <- range(f.truth.eval, f.hat)
  } else {
    ylims <- range(f.hat)
  }

  plot.args <- list(
    x = originalX[order(originalX)],
    y = f.hat[order(originalX)],
    ylim = ylims,
    xlab = xvar,
    ylab = sprintf("f(%s)", xvar),
    type = "n",
    # xlim = rev(range(l)),
    # las = 1,
    cex.lab = 1.5,
    cex.axis = 1.5,
    cex = 1.5,
    # bty = "n",
    # mai=c(1,1,0.1,0.2),
    # tcl = -0.5,
    # omi = c(0.2,1,0.2,0.2),
    family = "serif"
  )
  new.args <- list(...)
  if (length(new.args)) {
    new.plot.args <- new.args[names(new.args) %in% c(
      names(par()),
      names(formals(plot.default))
    )]
    plot.args[names(new.plot.args)] <- new.plot.args
  }
  do.call("plot", plot.args)
  abline(h = 0, lwd = 1, col = "gray")
  lines(originalX[order(originalX)], f.hat[order(originalX)], col = col[1], lwd = 3)
  if (rug) graphics::rug(originalX, side = 1)
  if (!missing(f.truth)) {
    lines(seqs[order(seqs)], f.truth.eval[order(seqs)], col = col[2], lwd = 3)
  }
  if (!missing(f.truth)) {
    legend(legend.position,
           c("Estimated", "Truth"),
           col = col, cex = 1, bty = "n", lwd = 3
    )
  }
}


# this extrapolates to the entire sample. not just the training set
plotInterADNI <- function(object, x, xvar, s,
                          design, # this contains user defined expand matrix
                          e, # this is E vector for whole sample
                          apoe = TRUE,
                          xlab = "supramarginal gyrus right", ylab = "Mini-Mental State Examination",
                          legend.position = "bottomleft", main = "", rug = TRUE,
                          color = sail:::cbbPalette[c(6,4,7)], legend = TRUE) {

  # cv_obj = cvfit; original_name = "X60"; sail_name = "X19"; xlab =  "supramarginal gyrus right";
  # lambda_type = "lambda.min";ylab = "Mini-Mental State Examination";
  # color = RColorBrewer::brewer.pal(9,"Set1"); legend = TRUE; ylim =  c(15,30)
  # ==================
  # browser()
  ind <- which(object$vnames %in% xvar)
  allCoefs <- coef(object, s = s)
  a0 <- allCoefs[1, ]
  betas <- as.matrix(allCoefs[object$main.effect.names[ind], , drop = FALSE])
  alphas <- as.matrix(allCoefs[object$interaction.names[ind], , drop = FALSE])
  betaE <- as.matrix(allCoefs["E", , drop = FALSE])
  betaAPOE <- as.matrix(allCoefs["APOE_bin", , drop = FALSE])
  betaAPOEinter <- as.matrix(allCoefs["APOE_bin:E", , drop = FALSE])
  # if you dont want to extrapolate, un-comment the following lines
  # design.mat.main <- object$design[, object$main.effect.names[ind], drop = FALSE]
  # design.mat.int <- object$design[, object$interaction.names[ind], drop = FALSE]

  design.mat.main <- design[, object$main.effect.names[ind], drop = FALSE]
  design.mat.int <- design[, object$main.effect.names[ind], drop = FALSE] * e
  apoee4 <- design[, "APOE_bin"]
  apoee4inter <- design[, "APOE_bin"] * e

  # originalE <- object$design[, "E", drop = FALSE] # this is the centered E
  # originalX <- x

  originalE <- e # this is the centered E
  originalX <- x


  # f.hat <- drop(a0 + design.mat %*% betas)
  f.hat <- drop(originalE * as.vector(betaE) + apoee4 * as.vector(betaAPOE) +
                  design.mat.main %*% betas + design.mat.int %*% alphas +
                  apoee4inter * as.vector(betaAPOEinter))
  # f.hat <- drop(originalE * as.vector(betaE)  + design.mat.int %*% alphas)
  # f.hat <- drop(design.mat.int %*% alphas)
  ylims <- range(f.hat)

  # dfs <- cv_obj$sail.fit$df
  # lin_pred <- coef(cv_obj, s = lambda_type)["(Intercept)",,drop=T] +
  #   # cv_obj$sail.fit$design[,paste("X97",seq_len(dfs), sep = "_")] %*%
  #   # coef(cv_obj, s = lambda_type)[paste("X97",seq_len(dfs), sep = "_"),,drop=F] +
  #   cv_obj$sail.fit$design[,paste(sail_name,seq_len(dfs), sep = "_")] %*%
  #   coef(cv_obj, s = lambda_type)[paste(sail_name,seq_len(dfs), sep = "_"),,drop=F] +
  #   cv_obj$sail.fit$design[,"X_E"] %*%
  #   coef(cv_obj, s = lambda_type)["X_E",,drop=F] +
  #   cv_obj$sail.fit$design[,paste0(paste(sail_name,seq_len(dfs), sep = "_"),":X_E")] %*%
  #   coef(cv_obj, s = lambda_type)[paste0(paste(sail_name,seq_len(dfs), sep = "_"),":X_E"),,drop=F] +
  #   cv_obj$sail.fit$design[,paste("X98",seq_len(dfs), sep = "_")] %*%
  #   coef(cv_obj, s = lambda_type)[paste("X98",seq_len(dfs), sep = "_"),,drop=F] +
  #   cv_obj$sail.fit$design[,paste0(paste("X98",seq_len(dfs), sep = "_"),":X_E")] %*%
  #   coef(cv_obj, s = lambda_type)[paste0(paste("X98",seq_len(dfs), sep = "_"),":X_E"),,drop=F]

  # all(rownames(coef(cv_obj, s = lambda_type))[-1] ==  colnames(cv_obj$sail.fit$design))
  # lin_pred <- cbind(1,cv_obj$sail.fit$design) %*% coef(cv_obj, s = lambda_type)


  control = drop(unique(originalE))[1]
  mci = drop(unique(originalE))[2]
  ad = drop(unique(originalE))[3]

  apoe_no <- drop(unique(apoee4))[1]
  apoe_yes <- drop(unique(apoee4))[2]

  cont_index0 <- which(apoee4==apoe_no & originalE==control)
  cont_index1 <- which(apoee4==apoe_yes & originalE==control)
  mci_index0 <- which(apoee4==apoe_no & originalE==mci)
  mci_index1 <- which(apoee4==apoe_yes & originalE==mci)
  ad_index0 <- which(apoee4==apoe_no & originalE==ad)
  ad_index1 <- which(apoee4==apoe_yes & originalE==ad)


  # browser()
  # 1=control, 2=MCI (Mild Cognitive Impairment) and 3=Alzeimer Disease
  # cont_index <- which(originalE==control)
  # mci_index <- which(originalE==mci)
  # ad_index <- which(originalE==ad)
  #
  # cont_pred <- f.hat[cont_index]
  # mci_pred <- f.hat[mci_index]
  # ad_pred <- f.hat[ad_index]

  cont_pred0 <- f.hat[cont_index0]
  cont_pred1 <- f.hat[cont_index1]
  mci_pred0 <- f.hat[mci_index0]
  mci_pred1 <- f.hat[mci_index1]
  ad_pred0 <- f.hat[ad_index0]
  ad_pred1 <- f.hat[ad_index1]

  min.length.top <- range(f.hat)[1] ; max.length.top <- range(f.hat)[2]
  par(mai=c(1,1,1,0.2))
  plot(originalX, f.hat,
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
       ylim = ylims)
  axis(1, labels = T, cex.axis = 2)

  if (apoe) {
    # points(X[exposed_index,original_name], e1, pch = 19, col = color[2], cex = 1.5)
    # points(X[unexposed_index,original_name], e0, pch = 19, col = color[1], cex = 1.5)
    lines(originalX[cont_index1][order(originalX[cont_index1])], cont_pred1[order(originalX[cont_index1])], col = color[1], lwd = 3)
    lines(originalX[mci_index1][order(originalX[mci_index1])], mci_pred1[order(originalX[mci_index1])], col = color[2], lwd = 3)
    lines(originalX[ad_index1][order(originalX[ad_index1])], ad_pred1[order(originalX[ad_index1])], col = color[3], lwd = 3)
  } else {
    lines(originalX[cont_index0][order(originalX[cont_index0])], cont_pred0[order(originalX[cont_index0])], col = color[1], lwd = 3)
    lines(originalX[mci_index0][order(originalX[mci_index0])], mci_pred0[order(originalX[mci_index0])], col = color[2], lwd = 3)
    lines(originalX[ad_index0][order(originalX[ad_index0])], ad_pred0[order(originalX[ad_index0])], col = color[3], lwd = 3)
  }

  # if (legend) legend("bottomright", c("APOE = 1","APOE = 0"), col = color[2:1], pch = 19, cex = 2, bty = "n")
  # text(3, main, cex = 2)
  if (legend) legend(legend.position, c("Control", "Mild Cognitive Impairment","Alzeimer Disease"),
                     col = color[1:3], pch = 19, cex = 2, bty = "n")

  if (rug) graphics::rug(originalX, side = 1)
}





registerDoMC(cores = 10)
dat$xtrain %>% dim
cvfit <- cv.sail(x = dat$xtrain, y = dat$ytrain, e = dat$etrain,
                 # basis = f.basis,
               expand = FALSE,
               center.x = F,
               center.e = T,
               group = attr(model_mat, "assign"),
               alpha = 0.1,
               maxit = 250,
               parallel = TRUE,
               nfolds = 10,
               # thresh = 1e-02,
               strong = TRUE,
               verbose = 2)
plot(cvfit)
predict(cvfit, type="non", s = "lambda.min")

yvalid_hat <- predict(fit, newx = dat$xtest, newe = dat$etest, s = fit$lambda.min)
(msevalid <- mean((dat$ytest - drop(yvalid_hat))^2))

(nzcoef <- predict(fit, s = fit$lambda.min, type = "nonzero"))


fitGL <- glinternet(X = dat_lasso$xtrain_lasso, Y = dat_lasso$ytrain,
                    numLevels = c(3, rep(1, ncol(dat_lasso$xtrain_lasso)-2), 2),
                    # numLevels = c(2, rep(1, ncol(dat_lasso$xtrain_lasso)-1)),
                    nLambda = 100, interactionCandidates = c(1),
                    verbose = T)

ytest_hat <- predict(fitGL, X = dat_lasso$xtest_lasso)
msetest <- colMeans((dat_lasso$ytest - ytest_hat)^2)
lambda.min.index <- as.numeric(which.min(msetest))
lambda.min <- fitGL$lambda[which.min(msetest)]

yvalid_hat <- predict(fitGL, X = dat_lasso$xvalid_lasso, lambda = lambda.min)
(msevalid <- mean((dat_lasso$yvalid - drop(yvalid_hat))^2))
tc <- coef(fitGL, lambdaIndex = lambda.min.index)

(mains <- c(colnames(dat_lasso$xtrain_lasso)[tc[[1]]$mainEffects$cont],
  colnames(dat_lasso$xtrain_lasso)[tc[[1]]$mainEffects$cat]))

(inters <- c(paste(colnames(dat_lasso$xtrain_lasso)[tc[[1]]$interactions$catcont[,1]],
                 colnames(dat_lasso$xtrain_lasso)[tc[[1]]$interactions$catcont[,2]], sep=":"),
             paste(colnames(dat_lasso$xtrain_lasso)[tc[[1]]$interactions$contcont[,1]],
                   colnames(dat_lasso$xtrain_lasso)[tc[[1]]$interactions$contcont[,2]], sep=":"),
             paste(colnames(dat_lasso$xtrain_lasso)[tc[[1]]$interactions$catcat[,1]],
                   colnames(dat_lasso$xtrain_lasso)[tc[[1]]$interactions$catcat[,2]], sep=":")))


fitGL <- glinternet.cv(X = dat_lasso$xtrain_lasso, Y = dat_lasso$ytrain,
                       family = "gaussian",
                       nFolds = 10,
                       # numLevels = c(2, rep(1, ncol(dat_lasso$xtrain_lasso)-2), 3),
                       numLevels = c(2, rep(1, ncol(dat_lasso$xtrain_lasso)-1)),
                       nLambda = 100, interactionCandidates = c(1),
                       verbose = T)
plot(fitGL)
yvalid_hat <- predict(fitGL, X = dat_lasso$xtest_lasso, lambda = "lambdaHat")
(msevalid <- mean((dat_lasso$ytest - drop(yvalid_hat))^2))


tc <- coef(fitGL, lambdaType = "lambdaHat")
(mains <- colnames(dat_lasso$xtrain_lasso)[tc[[1]]$mainEffects$cont])
(inters <- paste(colnames(dat_lasso$xtrain_lasso)[tc[[1]]$interactions$catcont[,1]],
                 colnames(dat_lasso$xtrain_lasso)[tc[[1]]$interactions$catcont[,2]], sep=":"))



fit <- glmnet(x = dat_lasso[["xtrain_lasso"]], y = dat_lasso[["ytrain"]],
              alpha = 1)

ytest_hat <- predict(fit, newx = dat_lasso[["xtest_lasso"]])
msetest <- colMeans((dat_lasso[["ytest"]] - ytest_hat)^2)
lambda.min.index <- as.numeric(which.min(msetest))
lambda.min <- fit$lambda[which.min(msetest)]

yvalid_hat <- predict(fit, newx = dat_lasso[["xvalid_lasso"]], s = lambda.min)
(msevalid <- mean((dat_lasso[["yvalid"]] - drop(yvalid_hat))^2))

nzcoef <- coef(fit, s = lambda.min)[nonzeroCoef(coef(fit, s = lambda.min)),,drop=F]


fit <- cv.glmnet(x = dat_lasso[["xtrain_lasso"]], y = dat_lasso[["ytrain"]],
              alpha = 1, nfolds = 10)
plot(fit)
yvalid_hat <- predict(fit, newx = dat_lasso[["xtest_lasso"]], s = fit$lambda.min)
(msevalid <- mean((dat_lasso[["ytest"]] - drop(yvalid_hat))^2))

nzcoef <- coef(fit, s = lambda.min)[nonzeroCoef(coef(fit, s = lambda.min)),,drop=F]








# simulator style ---------------------------------------------------------

rm(list=ls())
pacman::p_load(simulator) # this file was created under simulator version 0.2.1
devtools::load_all()
pacman::p_load(data.table)
pacman::p_load(magrittr)
pacman::p_load(genefilter)
pacman::p_load(tidyverse)
pacman::p_load(doParallel)
pacman::p_load(splines)
pacman::p_load(foreach)
pacman::p_load(doMC)
pacman::p_load(glmnet)
pacman::p_load_current_gh('sahirbhatnagar/glinternet')
pacman::p_load(LassoBacktracking)
pacman::p_load(SAM)
pacman::p_load(gamsel)
pacman::p_load_gh('asadharis/HierBasis')

source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_git_v2/sail/my_sims/model_functions.R")
source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_git_v2/sail/my_sims/method_functions_rda.R")
source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_git_v2/sail/my_sims/eval_functions.R")

amy_mat <- read.csv("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_git_v2/sail/rda/csf_amyloid_final.csv", stringsAsFactors = FALSE)
covr <- read.csv("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_git_v2/sail/rda/covariates.csv", stringsAsFactors = FALSE, sep = ";")
# used for plotting
DT <- dplyr::inner_join(amy_mat, covr, by = c("PTID" = "IID")) %>%
  select(-AV45_path_bl)
X <- DT %>% select(starts_with("X"), diag_3bl.x, APOE_bin) %>%
  mutate(diag_3bl.x = diag_3bl.x - 1) %>%
  as.matrix()
Xnorm <- sail:::standardize(X, center = TRUE, normalize = TRUE)$x
E <- Xnorm[, "diag_3bl.x"]

# we didnt fit spam or gamsel because it wouldnt let expansion of binary predictors
# sim <- new_simulation(name = "rda_ADNI_may_17_2018",
#                       label = "rda_ADNI_may_17_2018",
#                       dir = ".") %>%
#   generate_model(make_ADNI_data_split, seed = 12345,
#                  phenoVariable = "MMSCORE_bl",
#                  exposure = "diag_3bl.x",n_train_test = 200) %>%
#   simulate_from_model(nsim = 1, index = 1:2)  %>%
#   run_method(list(lassosplit,lassosplitadaptive, sailsplitweak, #lassoBTsplit, sailsplitweak, sailsplitadaptiveweak
#                   sailsplit, sailsplitadaptive),
#              parallel = list(socket_names = 20,
#                              libraries = c("LassoBacktracking", "glinternet","glmnet","splines",
#                                            "magrittr","sail","gamsel","SAM","HierBasis","simulator", "parallel")))

#
# draws(sim)@draws$r2.1$xtrain_lasso[, -which(colnames(draws(sim)@draws$r2.1$xtrain_lasso) %in% c("APOE_bin", "E"))]
# draws(sim)@draws$r2.1$group

sim <- new_simulation(name = "rda_ADNI_may_17_2018v2",
                      label = "rda_ADNI_may_17_2018v2",
                      dir = ".") %>%
  generate_model(make_ADNI_data_split, seed = 12345,
                 phenoVariable = "MMSCORE_bl",
                 exposure = "diag_3bl.x",n_train_test = 250) %>%
  simulate_from_model(nsim = 6, index = 1:35)  %>%
  run_method(list(lassosplit, sailsplit, GLinternetsplit,lassoBTsplit, Hiersplit),
             parallel = list(socket_names = 35,
             libraries = c("LassoBacktracking", "glinternet","glmnet","splines",
             "magrittr","sail","gamsel","SAM","HierBasis","simulator", "parallel")))
simulator::save_simulation(sim)

sim <- simulator::load_simulation("rda_ADNI_may_17_2018v2")

sim <- sim %>% evaluate(list(msevalid, nactive, r2))
sim %>% plot_eval(metric_name = "nactive")
sim %>% plot_eval(metric_name = "mse")
sim %>% plot_eval(metric_name = "r2")
DT <- as.data.frame(evals(sim)) %>% as.data.table()

draw_ind <- DT[Method=="sail"][which.min(mse)]$Draw
# DT[Method=="sail"][order(mse)]$mse
# draw_ind <- DT[Method=="sail"][which(mse >= mean(mse)-.0036 & mse <= mean(mse)+.002)]$Draw
dat <- draws(sim)@draws[[draw_ind]]
simulator::output(sim)[[2]]@out[[draw_ind]]
dat$xtrain

fit <- sail(x = dat$xtrain, y = dat$ytrain, e = dat$etrain,
            # basis = f.basis,
            expand = FALSE,
            center.x = F,
            center.e = T,
            group = dat$group,
            alpha = 0.1,
            maxit = 250,
            # thresh = 1e-02,
            strong = TRUE,
            verbose = 2)
plot(fit)
ytest_hat <- predict(fit, newx = dat$xtest, newe = dat$etest)
msetest <- colMeans((dat$ytest - ytest_hat)^2)
plot(log(fit$lambda), msetest)
lambda.min.index <- as.numeric(which.min(msetest))
lambda.min <- fit$lambda[which.min(msetest)]

yvalid_hat <- predict(fit, newx = dat$xvalid, newe = dat$evalid, s = lambda.min)
(msevalid <- mean((dat$yvalid - drop(yvalid_hat))^2))

(nzcoef <- predict(fit, s = lambda.min, type = "nonzero"))

X <- do.call(rbind, list(dat$xtrain_lasso, dat$xtest_lasso, dat$xvalid_lasso))
design <- do.call(rbind, list(dat$xtrain, dat$xtest, dat$xvalid))

# names taken from email from lai (search for ADNI in gmail)
plotMainADNI(fit, x = DT[as.numeric(as.character(rownames(design))), "X175"],
             xvar = paste0("bs(X175, 3)",1:3),design = design, s = lambda.min,
             ylab = "f(X175)", xlab = "Cuneus right")
plotMainADNI(fit, x = DT[as.numeric(as.character(rownames(design))), "X196"],
             xvar = paste0("bs(X196, 3)",1:3),design = design, s = lambda.min,
             ylab = "f(X196)", xlab = "Lateral occipitotemporal gyrus left")
plotMainADNI(fit, x = DT[as.numeric(as.character(rownames(design))), "X154"],
             xvar = paste0("bs(X154, 3)",1:3),design = design, s = lambda.min,
             ylab = "f(X154)", xlab =  "Middle occipital gyrus left")
plotMainADNI(fit, x = DT[as.numeric(as.character(rownames(design))), "X24"],
             xvar = paste0("bs(X24, 3)",1:3),design = design, s = lambda.min,
             ylab = "f(X24)", xlab =  "Middle occipital gyrus left")

par(mfrow=c(1,2), tcl=-0.5, family="serif",
    omi=c(0.2,0.2,0,0))
# mai = c(1,1,2,0))
# oma = c(1,1,2,1))
# plot_apoe_inter(cv_obj = cvfit, X = X, original_name = "X60", sail_name = "X19", xlab =  "supramarginal gyrus right", legend = T,
#                 apoe = FALSE)
# plot_apoe_inter(cv_obj = cvfit, X = X, original_name = "X60", sail_name = "X19", xlab =  "supramarginal gyrus right", legend = F,
#                 apoe = TRUE, ylab = "", main = "APOE e4 = 1")
plotInterADNI(fit, x = DT[as.numeric(as.character(rownames(design))), "X175"],
              xvar = paste0("bs(X175, 3)",1:3), design = design, s = lambda.min,
              e = E, apoe = FALSE, legend = T, legend.position = "bottomleft",
              ylab = "f(X175)", xlab = "Cuneus right", main = "APOE e4 = 0")
plotInterADNI(fit, x = DT[as.numeric(as.character(rownames(design))), "X175"],
              xvar = paste0("bs(X175, 3)",1:3), design = design, s = lambda.min,
              e = E, apoe = TRUE, legend = F,
              ylab = "f(X175)", xlab = "Cuneus right", main = "APOE e4 = 1")
dev.off()


par(mfrow=c(1,2), tcl=-0.5, family="serif",
    omi=c(0.2,0.2,0,0))
# mai = c(1,1,2,0))
# oma = c(1,1,2,1))
# plot_apoe_inter(cv_obj = cvfit, X = X, original_name = "X60", sail_name = "X19", xlab =  "supramarginal gyrus right", legend = T,
#                 apoe = FALSE)
# plot_apoe_inter(cv_obj = cvfit, X = X, original_name = "X60", sail_name = "X19", xlab =  "supramarginal gyrus right", legend = F,
#                 apoe = TRUE, ylab = "", main = "APOE e4 = 1")
plotInterADNI(fit, x = DT[as.numeric(as.character(rownames(design))), "X154"],
              xvar = paste0("bs(X154, 3)",1:3), design = design, s = lambda.min,
              e = E, apoe = FALSE, legend = F,
              ylab = "f(X154)", xlab = "Middle occipital gyrus left", main = "APOE e4 = 0")
plotInterADNI(fit, x = DT[as.numeric(as.character(rownames(design))), "X154"],
              xvar = paste0("bs(X154, 3)",1:3), design = design, s = lambda.min,
              e = E, apoe = TRUE, legend = T,
              ylab = "f(X154)", xlab = "Middle occipital gyrus left", main = "APOE e4 = 1")
dev.off()



DT[, mean(mse), by = Method]
DT[, median(mse), by = Method]
DT[, sd(mse), by = Method]
sim
simulator::output(sim)[[2]]@out$r1.1$active
simulator::output(sim)[[3]]@out$r1.1$active
simulator::output(sim)[[3]]@out$r2.1$active
# not-used ----------------------------------------------------------------


vars <- c("bs(X6, 3)1","bs(X6, 3)2","bs(X6, 3)3","bs(X154, 3)1","bs(X154, 3)2",
  "bs(X154, 3)3","bs(X175, 3)1","bs(X175, 3)2","bs(X175, 3)3","APOEbin",
  "bs(X6, 3)1:E","bs(X6, 3)2:E","bs(X6, 3)3:E","APOEbin:E","E")

(ints <- grep(":", vars, value=T))
(mains <- setdiff(vars, ints))
apoe <- if(any(mains=="APOEbin")) "APOEbin" else NULL
environ <- if(any(mains=="E")) "E" else NULL
mains_wo_apoe <- setdiff(mains, c("APOEbin","E"))

orig_mains <- unique(stringr::str_extract(mains_wo_apoe, "X\\d*"))
c(orig_mains, apoe, environ)



f.basis <- function(i) i
fit <- sail(x = X, y = Y, e = E, basis = f.basis, alpha = 0.5, thresh = 1e-02, strong = FALSE,
            verbose = 2)
plot(fit)

fit <- sail(
  x = model_mat,
  # x = as.matrix(DT[,brain_regions[1:10]]),
  y = Y, e = E,
  expand = FALSE,
  thresh = 5e-3,
  # center.e = FALSE,
  alpha = 0.2,
  # fdev = 1e-8,
  group = attr(model_mat, "assign"),
  verbose = 2)

plot(fit)
fit
as.matrix(coef(fit)[nonzero(coef(fit)),,])

as.matrix(fit$alpha[nonzero(fit$alpha),,])

help(sail)

fit
any(!fit$converged)
plot(fit)
fit$active
coef(fit)

tt <- KKT(b0 = fit$a0, betaE = fit$bE, beta = fit$beta, gamma = fit$gamma,
          alpha = fit$alpha, y = Y, phij = fit$Phi_j, xe_phij = fit$XE_Phi_j,
          e = E, df = fit$df,
          lambda = fit$lambda, lambda2 = fit$lambda2, group = fit$group,
          we = fit$we, wj = fit$wj, wje = fit$wje, thr = 1e-1, loss = "ls")

registerDoMC(cores = 8)
system.time(
  cvfit <- cv.sail(x = X, y = Y, e = E, df = 3, degree = 3, maxit = 1000,
                   alpha = 0.1,
                   nlambda = 100,
                   nfolds = 10,
                   verbose = TRUE,
                   thresh = 1e-4, parallel = TRUE)
)

plot(cvfit)
coef(cvfit, s = "lambda.min")[nonzero(coef(cvfit, s = "lambda.min")),,drop=F]
coef(cvfit, s = "lambda.1se")[nonzero(coef(cvfit, s = "lambda.1se")),,drop=F]
