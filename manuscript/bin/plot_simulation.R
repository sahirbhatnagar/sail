# file used to create main effect plots for CSDA submission
# Updated October 10, 2021

## ---- packages-simulation-plots ----

pacman::p_load(simulator) # this file was created under simulator version 0.2.1
pacman::p_load(splines)
pacman::p_load(magrittr)
# pacman::p_load(methods)
# pacman::p_load(glmnet)
# pacman::p_load(LassoBacktracking)
# pacman::p_load_current_gh('sahirbhatnagar/glinternet')
# pacman::p_load(gbm)
# # remotes::install_github('asadharis/HierBasis')
# library("HierBasis")
# devtools::load_all()
# pacman::p_load(SAM)
# pacman::p_load(gamsel)
# pacman::p_load(cowplot)
pacman::p_load(dplyr)
pacman::p_load(tidyr)
pacman::p_load(data.table)
pacman::p_load(ggplot2)
pacman::p_load(latex2exp)
pacman::p_load(lemon)
pacman::p_load(here)
# pacman::p_load(sail)
pacman::p_load(truncnorm)
pacman::p_load(foreach)
pacman::p_load(doMC)
pacman::p_load(latex2exp)
pacman::p_load(UpSetR)
pacman::p_load(fst)


## ---- main-effects-simulation ----

nonlinear <- TRUE

if (nonlinear) {
  # non linear
  f1 <- function(x) 5 * x
  f2 <- function(x) 3 * (2 * x - 1)^2
  f3 <- function(x) 4 * sin(2 * pi * x) / (2 - sin(2 * pi * x))
  f4 <- function(x) 6 * (0.1 * sin(2 * pi * x) + 0.2 * cos(2 * pi * x) +
                           0.3 * sin(2 * pi * x)^2 + 0.4 * cos(2 * pi * x)^3 +
                           0.5 * sin(2 * pi * x)^3)
} else {
  # linear (Scenario 2)
  f1 <- function(x) -1.5 * (x-2)
  f2 <- function(x) 1 * (x+1)
  f3 <- function(x) 1.5 * x
  f4 <- function(x) -2 * x
}

sim <- simulator::load_simulation("aug_14_2021", dir = here::here())

# used for fit objects
f.hat.fit <- function(object, xvar, x = truncnorm::rtruncnorm(1000, a = 0, b = 1)){
  # x is original X, not the expanded version
  # browser()
  # xvar = "X2"
  ind <- object$fit$group == which(object$fit$vnames == xvar)
  allCoefs <- object$beta
  a0 <- allCoefs[1,]

  # object$fit$alpha
  # df <- object$originalX[,object$active[!grepl(":|E",object$active)], drop=F] %>% as.data.frame()
  # df <- object$originalX[,object$causal[!grepl(":|E",object$causal)], drop=F] %>% as.data.frame()
  # df$E <- object$originale
  # df$Y <- object$originaly
  # library(visreg)
  #
  # fit <- lm(Y ~ bs(X1, degree = 5) + bs(X2, degree = 5) + (bs(X3, degree = 5) + bs(X4, degree = 5)) * E, data = df )
  # summary(fit)
  # termplot(fit, terms = 1)
  # tt <- visreg(fit = fit, "X2", alpha = 1)
  # tt$fit
  # par(mfrow=c(1,2))
  # dev.off()
  # plot(tt$fit$X2, tt$fit$visregFit, type = "l", col = sail:::cbbPalette[3], lwd = 3, ylim = c(-1,5))
  # curve(f2, add = TRUE, lwd = 3, col = sail:::cbbPalette[7])
  # visreg(fit = fit, "X4", by = "E", breaks = c(E=-1,0,1), overlay = TRUE)
  # summary(fit)
  # termplot(fit)

  betas <- as.matrix(allCoefs[object$fit$main.effect.names[ind],,drop = FALSE])
  # design.mat <- object$fit$design_not_centered[,object$fit$main.effect.names[ind],drop = FALSE]
  # originalX <- x[,unique(object$fit$group[ind])]
  originalX <- object$originalX[,xvar]
  # originalX <- truncnorm::rtruncnorm(1000, a = -1, b = 1)
  # all.equal(splines::bs(originalX, degree = 5, intercept = FALSE)[,1],
  # design.mat[,1])
  # fhat <- drop(design.mat %*% betas)
  expansion <- object$fit$basis(originalX)
  design.mat <- scale(expansion, center = xvar %in% c("X2","X4"), scale = FALSE)
  # if (xvar %in% c("X2","X4")) avx <- attr(design.mat, "scaled:center")
  # non_zero <- apply(design.mat==0, 1, function(x) sum(x)!=5)
  # design.mat <- design.mat[non_zero,]
  # crossprod(avx,betas)
  # originalX <- originalX[non_zero]

  # design.mat <- scale(splines::bs(originalX, degree = 5, intercept = FALSE), center = FALSE, scale = FALSE)
  # design.mat <- object$fit$design[,object$fit$main.effect.names[ind]]
  # fhat <- drop(design.mat %*% betas) + ifelse(xvar %in% c("X2","X4"), crossprod(attr(design.mat, "scaled:center"),betas),0)
  # fhat <- drop(design.mat %*% betas) + ifelse(xvar %in% c("X2","X4"), a0,0)
  fhat <- drop(design.mat %*% betas)
  # plot(originalX[order(originalX)], fhat[order(originalX)], type = "l", ylim = c(-5,5))
  # curve(f2-a0, add = FALSE, lwd = 3, col = sail:::cbbPalette[7], ylim = c(-1,4))
  # lines(originalX[order(originalX)], fhat[order(originalX)])
  return(list(X = originalX[order(originalX)], fX = fhat[order(originalX)]))
}

out <- simulator::output(sim, subset = 1, methods = "sail")
# save(out, file = here::here("manuscript/results/simulation_sail_scen1.rdata"), compress = "xz")
pdf(file=here::here("manuscript/bin/figure/","main-effects-simulation-1.pdf"), width=11, height=8)
par(mfrow = c(2,2))
for (xt in 1:4){
  # xt = 2
  xvar = paste0("X",xt)

  if (nonlinear) {

    main <- switch(xvar,
                   X1 = latex2exp::TeX("$f(x_1) = 5x_1$"),
                   X2 = latex2exp::TeX("$f(x_2) = 3 (2x_2 - 1)^2$"),
                   X3 = latex2exp::TeX("$f(x_3) = \\frac{4\\sin(2\\pi x_3)}{2 - \\sin(2\\pi x_3)}$"),
                   X4 = as.list(expression(paste("f(x", phantom()[{
                     paste("4")
                   }], "",") = 6(0.1sin(2",pi,"x", phantom()[{
                     paste("4")
                   }], "",") + 0.2cos(2",pi,"x", phantom()[{
                     paste("4")
                   }], "",") + 0.3sin(2",pi,"x", phantom()[{
                     paste("4")
                   }], "",")","",
                   phantom()^{paste("2")},"+ "),
                   paste("     ","0.4cos(2",pi,"x", phantom()[{
                     paste("4")
                   }], "",")","",
                   phantom()^{paste("3")} ,"+ 0.5sin(2",pi,"x", phantom()[{
                     paste("4")
                   }], "",")","",
                   phantom()^{paste("3")},")"))))
  } else {
    main <- switch(xvar,
                   X1 = latex2exp::TeX("$f(x_1) = 0.5 x_1$"),
                   X2 = latex2exp::TeX("$f(x_2) = x_2$"),
                   X3 = latex2exp::TeX("$f(x_3) = 1.5 x_3$"),
                   X4 = latex2exp::TeX("$f(x_4) = 2 x_4$"))
  }



  # fhats <- lapply(dat_list, function(i) f.hat(object = i$sail.fit, xvar = xvar, s = i[[lambda.type]])[["fX"]])
  # iterate over each simulation
  fhats <- lapply(out@out, function(i) f.hat.fit(object = i, xvar = xvar)[["fX"]])
  fhats_dat <- do.call(cbind,fhats)
  Xs <- lapply(out@out, function(i) f.hat.fit(object = i, xvar = xvar)[["X"]])
  Xs_dat <- do.call(cbind,Xs)

  # fhats <- f.hat.fit(object = fit, xvar = xvar, x = dd[["xtrain"]])[["fX"]]
  # fhats_dat <- do.call(cbind,list(fhats))

  # fhats_dat <- fhats[[xt]]
  # Xs_dat <- Xs[[xt]]
  #

  ff <- get(paste0("f",xt), mode="function")

  ylims <- range(fhats_dat,ff(seq(0,1, length.out = 100)), na.rm = TRUE)
  ylims[1] <- if (nonlinear) ylims[1] else ylims[1]-0.5
  ylims[2] <- if (nonlinear) ylims[2] else ylims[2] + 1

  # if(xvar=="X4") ylims <- c(-4,6.5) else if (xvar=="X2") ylims <- c(0, 5)

  # Xs <- lapply(dat_list, function(i) f.hat(object = i$sail.fit, xvar = xvar, s = i[[lambda.type]])[["X"]])
  # Xs <- f.hat.fit(object = fit, xvar = xvar, x = dd[["xtrain"]])[["X"]]
  # Xs_dat <- do.call(cbind,list(Xs))
  xlims <- range(Xs_dat)


  plot.args <- list(x = Xs_dat[,1],
                    y = fhats_dat[,1],
                    ylim = c(ylims[1], ylims[2]),
                    xlim = xlims,
                    xlab = TeX(sprintf("$x_{%d}$", xt)),
                    ylab = TeX(sprintf("$f(x_{%d})$", xt)),
                    type = "n",
                    # xlim = rev(range(l)),
                    # las = 1,
                    cex.lab = 1.5,
                    cex.axis = 1.5,
                    cex = 1.5,
                    cex.main = 1.5,
                    bty = "n",
                    main = if(xvar=="X4" & nonlinear) "" else main,
                    # mai=c(1.02,2 , 0.82, 0.42),
                    # tcl = -0.5,
                    # oma = c(0,2,0,0),
                    family = "serif")
  # plot.new()
  par(mar = c(5.1,5, 4.1, 2.1))
  do.call("plot", plot.args)
  abline(h = 0, lwd = 1, col = "gray")
  for (i in seq_len(ncol(Xs_dat))) {
    lines(Xs_dat[,i], fhats_dat[,i], col = sail:::cbbPalette[3], lwd = 1)
  }


  curve(ff, add = TRUE, lwd = 3, col = sail:::cbbPalette[7])
  if(xvar=="X4" & nonlinear) mtext(do.call(expression, main),side=3, line = c(1,-1) , cex = 1.3,
                                   family = "serif")
  if(xvar=="X1") legend("topleft", c("Truth", "Estimated"),
                        cex = 1.2, bty = "n", lwd = 2,
                        col = sail:::cbbPalette[c(7,3)])

}

dev.off()



## ---- persp-effects ----

nonlinear <- TRUE

if (nonlinear) {
  f3.persp = function(X, E) {
    # E * as.vector(DT$betaE) + DT$f3.f(X) +
    E * f3(X)
  }

  f4.persp = function(X, E) {
    # E * as.vector(DT$betaE) + DT$f4.f(X) +
    E * f4(X)
  }
} else {

  f3.persp = function(X, E) {
    # E * as.vector(DT$betaE) + DT$f3.f(X) +
    E * f3(X)
  }

  f4.persp = function(X, E) {
    # E * as.vector(DT$betaE) + DT$f4.f(X) +
    -1.5 * E * f4(X)
  }
}

# i = 4
# xv <- paste0("X",4)
# sail:::plotInter(object = dat_list[[120]]$fit,
#           x = dat_list[[120]]$x,
#           xvar = xv, s = dat_list[[120]][[lambda.type]],
#           f.truth = f4.persp)

# used for cvfit
f.hat.inter <- function(object, xvar, s, f.truth){

  ind <- object$group == which(object$vnames == xvar)
  allCoefs <- coef(object, s = s)
  a0 <- allCoefs[1,]

  betas <- as.matrix(allCoefs[object$main.effect.names[ind],,drop = FALSE])
  alphas <- as.matrix(allCoefs[object$interaction.names[ind],,drop = FALSE])
  betaE <- as.matrix(allCoefs["E",,drop = FALSE])

  design.mat.main <- object$design[,object$main.effect.names[ind],drop = FALSE]
  design.mat.int <- object$design[,object$interaction.names[ind],drop = FALSE]
  originalE <- object$design[,"E",drop = FALSE] # this is the centered E
  originalX <- object$x[,unique(object$group[ind])]

  # fhat <- drop(originalE %*% betaE + design.mat.main %*% betas + design.mat.int %*% alphas)
  fhat <- drop(design.mat.int %*% alphas)
  ftruth <- drop(f.truth(originalX,originalE))

  l2norm(fhat - ftruth)
  # dist(rbind(fhat,ftruth), method = "minkowski", p=3)
  # return(list(fhat = fhat, ftruth = ftruth, l2 = l2norm(fhat - ftruth)))
}

# used for fit (objects of class sail.fit)
f.hat.inter.fit <- function(object, xvar, s, f.truth, x){
# browser()
  ind <- object$group == which(object$vnames == xvar)
  allCoefs <- coef(object, s = s)
  a0 <- allCoefs[1,]

  betas <- as.matrix(allCoefs[object$main.effect.names[ind],,drop = FALSE])
  alphas <- as.matrix(allCoefs[object$interaction.names[ind],,drop = FALSE])
  betaE <- as.matrix(allCoefs["E",,drop = FALSE])

  design.mat.main <- object$design[,object$main.effect.names[ind],drop = FALSE]
  design.mat.int <- object$design[,object$interaction.names[ind],drop = FALSE]
  originalE <- object$design[,"E",drop = FALSE] # this is the centered E
  originalX <- x[,unique(object$group[ind])]

  # fhat <- drop(originalE %*% betaE + design.mat.main %*% betas + design.mat.int %*% alphas)
  fhat <- drop(design.mat.int %*% alphas)
  ftruth <- drop(f.truth(originalX,originalE))

  sail:::l2norm(fhat - ftruth)
  # dist(rbind(fhat,ftruth), method = "minkowski", p=3)
  # return(list(fhat = fhat, ftruth = ftruth, l2 = l2norm(fhat - ftruth)))
}



## ---- X3-interactions -------------

# saveRDS(out@out, file = here::here("manuscript/results/simulation_output_sail_scenario1.rds"))
lambda.type <- "lambda.min"
# out@out$r1.1[[lambda.type]]
# out@out$r1.1[["originalX"]]
resX3 <- sapply(out@out, function(i) f.hat.inter.fit(object = i$fit,
                                                      xvar = "X3",
                                                      s = i[[lambda.type]],
                                                      f.truth = f3.persp,
                                                      x = i[["originalX"]]))
# resX3
indX3 <- match(quantile(resX3, probs = c(0.25,0.50,0.75), type = 1), resX3)
# resX3[indX3]
# plot(dat_list[[indX3[1]]]$fit)

i = 3
xv <- paste0("X",i)
pdf(file=here::here("manuscript/source/figure/","sail_intertruth_X3_paramIndex1_200sims.pdf"), width=11, height=8)
sail:::plotInter(object = out@out[[indX3[1]]]$fit,
                 xvar = xv,
                 x = out@out[[indX3[1]]][["originalX"]],
                 s = out@out[[indX3[1]]][[lambda.type]],
                 truthonly = TRUE,
                 f.truth = f3.persp)
dev.off()

pdf(file=here::here("manuscript/source/figure/","sail_inter25_X3_paramIndex1_200sims.pdf"),
    width=11,height=8)
sail:::plotInter(object = out@out[[indX3[1]]]$fit,
          xvar = xv,
          x = out@out[[indX3[1]]][["originalX"]],
          s = out@out[[indX3[1]]][[lambda.type]],
          title_z = "Estimated: 25th Percentile")
dev.off()

pdf(file=here::here("manuscript/source/figure/","sail_inter50_X3_paramIndex1_200sims.pdf"),
    width=11,height=8)
sail:::plotInter(object = out@out[[indX3[2]]]$fit,
                 xvar = xv,
                 x = out@out[[indX3[2]]][["originalX"]],
                 s = out@out[[indX3[2]]][[lambda.type]],
                 title_z = "Estimated: 50th Percentile")
dev.off()


pdf(file=here::here("manuscript/source/figure/","sail_inter75_X3_paramIndex1_200sims.pdf"),
    width=11,height=8)
# png(file="/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/mcgillsims/figures/gendata2_inter75_n200_p1000_SNR2_betaE1_df5_degree3_alpha2_1a_500sims.png",
#     width=11,height=8, units = "in", res = 125)
sail:::plotInter(object = out@out[[indX3[3]]]$fit,
                 xvar = xv,
                 x = out@out[[indX3[3]]][["originalX"]],
                 s = out@out[[indX3[3]]][[lambda.type]],
                 title_z = "Estimated: 75th Percentile")
dev.off()


## ---- X4-interactions ---------------------------------------------------------


resX4 <- sapply(out@out, function(i) f.hat.inter.fit(object = i$fit,
                                                      xvar = "X4",
                                                      s = i[[lambda.type]],
                                                      f.truth = f4.persp,
                                                      x = i[["originalX"]]))
# resX4
indX4 <- match(quantile(resX4, probs = c(0.25,0.5,0.75), type = 1), resX4)
# resX4[indX4]
# plot(dat_list[[indX4[1]]]$fit)
# coef(dat_list[[indX4[2]]]$fit, s = lambda.type)[nonzero(coef(dat_list[[indX4[2]]]$fit, s = lambda.type)),,drop=F]

# i = 4
# xv <- paste0("X",i)
# plotInter(object = dat_list[[indX4[1]]]$sail.fit,
#           xvar = xv,
#           s = dat_list[[indX4[1]]][[lambda.type]],
#           f.truth = f4.persp,
#           simulation = T,
#           npoints = 40)

# devtools::load_all("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/")

i = 4
xv <- paste0("X",i)
pdf(file=here::here("manuscript/source/figure/","sail_intertruth_X4_paramIndex1_200sims.pdf"),
    width=11, height=8)

sail:::plotInter(object = out@out[[indX4[1]]]$fit,
                 xvar = xv,
                 x = out@out[[indX4[1]]][["originalX"]],
                 s = out@out[[indX4[1]]][[lambda.type]],
                 truthonly = TRUE,
                 f.truth = f4.persp)
dev.off()

pdf(file=here::here("manuscript/source/figure/","sail_inter25_X4_paramIndex1_200sims.pdf"),
    width=11,height=8)
sail:::plotInter(object = out@out[[indX4[1]]]$fit,
                 xvar = xv,
                 x = out@out[[indX4[1]]][["originalX"]],
                 s = out@out[[indX4[1]]][[lambda.type]],
                 title_z = "Estimated: 25th Percentile")
dev.off()

pdf(file=here::here("manuscript/source/figure/","sail_inter50_X4_paramIndex1_200sims.pdf"),
    width=11,height=8)
sail:::plotInter(object = out@out[[indX4[2]]]$fit,
                 xvar = xv,
                 x = out@out[[indX4[2]]][["originalX"]],
                 s = out@out[[indX4[2]]][[lambda.type]],
                 title_z = "Estimated: 50th Percentile")
dev.off()


pdf(file=here::here("manuscript/source/figure/","sail_inter75_X4_paramIndex1_200sims.pdf"),
    width=11, height=8)
sail:::plotInter(object = out@out[[indX4[3]]]$fit,
                 xvar = xv,
                 x = out@out[[indX4[3]]][["originalX"]],
                 s = out@out[[indX4[3]]][[lambda.type]],
                 title_z = "Estimated: 75th Percentile")
dev.off()


## ---- upset-plot ---------------------



# out@out[[1]]$fit$lambda

upset_data <- lapply(1:5, function(mods) {
  simulator::output(sim, subset = mods, methods = "sail")
})

# had to do this manually since the for loop or lapply isnt working
mods <- 5
out <- upset_data[[mods]]

vnames <- c(out@out[[1]]$fit$vnames, "E", paste0(out@out[[1]]$fit$vnames, ":E"))
  show <- c("X1","X2","X3","X4","E","X3:E","X4:E") # this will be the same as truth for hierarchy=TRUE
  truth <- c("X1","X2","X3","X4","E","X3:E","X4:E") # this will be the same as truth for hierarchy=TRUE
  # truth <- c("X1","X2","E","X3:E","X4:E") # this is for hierarchy=FALSE
  negatives <- setdiff(vnames, truth)
  lambda.type <- "lambda.min"

  # dat_list[[1]]$fit$active[[which(dat_list[[1]][[lambda.type]]==dat_list[[1]]$fit$lambda )]]
  active_set <- do.call(rbind,lapply(out@out, function(i){
    1 * (show %in% i$fit$active[[which(i[[lambda.type]]==i$fit$lambda)]])
  })) %>% as.data.frame()

  colnames(active_set) <- show

  # n_active <- sapply(dat_list, function(i){
  #   length(i$fit$active[[which(i[[lambda.type]]==i$fit$lambda)]])
  # })
  #
  # active_set$`Number Active` <- n_active
  #
  # # r2 <- sapply(dat_list, function(i){
  # #   cor(predict(i, s = i[[lambda.type]]),
  # #       i$fit$y)^2})
  #
  # # active_set$`R-Squared` <- r2
  #
  # CVMSE <- sapply(dat_list, function(i){
  #   # i$cvm[which(i[[lambda.type]]==i$lambda)]
  #   i$msevalid
  # })
  #
  # active_set$`Test Set MSE` <- CVMSE

  # correct_sparsity <- sapply(dat_list, function(i) {
  #   correct_nonzeros <- sum(i$sail.fit$active[[which(i[[lambda.type]]==i$lambda)]] %in% truth)
  #   correct_zeros <- length(setdiff(negatives,i$sail.fit$active[[which(i[[lambda.type]]==i$lambda)]]))
  #   #correct sparsity
  #   (1 / length(vnames)) * (correct_nonzeros + correct_zeros)
  # })
  #
  # active_set$`Correct Sparsity` <- correct_sparsity


  # tpr <- sapply(dat_list, function(i) {
  #       length(intersect(i$sail.fit$active[[which(i[[lambda.type]]==i$lambda)]],
  #                                      truth))/length(truth)
  #                   })
  #
  # active_set$`True_Positive_Rate` <- tpr
  #
  # fpr <- sapply(dat_list, function(i) {
  #   active <- i$sail.fit$active[[which(i[[lambda.type]]==i$lambda)]]
  #   FPR <- sum(active %ni% truth) / length(negatives)
  #   FPR
  # })
  #
  # active_set$`False_Positive_Rate` <- fpr


  # Myfunc <- function(row, release, rating) {
  #   data <- (row["ReleaseDate"] %in% release) & (row["AvgRating"] > rating)
  # }

  # head(active_set)
  # png("mcgillsims/figures/upset_selection.png", width = 11, height = 8, units = "in", res = 100)

  pdf(here::here("manuscript/bin/figure",sprintf("upset_selection_sail_paramIndex%i.pdf",mods)), width = 10, height = 8)
  upset(active_set,
        sets = rev(show),
        # sets =
        point.size = 3.5, line.size = 2,
        mainbar.y.label = "Frequency", sets.x.label = "Selection Frequency",
        text.scale = 2,
        order.by = "freq",
        keep.order = TRUE#,
        # boxplot.summary = c("Number Active"),
        # boxplot.summary = c("Test Set MSE","Number Active"),
        # queries = list(
        # list(query = intersects, params = list(truth), active = T, color = "#D55E00"))#,
        #     attribute.plots = list(gridrows = 45,
        #                            plots = list(list(plot = scatter_plot,
        # x = "False_Positive_Rate", y = "True_Positive_Rate", queries = F)))
  )
  dev.off()

  # head(active_set)
  # upset(movies, main.bar.color = "black", queries = list(list(query = intersects,
  #     params = list("Drama"), color = "red", active = F), list(query = intersects,
  #     params = list("Action", "Drama"), active = T), list(query = intersects,
  #     params = list("Drama", "Comedy", "Action"), color = "orange", active = T)),
  #     attribute.plots = list(gridrows = 45, plots = list(list(plot = scatter_plot,
  #     x = "ReleaseDate", y = "AvgRating", queries = T), list(plot = scatter_plot,
  #     x = "AvgRating", y = "Watches", queries = F)), ncols = 2), query.legend = "bottom")



## ---- not-used-underneath ----



# sim <- sim %>% evaluate(list(msevalid, tpr, fpr, nactive, r2))

simulator::plot_eval(sim,"r2")
simulator::plot_eval(sim,"tpr")
simulator::plot_eval(sim,"fpr")
simulator::tabulate_eval(sim, "fpr")

devtools::load_all()
tt <- simulator::draws(sim)@draws
# tt$r1.1
# pt <- simulator::output(sim)@out
# pt$r1.1 %>% names

draw <- tt$r1.1


sail:::lspath

help(sail)
library(tictoc)
tic()
fit <- sail(x = draw[["xtrain"]], y = draw[["ytrain"]], e = draw[["etrain"]],
            basis = function(i) splines::bs(i, degree = 5), group.penalty = "gglasso", verbose = 2,
            nlambda = 25, center.x = TRUE, center.e = TRUE, dfmax = 12, maxit = 50)
toc()
fit
ytest_hat <- predict(fit, newx = draw[["xtest"]], newe = draw[["etest"]])

msetest <- colMeans((draw[["ytest"]] - ytest_hat)^2)
sd(msetest)

lambda.min.index <- as.numeric(which.min(msetest))
lambda.min <- fit$lambda[which.min(msetest)]
fit
coef(fit)
# lambda.min <- 0.0002229

plot(-log(fit$lambda), msetest, type = "l")
abline(v = -log(lambda.min))

yvalid_hat <- predict(fit, newx = draw[["xvalid"]], newe = draw[["evalid"]], s = lambda.min)
msevalid <- mean((draw[["yvalid"]] - drop(yvalid_hat))^2)

nzcoef <- predict(fit, s = lambda.min, type = "nonzero")



# pacman::p_load(sail)
# rm(list=ls())
# dev.off()
# devtools::load_all("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/")


# Gendata2-Hierarchy=TRUE ----------------------------------------------------------------

lambda.type <- "lambda.min"
nonlinear <- TRUE

# Gendata2-Hierarchy=TRUE main effects ----------------------------------------------------------------

# this is used for cvfit objects
# f.hat <- function(object, xvar, s){
#   # browser()
#   ind <- object$group == which(object$vnames == xvar)
#   allCoefs <- coef(object, s = s)
#   a0 <- allCoefs[1,]
#   betas <- as.matrix(allCoefs[object$main.effect.names[ind],,drop = FALSE])
#   design.mat <- object$design[,object$main.effect.names[ind],drop = FALSE]
#   originalX <- object$x[,unique(object$group[ind])]
#
#   # fhat <- drop(a0 + design.mat %*% betas)
#   fhat <- drop(design.mat %*% betas)
#   return(list(X = originalX[order(originalX)], fX = fhat[order(originalX)]))
# }





if (nonlinear) {
  # non linear
  f1 <- function(x) 5 * x
  f2 <- function(x) 3 * (2 * x - 1)^2
  f3 <- function(x) 4 * sin(2 * pi * x) / (2 - sin(2 * pi * x))
  f4 <- function(x) 6 * (0.1 * sin(2 * pi * x) + 0.2 * cos(2 * pi * x) +
                           0.3 * sin(2 * pi * x)^2 + 0.4 * cos(2 * pi * x)^3 +
                           0.5 * sin(2 * pi * x)^3)
} else {
  # linear (Scenario 2)
  f1 <- function(x) -1.5 * (x-2)
  f2 <- function(x) 1 * (x+1)
  f3 <- function(x) 1.5 * x
  f4 <- function(x) -2 * x
}

f.hat.fit <- function(object, xvar, x = truncnorm::rtruncnorm(1000, a = 0.0, b = 1), lambda.min){
  # x is original X, not the expanded version
  # browser()
  ind <- object$group == which(object$vnames == xvar)
  allCoefs <- coef(object, s = lambda.min)
  a0 <- allCoefs[1,]

  betas <- as.matrix(allCoefs[object$main.effect.names[ind],,drop = FALSE])
  # design.mat <- object$fit$design_not_centered[,object$fit$main.effect.names[ind],drop = FALSE]
  # originalX <- x[,unique(object$fit$group[ind])]
  originalX <- x
  # all.equal(splines::bs(originalX, degree = 5, intercept = FALSE)[,1],
  # design.mat[,1])
  # fhat <- drop(a0 + design.mat %*% betas)
  cent <- xvar %in% c("X2","X4")
  design.mat <- scale(splines::bs(originalX, degree = 5, intercept = FALSE), center = cent, scale = FALSE)
  # design.mat <- scale(splines::bs(originalX, degree = 5, intercept = FALSE), center = FALSE, scale = FALSE)
  fhat <- drop(design.mat %*% betas)
  # plot(originalX[order(originalX)], fhat[order(originalX)], type = "l", ylim = c(-5,5))
  # curve(f2, add = TRUE, lwd = 3, col = sail:::cbbPalette[7])
  # lines(originalX[order(originalX)], fhat[order(originalX)])
  return(list(X = originalX[order(originalX)], fX = fhat[order(originalX)]))
}

dev.off()

# pdf(file=here::here("manuscript/source/figure/","sail_main_eff_paramIndex1_200sims.pdf",
#     width=11,height=8)
par(mfrow = c(2,2))
newx <- truncnorm::rtruncnorm(1000, a = 0.0, b = 1)
for (xt in 1:4){
  # xt = 4
  xvar = paste0("X",xt)

  if (nonlinear) {

    main <- switch(xvar,
                   X1 = latex2exp::TeX("$f(x_1) = 5x_1$"),
                   X2 = latex2exp::TeX("$f(x_2) = 3 (2x_2 - 1)^2$"),
                   X3 = latex2exp::TeX("$f(x_3) = \\frac{4\\sin(2\\pi x_3)}{2 - \\sin(2\\pi x_3)}$"),
                   X4 = as.list(expression(paste("f(x", phantom()[{
                     paste("4")
                   }], "",") = 6(0.1sin(2",pi,"x", phantom()[{
                     paste("4")
                   }], "",") + 0.2cos(2",pi,"x", phantom()[{
                     paste("4")
                   }], "",") + 0.3sin(2",pi,"x", phantom()[{
                     paste("4")
                   }], "",")","",
                   phantom()^{paste("2")},"+ "),
                   paste("     ","0.4cos(2",pi,"x", phantom()[{
                     paste("4")
                   }], "",")","",
                   phantom()^{paste("3")} ,"+ 0.5sin(2",pi,"x", phantom()[{
                     paste("4")
                   }], "",")","",
                   phantom()^{paste("3")},")"))))
  } else {
    main <- switch(xvar,
                   X1 = latex2exp::TeX("$f(x_1) = 0.5 x_1$"),
                   X2 = latex2exp::TeX("$f(x_2) = x_2$"),
                   X3 = latex2exp::TeX("$f(x_3) = 1.5 x_3$"),
                   X4 = latex2exp::TeX("$f(x_4) = 2 x_4$"))
  }



  # fhats <- lapply(dat_list, function(i) f.hat(object = i$sail.fit, xvar = xvar, s = i[[lambda.type]])[["fX"]])
  # fhats <- lapply(dat_list, function(i) f.hat.fit(object = i$fit, xvar = xvar, s = i$lambda.min, x = i$x)[["fX"]])
  # fhats_dat <- do.call(cbind,fhats)
  # Xs <- lapply(dat_list, function(i) f.hat.fit(object = i$fit, xvar = xvar, s = i$lambda.min, x = i$x)[["X"]])
  # Xs_dat <- do.call(cbind,Xs)

  fhats <- f.hat.fit(object = fit, xvar = xvar, lambda.min = lambda.min, x = newx)[["fX"]]
  fhats_dat <- do.call(cbind,list(fhats))
  ylims <- range(fhats_dat)
  ylims[1] <- if (nonlinear) ylims[1] else ylims[1]-0.5
  ylims[2] <- if (nonlinear) ylims[2] + 4.5 else ylims[2] + 1

  # Xs <- lapply(dat_list, function(i) f.hat(object = i$sail.fit, xvar = xvar, s = i[[lambda.type]])[["X"]])
  Xs <- f.hat.fit(object = fit, xvar = xvar, lambda.min = lambda.min, x = newx)[["X"]]
  Xs_dat <- do.call(cbind,list(Xs))
  xlims <- range(Xs_dat)


  plot.args <- list(x = Xs,
                    y = fhats,
                    ylim = c(ylims[1], ylims[2]),
                    xlim = xlims,
                    xlab = TeX(sprintf("$x_{%d}$", xt)),
                    ylab = TeX(sprintf("$f(x_{%d})$", xt)),
                    type = "n",
                    # xlim = rev(range(l)),
                    # las = 1,
                    cex.lab = 1.5,
                    cex.axis = 1.5,
                    cex = 1.5,
                    cex.main = 1.5,
                    bty = "n",
                    main = if(xvar=="X4" & nonlinear) "" else main,
                    # mai=c(1.02,2 , 0.82, 0.42),
                    # tcl = -0.5,
                    # oma = c(0,2,0,0),
                    family = "serif")
  # plot.new()
  par(mar = c(5.1,5, 4.1, 2.1))
  do.call("plot", plot.args)
  abline(h = 0, lwd = 1, col = "gray")
  lines(Xs, fhats)
  # for (i in seq_len(ncol(Xs_dat))) {
  #   lines(Xs_dat[,i], fhats_dat[,i], col = sail:::cbbPalette[3], lwd = 1)
  # }

  ff <- get(paste0("f",xt), mode="function")
  curve(ff, add = TRUE, lwd = 3, col = sail:::cbbPalette[7])
  if(xvar=="X4" & nonlinear) mtext(do.call(expression, main),side=3, line = c(1,-1) , cex = 1.3,
                                   family = "serif")
  if(xvar=="X1") legend("topleft", c("Truth", "Estimated"),
                        cex = 1.2, bty = "n", lwd = 2,
                        col = sail:::cbbPalette[c(7,3)])

}
dev.off()

