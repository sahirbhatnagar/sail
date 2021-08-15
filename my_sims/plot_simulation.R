pacman::p_load(simulator) # this file was created under simulator version 0.2.1
pacman::p_load(splines)
pacman::p_load(magrittr)
pacman::p_load(methods)
pacman::p_load(glmnet)
pacman::p_load(LassoBacktracking)
pacman::p_load_current_gh('sahirbhatnagar/glinternet')
pacman::p_load(gbm)
# remotes::install_github('asadharis/HierBasis')
library("HierBasis")
devtools::load_all()
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
pacman::p_load(foreach)
pacman::p_load(doMC)
pacman::p_load(latex2exp)
pacman::p_load(UpSetR)
pacman::p_load(truncnorm)

nonlinear = TRUE

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

sim <- simulator::load_simulation("aug_12_2021")
sim

# used for fit objects
f.hat.fit <- function(object, xvar, x = truncnorm::rtruncnorm(1000, a = 0, b = 1)){
  # x is original X, not the expanded version
  browser()
  ind <- object$fit$group == which(object$fit$vnames == xvar)
  allCoefs <- object$beta
  # a0 <- allCoefs[1,]

  betas <- as.matrix(allCoefs[object$fit$main.effect.names[ind],,drop = FALSE])
  # design.mat <- object$fit$design_not_centered[,object$fit$main.effect.names[ind],drop = FALSE]
  # originalX <- x[,unique(object$fit$group[ind])]
  originalX <- x
  # all.equal(splines::bs(originalX, degree = 5, intercept = FALSE)[,1],
  # design.mat[,1])
  # fhat <- drop(a0 + design.mat %*% betas)
  design.mat <- scale(splines::bs(originalX, degree = 5, intercept = FALSE), center = TRUE, scale = FALSE)
  design.mat <- scale(splines::bs(originalX, degree = 5, intercept = FALSE), center = FALSE, scale = FALSE)
  fhat <- drop(design.mat %*% betas)
  plot(originalX[order(originalX)], fhat[order(originalX)], type = "l", ylim = c(-5,5))
  curve(f4, add = TRUE, lwd = 3, col = sail:::cbbPalette[7])
lines(originalX[order(originalX)], fhat[order(originalX)])
  return(list(X = originalX[order(originalX)], fX = fhat[order(originalX)]))
}

out <- simulator::output(sim, index = 2)@out

f.hat.fit(object = out[[1]], xvar = "X4")

out[[1]]
fit <- out[[1]]
dd <- draw@draws[[1]]


# ran this once takes a while
# fhats <- lapply(1:4, function(k) {
#
#   xvar <- paste0("X", k)
#
#   do.call(cbind,
#           lapply(1:35, function(j) {
#             out <- simulator::output(sim, index = j)@out
#             draw <- simulator::draws(sim, index = j)@draws
#
#             do.call(cbind,
#                     lapply(1:6, function(i) f.hat.fit(object = out[[i]],
#                                                       xvar = xvar,
#                                                       x = draw[[i]][["xtrain"]])[["fX"]])
#             )
#           }
#           )
#   )
# }
# )
#
# saveRDS(fhats, file = "my_sims/simulation_results/fhats.rds")
#
# Xs <- lapply(1:4, function(k) {
#
#   xvar <- paste0("X", k)
#
#   do.call(cbind,
#           lapply(1:35, function(j) {
#             out <- simulator::output(sim, index = j)@out
#             draw <- simulator::draws(sim, index = j)@draws
#
#             do.call(cbind,
#                     lapply(1:6, function(i) f.hat.fit(object = out[[i]],
#                                                       xvar = xvar,
#                                                       x = draw[[i]][["xtrain"]])[["X"]])
#             )
#           }
#           )
#   )
# }
# )
#
#
# saveRDS(Xs, file = "my_sims/simulation_results/Xs.rds")

fhats <- readRDS(here::here("my_sims/simulation_results/fhats.rds"))
Xs <- readRDS(here::here("my_sims/simulation_results/Xs.rds"))

par(mfrow = c(2,2))
for (xt in 1:4){
  xt = 2
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

  # fhats <- f.hat.fit(object = fit, xvar = xvar, x = dd[["xtrain"]])[["fX"]]
  # fhats_dat <- do.call(cbind,list(fhats))

  fhats_dat <- fhats[[xt]]
  Xs_dat <- Xs[[xt]]

  ylims <- range(fhats_dat)+3
  # ylims[1] <- if (nonlinear) ylims[1] else ylims[1]-0.5
  # ylims[2] <- if (nonlinear) ylims[2] + 4.5 else ylims[2] + 1

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
    lines(Xs_dat[,i], fhats_dat[,i]+3, col = sail:::cbbPalette[3], lwd = 1)
  }

  ff <- get(paste0("f",xt), mode="function")
  curve(ff, add = TRUE, lwd = 3, col = sail:::cbbPalette[7])
  if(xvar=="X4" & nonlinear) mtext(do.call(expression, main),side=3, line = c(1,-1) , cex = 1.3,
                                   family = "serif")
  if(xvar=="X1") legend("topleft", c("Truth", "Estimated"),
                        cex = 1.2, bty = "n", lwd = 2,
                        col = sail:::cbbPalette[c(7,3)])

}









fit$fit$design






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

# pdf(file="/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_git_v2/sail/my_sims/figures/sail_main_eff_paramIndex1_200sims.pdf",
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

