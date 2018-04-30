pacman::p_load(splines)
pacman::p_load(magrittr)
pacman::p_load(foreach)
pacman::p_load(methods)
pacman::p_load(doMC)
pacman::p_load(latex2exp)
pacman::p_load(UpSetR)
pacman::p_load(sail)
# rm(list=ls())
# dev.off()
# devtools::load_all("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/")


# Gendata2-Hierarchy=TRUE ----------------------------------------------------------------

# files = list.files(path = '/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/mcgillsims/gendata2_p1000_1a',
                   # pattern = '*.rds', full.names = TRUE)
# files = list.files(path = '/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/mcgillsims/gendata2_p1000_weak_hier',
#                    pattern = '*.rds', full.names = TRUE)
# files <- list.files(path = '/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/mcgillsims/gendata2/',
#                    pattern = '*.rds', full.names = TRUE)
# files <- list.files(path = '/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/mcgillsims/gendata2_p1000_1c_2_3_6',
#                     pattern = 'degree1_alpha05_6_', full.names = TRUE)
files <- list.files(path = '/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/mcgillsims/thesis_p1000_1a',
                    pattern = '^fit', full.names = TRUE)
dat_list <- lapply(files, function (x) readRDS(x))
lambda.type <- "lambda.min"
nonlinear <- TRUE

# Gendata2-Hierarchy=TRUE main effects ----------------------------------------------------------------

# this is used for cvfit objects
f.hat <- function(object, xvar, s){
# browser()
  ind <- object$group == which(object$vnames == xvar)
  allCoefs <- coef(object, s = s)
  a0 <- allCoefs[1,]
  betas <- as.matrix(allCoefs[object$main.effect.names[ind],,drop = FALSE])
  design.mat <- object$design[,object$main.effect.names[ind],drop = FALSE]
  originalX <- object$x[,unique(object$group[ind])]

  # fhat <- drop(a0 + design.mat %*% betas)
  fhat <- drop(design.mat %*% betas)
  return(list(X = originalX[order(originalX)], fX = fhat[order(originalX)]))
}

# used for fit objects
f.hat.fit <- function(object, xvar, s, x){
  # browser()
  ind <- object$group == which(object$vnames == xvar)
  allCoefs <- coef(object, s = s)
  a0 <- allCoefs[1,]
  betas <- as.matrix(allCoefs[object$main.effect.names[ind],,drop = FALSE])
  design.mat <- object$design[,object$main.effect.names[ind],drop = FALSE]
  originalX <- x[,unique(object$group[ind])]

  # fhat <- drop(a0 + design.mat %*% betas)
  fhat <- drop(design.mat %*% betas)
  return(list(X = originalX[order(originalX)], fX = fhat[order(originalX)]))
}

# i = dat_list[[1]]
# plot(i$fit)
# f.hat.fit(object = i$fit, xvar = "X1", s = i$lambda.min, x = i$x)[["fX"]]


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

dev.off()

pdf(file="/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/mcgillsims/figures/gendata2_main_n200_p1000_SNR2_betaE1_df5_degree3_alpha2_1a_500sims.pdf",
    width=11,height=8)
par(mfrow = c(2,2))
for (xt in 1:4){
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
  fhats <- lapply(dat_list, function(i) f.hat.fit(object = i$fit, xvar = xvar, s = i$lambda.min, x = i$x)[["fX"]])
  fhats_dat <- do.call(cbind,fhats)
  ylims <- range(fhats_dat)
  ylims[1] <- if (nonlinear) ylims[1] else ylims[1]-0.5
  ylims[2] <- if (nonlinear) ylims[2] + 2.5 else ylims[2] + 1

  # Xs <- lapply(dat_list, function(i) f.hat(object = i$sail.fit, xvar = xvar, s = i[[lambda.type]])[["X"]])
  Xs <- lapply(dat_list, function(i) f.hat.fit(object = i$fit, xvar = xvar, s = i$lambda.min, x = i$x)[["X"]])
  Xs_dat <- do.call(cbind,Xs)
  xlims <- range(Xs_dat)


  plot.args <- list(x = Xs[[1]],
                    y = fhats[[1]],
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

  ff <- get(paste0("f",xt), mode="function")
  curve(ff, add = TRUE, lwd = 3, col = sail:::cbbPalette[7])
  if(xvar=="X4" & nonlinear) mtext(do.call(expression, main),side=3, line = c(1,-1) , cex = 1.3,
                       family = "serif")

}
dev.off()


# Gendata2-Hierarchy=TRUE Interactions ---------------------------------------------------

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

i = 4
xv <- paste0("X",4)
plotInter(object = dat_list[[120]]$sail.fit, xvar = xv, s = dat_list[[120]][[lambda.type]],
          f.truth = f4.persp,
          simulation = T,
          npoints = 40)


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

resX3 <- sapply(dat_list, function(i) f.hat.inter(object = i$sail.fit,
                                                  xvar = "X3",
                                                  s = i[[lambda.type]],
                                                  f.truth = f3.persp))
resX3
indX3 <- match(quantile(resX3, probs = c(0.25,0.50,0.75), type = 1), resX3)
resX3[indX3]
plot(dat_list[[indX3[1]]])
dev.off()


devtools::load_all("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/")

i = 3
xv <- paste0("X",i)
pdf(file="/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/mcgillsims/figures/gendata2_inter_truth_n200_p1000_SNR2_betaE1_df5_degree3_alpha2_1a_500sims.pdf",
    width=11,height=8)
plotInter(object = dat_list[[indX3[1]]]$sail.fit, xvar = xv, s = dat_list[[indX3[1]]][[lambda.type]],
          truthonly = TRUE,
          f.truth = f3.persp,
          simulation = T,
          npoints = 30)
dev.off()

pdf(file="/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/mcgillsims/figures/gendata2_inter25_n200_p1000_SNR2_betaE1_df5_degree3_alpha2_1a_500sims.pdf",
    width=11,height=8)
plotInter(object = dat_list[[indX3[1]]]$sail.fit, xvar = xv, s = dat_list[[indX3[1]]][[lambda.type]],
          simulation = T,
          npoints = 30, title_z = "Estimated: 25th Percentile")
dev.off()

pdf(file="/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/mcgillsims/figures/gendata2_inter50_n200_p1000_SNR2_betaE1_df5_degree3_alpha2_1a_500sims.pdf",
    width=11,height=8)
plotInter(object = dat_list[[indX3[2]]]$sail.fit, xvar = xv, s = dat_list[[indX3[2]]][[lambda.type]],
          simulation = T,
          npoints = 30, title_z = "Estimated: 50th Percentile")
dev.off()


pdf(file="/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/mcgillsims/figures/gendata2_inter75_n200_p1000_SNR2_betaE1_df5_degree3_alpha2_1a_500sims.pdf",
    width=11,height=8)
# png(file="/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/mcgillsims/figures/gendata2_inter75_n200_p1000_SNR2_betaE1_df5_degree3_alpha2_1a_500sims.png",
#     width=11,height=8, units = "in", res = 125)
plotInter(object = dat_list[[indX3[3]]]$sail.fit, xvar = xv, s = dat_list[[indX3[3]]][[lambda.type]],
          simulation = T,
          npoints = 30, title_z = "Estimated: 75th Percentile")
dev.off()


resX4 <- sapply(dat_list, function(i) f.hat.inter(object = i$sail.fit,
                                                  xvar = "X4",
                                                  s = i[[lambda.type]],
                                                  f.truth = f4.persp))
resX4
indX4 <- match(quantile(resX4, probs = c(0.25,0.5,0.75), type = 1), resX4)
resX4[indX4]
plot(dat_list[[indX4[1]]])
coef(dat_list[[indX4[2]]], s = lambda.type)[nonzero(coef(dat_list[[indX4[2]]], s = lambda.type)),,drop=F]
dev.off()
i = 4
xv <- paste0("X",i)
plotInter(object = dat_list[[indX4[1]]]$sail.fit, xvar = xv, s = dat_list[[indX4[1]]][[lambda.type]],
          f.truth = f4.persp,
          simulation = T,
          npoints = 40)

devtools::load_all("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/")

i = 4
xv <- paste0("X",i)
pdf(file="/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/mcgillsims/figures/gendata2_inter_truth_X4_n200_p1000_SNR2_betaE1_df5_degree3_alpha2_1a_500sims.pdf",
    width=11,height=8)
plotInter(object = dat_list[[indX4[1]]]$sail.fit, xvar = xv, s = dat_list[[indX4[1]]][[lambda.type]],
          truthonly = TRUE,
          f.truth = f4.persp,
          simulation = T,
          npoints = 30)
dev.off()

pdf(file="/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/mcgillsims/figures/gendata2_inter25_X4_n200_p1000_SNR2_betaE1_df5_degree3_alpha2_1a_500sims.pdf",
    width=11,height=8)
plotInter(object = dat_list[[indX4[1]]]$sail.fit, xvar = xv, s = dat_list[[indX4[1]]][[lambda.type]],
          simulation = T,
          npoints = 30, title_z = "Estimated: 25th Percentile")
dev.off()

pdf(file="/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/mcgillsims/figures/gendata2_inter50_X4_n200_p1000_SNR2_betaE1_df5_degree3_alpha2_1a_500sims.pdf",
    width=11,height=8)
plotInter(object = dat_list[[indX4[2]]]$sail.fit, xvar = xv, s = dat_list[[indX4[2]]][[lambda.type]],
          simulation = T,
          npoints = 30, title_z = "Estimated: 50th Percentile")
dev.off()


pdf(file="/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/mcgillsims/figures/gendata2_inter75_X4_n200_p1000_SNR2_betaE1_df5_degree3_alpha2_1a_500sims.pdf",
    width=11,height=8)
plotInter(object = dat_list[[indX4[3]]]$sail.fit, xvar = xv, s = dat_list[[indX4[3]]][[lambda.type]],
          simulation = T,
          npoints = 30, title_z = "Estimated: 75th Percentile")
dev.off()




# Gendata2-Hierarchy=TRUE Selection Rates ---------------------------------

vnames <- c(dat_list[[1]]$sail.fit$vnames, "E", paste0(dat_list[[1]]$sail.fit$vnames, ":E"))
show <- c("X1","X2","X3","X4","E","X3:E","X4:E") # this will be the same as truth for hierarchy=TRUE
truth <- c("X1","X2","X3","X4","E","X3:E","X4:E") # this will be the same as truth for hierarchy=TRUE
# truth <- c("X1","X2","E","X3:E","X4:E") # this is for hierarchy=FALSE
negatives <- setdiff(vnames, truth)

active_set <- do.call(rbind,lapply(dat_list, function(i){
  1 * (show %in% i$sail.fit$active[[which(i[[lambda.type]]==i$lambda)]])
})) %>% as.data.frame()

colnames(active_set) <- show

n_active <- sapply(dat_list, function(i){
  length(i$sail.fit$active[[which(i[[lambda.type]]==i$lambda)]])
})

active_set$`Number Active` <- n_active

r2 <- sapply(dat_list, function(i){
  cor(predict(i, s = i[[lambda.type]]),
      i$sail.fit$y)^2})

active_set$`R-Squared` <- r2

CVMSE <- sapply(dat_list, function(i){
  i$cvm[which(i[[lambda.type]]==i$lambda)]
  })

active_set$`10 Fold-CV MSE` <- CVMSE

correct_sparsity <- sapply(dat_list, function(i) {
  correct_nonzeros <- sum(i$sail.fit$active[[which(i[[lambda.type]]==i$lambda)]] %in% truth)
  correct_zeros <- length(setdiff(negatives,i$sail.fit$active[[which(i[[lambda.type]]==i$lambda)]]))
  #correct sparsity
  (1 / length(vnames)) * (correct_nonzeros + correct_zeros)
})

active_set$`Correct Sparsity` <- correct_sparsity


tpr <- sapply(dat_list, function(i) {
      length(intersect(i$sail.fit$active[[which(i[[lambda.type]]==i$lambda)]],
                                     truth))/length(truth)
                  })

active_set$`True_Positive_Rate` <- tpr

fpr <- sapply(dat_list, function(i) {
  active <- i$sail.fit$active[[which(i[[lambda.type]]==i$lambda)]]
  FPR <- sum(active %ni% truth) / length(negatives)
  FPR
})

active_set$`False_Positive_Rate` <- fpr


# Myfunc <- function(row, release, rating) {
#   data <- (row["ReleaseDate"] %in% release) & (row["AvgRating"] > rating)
# }

head(active_set)
png("mcgillsims/figures/upset_selection.png", width = 11, height = 8, units = "in", res = 100)

pdf("mcgillsims/figures/upset_selection_truth_linear_sail_degree1.pdf", width = 10, height = 8)
upset(active_set,
      sets = rev(show),
      # sets =
      point.size = 3.5, line.size = 2,
      mainbar.y.label = "Frequency", sets.x.label = "Selection Frequency",
      text.scale = 2, order.by = "freq", keep.order = TRUE,
      boxplot.summary = c("Number Active"),
      # boxplot.summary = c("Correct Sparsity","Number Active")#,
      queries = list(
      list(query = intersects, params = list(truth), active = T, color = "#D55E00"))#,
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





# Results from simulator package -----------------------------------------
# comparing methods
rm(list=ls())
pacman::p_load(cowplot)
pacman::p_load(tidyverse)
pacman::p_load(data.table)

df <- do.call(rbind, lapply(list.files(path = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/mcgillsims/sim_results/",
                                       pattern = "*.rds", full.names = TRUE), readRDS))

head(df)

df <- df %>% separate(Model, into = c("simnames","betaE","corr","lambda.type","n","p","parameterIndex","SNR_2"),
                sep = "/")

df2 <- readRDS("/mnt/GREENWOOD_SCRATCH/sahir.bhatnagar/sail_simulations_other/res_lassoBT_glinternet42e1653665d9.rds") %>%
  separate(Model, into = c("simnames","betaE","corr","lambda.type","n","p","parameterIndex","SNR_2"),
           sep = "/")

colnames(df) %in% colnames(df2)

DT <- data.table::rbindlist(list(df, df2), use.names = T, fill = T)

# DT[Method=="GLinternet", cvmse:= mse]
# DT[Method=="GLinternet", cvmse:= ifelse(cvmse<mse, mse, cvmse), by = parameterIndex]

trop = RSkittleBrewer::RSkittleBrewer("trop")

DT[parameterIndex=="parameterIndex_1", table(Method)]


cbbPalette <- function() {
  c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
}
trop <- RSkittleBrewer::RSkittleBrewer("trop")
gg_sy <- theme(legend.position = "bottom", axis.text = element_text(size = 20),
               axis.title = element_text(size = 20), legend.text = element_text(size = 20),
               legend.title = element_text(size = 20),plot.title = element_text(size = 20) )




# Scenario 1 --------------------------------------------------------------

(p1_cvmse <- ggplot(DT[parameterIndex=="parameterIndex_1"], aes(Method, cvmse, fill = Method)) +
    ggplot2::geom_boxplot() +
    scale_fill_manual(values=cbbPalette()[c(2,3,4,7)]) +
    ggplot2::labs(y = "Minimum 10-Fold CV MSE", title = "") + xlab("") + gg_sy +
  theme(legend.position = "none"))

save_plot("mcgillsims/figures/p1_cvmse.pdf", p1_cvmse,
          base_height = 7, base_width = 9)

(p1_time <- ggplot(DT[parameterIndex=="parameterIndex_1"], aes(Method, time/60, fill = Method)) +
    ggplot2::geom_boxplot() +
    scale_fill_manual(values=cbbPalette()[c(2,3,4,7)]) +
    ylim(c(0,6)) +
    # scale_y_continuous(breaks = seq(0, 6, 1)) +
    ggplot2::labs(y = "Run Time (min) for Entire Solution Path and 10-Fold CV", title = "") + xlab("") + gg_sy +
    theme(legend.position = "none"))

save_plot("mcgillsims/figures/p1_time.pdf", p1_time,
          base_height = 8, base_width = 9)

(p1_fprtpr <- ggplot(DT[parameterIndex=="parameterIndex_1"], aes(fpr, tpr,group=Method, color = Method)) +
    ggplot2::geom_point(size=2) +
    # geom_jitter(size=3.5)+
    scale_color_manual(values=cbbPalette()[c(2,3,4,7)]) + gg_sy + coord_equal(ratio=.07)+
    ggplot2::labs(y = "True Positive Rate", title = "") +
    xlab("False Positive Rate"))
save_plot("mcgillsims/figures/p1_fprtpr.pdf", p1_fprtpr,
          base_height = 7, base_width = 9)


(p1_gl <- ggplot(DT[Method=="GLinternet" & parameterIndex == "parameterIndex_1"], aes(x=mse, y=cvmse, colour = mse > cvmse)) +
    geom_point(size=2.8) + gg_sy + geom_abline(intercept = 0, slope = 1) +
    ylim(c(0,43)) + xlim(c(0,45)) +
    coord_equal(ratio=1) +
    xlab("") +
    ylab("10-Fold CV MSE") + labs(title = "GLinternet") +
    # geom_point(aes(x = C1, y = C2, colour = C1 >0)) +
    scale_colour_manual(name = 'PC1 > 0', values = setNames(c(cbbPalette()[c(7)],"black"),c(T, F))) +
    theme(legend.position = "none")) #+ annotate("text", label = sprintf("%0.2f %% of simulations "), x = 30, y = 5, size = 8, colour = "black"))
save_plot("mcgillsims/figures/cvmse_mse_GLinternet.pdf", p1_gl,
          base_height = 7, base_width = 7)


(p1_gl <- ggplot(DT[Method=="lassoBT" & parameterIndex == "parameterIndex_1"], aes(x=mse, y=cvmse, colour = mse > cvmse)) +
  geom_point(size=2.8) + gg_sy + geom_abline(intercept = 0, slope = 1) +
    ylim(c(0,43)) + xlim(c(0,40)) +
  coord_equal(ratio=1) +
    xlab("Training Set MSE") + ylab("10-Fold CV MSE") + labs(title = "lassoBT") +
  # geom_point(aes(x = C1, y = C2, colour = C1 >0)) +
  scale_colour_manual(name = 'PC1 > 0', values = setNames(c(cbbPalette()[c(7)],"black"),c(T, F))) +
  theme(legend.position = "none")) #+ annotate("text", label = sprintf("%0.2f %% of simulations "), x = 30, y = 5, size = 8, colour = "black"))
save_plot("mcgillsims/figures/cvmse_mse_lassoBT.pdf", p1_gl,
          base_height = 7, base_width = 7)


(p1_gl <- ggplot(DT[Method=="sail" & parameterIndex == "parameterIndex_1"], aes(x=mse, y=cvmse, colour = mse > cvmse)) +
    geom_point(size=2.8) + gg_sy + geom_abline(intercept = 0, slope = 1) +
    ylim(c(0,43)) + xlim(c(0,40)) +
    coord_equal(ratio=1) +
    xlab("Training Set MSE") + ylab("10-Fold CV MSE") + labs(title = "sail") +
    # geom_point(aes(x = C1, y = C2, colour = C1 >0)) +
    scale_colour_manual(name = 'PC1 > 0', values = setNames(c(cbbPalette()[c(7)],"black"),c(T, F))) +
    theme(legend.position = "none")) #+ annotate("text", label = sprintf("%0.2f %% of simulations "), x = 30, y = 5, size = 8, colour = "black"))
save_plot("mcgillsims/figures/cvmse_mse_sail.pdf", p1_gl,
          base_height = 7, base_width = 7)

(p1_gl <- ggplot(DT[Method=="lasso" & parameterIndex == "parameterIndex_1"], aes(x=mse, y=cvmse, colour = mse > cvmse)) +
    geom_point(size=2.8) + gg_sy + geom_abline(intercept = 0, slope = 1) +
    ylim(c(0,43)) + xlim(c(0,40)) +
    coord_equal(ratio=1) +
    xlab("") +
    ylab("10-Fold CV MSE") + labs(title = "lasso") +
    # geom_point(aes(x = C1, y = C2, colour = C1 >0)) +
    scale_colour_manual(name = 'PC1 > 0', values = setNames(c(cbbPalette()[c(7)],"black"),c(T, F))) +
    theme(legend.position = "none")) #+ annotate("text", label = sprintf("%0.2f %% of simulations "), x = 30, y = 5, size = 8, colour = "black"))
save_plot("mcgillsims/figures/cvmse_mse_lasso.pdf", p1_gl,
          base_height = 7, base_width = 7)








DT[Method=="GLinternet" & parameterIndex == "parameterIndex_1"][mse>cvmse]

abline(a=0, b=1)

gl <- melt(DT[Method=="GLinternet" & parameterIndex == "parameterIndex_1", c("Draw","mse","cvmse")], id.vars = "Draw")
gl[variable=="mse", variable:="MSE"]
gl[variable=="cvmse", variable:="10-Fold CV MSE"]

(p1_gl <- ggplot(gl, aes(variable, value, fill = variable)) +
    ggplot2::geom_boxplot() +
    scale_fill_manual(values=cbbPalette()[c(8,4)]) +
    ggplot2::labs(y = "", title = "GLinternet: MSE vs. 10-Fold CV MSE") + xlab("") + gg_sy +
  theme(legend.position = "none"))
save_plot("mcgillsims/figures/glint_cv_mse_scenario1.pdf", p1_gl,
          base_height = 6, base_width = 8)




# Scenario 5 --------------------------------------------------------------

(p5_cvmse <- ggplot(DT[parameterIndex=="parameterIndex_5"], aes(Method, cvmse, fill = Method)) +
   ggplot2::geom_boxplot() +
   scale_fill_manual(values=cbbPalette()[c(2,3,4,7)]) +
   ggplot2::labs(y = "Minimum 10-Fold CV MSE", title = "") + xlab("") + gg_sy +
   theme(legend.position = "none"))

save_plot("mcgillsims/figures/p5_cvmse.pdf", p5_cvmse,
          base_height = 7, base_width = 9)

(p5_fprtpr <- ggplot(DT[parameterIndex=="parameterIndex_5"], aes(fpr, tpr,group=Method, color = Method)) +
    # ggplot2::geom_point(size=2) +
    geom_jitter(size=3.5)+
    scale_color_manual(values=cbbPalette()[c(2,3,4,7)]) + gg_sy + coord_equal(ratio=.07)+
    ggplot2::labs(y = "True Positive Rate", title = "") +
    xlab("False Positive Rate"))
save_plot("mcgillsims/figures/p5_fprtpr.pdf", p5_fprtpr,
          base_height = 7, base_width = 9)



gl <- melt(DT[Method=="GLinternet" & parameterIndex == "parameterIndex_5", c("Draw","mse","cvmse")], id.vars = "Draw")
gl[variable=="mse", variable:="MSE"]
gl[variable=="cvmse", variable:="10-Fold CV MSE"]

(p1_gl <- ggplot(gl, aes(variable, value, fill = variable)) +
    ggplot2::geom_boxplot() +
    scale_fill_manual(values=cbbPalette()[c(8,4)]) +
    ggplot2::labs(y = "", title = "GLinternet: MSE vs. 10-Fold CV MSE") + xlab("") + gg_sy +
    theme(legend.position = "none"))
save_plot("mcgillsims/figures/glint_cv_mse_scenario1.pdf", p1_gl,
          base_height = 6, base_width = 8)




# Scenario 4 --------------------------------------------------------------

(p4_cvmse <- ggplot(DT[parameterIndex=="parameterIndex_4"], aes(Method, cvmse, fill = Method)) +
   ggplot2::geom_boxplot() +
   scale_fill_manual(values=cbbPalette()[c(2,3,4,7)]) +
   ggplot2::labs(y = "Minimum 10-Fold CV MSE", title = "") + xlab("") + gg_sy +
   theme(legend.position = "none"))

save_plot("mcgillsims/figures/p4_cvmse.pdf", p4_cvmse,
          base_height = 7, base_width = 9)

(p4_fprtpr <- ggplot(DT[parameterIndex=="parameterIndex_4"], aes(fpr, tpr,group=Method, color = Method)) +
    # ggplot2::geom_point(size=2) +
    geom_jitter(size=3.5)+
    scale_color_manual(values=cbbPalette()[c(2,3,4,7)]) + gg_sy + coord_equal(ratio=.07)+
    ggplot2::labs(y = "True Positive Rate", title = "") +
    xlab("False Positive Rate"))
save_plot("mcgillsims/figures/p4_fprtpr.pdf", p5_fprtpr,
          base_height = 7, base_width = 9)



gl <- melt(DT[Method=="GLinternet" & parameterIndex == "parameterIndex_5", c("Draw","mse","cvmse")], id.vars = "Draw")
gl[variable=="mse", variable:="MSE"]
gl[variable=="cvmse", variable:="10-Fold CV MSE"]

(p1_gl <- ggplot(gl, aes(variable, value, fill = variable)) +
    ggplot2::geom_boxplot() +
    scale_fill_manual(values=cbbPalette()[c(8,4)]) +
    ggplot2::labs(y = "", title = "GLinternet: MSE vs. 10-Fold CV MSE") + xlab("") + gg_sy +
    theme(legend.position = "none"))
save_plot("mcgillsims/figures/glint_cv_mse_scenario1.pdf", p1_gl,
          base_height = 6, base_width = 8)





