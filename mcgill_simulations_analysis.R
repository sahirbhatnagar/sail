pacman::p_load(splines)
pacman::p_load(magrittr)
pacman::p_load(foreach)
pacman::p_load(methods)
pacman::p_load(doMC)

# rm(list=ls())
# dev.off()
devtools::load_all()


# Gendata2 ----------------------------------------------------------------

files = list.files(path = '/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/mcgillsims',
                   pattern = '*.rds', full.names = TRUE)
dat_list = lapply(files, function (x) readRDS(x))
f.hat <- function(object, xvar, s){

  ind <- object$group == which(object$vnames == xvar)
  allCoefs <- coef(object, s = s)
  a0 <- allCoefs[1,]
  betas <- as.matrix(allCoefs[object$main.effect.names[ind],,drop = FALSE])
  design.mat <- object$design[,object$main.effect.names[ind],drop = FALSE]
  originalX <- object$x[,unique(object$group[ind])]

  fhat <- drop(a0 + design.mat %*% betas)
  # fhat <- drop(design.mat %*% betas)
  return(list(X = originalX[order(originalX)], fX = fhat[order(originalX)]))
}
f1 <- function(t) 5 * t
f2 <- function(t) 4.5 * (2 * t - 1) ^ 2
f3 <- function(t) 4 * sin(2 * pi * t) / (2 - sin(2 * pi * t))
f4 <- function(t) 6 * (0.1 * sin(2 * pi * t) + 0.2 * cos(2 * pi * t) +
                         0.3 * sin(2 * pi * t) ^ 2 + 0.4 * cos(2 * pi * t) ^ 3 +
                         0.5 * sin(2 * pi * t) ^ 3)


xvar = "X4"

fhats <- lapply(dat_list, function(i) f.hat(object = i$sail.fit, xvar = xvar, s = i$lambda.1se)[["fX"]])
fhats_dat <- do.call(cbind,fhats)
ylims <- range(fhats_dat)

Xs <- lapply(dat_list, function(i) f.hat(object = i$sail.fit, xvar = xvar, s = i$lambda.1se)[["X"]])
Xs_dat <- do.call(cbind,Xs)
xlims <- range(Xs_dat)


plot.args <- list(x = Xs[[1]],
                  y = fhats[[1]],
                  ylim = c(ylims[1], ylims[2]+2),
                  xlim = xlims,
                  xlab = xvar,
                  ylab = sprintf("f(%s)",xvar),
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
                  family = "serif")
do.call("plot", plot.args)
abline(h = 0, lwd = 1, col = "gray")
for (i in seq_len(ncol(Xs_dat))) {
  lines(Xs_dat[,i], fhats_dat[,i], col = cbbPalette()[3], lwd = 1)
}
curve(f4, add = TRUE)



# Gendata ----------------------------------------------------------------

files = list.files(
  path = '/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/mcgillsims/gendata/',
  pattern = '*.rds', full.names = TRUE)
dat_list = lapply(files, function (x) readRDS(x))
f.hat <- function(object, xvar, s){

  ind <- object$group == which(object$vnames == xvar)
  allCoefs <- coef(object, s = s)
  a0 <- allCoefs[1,]
  betas <- as.matrix(allCoefs[object$main.effect.names[ind],,drop = FALSE])
  design.mat <- object$design[,object$main.effect.names[ind],drop = FALSE]
  originalX <- object$x[,unique(object$group[ind])]

  fhat <- drop(a0 + design.mat %*% betas)
  # fhat <- drop(design.mat %*% betas)
  return(list(X = originalX[order(originalX)], fX = fhat[order(originalX)]))
}
f1 <- function(t) 5 * t
f2 <- function(t) 4.5 * (2 * t - 1) ^ 2
f3 <- function(t) 4 * sin(2 * pi * t) / (2 - sin(2 * pi * t))
f4 <- function(t) 6 * (0.1 * sin(2 * pi * t) + 0.2 * cos(2 * pi * t) +
                         0.3 * sin(2 * pi * t) ^ 2 + 0.4 * cos(2 * pi * t) ^ 3 +
                         0.5 * sin(2 * pi * t) ^ 3)


xvar = "V4"

fhats <- lapply(dat_list, function(i) f.hat(object = i$sail.fit, xvar = xvar, s = i$lambda.1se)[["fX"]])
fhats_dat <- do.call(cbind,fhats)
ylims <- range(fhats_dat)

Xs <- lapply(dat_list, function(i) f.hat(object = i$sail.fit, xvar = xvar, s = i$lambda.1se)[["X"]])
Xs_dat <- do.call(cbind,Xs)
xlims <- range(Xs_dat)


plot.args <- list(x = Xs[[1]],
                  y = fhats[[1]],
                  ylim = c(ylims[1], ylims[2]+0),
                  xlim = xlims,
                  xlab = xvar,
                  ylab = sprintf("f(%s)",xvar),
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
                  family = "serif")
do.call("plot", plot.args)
abline(h = 0, lwd = 1, col = "gray")
for (i in seq_len(ncol(Xs_dat))) {
  lines(Xs_dat[,i], fhats_dat[,i], col = cbbPalette()[3], lwd = 1)
}
curve(f4, add = TRUE)
