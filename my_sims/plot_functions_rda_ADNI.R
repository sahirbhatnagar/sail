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
                          stratify = FALSE, # stratify by APOE?
                          apoe = TRUE,
                          ...,
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
  f.hat <- if (!stratify) {
    drop(originalE * as.vector(betaE) +
           design.mat.main %*% betas + design.mat.int %*% alphas)
  } else {
    drop(originalE * as.vector(betaE) + apoee4 * as.vector(betaAPOE) +
           design.mat.main %*% betas + design.mat.int %*% alphas +
           apoee4inter * as.vector(betaAPOEinter))
  }
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
  cont_index <- which(originalE==control)
  mci_index <- which(originalE==mci)
  ad_index <- which(originalE==ad)

  cont_pred <- f.hat[cont_index]
  mci_pred <- f.hat[mci_index]
  ad_pred <- f.hat[ad_index]

  cont_pred0 <- f.hat[cont_index0]
  cont_pred1 <- f.hat[cont_index1]
  mci_pred0 <- f.hat[mci_index0]
  mci_pred1 <- f.hat[mci_index1]
  ad_pred0 <- f.hat[ad_index0]
  ad_pred1 <- f.hat[ad_index1]

  min.length.top <- range(f.hat)[1] ; max.length.top <- range(f.hat)[2]
  # par(mai=c(1,1,1,0.2))
  plot(originalX, f.hat,
       pch = 19,
       ylab = ylab,
       xlab = xlab,
       col = color[1],
       bty="n",
       xaxt="n",
       type = "n",
       # cex.lab = 2,
       # cex.axis = 2,
       # cex = 2,
       main = main,
       # cex.main = 2.5,
       # ylim = c(min.length.top-3, max.length.top+3),
       # ylim = ylims,
       ...)
  # axis(1, labels = T, cex.axis = 2)
  axis(1, labels = T)

  if (!stratify) {
    lines(originalX[cont_index][order(originalX[cont_index])], cont_pred[order(originalX[cont_index])], col = color[1], lwd = 3)
    lines(originalX[mci_index][order(originalX[mci_index])], mci_pred[order(originalX[mci_index])], col = color[2], lwd = 3)
    lines(originalX[ad_index][order(originalX[ad_index])], ad_pred[order(originalX[ad_index])], col = color[3], lwd = 3)
  } else {
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
  }

  # if (legend) legend("bottomright", c("APOE = 1","APOE = 0"), col = color[2:1], pch = 19, cex = 2, bty = "n")
  # text(3, main, cex = 2)
  if (legend) legend(legend.position, c("Control", "Mild Cognitive Impairment","Alzeimer Disease"),
                     col = color[1:3], pch = 19, cex = 2, bty = "n")

  if (rug) graphics::rug(originalX, side = 1)
}

