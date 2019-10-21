# this extrapolates to the entire sample. not just the training set
# plotMainADNI <- function(object, x, design, xvar, s, f.truth,
#                          legend.position = "bottomleft", rug = TRUE, col=sail:::cbbPalette[c(6,4,7)],...) {
#
#   # browser()
#
#   if (length(s) > 1) {
#     s <- s[[1]]
#     warning("More than 1 s value provided. Only first element will be used for the estimated coefficients.")
#   }
#
#   ind <- which(object$vnames %in% xvar)
#   allCoefs <- coef(object, s = s)
#   a0 <- allCoefs[1, ]
#   #betas <- as.matrix(allCoefs[object$main.effect.names[c((3*ind-2),(3*ind-1),(3*ind))], , drop = FALSE])
#   #design.mat <- design[, xvar, drop = FALSE]
#   originalX <- x
#
#   # f.hat <- drop(a0 + design.mat %*% betas)
#   f.hat <- drop(design.mat %*% betas)
#   if (!missing(f.truth)) {
#     seqs <- seq(range(originalX)[1], range(originalX)[2], length.out = 100)
#     f.truth.eval <- f.truth(seqs)
#     ylims <- range(f.truth.eval, f.hat)
#   } else {
#     ylims <- range(f.hat)
#   }
#
#   plot.args <- list(
#     x = originalX[order(originalX)],
#     y = f.hat[order(originalX)],
#     ylim = ylims,
#     xlab = xvar,
#     ylab = sprintf("f(%s)", xvar),
#     type = "n",
#     # xlim = rev(range(l)),
#     # las = 1,
#     cex.lab = 1.5,
#     cex.axis = 1.5,
#     cex = 1.5,
#     # bty = "n",
#     # mai=c(1,1,0.1,0.2),
#     # tcl = -0.5,
#     # omi = c(0.2,1,0.2,0.2),
#     family = "serif"
#   )
#   new.args <- list(...)
#   if (length(new.args)) {
#     new.plot.args <- new.args[names(new.args) %in% c(
#       names(par()),
#       names(formals(plot.default))
#     )]
#     plot.args[names(new.plot.args)] <- new.plot.args
#   }
#   do.call("plot", plot.args)
#   abline(h = 0, lwd = 1, col = "gray")
#   lines(originalX[order(originalX)], f.hat[order(originalX)], col = col[1], lwd = 3)
#   if (rug) graphics::rug(originalX, side = 1)
#   if (!missing(f.truth)) {
#     lines(seqs[order(seqs)], f.truth.eval[order(seqs)], col = col[2], lwd = 3)
#   }
#   if (!missing(f.truth)) {
#     legend(legend.position,
#            c("Estimated", "Truth"),
#            col = col, cex = 1, bty = "n", lwd = 3
#     )
#   }
# }








# this extrapolates to the entire sample. not just the training set
# plotInterADNI <- function(object, x, xvar, s,
#                           design, # this contains user defined expand matrix
#                           e, # this is E vector for whole sample
#                           xlab , ylab ,
#                           legend.position = "bottomleft", main = "", rug = TRUE,
#                           color = sail:::cbbPalette[c(6,4,7)], legend = TRUE) {
#
#   # cv_obj = cvfit; original_name = "X60"; sail_name = "X19"; xlab =  "supramarginal gyrus right";
#   # lambda_type = "lambda.min";ylab = "Mini-Mental State Examination";
#   # color = RColorBrewer::brewer.pal(9,"Set1"); legend = TRUE; ylim =  c(15,30)
#   # ==================
#   # browser()
#   ind <- which(object$vnames %in% xvar)
#   allCoefs <- coef(object, s = s)
#   a0 <- allCoefs[1, ]
#   betas <- as.matrix(allCoefs[object$main.effect.names[ind], , drop = FALSE])
#   alphas <- as.matrix(allCoefs[object$interaction.names[ind], , drop = FALSE])
#   betaE <- as.matrix(allCoefs["E", , drop = FALSE])
#
#   # if you dont want to extrapolate, un-comment the following lines
#    #design.mat.main <- object$design[, object$main.effect.names[ind], drop = FALSE]
#   #design.mat.int <- object$design[, object$interaction.names[ind], drop = FALSE]
#
#   design.mat.main <- design[, xvar, drop = FALSE]
#   design.mat.int <- design[, xvar, drop = FALSE] * e
#
#
#   # originalE <- object$design[, "E", drop = FALSE] # this is the centered E
#   # originalX <- x
#
#   originalE <- e # this is the centered E
#   originalX <- x
#
#
#   # f.hat <- drop(a0 + design.mat %*% betas)
#   f.hat <- drop(originalE * as.vector(betaE) +
#                   design.mat.main %*% betas + design.mat.int %*% alphas)
#
#   # f.hat <- drop(originalE * as.vector(betaE)  + design.mat.int %*% alphas)
#   # f.hat <- drop(design.mat.int %*% alphas)
#   ylims <- range(f.hat)
#
#
#
#
#   intervention = drop(unique(originalE))[1]
#   control = drop(unique(originalE))[2]
#
#
#
#   # 0=control, 1=Intervention
#   cont_index <- which(originalE==control)
#   int_index <- which(originalE==intervention)
#
#
#   cont_pred <- f.hat[cont_index]
#   int_pred <- f.hat[int_index]
#
#
#
#
#   min.length.top <- range(f.hat)[1] ; max.length.top <- range(f.hat)[2]
#   par(mai=c(1,1,1,0.2))
#   plot(originalX, f.hat,
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
#        main = main,
#        cex.main = 2.5,
#        # ylim = c(min.length.top-3, max.length.top+3),
#        ylim = ylims)
#   axis(1, labels = T, cex.axis = 2)
#
#
#   lines(originalX[cont_index][order(originalX[cont_index])], cont_pred[order(originalX[cont_index])], col = color[1], lwd = 3)
#   lines(originalX[int_index][order(originalX[int_index])], int_pred[order(originalX[int_index])], col = color[2], lwd = 3)
#   legend("bottomleft",c("Intervention","Control"),col=c(color[2],color[1]),lwd=c(3,3))
#
#
# }




plotInterPRS <- function(object,
                         originalDataNotCentered,
                         xvar,
                         s,
                         design, # this contains the design matrix you want to extrapolate to. if missing, then it will use design from object
                         evar, # this is E vector for whole sample
                         xlab = xvar,
                         degree,
                         ylab = "Marginal Risk",
                         legend.position = "bottomleft",
                         main = "",
                         rug = TRUE,
                         color = sail:::cbbPalette[c(6,4,7)],
                         legend = TRUE) {

  # cv_obj = cvfit; original_name = "X60"; sail_name = "X19"; xlab =  "supramarginal gyrus right";
  # lambda_type = "lambda.min";ylab = "Mini-Mental State Examination";
  # color = RColorBrewer::brewer.pal(9,"Set1"); legend = TRUE; ylim =  c(15,30)
  # ==================
  # browser()
  ind <- which(object$main.effect.names == paste(xvar,1:degree, sep = "_"))
  allCoefs <- coef(object, s = s)
  a0 <- allCoefs[1, ]
  betas <- as.matrix(allCoefs[object$main.effect.names[ind], , drop = FALSE])
  alphas <- as.matrix(allCoefs[object$interaction.names[ind], , drop = FALSE])
  betaE <- as.matrix(allCoefs["E", , drop = FALSE])

  if(missing(design)) {
    design.mat.main <- object$design[, object$main.effect.names[ind], drop = FALSE]
    design.mat.int <- object$design[, object$interaction.names[ind], drop = FALSE]
    design.mat.e <- object$design[,"E"] # this is the non-centered E
  } else {
    design.mat.main <- design[, object$main.effect.names[ind], drop = FALSE]
    design.mat.int <- design[, object$interaction.names[ind], drop = FALSE]
    design.mat.e <- design[,"E"] # this is the non-centered E
  }

  # design.mat.main <- design[, xvar, drop = FALSE]
  # design.mat.int <- design[, xvar, drop = FALSE] * e

  # originalE <- object$design[, "E", drop = FALSE] # this is the centered E
  # originalX <- x
  originalX <- originalDataNotCentered[,xvar]


  # f.hat <- drop(a0 + design.mat %*% betas)
  f.hat <- drop(design.mat.e * as.vector(betaE) +
                  design.mat.main %*% betas + design.mat.int %*% alphas)

  # f.hat <- drop(originalE * as.vector(betaE)  + design.mat.int %*% alphas)
  # f.hat <- drop(design.mat.int %*% alphas)
  ylims <- range(f.hat)

  intervention <- drop(unique(design.mat.e))[1]
  control <- drop(unique(design.mat.e))[2]

  # 0=control, 1=Intervention
  cont_index <- which(design.mat.e==control)
  int_index <- which(design.mat.e==intervention)

  cont_pred <- f.hat[cont_index]
  int_pred <- f.hat[int_index]

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


  lines(originalX[cont_index][order(originalX[cont_index])], cont_pred[order(originalX[cont_index])], col = color[1], lwd = 3)
  lines(originalX[int_index][order(originalX[int_index])], int_pred[order(originalX[int_index])], col = color[2], lwd = 3)

  if (legend) {
    legend(legend.position,c("Intervention","Control"),col=c(color[2],color[1]),lwd=c(3,3))
  }

  return(list(control_x = originalX[cont_index][order(originalX[cont_index])],
              control_y = cont_pred[order(originalX[cont_index])],
              intervention_x = originalX[int_index][order(originalX[int_index])],
              intervention_y = int_pred[order(originalX[int_index])]))

}















