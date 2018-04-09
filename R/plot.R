######################################
#' R Source code file for plotting functions
#' plotSailCoef is called by plot.sail and is not exported
#' plotMain and plotInter are exported
#' Author: Sahir Bhatnagar
#' Created: 2016
#' Updated: April 9, 2018
#####################################


#' Plot the coefficient plot produced by sail
#'
#' @description Plot the coefficient plot produced by sail
plotSailCoef <- function(coefs, lambda, group, df, dev, vnames, environ,
                         alpha = 1, legend.loc, label = FALSE, log.l = TRUE,
                         norm = FALSE, ...) {

  # browser()
  if (alpha < 0 | alpha > 1) {
    warning("alpha must be in the range [0,1]. Setting alpha = 1")
    alpha <- 1
  }

  if (norm) { # not implemented for now
    # Y <- predict(x, type = "norm")
    # index <- Y
    # approx.f = 1
    # if (any(x$group == 0))
    #   Y <- Y[-1, ]
    # nonzero <- which(apply(abs(Y), 1, sum) != 0)
    # Y <- Y[nonzero, ]
    # g <- 1:nrow(Y)
  } else {
    if (length(dim(coefs)) == 3) {
      beta <- matrix(coefs[, -1, , drop = FALSE], ncol = dim(coefs)[3])
    } else {
      beta <- coefs
    }
    penalized <- which(group != 0)
    nonzero <- which(apply(abs(beta), 1, sum) != 0)
    ind <- intersect(penalized, nonzero)
    Y <- as.matrix(beta[ind, , drop = FALSE])
    g <- as.numeric(as.factor(group[ind]))
  }
  p <- nrow(Y)
  l <- lambda
  n.g <- max(g)
  if (log.l) {
    l <- log(l)
    index <- l
    approx.f <- 0
    xlab <- expression(log(lambda))
  } else {
    xlab <- expression(lambda)
    index <- lambda
    approx.f <- 0
  }

  ylims <- if (!missing(environ)) range(Y, environ) else range(Y)
  plot.args <- list(
    x = l, y = 1:length(l), ylim = ylims,
    xlab = xlab, ylab = "", type = "n",
    xlim = rev(range(l)),
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
  if (plot.args$ylab == "") {
    ylab <- if (norm) {
      expression("||" * hat(theta) * "||")
    } else {
      expression(hat(theta))
    }
    mtext(ylab, 2, 3.5, las = 1, adj = 0, cex = 2)
  }
  abline(h = 0, lwd = 0.8, col = "gray")
  cols <- hcl(
    h = seq(15, 375, len = max(4, n.g + 1)), l = 60,
    c = 150, alpha = alpha
  )
  cols <- if (n.g == 2) cols[c(1, 3)] else cols[1:n.g]
  line.args <- list(
    col = cols, lwd = 1 + 2 * exp(-p / 20),
    lty = 1, pch = ""
  )
  if (length(new.args)) {
    line.args[names(new.args)] <- new.args
  }
  line.args$x <- l
  line.args$y <- t(Y)
  line.args$col <- line.args$col[g]
  line.args$lty <- rep(line.args$lty, length.out = max(g))
  line.args$lty <- line.args$lty[g]
  do.call("matlines", line.args)
  if (!missing(environ)) lines(l, environ, lwd = line.args$lwd)
  if (!missing(legend.loc)) {
    legend.args <- list(
      col = cols, lwd = line.args$lwd,
      lty = line.args$lty, legend = vnames
    )
    if (length(new.args)) {
      new.legend.args <- new.args[names(new.args) %in%
        names(formals(legend))]
      legend.args[names(new.legend.args)] <- new.legend.args
    }
    legend.args$x <- legend.loc
    do.call("legend", legend.args)
  }
  if (label) {
    ypos <- Y[, ncol(Y)]
    text(min(l), ypos, names(ypos), xpd = NA, adj = c(
      0,
      NA
    ))
  }

  atdf <- pretty(index)
  prettydf <- stats::approx(
    x = index, y = df, xout = atdf, rule = 2,
    method = "constant", f = approx.f
  )$y
  axis(3,
    at = atdf, labels = prettydf, mgp = c(3, .3, 0),
    tcl = NA,
    cex.axis = 1.2
  )
}





#' @title Plot Estimated Component Smooth Functions for Main Effects
#' @description Takes a fitted sail object produced by \code{sail()} or
#'   \code{cv.sail()$sail.fit} and plots the component smooth function for a
#'   pre-specified variable at a given value of lambda and on the scale of the
#'   linear predictor. Currently only implemented for \code{type="gaussian"}
#' @param object a fitted \code{sail} object as produced by \code{sail()} or
#'   \code{cv.sail()$sail.fit}
#' @param x original data supplied to the original call to \code{\link{sail}}
#' @param xvar a character corresponding to the predictor to be plotted. Only
#'   one variable name should be supplied, if more than one is supplied, only
#'   the first element will be plotted. This variable name must be in
#'   \code{colnames(x)}.
#' @param s a single value of the penalty parameter \code{lambda} at which
#'   coefficients will be extracted via the \code{coef} method for objects of
#'   class \code{"sail"}. If more than one is supplied, only the first one will
#'   be used.
#' @param f.truth true function. Only used for simulation purposes when the
#'   truth is known. The function takes as a input a numeric vector
#'   corresponding the \code{xvar} column in \code{x} of length \code{nrow(x)}.
#'   A second line will be plotted for the truth and a legend is added to the
#'   plot.
#' @param col color of the line. The first element corresponds to the color used
#'   for the estimated function and the second element is for the true function
#'   (if \code{f.truth} is specified). Default: c("#D55E00", "#009E73")
#' @param legend.position position of the legend. Only used when \code{f.truth}
#'   is specified. Default: 'bottomleft'. Can be a single keyword from the list
#'   "bottomright", "bottom", "bottomleft", "left", "topleft", "top",
#'   "topright", "right" and "center". This places the legend on the inside of
#'   the plot frame at the given location. Partial argument matching is used.
#' @param rug adds a rug representation (1-d plot) of the data to the plot, logical. Default: TRUE.
#' @param ... other graphical paramters passed to \code{plot}.
#' @return A plot is produced and nothing is returned
#' @details The linear predictor \eqn{\beta_0 + basis(xvar) * \beta_xvar} is
#'   plotted against \code{xvar}, where \code{basis} is the expansion provided
#'   in the original call to \code{sail}.
#' @examples
#' \dontrun{
#' if(interactive()){
#' # Parallel
#' library(doMC)
#' registerDoMC(cores = 4)
#' data("sailsim")
#' f.basis <- function(i) splines::bs(i, degree = 5)
#' cvfit <- cv.sail(x = sailsim$x, y = sailsim$y, e = sailsim$e,
#'                  basis = f.basis, nfolds = 10, parallel = TRUE)
#' # plot cv-error curve
#' plot(cvfit)
#' # non-zero estimated coefficients at lambda.min
#' predict(cvfit, type = "nonzero", s="lambda.min")
#' # plot main effect for X4 with a line for the truth also
#' plotMain(cvfit$sail.fit, x = sailsim$x, xvar = "X4",
#'          s = cvfit$lambda.min, f.truth = sailsim$f4)
#'  }
#' }
#' @seealso \code{\link{coef.sail}} \code{\link{predict.sail}}, \code{\link[graphics]{rug}}
#' @rdname plotMain
#' @export
plotMain <- function(object, x, xvar, s, f.truth, col = c("#D55E00", "#009E73"),
                     legend.position = "bottomleft", rug = TRUE, ...) {

  # browser()
  if (length(xvar) > 1) {
    xvar <- xvar[[1]]
    warning("More than 1 xvar provided. Only first element will be plotted.")
  }

  if (length(s) > 1) {
    s <- s[[1]]
    warning("More than 1 s value provided. Only first element will be used for the estimated coefficients.")
  }

  ind <- object$group == which(object$vnames == xvar)
  allCoefs <- coef(object, s = s)
  a0 <- allCoefs[1, ]
  betas <- as.matrix(allCoefs[object$main.effect.names[ind], , drop = FALSE])
  design.mat <- object$design[, object$main.effect.names[ind], drop = FALSE]
  originalX <- x[, unique(object$group[ind])]

  # f.hat <- drop(a0 + design.mat %*% betas)
  f.hat <- drop(design.mat %*% betas)
  if (!missing(f.truth)) {
    seqs <- seq(range(originalX)[1],range(originalX)[2], length.out = 100)
    f.truth.eval <- f.truth(seqs)
    ylims <- range(f.truth.eval, f.hat)
  } else { ylims <- range(f.hat) }

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



#' Plot Interaction Effects from sail object
#' @param object sail object
#'
#' @export
plotInter <- function(object, x, xvar, s, f.truth, simulation = TRUE, truthonly = FALSE,
                      npoints = 30, col = c("#56B4E9", "#D55E00"), title_z,
                      legend.position = "bottomleft", ...) {

  browser()
  if (length(xvar) > 1) {
    xvar <- xvar[[1]]
    warning("More than 1 xvar provided. Only first element will be plotted.")
  }

  if (length(s) > 1) {
    s <- s[[1]]
    warning("More than 1 s value provided. Only first element will be used for the estimated coefficients.")
  }

  ind <- object$group == which(object$vnames == xvar)
  allCoefs <- coef(object, s = s)
  a0 <- allCoefs[1, ]

  betas <- as.matrix(allCoefs[object$main.effect.names[ind], , drop = FALSE])
  alphas <- as.matrix(allCoefs[object$interaction.names[ind], , drop = FALSE])
  betaE <- as.matrix(allCoefs["E", , drop = FALSE])

  design.mat.main <- object$design[, object$main.effect.names[ind], drop = FALSE]
  design.mat.int <- object$design[, object$interaction.names[ind], drop = FALSE]
  originalE <- object$design[, "E", drop = FALSE] # this is the centered E
  originalX <- x[, unique(object$group[ind])]

  # f.hat <- drop(a0 + design.mat %*% betas)
  # f.hat <- drop(originalE %*% betaE + design.mat.main %*% betas + design.mat.int %*% alphas)

  # all.equal((e * standardize(splines::bs(x, df = object$df, degree = object$degree))$x),
  #           sweep(standardize(splines::bs(x, df = object$df, degree = object$degree))$x, 1, e, FUN = "*"))

  x <- seq(range(originalX)[1], range(originalX)[2], length.out = npoints)
  e <- seq(range(originalE)[1], range(originalE)[2], length.out = npoints)

  # show interaction effect only for simulation
  if (simulation) {
    f.est <- function(X, E) {
      # E * as.vector(betaE) + standardize(splines::bs(X, df = object$df, degree = object$degree))$x %*% betas +
      (E * standardize(splines::bs(X, df = object$df, degree = object$degree))$x) %*% alphas
    }
  } else {
    f.est <- function(X, E) {
      E * as.vector(betaE) + standardize(splines::bs(X, df = object$df, degree = object$degree))$x %*% betas +
        (E * standardize(splines::bs(X, df = object$df, degree = object$degree))$x) %*% alphas
    }
  }

  # f.truth <- function(x, e) { e * DT$f4(x)  }
  z.est <- outer(x, e, f.est)

  if (!missing(f.truth)) {
    z.truth <- outer(x, e, f.truth)
    z_range <- c(min(z.est, z.truth), max(z.est, z.truth))
  } else {
    z_range <- c(min(z.est), max(z.est))
  }

  op <- par(bg = "white")

  # c(bottom, left, top, right)
  if (truthonly) {
    par(mfrow = c(1, 1), tcl = -0.5, family = "serif", omi = c(0.2, 0.2, 0, 0))
    par(mai = c(0., 0.2, 0.4, 0.))
    graphics::persp(x, e, z.truth,
      zlim = z_range,
      theta = 30, phi = 30,
      ltheta = 120, expand = 0.5,
      r = 2, shade = 0.3, axes = TRUE, scale = TRUE, box = T,
      nticks = 5,
      # ticktype="detailed",
      col = col[2],
      cex.lab = 3,
      cex.main = 3,
      xlab = sprintf("f(%s)", xvar),
      ylab = "X_E",
      zlab = "Y", main = "Truth"
    )
  } else if (!missing(f.truth)) {
    par(mfrow = c(1, 2), tcl = -0.5, family = "serif", omi = c(0.2, 0.2, 0, 0))
    par(mai = c(0., 0.8, 0.6, 0.))
    graphics::persp(x, e, z.truth,
      zlim = z_range,
      theta = 30, phi = 30,
      ltheta = 120, expand = 0.5,
      r = 2, shade = 0.3, axes = TRUE, scale = TRUE, box = T,
      nticks = 5,
      # ticktype="detailed",
      col = col[2],
      cex.lab = 3,
      cex.main = 3,
      xlab = sprintf("f(%s)", xvar),
      ylab = "X_E",
      zlab = "Y", main = "Truth"
    )
    graphics::persp(x, e, z.est,
      theta = 30, phi = 30,
      ltheta = 120, expand = 0.5,
      r = 2, shade = 0.3, axes = TRUE, scale = TRUE, box = T,
      nticks = 5,
      zlim = z_range,
      cex.lab = 3,
      cex.main = 3,
      # ticktype="detailed",
      col = col[1],
      xlab = sprintf("f(%s)", xvar),
      ylab = "X_E",
      zlab = "Y", main = title_z
    )
  } else {
    par(mfrow = c(1, 1), tcl = -0.5, family = "serif", omi = c(0.2, 0.2, 0, 0))
    par(mai = c(0., 0.2, 0.4, 0.))
    graphics::persp(x, e, z.est,
      cex.lab = 3,
      cex.main = 3,
      theta = 30, phi = 30,
      ltheta = 120, expand = 0.5,
      r = 2, shade = 0.3, axes = TRUE, scale = TRUE, box = T,
      nticks = 5,
      zlim = z_range,
      col = col[1],
      # xlab=sprintf("f(x_%s)",as.numeric(gsub("X","",4))),
      xlab = sprintf("f(%s)", xvar),
      ylab = "X_E",
      zlab = "Y",
      # main=sprintf("Estimated Interaction Effect for %s",xvar)
      main = title_z
    )
  }
}
