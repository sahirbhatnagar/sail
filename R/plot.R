#' Plot the cross-validation curve produced by cv.sail
#'
#' @description Plots the cross-validation curve, and upper and lower standard
#'   deviation curves, as a function of the \eqn{\lambda_\beta} and
#'   \eqn{\lambda_\gamma} values used. Using \code{ggplot2} facet plots, each
#'   facet represents a unique value for \eqn{\lambda_\gamma}, and the x-axis is
#'   the sequence of corresponding \eqn{\lambda_\beta}
#' @param x fitted \code{cv.sail} object
#' @details A plot is produced, and nothing is returned. A colored vertical line
#'   is drawn at the pair of tuning parameters that leads to the minimum CV
#'   error and another is drawn at the 1 standard error rule pair of tuning
#'   parameters
#' @seealso \code{\link{sail}} and \code{\link{cv.sail}}
#' @author
#' Sahir Bhatnagar
#'
#' Maintainer: Sahir Bhatnagar \email{sahir.bhatnagar@@mail.mcgill.ca}
#' @import data.table
#' @export

plot.cv.sail_old <- function(x) {
  pacman::p_load(ggplot2)
  pacman::p_load(data.table)
  pacman::p_load(latex2exp)

  # x = cvfit
  # ====

  # browser()
  cvobj <- x

  d <- data.frame(cvobj$df,
    lambda.min.beta = cvobj$lambda.min.beta,
    lambda.1se.beta = cvobj$lambda.1se.beta,
    row.names = rownames(cvobj$df)
  )


  # needed to get colored lines
  d2 <- data.table::melt(d[which(rownames(d) %in% c(cvobj$lambda.min.name, cvobj$lambda.1se.name)), , drop = F],
    measure.vars = c("lambda.min.beta", "lambda.1se.beta")
  )

  # d2 <- as.data.table(d2)
  d2 <- transform(d2, variable = gsub(".beta", "", variable))

  appender <- function(string) latex2exp::TeX(paste("$\\log(\\lambda_{\\gamma}) = $", string))

  p <- ggplot2::ggplot(
    d,
    ggplot2::aes(log(lambda.beta),
      ymin = lower,
      ymax = upper
    )
  )

  l <- ggplot2::ggplot_build(p)
  p + ggplot2::geom_errorbar(color = "grey", width = 0.5) +
    ggplot2::geom_point(aes(x = log(lambda.beta), y = mse), colour = "red") +
    # theme_bw() +
    # ylim(c(min(d$lower) - 10 , max(d$upper) + 500)) +
    ggplot2::facet_wrap(~log.gamma,
      scales = "fixed",
      # switch = "x",
      labeller = ggplot2::as_labeller(appender, default = ggplot2::label_parsed)
    ) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_text(size = ggplot2::rel(1.3)),
      legend.position = "bottom"
    ) +
    ggplot2::xlab(TeX("$\\log(\\lambda_{\\beta})$")) +
    ggplot2::geom_vline(
      data = d2[(d2$lambda.beta == d2$value & d2$variable == "lambda.1se"), ],
      aes(xintercept = log(value), colour = variable), size = 0.7, linetype = 1
    ) +
    geom_vline(
      data = d2[(d2$lambda.beta == d2$value & d2$variable == "lambda.min"), ],
      aes(xintercept = log(value), colour = variable), size = 0.7, linetype = 1
    ) +
    ggplot2::scale_color_discrete(name = "") +
    ggplot2::geom_text(aes(label = nz.main, x = log(lambda.beta), y = Inf, vjust = 1)) +
    ggplot2::geom_text(aes(
      label = nz.interaction, x = log(lambda.beta), y = Inf,
      vjust = 2
    )) +
    ggplot2::ylab(c("5 fold CV MSE")) #+
  # coord_cartesian(ylim = c(l$panel$ranges[[1]]$y.range[1], l$panel$ranges[[1]]$y.range[2]*1.1))
}


#' @export
plot.cv.sail <- function(x, sign.lambda = 1, ...) {
  cvobj <- x
  xlab <- "log(Lambda)"
  if (sign.lambda < 0) xlab <- paste("-", xlab, sep = "")
  plot.args <- list(x = sign.lambda * log(cvobj$lambda), y = cvobj$cvm, ylim = range(cvobj$cvup, cvobj$cvlo), xlab = xlab, ylab = cvobj$name, type = "n")
  new.args <- list(...)
  if (length(new.args)) plot.args[names(new.args)] <- new.args
  do.call("plot", plot.args)
  error.bars(sign.lambda * log(cvobj$lambda), cvobj$cvup, cvobj$cvlo, width = 0.01, col = "darkgrey")
  points(sign.lambda * log(cvobj$lambda), cvobj$cvm, pch = 20, col = "red")
  axis(side = 3, at = sign.lambda * log(cvobj$lambda), labels = paste(cvobj$nz), tick = FALSE, line = 0)
  abline(v = sign.lambda * log(cvobj$lambda.min), lty = 3)
  abline(v = sign.lambda * log(cvobj$lambda.1se), lty = 3)
  invisible()
}




# plotCoefSail <- function(beta, norm, lambda, df, dev, label = FALSE,
#                          xvar = c("norm", "lambda", "dev"),
#                          xlab = iname, ylab = "Coefficients", ...) {
#   which = nonzero(beta)
#   nwhich = length(which)
#   switch(nwhich + 1, `0` = {
#     warning("No plot produced since all coefficients zero")
#     return()
#   }, `1` = warning("1 or less nonzero coefficients; glmnet plot is not meaningful"))
#   beta = as.matrix(beta[which, , drop = FALSE])
#   xvar = match.arg(xvar)
#   switch(xvar, norm = {
#     index = if (missing(norm)) apply(abs(beta), 2, sum) else norm
#     iname = "L1 Norm"
#     approx.f = 1
#   }, lambda = {
#     index = log(lambda)
#     iname = "Log Lambda"
#     approx.f = 0
#   }, dev = {
#     index = dev
#     iname = "Fraction Deviance Explained"
#     approx.f = 1
#   })
#   dotlist = list(...)
#   type = dotlist$type
#   if (is.null(type))
#     matplot(index, t(beta), lty = 1, xlab = xlab, ylab = ylab,
#             type = "l", ...)
#   else matplot(index, t(beta), lty = 1, xlab = xlab, ylab = ylab,
#                ...)
#   atdf = pretty(index)
#   prettydf = approx(x = index, y = df, xout = atdf, rule = 2,
#                     method = "constant", f = approx.f)$y
#   axis(3, at = atdf, labels = prettydf, tcl = NA)
#   if (label) {
#     nnz = length(which)
#     xpos = max(index)
#     pos = 4
#     if (xvar == "lambda") {
#       xpos = min(index)
#       pos = 2
#     }
#     xpos = rep(xpos, nnz)
#     ypos = beta[, ncol(beta)]
#     text(xpos, ypos, paste(which), cex = 0.5, pos = pos)
#   }
# }

# myplotgrp(fit, log.l = T, ylab = "main")

#' Plot the coefficient plot produced by sail
#'
#' @description Plot the coefficient plot produced by sail
#' @export

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





#' Plot Main Effects from sail object
#' @param object sail object
#' @export
plotMain <- function(object, xvar, s, f.truth, col = c("#D55E00", "#009E73"), legend.position = "bottomleft", ...) {

  # browser()
  ind <- object$group == which(object$vnames == xvar)
  allCoefs <- coef(object, s = s)
  a0 <- allCoefs[1, ]
  betas <- as.matrix(allCoefs[object$main.effect.names[ind], , drop = FALSE])
  design.mat <- object$design[, object$main.effect.names[ind], drop = FALSE]
  originalX <- object$x[, unique(object$group[ind])]

  f.hat <- drop(a0 + design.mat %*% betas)
  # f.hat <- drop(design.mat %*% betas)
  ylims <- if (!missing(f.truth)) range(f.truth, f.hat) else range(f.hat)

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
  rug(originalX, side = 1)
  if (!missing(f.truth)) {
    lines(originalX[order(originalX)], f.truth[order(originalX)], col = col[2], lwd = 3)
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
#' @export
plotInter <- function(object, xvar, s, f.truth, simulation = TRUE, truthonly = FALSE,
                      npoints = 30, col = c("#56B4E9", "#D55E00"), title_z,
                      legend.position = "bottomleft", ...) {

  # browser()

  ind <- object$group == which(object$vnames == xvar)
  allCoefs <- coef(object, s = s)
  a0 <- allCoefs[1, ]

  betas <- as.matrix(allCoefs[object$main.effect.names[ind], , drop = FALSE])
  alphas <- as.matrix(allCoefs[object$interaction.names[ind], , drop = FALSE])
  betaE <- as.matrix(allCoefs["E", , drop = FALSE])

  design.mat.main <- object$design[, object$main.effect.names[ind], drop = FALSE]
  design.mat.int <- object$design[, object$interaction.names[ind], drop = FALSE]
  originalE <- object$design[, "E", drop = FALSE] # this is the centered E
  originalX <- object$x[, unique(object$group[ind])]

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
    persp(x, e, z.truth,
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
    persp(x, e, z.truth,
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
    persp(x, e, z.est,
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
    persp(x, e, z.est,
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
