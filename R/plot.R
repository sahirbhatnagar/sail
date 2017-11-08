#' Plot the cross-validation curve produced by cv.shim
#'
#' @description Plots the cross-validation curve, and upper and lower standard
#'   deviation curves, as a function of the \eqn{\lambda_\beta} and
#'   \eqn{\lambda_\gamma} values used. Using \code{ggplot2} facet plots, each
#'   facet represents a unique value for \eqn{\lambda_\gamma}, and the x-axis is
#'   the sequence of corresponding \eqn{\lambda_\beta}
#' @param x fitted \code{cv.shim} object
#' @details A plot is produced, and nothing is returned. A colored vertical line
#'   is drawn at the pair of tuning parameters that leads to the minimum CV
#'   error and another is drawn at the 1 standard error rule pair of tuning
#'   parameters
#' @seealso \code{\link{shim}} and \code{\link{cv.shim}}
#' @author
#' Sahir Bhatnagar
#'
#' Maintainer: Sahir Bhatnagar \email{sahir.bhatnagar@@mail.mcgill.ca}
#' @import data.table
#' @export

plot.cv.funshim <- function(x) {

  pacman::p_load(ggplot2)
  pacman::p_load(data.table)
  pacman::p_load(latex2exp)

  # x = cvfit
  #====
  cvobj <- x

  d <- data.table::as.data.table(transform(cvobj$df,
                               lambda.min.beta = cvobj$lambda.min.beta,
                               lambda.1se.beta = cvobj$lambda.1se.beta),
                     keep.rownames = TRUE)


  # needed to get colored lines
  d2 <- data.table::melt(d[rn %in% c(cvobj$lambda.min.name, cvobj$lambda.1se.name)],
                         measure.vars = c("lambda.min.beta","lambda.1se.beta"))

  d2[,variable := gsub(".beta", "",variable)]

  appender <- function(string) latex2exp::TeX(paste("$\\log(\\lambda_{\\gamma}) = $",string))

  p <- ggplot2::ggplot(d,
              ggplot2::aes(log(lambda.beta),
                  ymin = lower,
                  ymax = upper))

  l <- ggplot2::ggplot_build(p)
  p + ggplot2::geom_errorbar(color = "grey", width = 0.5) +
    geom_point(aes(x = log(lambda.beta), y = mse), colour = "red") +
    # theme_bw() +
    # ylim(c(min(d$lower) - 10 , max(d$upper) + 500)) +
    facet_wrap(~log.gamma, scales = "fixed",
               #switch = "x",
               labeller = as_labeller(appender, default = label_parsed)) +
    theme(strip.background = element_blank(),
          strip.text.x = element_text(size = rel(1.3)),
          legend.position = "bottom") +
    xlab(TeX("$\\log(\\lambda_{\\beta})$")) +
    geom_vline(data = d2[lambda.beta == value & variable == "lambda.1se"],
               aes(xintercept = log(value), colour = variable), size = 0.7, linetype = 1) +
    geom_vline(data = d2[lambda.beta == value & variable == "lambda.min"],
               aes(xintercept = log(value), colour = variable),size = 0.7, linetype = 1) +
    scale_color_discrete(name="") +
    geom_text(aes(label = nz.main, x = log(lambda.beta), y = Inf, vjust = 1)) +
    geom_text(aes(label = nz.interaction, x = log(lambda.beta), y = Inf,
                  vjust = 2)) +
    ylab(c("5 fold CV MSE")) #+
    # coord_cartesian(ylim = c(l$panel$ranges[[1]]$y.range[1], l$panel$ranges[[1]]$y.range[2]*1.1))
}


#' Plot the coefficient plot produced by shim
#'
#' @description Plot the coefficient plot produced by shim
#'
#' @export

plotCoefShim <- function(beta, norm, lambda, df, dev, label = FALSE,
                         xvar = c("norm", "lambda", "dev"),
                         xlab = iname, ylab = "Coefficients", ...) {
  which = nonzero(beta)
  nwhich = length(which)
  switch(nwhich + 1, `0` = {
    warning("No plot produced since all coefficients zero")
    return()
  }, `1` = warning("1 or less nonzero coefficients; glmnet plot is not meaningful"))
  beta = as.matrix(beta[which, , drop = FALSE])
  xvar = match.arg(xvar)
  switch(xvar, norm = {
    index = if (missing(norm)) apply(abs(beta), 2, sum) else norm
    iname = "L1 Norm"
    approx.f = 1
  }, lambda = {
    index = log(lambda)
    iname = "Log Lambda"
    approx.f = 0
  }, dev = {
    index = dev
    iname = "Fraction Deviance Explained"
    approx.f = 1
  })
  dotlist = list(...)
  type = dotlist$type
  if (is.null(type))
    matplot(index, t(beta), lty = 1, xlab = xlab, ylab = ylab,
            type = "l", ...)
  else matplot(index, t(beta), lty = 1, xlab = xlab, ylab = ylab,
               ...)
  atdf = pretty(index)
  prettydf = approx(x = index, y = df, xout = atdf, rule = 2,
                    method = "constant", f = approx.f)$y
  axis(3, at = atdf, labels = prettydf, tcl = NA)
  if (label) {
    nnz = length(which)
    xpos = max(index)
    pos = 4
    if (xvar == "lambda") {
      xpos = min(index)
      pos = 2
    }
    xpos = rep(xpos, nnz)
    ypos = beta[, ncol(beta)]
    text(xpos, ypos, paste(which), cex = 0.5, pos = pos)
  }
}





plot.grpreg <- function(x, alpha=1, legend.loc, label=FALSE, log.l=FALSE, norm=FALSE, ...) {

  browser()

    beta <- as(x$beta, "matrix")

    penalized <- which(x$group!=0)
    nonzero <- which(apply(abs(beta),1,sum)!=0)
    ind <- intersect(penalized, nonzero)
    Y <- beta[ind,,drop=FALSE]
    g <- as.numeric(as.factor(x$group[ind]))

  p <- nrow(Y)
  l <- x$lambda.beta
  n.g <- max(g)

  if (log.l) {
    l <- log(l)
    xlab <- expression(log(lambda))
  } else xlab <- expression(lambda)

  plot.args <- list(x=l, y=1:length(l), ylim=range(Y), xlab=xlab, ylab="", type="n", xlim=rev(range(l)), las=1, bty="n")
  new.args <- list(...)
  if (length(new.args)) {
    new.plot.args <- new.args[names(new.args) %in% c(names(par()), names(formals(plot.default)))]
    plot.args[names(new.plot.args)] <- new.plot.args
  }
  do.call("plot", plot.args)
  if (plot.args$ylab=="") {
    ylab <- if (norm) expression("||"*hat(beta)*"||") else expression(hat(beta))
    mtext(ylab, 2, 3.5, las=1, adj=0)
  }
  abline(h=0,lwd=0.5,col="gray")

  cols <- hcl(h=seq(15,375,len=max(4,n.g+1)),l=60,c=150,alpha=alpha)
  cols <- if (n.g==2) cols[c(1,3)] else cols[1:n.g]
  line.args <- list(col=cols, lwd=1+2*exp(-p/20), lty=1, pch="")
  if (length(new.args)) line.args[names(new.args)] <- new.args
  line.args$x <- l
  line.args$y <- t(Y)
  line.args$col <- line.args$col[g]
  line.args$lty <- rep(line.args$lty, length.out=max(g))
  line.args$lty <- line.args$lty[g]
  do.call("matlines",line.args)

  if(!missing(legend.loc)) {
    legend.args <- list(col=cols, lwd=line.args$lwd, lty=line.args$lty, legend=names(x$group.multiplier))
    if (length(new.args)) {
      new.legend.args <- new.args[names(new.args) %in% names(formals(legend))]
      legend.args[names(new.legend.args)] <- new.legend.args
    }
    legend.args$x <- legend.loc
    do.call("legend",legend.args)
  }
  if (label) {
    ypos <- Y[,ncol(Y)]
    text(-0.001, ypos, names(ypos), xpd=NA, adj=c(0,NA))
  }
}
