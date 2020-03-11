######################################
#' R Source code file for predict, coef, plot and print methods for the sail package
#' Author: Sahir Bhatnagar
#' Created: 2016
#' Updated: April 6, 2018
#####################################


#' @title Make predictions from a sail object
#' @description Similar to other predict methods, this functions predicts fitted
#'   values, logits, coefficients and more from a fitted \code{sail} object.
#' @param object Fitted \code{sail} model object
#' @param newx matrix of new values for \code{x} at which predictions are to be
#'   made. Do not include the intercept (this function takes care of that). Must
#'   be a matrix. This argument is not used for
#'   \code{type=c("coefficients","nonzero")}. This matrix will be passed to
#'   \code{\link{design_sail}} to create the design matrix necessary for
#'   predictions. This matrix must have the same number of columns originally
#'   supplied to the \code{sail} fitting function.
#' @param newe vector of new values for the exposure variable \code{e}. This is
#'   passed to the \code{design_sail} function.
#' @param s Value(s) of the penalty parameter \code{lambda} at which predictions
#'   are required. Default is the entire sequence used to create the model.
#' @param type Type of prediction required. Type \code{"link"} gives the linear
#'   predictors for \code{"binomial"} (not implemented yet); for
#'   \code{"gaussian"} models it gives the fitted values. Type \code{"response"}
#'   gives the fitted probabilities for \code{"binomial"} (not implemented yet),
#'   for \code{"gaussian"} type \code{"response"} is equivalent to type
#'   \code{"link"}. Type \code{"coefficients"} computes the coefficients at the
#'   requested values for \code{s}.  Note that for \code{"binomial"} models,
#'   results are returned only for the class corresponding to the second level
#'   of the factor response (not implemented yet). Type \code{"class"} applies
#'   only to \code{"binomial"} models, and produces the class label
#'   corresponding to the maximum probability (not implemented yet). Type
#'   \code{"nonzero"} returns a list of the the nonzero coefficients for each
#'   value of \code{s}. Default: "link"
#' @param ... currently ignored
#' @return The object returned depends on type.
#' @details \code{s} is the new vector at which predictions are requested. If
#'   \code{s} is not in the lambda sequence used for fitting the model, the
#'   predict function will use linear interpolation to make predictions. The new
#'   values are interpolated using a fraction of predicted values from both left
#'   and right lambda indices. \code{coef(...)} is equivalent to
#'   \code{predict(sail.object, type="coefficients",...)}
#' @examples
#' f.basis <- function(i) splines::bs(i, degree = 3)
#' fit <- sail(x = sailsim$x, y = sailsim$y, e = sailsim$e,
#'             basis = f.basis, dfmax = 5, nlambda = 50)
#' predict(fit) # predicted response for whole solution path
#' predict(fit, s = 0.45) # predicted response for a single lambda value
#' predict(fit, s = c(2.15, 0.32, 0.40), type="nonzero") # nonzero coefficients
#' @seealso \code{\link{predict.cv.sail}}
#' @note When the coef method is called, the alpha values, which represent the
#'   interaction term are returned. This alpha is the product of beta_e,gamma_j
#'   and theta_j
#' @rdname predict.sail
#' @export
predict.sail <- function(object, newx, newe, s = NULL,
                         type = c(
                           "link", "response", "coefficients",
                           "nonzero", "class"
                         ), ...) {
  type <- match.arg(type)

  if (missing(newx)) {
    if (!match(type, c("coefficients", "nonzero"), FALSE)) {
      newx <- object$design # this would already have gone through the standardize function
    }
  } else if (!missing(newx) & missing(newe)) {
    stop("newe is missing. please supply the vector of the environment variable.")
  } else if (!missing(newx) & !missing(newe)) {
    newx <- design_sail(
      x = newx, e = newe, expand = object$expand, group = object$group, basis = object$basis, nvars = object$nvars,
      vnames = object$vnames, center.x = object$center.x, center.e = object$center.e
    )$design
  }

  a0 <- t(as.matrix(object$a0))
  rownames(a0) <- "(Intercept)"
  # this has only values for which lambda did converge
  nbeta <- rbind(a0, object$beta, E = object$bE, object$alpha)

  if (!is.null(s)) {
    vnames <- dimnames(nbeta)[[1]]
    dimnames(nbeta) <- list(NULL, NULL)
    lambda <- object$lambda
    lamlist <- lambda.interp(lambda, s)
    if (length(s) == 1) {
      nbeta <- nbeta[, lamlist$left, drop = FALSE] * lamlist$frac +
        nbeta[, lamlist$right, drop = FALSE] * (1 - lamlist$frac)
    } else {
      nbeta <- nbeta[, lamlist$left, drop = FALSE] %*% diag(lamlist$frac) +
        nbeta[, lamlist$right, drop = FALSE] %*% diag(1 - lamlist$frac)
    }
    dimnames(nbeta) <- list(vnames, paste(seq(along = s)))
  }

  if (type == "coefficients") return(nbeta)


  if (type == "nonzero") {
    nbeta.mat <- as.matrix(nbeta)
    if (length(s) == 1) {
      return(nbeta.mat[nonzero(nbeta.mat, bystep = TRUE)[[1]], , drop = FALSE])
    } else {
      nzs <- nonzero(nbeta.mat, bystep = TRUE)
      return(lapply(seq_along(nzs), function(i) nbeta.mat[nzs[[i]], i, drop = FALSE]))
    }
  }

  # this is used by the cv.lspath function to calculate predicted values
  # which will subsequently be used for calculating MSE for each fold
  if (type == "link") {
    nfit <- as.matrix(methods::cbind2(1, newx) %*% nbeta)
    return(nfit)
  }
}

#' @inheritParams predict.sail
#' @rdname predict.sail
#' @export
coef.sail <- function(object, s = NULL, ...) {
  stats::predict(object, s = s, type = "coefficients", ...)
}



#' @title Make predictions from a \code{cv.sail} object
#' @description This function makes predictions from a cross-validated sail
#'   model, using the stored "sail.fit" object, and the optimal value chosen for
#'   lambda.
#' @param object fitted \code{cv.sail} object
#' @inheritParams predict.sail
#' @param s Value(s) of the penalty parameter \code{lambda} at which predictions
#'   are required. Default is the value \code{s="lambda.1se"} stored on the CV
#'   \code{object}. Alternatively \code{s="lambda.min"} can be used. If \code{s}
#'   is numeric, it is taken as the value(s) of \code{lambda} to be used.
#' @param ... other arguments passed to \code{\link{predict.sail}}
#' @return The object returned depends the ... argument which is passed on to
#'   the predict method for \code{sail} objects.
#' @details This function makes it easier to use the results of cross-validation
#'   to make a prediction.
#' @examples
#' if(interactive()){
#' data("sailsim")
#' library(doParallel)
#' registerDoParallel(cores = 2)
#' f.basis <- function(i) splines::bs(i, degree = 3)
#' cvfit <- cv.sail(x = sailsim$x, y = sailsim$y, e = sailsim$e,
#'                   basis = f.basis, nfolds = 3, dfmax = 5, parallel = TRUE)
#' predict(cvfit) # predict at "lambda.1se"
#' predict(cvfit, s = "lambda.min") # predict at "lambda.min"
#' predict(cvfit, s = 0.5) # predict at specific value of lambda
#' predict(cvfit, type = "nonzero") # non-zero coefficients at lambda.1se
#'  }
#' @rdname predict.cv.sail
#' @seealso \code{\link{predict.sail}}
#' @export
predict.cv.sail <- function(object, newx, newe, s = c("lambda.1se", "lambda.min"), ...) {
  if (is.numeric(s)) {
    lambda <- s
  } else
  if (is.character(s)) {
    s <- match.arg(s)
    lambda <- object[[s]]
  }
  else {
    stop("Invalid form for s")
  }
  stats::predict(object$sail.fit, newx, newe, s = lambda, ...)
}


#' @inheritParams predict.cv.sail
#' @rdname predict.cv.sail
#' @export
coef.cv.sail <- function(object, s = c( "lambda.1se","lambda.min"), ...) {
  if (is.numeric(s)) {
    lambda <- s
  } else
  if (is.character(s)) {
    s <- match.arg(s)
    lambda <- object[[s]]
  }
  else {
    stop("Invalid form for s")
  }
  stats::coef(object$sail.fit, s = lambda, ...)
}


#' @title Print Method for \code{sail} object
#' @description Print a summary of the \code{sail} path at each step along the
#'   path.
#' @param x fitted \code{sail} object
#' @param digits significant digits in printout. Default: \code{max(3,
#'   getOption("digits") - 3)}
#' @param ... additional print arguments
#' @return OUTPUT_DESCRIPTION
#' @details The call that produced the object \code{x} is printed, followed by a
#'   five-column matrix with columns \code{df_main}, \code{df_interaction},
#'   \code{df_environment}, \code{\%Dev} and \code{Lambda}. The \code{df_}
#'   columns are the corresponding number of nonzero coefficients for main
#'   effects, interactions and exposure, respectively. \code{\%dev} is the
#'   percent deviance explained (relative to the null deviance). For
#'   \code{type="gaussian"} this is the r-squared.
#' @examples
#' if(interactive()){
#' data("sailsim")
#' f.basis <- function(i) splines::bs(i, degree = 3)
#' fit <- sail(x = sailsim$x, y = sailsim$y, e = sailsim$e,
#'             basis = f.basis, dfmax = 5, nlambda = 50)
#' fit
#'  }
#' @rdname print.sail
#' @seealso \code{\link{sail}}, \code{\link{cv.sail}}
#' @export
print.sail <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall: ", deparse(x$call), "\n\n")
  print(cbind(
    df_main = x$dfbeta,
    df_interaction = x$dfalpha,
    df_environment = x$dfenviron,
    `%Dev` = signif(x$dev.ratio, digits),
    Lambda = signif(x$lambda, digits)
  ))
}


#' @title Plot Method for \code{sail} object
#' @description Produces a coefficient profile plot of the coefficient paths for
#'   a fitted \code{sail} object. Both main effects and interactions (if
#'   present) are plotted.
#' @param x fitted \code{sail} object
#' @param type which type of predictors should be plotted. \code{type="both"}
#'   will plot the solution path for main effects and interactions,
#'   \code{type="main"} will only plot solution path of main effects (this also
#'   includes the exposure variable) and \code{type="interaction"} will only
#'   plot solution path for interaction effects Default: c("both", "main",
#'   "interaction"). Default: \code{type="both"}.
#' @param ... other graphical paramters passed to \code{plot}
#' @details A coefficient profile plot is produced
#' @return A plot is produced and nothing is returned
#' @examples
#' if(interactive()){
#' data("sailsim")
#' f.basis <- function(i) splines::bs(i, degree = 3)
#' fit <- sail(x = sailsim$x, y = sailsim$y, e = sailsim$e,
#'             basis = f.basis, dfmax = 5)
#' plot(fit)
#'  }
#' @rdname plot.sail
#' @seealso \code{\link{sail}}, \code{\link{cv.sail}}
#' @export
plot.sail <- function(x, type = c("both", "main", "interaction"), ...) {
  op <- graphics::par(no.readonly = TRUE)

  type <- match.arg(type)

  if (type != "main") {
    if (all(x$dfalpha == 0)) {
      warning("All interactions were estimated to be 0.\nPlotting solution path for main effects only")
      type <- "main"
    }
  }


  if (type == "main") {
    graphics::par(mar = 0.1 + c(4, 5, 2.5, 1))
    plotSailCoef(
      coefs = x$beta,
      environ = x$bE,
      lambda = x$lambda,
      df = x$dfbeta + x$dfenviron,
      group = x$group,
      dev = x$dev.ratio,
      vnames = x$vnames,
      ylab = "Main effects",
      ...
    )
  }

  if (type == "interaction") {
    graphics::par(mar = 0.1 + c(4, 5, 2.5, 1))
    plotSailCoef(
      coefs = x$alpha,
      lambda = x$lambda,
      df = x$dfalpha,
      group = x$group,
      dev = x$dev.ratio,
      vnames = x$vnames,
      ylab = "Interactions",
      ...
    )
  }


  if (type == "both") {
    # op <- graphics::par(
    #   mfrow = c(2, 1),
    #   mar = 0.1 + c(4.2, 4.0, 1, 1),
    #   oma = c(0, 1, 1, 0),
    #   cex.lab = 1.2, font.lab = 1.2, cex.axis = 1.2,
    #   cex.main = 1.2
    # )

    par(mfrow=c(2,1), tcl=-0.5, family="serif", omi=c(0.2,0.2,0,0),
        cex.lab = 1.2, font.lab = 1.2, cex.axis = 1.2,
        cex.main = 1.2)

    par(mar=c(2,4,2,3.2))
    plotSailCoef(
      coefs = x$beta,
      environ = x$bE,
      lambda = x$lambda,
      df = x$dfbeta + x$dfenviron,
      group = x$group,
      dev = x$dev.ratio,
      vnames = x$vnames,
      ylab = "Main effects",
      xlab = "", xaxt = "n",
      ...
    )

    par(mar=c(4,4,0,3.2))
    plotSailCoef(
      coefs = x$alpha,
      lambda = x$lambda,
      df = x$dfalpha,
      group = x$group,
      dev = x$dev.ratio,
      vnames = x$vnames,
      ylab = "Interactions",
      ...
    )

  }

  graphics::par(op)
}




#' @title Plot the cross-validation curve produced by \code{cv.sail}
#' @description Plots the cross-validation curve, and upper and lower standard
#'   deviation curves, as a function of the \code{lambda} values used.
#' @param x fitted \code{cv.sail} object
#' @param sign.lambda Either plot against \code{log(lambda)} (default) or its
#'   negative if \code{sign.lambda=-1}.
#' @param ... Other graphical parameters to plot
#' @return A plot is produced and nothing is returned
#' @details This is a port of \code{plot.cv.glmnet}
#' @examples
#' if(interactive()){
#' data("sailsim")
#' library(doParallel)
#' registerDoParallel(cores = 2)
#' f.basis <- function(i) splines::bs(i, degree = 3)
#' cvfit <- cv.sail(x = sailsim$x, y = sailsim$y, e = sailsim$e,
#'                   basis = f.basis, nfolds = 3, dfmax = 5, parallel = TRUE)
#' plot(cvfit)
#'  }
#' @rdname plot.cv.sail
#' @seealso \code{\link{sail}}, \code{\link{cv.sail}}
#' @references Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010).
#'   Regularization Paths for Generalized Linear Models via Coordinate Descent.
#'   Journal of Statistical Software, 33(1), 1-22.
#'   \url{http://www.jstatsoft.org/v33/i01/}.
#' @export
plot.cv.sail <- function(x, sign.lambda = 1, ...) {
  cvobj <- x
  xlab <- "log(Lambda)"
  if (sign.lambda < 0) xlab <- paste("-", xlab, sep = "")
  plot.args <- list(
    x = sign.lambda * log(cvobj$lambda), y = cvobj$cvm,
    ylim = range(cvobj$cvup, cvobj$cvlo), xlab = xlab, ylab = cvobj$name, type = "n"
  )
  new.args <- list(...)
  if (length(new.args)) plot.args[names(new.args)] <- new.args
  do.call("plot", plot.args)
  error.bars(sign.lambda * log(cvobj$lambda), cvobj$cvup, cvobj$cvlo, width = 0.01, col = "darkgrey")
  graphics::points(sign.lambda * log(cvobj$lambda), cvobj$cvm, pch = 20, col = "red")
  graphics::axis(side = 3, at = sign.lambda * log(cvobj$lambda), labels = paste(cvobj$nz), tick = FALSE, line = 0)
  graphics::abline(v = sign.lambda * log(cvobj$lambda.min), lty = 3)
  graphics::abline(v = sign.lambda * log(cvobj$lambda.1se), lty = 3)
  invisible()
}
