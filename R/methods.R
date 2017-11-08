#' Make predictions from a shim object
#'
#' @description this function only works for tuning parameter values defined by
#'   the shim_multiple_faster function. The interpolation feature is not working
#'   yet
#' @param s index of tuning parameter. Must be a character and an element of
#'   "s1","s2",...."s100", where "s100" is the index of the last pair of tuning
#'   parameters. Default is \code{NULL}
#' @param object Fitted shim model object
#' @param type Type of prediction required. Type "link" gives the fitted values.
#'   Type "response" for "gaussian" type "response" is equivalent to type
#'   "link". Type "coefficients" computes the coefficients at the requested
#'   values for s. Type "nonzero" returns a list of the indices of the nonzero
#'   coefficients for each value of s.
#' @export

predict.funshim <- function(object, newx, s = NULL,
                            type = c("link", "response", "coefficients",
                                     "nonzero", "class")) {

  # object = fit
  # type = "coefficients"
  #==================

  type = match.arg(type)

  if (missing(newx)) {
    if (!match(type, c("coefficients", "nonzero"), FALSE))
      stop("You need to supply a value for 'newx'")
  }

  a0 = t(as.matrix(object$b0))
  rownames(a0) = "(Intercept)"
  # this includes tuning parameters pairs that didnt converge
  nbeta = rbind(a0, object$beta, object$alpha)
  nbeta@Dimnames <- list(X = c("(Intercept)", object$main.effect.names,
                               object$interaction.names),
                         Y = paste0("s",seq_len(object$nlambda)))

  # this is the default returned by coef.shim i.e. any object of class shim
  # it will return all tuning parameters (including those that didnt converge)
  if (type == "coefficients" && is.null(s)) {
    return(nbeta)
  }

  if (type == "coefficients" && !is.null(s)) {
    return(nbeta[ , s, drop = F])
  }

  if (type == "nonzero") {
    nbeta = rbind(a0, object$beta, object$alpha)
    return(list(main = nonzero(nbeta[object$main.effect.names, , drop = FALSE], bystep = TRUE),
                interaction = nonzero(nbeta[object$interaction.names, , drop = FALSE], bystep = TRUE)))
  }

  if (inherits(newx, "sparseMatrix")) {
    newx = as(newx, "dgCMatrix")
  }

  # this is used by the cv_lspath function to calculate predicted values
  # which will subsequently be used for calculating MSE for each fold
  if (type == "link") {

    nfit = as.matrix(cbind2(1, newx) %*% nbeta)

    return(nfit)
  }

}

#' Get coefficients from a "shim" object
#'
#' @rdname predict.shim
#' @export

coef.funshim <- function(object, s = NULL) {
  predict(object, s = s, type = "coefficients")
}


#' Make predictions from a "cv.shim" object
#'
#' @param object object of class cv.shim from cv.shim function
#' @param s Value(s) of the penalty parameter lambda at which predictions are
#'   required. Default is the value \code{s="lambda.1se"} stored on the cv.shim
#'   object. Alternatively \code{s="lambda.min"} can be used.
#' @export

coef.cv.funshim <- function(object, s = c("lambda.1se", "lambda.min"), ...) {

  if (is.numeric(s) || s %ni% c("lambda.1se", "lambda.min")) stop("s must be in lambda.1se or lambda.min")

  s <- match.arg(s)

  lambda <- switch(s,
                   lambda.min = object$lambda.min.name,
                   lambda.1se = object$lambda.1se.name
  )

  coef(object$funshim.fit, s = lambda, ...)
}


#' Print Method for shim function
#'
#' @description print method for shim function
#' @export

print.shim <- function (x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall: ", deparse(x$call), "\n\n")
  print(cbind(DfBeta = x$dfbeta, DfAlpha = x$dfalpha,
              `%Dev` = signif(x$dev.ratio, digits),
              LambdaBeta = signif(x$lambda.beta, digits),
              LambdaGamma = signif(x$lambda.gamma, digits)))
}


print.funshim <- function (x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall: ", deparse(x$call), "\n\n")
  print(cbind(DfBeta = x$dfbeta, DfAlpha = x$dfalpha,
              `%Dev` = signif(x$dev.ratio, digits),
              LambdaBeta = signif(x$lambda.beta, digits),
              LambdaGamma = signif(x$lambda.gamma, digits)))
}


#' Plot Method for shim function
#'
#' @description plot method for shim function
#' @export

plot.shim <- function(x, xvar = c("norm", "lambda", "dev"), label = T,
                      ...) {
  xvar = match.arg(xvar)
  plotCoefShim(x$beta,
               lambda = x$lambda.beta,
               df = x$dfbeta,
               dev = x$dev.ratio,
               label = label,
               xvar = xvar, ...)
}



plot.funshim <- function(x, xvar = c("norm", "lambda", "dev"), label = T,
                      ...) {
  # xvar = match.arg(xvar)
  # plotCoefShim(x$beta,
  #              lambda = x$lambda.beta,
  #              df = x$dfbeta,
  #              dev = x$dev.ratio,
  #              label = label,
  #              xvar = xvar, ...)

  plot.grpreg(x = x, ...)

}
