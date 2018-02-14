#' Make predictions from a sail object
#'
#' @description this function only works for tuning parameter values defined by
#'   the sail_multiple_faster function. The interpolation feature is not working
#'   yet
#' @param s index of tuning parameter. Must be a character and an element of
#'   "s1","s2",...."s100", where "s100" is the index of the last pair of tuning
#'   parameters. Default is \code{NULL}
#' @param object Fitted sail model object
#' @param type Type of prediction required. Type "link" gives the fitted values.
#'   Type "response" for "gaussian" type "response" is equivalent to type
#'   "link". Type "coefficients" computes the coefficients at the requested
#'   values for s. Type "nonzero" returns a list of the indices of the nonzero
#'   coefficients for each value of s.
#' @export

predict.sail <- function(object, newx, s = NULL,
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

  # this is the default returned by coef.sail i.e. any object of class sail
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

#' Get coefficients from a "sail" object
#'
#' @rdname predict.sail
#' @export

coef.sail <- function(object, s = NULL) {
  predict(object, s = s, type = "coefficients")
}


#' Make predictions from a "cv.sail" object
#'
#' @param object object of class cv.sail from cv.sail function
#' @param s Value(s) of the penalty parameter lambda at which predictions are
#'   required. Default is the value \code{s="lambda.1se"} stored on the cv.sail
#'   object. Alternatively \code{s="lambda.min"} can be used.
#' @export

coef.cv.sail <- function(object, s = c("lambda.1se", "lambda.min"), ...) {

  if (is.numeric(s) || s %ni% c("lambda.1se", "lambda.min")) stop("s must be in lambda.1se or lambda.min")

  s <- match.arg(s)

  lambda <- switch(s,
                   lambda.min = object$lambda.min.name,
                   lambda.1se = object$lambda.1se.name
  )

  coef(object$sail.fit, s = lambda, ...)
}


#' Print Method for sail function
#'
#' @description print method for sail function
#' @export

print.sail <- function (x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall: ", deparse(x$call), "\n\n")
  print(cbind(df_main = x$dfbeta, df_interaction = x$dfalpha,
              df_environment = x$dfenviron,
              `%Dev` = signif(x$dev.ratio, digits),
              Lambda = signif(x$lambda, digits)))
}




#' Plot Method for sail function
#'
#' @description plot method for sail function
#' @export

plot.sail <- function(x, xvar = c("norm", "lambda", "dev"), label = T,
                      ...) {
  xvar = match.arg(xvar)
  plotCoefSail(x$beta,
               lambda = x$lambda.beta,
               df = x$dfbeta,
               dev = x$dev.ratio,
               label = label,
               xvar = xvar, ...)
}



plot.sail <- function(x, xvar = c("norm", "lambda", "dev"), label = T,
                      ...) {
  # xvar = match.arg(xvar)
  # plotCoefSail(x$beta,
  #              lambda = x$lambda.beta,
  #              df = x$dfbeta,
  #              dev = x$dev.ratio,
  #              label = label,
  #              xvar = xvar, ...)

  plot.grpreg(x = x, ...)

}
