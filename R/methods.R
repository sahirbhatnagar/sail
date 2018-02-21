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

predict.sail <- function(object, newx, newe, s = NULL,
                            type = c("link", "response", "coefficients",
                                     "nonzero", "class"), ...) {

  # object = fit
  # type = "coefficients"
  #==================

  type = match.arg(type)

  if (missing(newx)) {
    if (!match(type, c("coefficients", "nonzero"), FALSE))
      newx <- object$design
  } else if (!missing(newx) & missing(newe)) {
    stop("newe is missing. please supply the vector of the environment variable.")
  } else if (!missing(newx) & !missing(newe)) {
    newx <- design_sail(x = newx, e = newe, nvars = object$nvars,
                        vnames = object$vnames, df = object$df, degree = object$degree)$design
  }

  a0 <- t(as.matrix(object$a0))
  rownames(a0) = "(Intercept)"
  # this has only values for which lambda did converge
  nbeta <- rbind(a0, object$beta, E = object$bE, object$alpha)

  if (!is.null(s)) {
    vnames <- dimnames(nbeta)[[1]]
    dimnames(nbeta) <- list(NULL, NULL)
    lambda <- object$lambda
    lamlist <- lambda.interp(lambda, s)
    if(length(s) == 1) {
      nbeta = nbeta[, lamlist$left, drop=FALSE] * lamlist$frac +
        nbeta[, lamlist$right, drop=FALSE] * (1 - lamlist$frac)
    } else {
      nbeta = nbeta[, lamlist$left, drop=FALSE] %*% diag(lamlist$frac) +
        nbeta[, lamlist$right, drop=FALSE] %*% diag(1 - lamlist$frac)
    }
    dimnames(nbeta) <- list(vnames, paste(seq(along = s)))
  }

  if(type == "coefficients") return(nbeta)

  if(type == "nonzero") return(nonzero(nbeta[-1,,drop=FALSE],bystep=TRUE))


  # this is used by the cv_lspath function to calculate predicted values
  # which will subsequently be used for calculating MSE for each fold
  if (type == "link") {
    nfit <- as.matrix(methods::cbind2(1, newx) %*% nbeta)
    return(nfit)
  }

}


predict.cv.sail <- function(object,newx, newe, s=c("lambda.1se","lambda.min"),...){
  if(is.numeric(s)) lambda=s
  else
    if(is.character(s)){
      s=match.arg(s)
      lambda=object[[s]]
    }
  else stop("Invalid form for s")
  predict(object$sail.fit, newx, newe, s=lambda, ...)
}


#' Get coefficients from a "sail" object
#'
#' @rdname predict.sail
#' @export

coef.sail=function(object,s=NULL,exact=FALSE,...)
  predict(object,s=s,type="coefficients",exact=exact,...)


#' Make predictions from a "cv.sail" object
#'
#' @param object object of class cv.sail from cv.sail function
#' @param s Value(s) of the penalty parameter lambda at which predictions are
#'   required. Default is the value \code{s="lambda.1se"} stored on the cv.sail
#'   object. Alternatively \code{s="lambda.min"} can be used.
#' @export

coef.cv.sail=function(object,s=c("lambda.1se","lambda.min"),...){
  if(is.numeric(s))lambda=s
  else
    if(is.character(s)){
      s=match.arg(s)
      lambda=object[[s]]
    }
  else stop("Invalid form for s")
  coef(object$sail.fit,s=lambda,...)
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

plot.sail <- function(x, type = c("both","main","interaction"), ...) {

  op <- par(no.readonly=TRUE)

  type <- match.arg(type)

  if (type != "main") {
    if (all(x$dfalpha == 0)) {
      warning("All interactions were estimated to be 0.\nPlotting solution path for main effects only")
      type <- "main"
    }
  }


  if (type == "main") {

    par(mar = 0.1 + c(4, 5, 2.5, 1))
    plotSailCoef(coefs = x$beta,
                 environ = x$bE,
                 lambda = x$lambda,
                 df = x$dfbeta + x$dfenviron,
                 group = x$group,
                 dev = x$dev.ratio,
                 vnames = x$vnames,
                 ylab = "Main effects",
                 ...)

  }

  if (type == "interaction") {

    par(mar = 0.1 + c(4, 5, 2.5, 1))
    plotSailCoef(coefs = x$alpha,
                 lambda = x$lambda,
                 df = x$dfalpha,
                 group = x$group,
                 dev = x$dev.ratio,
                 vnames = x$vnames,
                 ylab = "Interactions",
                 ...)

  }


  if (type == "both") {

     op <- par(mfrow = c(2, 1),
              mar = 0.1 + c(4.2, 4.0, 1, 1),
              oma = c(0, 1, 1, 0),
              cex.lab = 1.2, font.lab = 1.2, cex.axis = 1.2,
              cex.main = 1.2)


    plotSailCoef(coefs = x$beta,
                 environ = x$bE,
                 lambda = x$lambda,
                 df = x$dfbeta + x$dfenviron,
                 group = x$group,
                 dev = x$dev.ratio,
                 vnames = x$vnames,
                 ylab = "Main effects",
                 xlab = "", xaxt='n',
                 ...)

    plotSailCoef(coefs = x$alpha,
                 lambda = x$lambda,
                 df = x$dfalpha,
                 group = x$group,
                 dev = x$dev.ratio,
                 vnames = x$vnames,
                 ylab = "Interactions",
                 ...)

    par(op)

  }

  par(op)

}



























