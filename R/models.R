#' Cross-validation for sail
#'
#' @description Does k-fold cross-validation for sail and determines the optimal
#'   pair of tuning parameters (\eqn{\lambda_\beta} and \eqn{\lambda_\gamma})
#'
#' @inheritParams sail
#' @param parallel If \code{TRUE}, use parallel \code{foreach} to fit each fold.
#'   Must register parallel before hand using the
#'   \code{\link[doMC]{registerDoMC}} function from the \code{doMC} package. See
#'   the example below for details.
#' @param type.measure loss to use for cross-validation. Currently only 1
#'   option. The default is \code{type.measure="mse"}, which uses squared-error
#'   for gaussian models
#' @param nfolds number of folds - default is 10. Although nfolds can be as
#'   large as the sample size (leave-one-out CV), it is not recommended for
#'   large datasets. Smallest value allowable is \code{nfolds=3}
#' @details The function runs sail nfolds+1 times; the first to get the tuning
#'   parameter sequences, and then the remainder to compute the fit with each of
#'   the folds omitted. The error is accumulated, and the average error and
#'   standard deviation over the folds is computed. Note also that the results
#'   of cv.sail are random, since the folds are selected at random using the
#'   \code{\link{createfolds}} function. Users can reduce this randomness by
#'   running cv.sail many times, and averaging the error curves.
#' @author Sahir Bhatnagar
#'
#'   Maintainer: Sahir Bhatnagar \email{sahir.bhatnagar@@mail.mcgill.ca}
#' @export

cv.sail <- function(x, y, e, ...,
                    weights,
                    lambda = NULL,
                    type.measure = c("mse", "deviance", "class", "auc", "mae"),
                    nfolds = 10, foldid, grouped = TRUE, keep = FALSE, parallel = FALSE) {
  if (missing(type.measure)) type.measure <- "default" else type.measure <- match.arg(type.measure)
  if (!is.null(lambda) && length(lambda) < 2) {
    stop("Need more than one value of lambda for cv.sail")
  }
  N <- nrow(x)
  if (missing(weights)) {
    weights <- rep(1, N)
  } else {
    weights <- as.double(weights)
  }
  y <- drop(y)
  sail.call <- match.call(expand.dots = TRUE)
  which <- match(c(
    "type.measure", "nfolds", "foldid", "grouped",
    "keep"
  ), names(sail.call), F)
  if (any(which)) {
    sail.call <- sail.call[-which]
  }
  sail.call[[1]] <- as.name("sail")
  sail.object <- sail(
    x = x, y = y, e = e, ...,
    weights = weights, lambda = lambda
  )

  sail.object$call <- sail.call

  ### Next line is commented out so each call generates its own lambda sequence
  # lambda <- sail.object$lambda

  # nz = sapply(predict(sail.object, type = "nonzero"), length)
  nz <- sapply(sail.object$active, length)
  if (missing(foldid)) foldid <- sample(rep(seq(nfolds), length = N)) else nfolds <- max(foldid)
  if (nfolds < 3) {
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  }
  outlist <- as.list(seq(nfolds))
  if (parallel) {
    outlist <- foreach(i = seq(nfolds), .packages = c("sail")) %dopar% {
      which <- foldid == i

      if (length(dim(y)) > 1) {
        y_sub <- y[!which, ]
      } else {
        y_sub <- y[!which]
      }

      if (length(dim(e)) > 1) {
        e_sub <- e[!which, ]
      } else {
        e_sub <- e[!which]
      }

      sail(
        x = x[!which, , drop = FALSE], y = y_sub, e = e_sub, ...,
        lambda = lambda, weights = weights[!which]
      )
    }
  } else {
    for (i in seq(nfolds)) {
      which <- foldid == i

      if (is.matrix(y)) {
        y_sub <- y[!which, ]
      } else {
        y_sub <- y[!which]
      }

      if (length(dim(e)) > 1) {
        e_sub <- e[!which, ]
      } else {
        e_sub <- e[!which]
      }

      outlist[[i]] <- sail(
        x = x[!which, , drop = FALSE], y = y_sub, e = e_sub, ...,
        lambda = lambda, weights = weights[!which]
      )
    }
  }

  fun <- paste("cv", class(sail.object)[[1]], sep = ".")
  lambda <- sail.object$lambda
  cvstuff <- do.call(fun, list(
    outlist, lambda, x, y, e, weights,
    foldid, type.measure, grouped, keep
  ))
  cvm <- cvstuff$cvm
  cvsd <- cvstuff$cvsd
  nas <- is.na(cvsd)
  if (any(nas)) {
    lambda <- lambda[!nas]
    cvm <- cvm[!nas]
    cvsd <- cvsd[!nas]
    nz <- nz[!nas]
  }
  cvname <- cvstuff$name
  out <- list(
    lambda = lambda, cvm = cvm, cvsd = cvsd, cvup = cvm +
      cvsd, cvlo = cvm - cvsd, nzero = nz, name = cvname,
    sail.fit = sail.object
  )
  if (keep) {
    out <- c(out, list(fit.preval = cvstuff$fit.preval, foldid = foldid))
  }
  lamin <- if (cvname == "AUC") {
    getmin(lambda, -cvm, cvsd)
  } else {
    getmin(lambda, cvm, cvsd)
  }
  obj <- c(out, as.list(lamin))
  class(obj) <- "cv.sail"
  obj
}
