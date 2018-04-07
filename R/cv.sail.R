#' Cross-validation for sail
#'
#' @description Does k-fold cross-validation for sail and determines the optimal
#'   tuning parameter \eqn{\lambda}.
#' @inheritParams sail
#' @param ... other arguments that can be passed to \code{\link{sail}}
#' @param lambda Optional user-supplied lambda sequence; default is NULL, and
#'   \code{\link{sail}} chooses its own sequence
#' @param type.measure loss to use for cross-validation. Currently only 3
#'   options are implemented. The default is \code{type.measure="deviance"},
#'   which uses squared-error for gaussian models (and is equivalent to
#'   \code{type.measure="mse"}) there). \code{type.measure="mae"} (mean absolute
#'   error) can also be used which measures the absolute deviation from the
#'   fitted mean to the response (\eqn{|y-\hat{y}|}).
#' @param nfolds number of folds. Although \code{nfolds} can be as large as the
#'   sample size (leave-one-out CV), it is not recommended for large datasets.
#'   Smallest value allowable is \code{nfolds=3}. Default: 10
#' @param foldid an optional vector of values between 1 and \code{nfold}
#'   identifying what fold each observation is in. If supplied,\code{nfold} can
#'   be missing. Often used when wanting to tune the second tuning parameter
#'   (\eqn{\alpha}) as well (see details).
#' @param grouped This is an experimental argument, with default \code{TRUE},
#'   and can be ignored by most users. This refers to computing \code{nfolds}
#'   separate statistics, and then using their mean and estimated standard error
#'   to describe the CV curve. If \code{grouped=FALSE}, an error matrix is built
#'   up at the observation level from the predictions from the \code{nfold}
#'   fits, and then summarized (does not apply to \code{type.measure="auc"}).
#'   Default: TRUE.
#' @param keep If \code{keep=TRUE}, a \emph{prevalidated} array is returned
#'   containing fitted values for each observation and each value of
#'   \code{lambda}. This means these fits are computed with this observation and
#'   the rest of its fold omitted. The \code{folid} vector is also returned.
#'   Default: FALSE
#' @param parallel If \code{TRUE}, use parallel \code{foreach} to fit each fold.
#'   Must register parallel before hand using the
#'   \code{\link[doMC]{registerDoMC}} function from the \code{doMC} package. See
#'   the example below for details. Default: FALSE
#' @return an object of class \code{"cv.sail"} is returned, which is a list with
#'   the ingredients of the cross-validation fit. \describe{ \item{lambda}{the
#'   values of converged \code{lambda} used in the fits.} \item{cvm}{The mean
#'   cross-validated error - a vector of length \code{length(lambda)}.}
#'   \item{cvsd}{estimate of standard error of \code{cvm}.} \item{cvup}{upper
#'   curve = \code{cvm+cvsd}.} \item{cvlo}{lower curve = \code{cvm-cvsd}.}
#'   \item{nzero}{number of non-zero coefficients at each \code{lambda}. This is
#'   the sum of the total non-zero main effects and interactions. Note that when
#'   \code{expand=TRUE}, we only count a variable once in the calculation of
#'   \code{nzero}, i.e., if a variable is expanded to three columns, then this
#'   is only counted once even though all three coefficients are estimated to be
#'   non-zero} \item{name}{a text string indicating type of measure (for
#'   plotting purposes).} \item{sail.fit}{a fitted \code{sail} object for the full
#'   data.} \item{lambda.min}{value of \code{lambda} that gives minimum
#'   \code{cvm}.} \item{lambda.1se}{largest value of \code{lambda} such that
#'   error is within 1 standard error of the minimum.} \item{fit.preval}{if
#'   \code{keep=TRUE}, this is the array of prevalidated fits. Some entries can
#'   be \code{NA}, if that and subsequent values of \code{lambda} are not
#'   reached for that fold} \item{foldid}{if \code{keep=TRUE}, the fold
#'   assignments used}}
#' @details The function runs \code{\link{sail}} \code{nfolds}+1 times; the
#'   first to get the \code{lambda} sequence, and then the remainder to compute
#'   the fit with each of the folds omitted. Note that a new lambda sequence is
#'   computed for each of the folds and then we use the \code{predict} method to
#'   get the solution path at each value of the original lambda sequence. The
#'   error is accumulated, and the average error and standard deviation over the
#'   folds is computed. Note that \code{cv.sail} does NOT search for values for
#'   \code{alpha}. A specific value should be supplied, else \code{alpha=0.5} is
#'   assumed by default. If users would like to cross-validate \code{alpha} as
#'   well, they should call \code{cv.sail} with a pre-computed vector
#'   \code{foldid}, and then use this same fold vector in separate calls to
#'   \code{cv.sail} with different values of \code{alpha}. Note also that the
#'   results of \code{cv.sail} are random, since the folds are selected at
#'   random. Users can reduce this randomness by running \code{cv.sail} many
#'   times, and averaging the error curves.
#' @note The skeleton of this function and the documentation were taken straight
#'   from the \code{glmnet} package. See references for details.
#' @references Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010).
#'   Regularization Paths for Generalized Linear Models via Coordinate Descent.
#'   Journal of Statistical Software, 33(1), 1-22.
#'   \url{http://www.jstatsoft.org/v33/i01/}.
#' @references Bhatnagar SR, Yang Y, Greenwood CMT. Sparse additive interaction
#'   models with the strong heredity property (2018+). Preprint.
#' @seealso \code{\link[splines]{bs}} \code{\link{sail}}
#' @examples
#' \dontrun{
#' if(interactive()){
#' f.basis <- function(i) splines::bs(i, degree = 5)
#' data("sailsim")
#' cvfit <- cv.sail(x = sailsim$x, y = sailsim$y, e = sailsim$e,
#'                  basis = f.basis, nfolds = 10)
#'
#' # Parallel
#' library(doMC)
#' registerDoMC(cores = 4)
#' cvfit <- cv.sail(x = sailsim$x, y = sailsim$y, e = sailsim$e,
#'                  parallel = TRUE, nlambda = 100, nfolds = 10)
#' # plot cross validated curve
#' plot(cvfit)
#' # plot solution path
#' plot(cvfit$sail.fit)
#'
#' # solution at lambda.min
#' coef(cvfit, s = "lambda.min")
#' # solution at lambda.1se
#' coef(cvfit, s = "lambda.1se")
#' # non-zero coefficients
#' coef(cvfit, s = "lambda.min")[nonzero(coef(cvfit, s = "lambda.min")),,drop=FALSE]
#'
#' # predicted response
#' predict(cvfit, s = "lambda.min")
#' predict(cvfit, s = "lambda.1se")
#' # predict response at any value for lambda
#' predict(cvfit, s = 2.1)
#'
#' # predict response for new data set
#' newx <- sailsim$x * 1.10
#' newe <- sailsim$e * 2
#' predict(cvfit, newx = newx, newe = newe, s = "lambda.min")
#'  }
#' }
#' @rdname cv.sail
#' @export
cv.sail <- function(x, y, e, ...,
                    weights,
                    lambda = NULL,
                    type.measure = c("mse", "deviance", "class", "auc", "mae"),
                    nfolds = 10, foldid, grouped = TRUE, keep = FALSE, parallel = FALSE) {

  if (!requireNamespace("foreach", quietly = TRUE)) {
    stop("Package \"foreach\" needed for this function to work in parallel. Please install it.",
         call. = FALSE)
  }

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
  sail.object <- sail::sail(
    x = x, y = y, e = e, ...,
    weights = weights, lambda = lambda
  )

  sail.object$call <- sail.call

  ### Next line is commented out so each call generates its own lambda sequence
  # lambda <- sail.object$lambda

  # nz = sapply(predict(sail.object, type = "nonzero"), length)
  nz <- sapply(sail.object$active, length)
  # if (missing(foldid)) foldid <- sample(rep(seq(nfolds), length = N)) else nfolds <- max(foldid)
  if (missing(foldid)) foldid <- createfolds(y = y, k = nfolds) else nfolds <- max(foldid)
  if (nfolds < 3) {
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  }
  outlist <- as.list(seq(nfolds))
  if (parallel) {
    outlist <- foreach::foreach(i = seq(nfolds), .packages = c("sail")) %dopar% {
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

      sail::sail(
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

      outlist[[i]] <- sail::sail(
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


#' Create CV Folds
#'
#' @description \code{createfolds} splits the data into \code{k} groups. Taken
#'   from the \code{caret} package (see references for details)
#' @param y vector of response
#' @param an integer for the number of folds.
#' @return A vector of CV fold ID's for each observation in \code{y}
#' @details For numeric y, the sample is split into groups sections based on
#'   percentiles and sampling is done within these subgroups
#' @references \url{http://topepo.github.io/caret/splitting.html}
createfolds <- function(y, k = 10, list = FALSE, returnTrain = FALSE) {
  if (class(y)[1] == "Surv") {
    y <- y[, "time"]
  }
  if (is.numeric(y)) {
    cuts <- floor(length(y) / k)
    if (cuts < 2) {
      cuts <- 2
    }
    if (cuts > 5) {
      cuts <- 5
    }
    breaks <- unique(stats::quantile(y, probs = seq(0, 1, length = cuts)))
    y <- cut(y, breaks, include.lowest = TRUE)
  }
  if (k < length(y)) {
    y <- factor(as.character(y))
    numInClass <- table(y)
    foldVector <- vector(mode = "integer", length(y))
    for (i in 1:length(numInClass)) {
      min_reps <- numInClass[i] %/% k
      if (min_reps > 0) {
        spares <- numInClass[i] %% k
        seqVector <- rep(1:k, min_reps)
        if (spares > 0) {
          seqVector <- c(seqVector, sample(1:k, spares))
        }
        foldVector[which(y == names(numInClass)[i])] <- sample(seqVector)
      }
      else {
        foldVector[which(y == names(numInClass)[i])] <- sample(1:k,
                                                               size = numInClass[i]
        )
      }
    }
  }
  else {
    foldVector <- seq(along = y)
  }
  if (list) {
    out <- split(seq(along = y), foldVector)
    names(out) <- paste("Fold", gsub(" ", "0", format(seq(along = out))),
                        sep = ""
    )
    if (returnTrain) {
      out <- lapply(out, function(data, y) y[-data], y = seq(along = y))
    }
  }
  else {
    out <- foldVector
  }
  out
}


#' Compute cross validation error
#'
#' @description functions used to calculate cross validation error and used by
#'   the \code{\link{cv.sail}} function
#'
#' @param outlist list of cross validated fitted models. List is of length equal
#'   to \code{nfolds} argument in \code{\link{cv.sail}} function
#' @param foldid numeric vector indicating which fold each observation belongs
#'   to
#' @inheritParams sail
#' @rdname cv.lspath
#' @seealso \code{\link{cv.sail}}
#' @details The output of the \code{cv.lspath} function only returns values for
#'   those tuning paramters that DID converge
#' @details Many have been taken verbatim from the \code{glmnet} package.
#' @references Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010).
#'   Regularization Paths for Generalized Linear Models via Coordinate Descent.
#'   Journal of Statistical Software, 33(1), 1-22.
#'   \url{http://www.jstatsoft.org/v33/i01/}.
cv.lspath <- function(outlist, lambda, x, y, e, weights,
                      foldid, type.measure, grouped, keep = FALSE) {
  typenames <- c(
    deviance = "Mean-Squared Error", mse = "Mean-Squared Error",
    mae = "Mean Absolute Error"
  )
  if (type.measure == "default") {
    type.measure <- "mse"
  }
  if (!match(type.measure, c("mse", "mae", "deviance"), FALSE)) {
    warning("Only 'mse', 'deviance' or 'mae'  available for Gaussian models; 'mse' used")
    type.measure <- "mse"
  }

  mlami <- max(sapply(outlist, function(obj) min(obj$lambda)))
  which_lam <- lambda >= mlami
  predmat <- matrix(NA, length(y), length(lambda))
  nfolds <- max(foldid)
  nlams <- double(nfolds)
  for (i in seq(nfolds)) {
    which <- foldid == i
    fitobj <- outlist[[i]]
    # fitobj$offset = FALSE
    preds <- predict(fitobj, newx = x[which, , drop = FALSE], newe = e[which], s = lambda[which_lam])
    nlami <- sum(which_lam)
    predmat[which, seq(nlami)] <- preds
    nlams[i] <- nlami
  }
  N <- length(y) - apply(is.na(predmat), 2, sum)
  cvraw <- switch(type.measure, mse = (y - predmat)^2, deviance = (y - predmat)^2, mae = abs(y - predmat))
  if ((length(y) / nfolds < 3) && grouped) {
    warning("Option grouped=FALSE enforced in cv.sail, since < 3 observations per fold",
            call. = FALSE
    )
    grouped <- FALSE
  }
  if (grouped) {
    cvob <- cvcompute(cvraw, weights, foldid, nlams)
    cvraw <- cvob$cvraw
    weights <- cvob$weights
    N <- cvob$N
  }
  cvm <- apply(cvraw, 2, stats::weighted.mean, w = weights, na.rm = TRUE)
  cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, stats::weighted.mean,
                     w = weights, na.rm = TRUE
  ) / (N - 1))
  out <- list(cvm = cvm, cvsd = cvsd, name = typenames[type.measure])
  if (keep) {
    out$fit.preval <- predmat
  }
  out
}

#' @describeIn cv.lspath Computations for crossvalidation error
#' @note taken verbatim from glmnet
cvcompute <- function(mat, weights, foldid, nlams) {
  ### Computes the weighted mean and SD within folds, and hence the se of the mean
  wisum <- tapply(weights, foldid, sum)
  nfolds <- max(foldid)
  outmat <- matrix(NA, nfolds, ncol(mat))
  good <- matrix(0, nfolds, ncol(mat))
  mat[is.infinite(mat)] <- NA # just in case some infinities crept in
  for (i in seq(nfolds)) {
    mati <- mat[foldid == i, , drop = FALSE]
    wi <- weights[foldid == i]
    outmat[i, ] <- apply(mati, 2, stats::weighted.mean, w = wi, na.rm = TRUE)
    good[i, seq(nlams[i])] <- 1
  }
  N <- apply(good, 2, sum)
  list(cvraw = outmat, weights = wisum, N = N)
}

#' @describeIn cv.lspath get lambda.min and lambda.1se
#' @note taken verbatim from glmnet
getmin <- function(lambda, cvm, cvsd) {
  cvmin <- min(cvm, na.rm = TRUE)
  idmin <- cvm <= cvmin
  lambda.min <- max(lambda[idmin], na.rm = TRUE)
  idmin <- match(lambda.min, lambda)
  semin <- (cvm + cvsd)[idmin]
  idmin <- cvm <= semin
  lambda.1se <- max(lambda[idmin], na.rm = TRUE)
  list(lambda.min = lambda.min, lambda.1se = lambda.1se)
}

#' @describeIn cv.lspath Interpolation function.
#' @note taken verbatim from glmnet
lambda.interp=function(lambda,s){
  ### lambda is the index sequence that is produced by the model
  ### s is the new vector at which evaluations are required.
  ### the value is a vector of left and right indices, and a vector of fractions.
  ### the new values are interpolated bewteen the two using the fraction
  ### Note: lambda decreases. you take:
  ### sfrac*left+(1-sfrac*right)

  if(length(lambda)==1){# degenerate case of only one lambda
    nums=length(s)
    left=rep(1,nums)
    right=left
    sfrac=rep(1,nums)
  }
  else{
    s[s > max(lambda)] = max(lambda)
    s[s < min(lambda)] = min(lambda)
    k=length(lambda)
    sfrac <- (lambda[1]-s)/(lambda[1] - lambda[k])
    lambda <- (lambda[1] - lambda)/(lambda[1] - lambda[k])
    coord <- approx(lambda, seq(lambda), sfrac)$y
    left <- floor(coord)
    right <- ceiling(coord)
    sfrac=(sfrac-lambda[right])/(lambda[left] - lambda[right])
    sfrac[left==right]=1
    sfrac[abs(lambda[left]-lambda[right])<.Machine$double.eps]=1

  }
  list(left=left,right=right,frac=sfrac)
}


