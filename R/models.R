#' Fit Strong Heredity Interaction Model
#'
#' @description function to fit the Strong Heredity Interaction Model for a
#'   sequence of tuning parameters. This is a penalized regression method that
#'   ensures the interaction term is non-zero only if its corresponding
#'   main-effects are non-zero.
#' @param x Design matrix of dimension \code{n x q}, where \code{n} is the
#'   number of subjects and q is the total number of variables; each row is an
#'   observation vector. This must include all main effects and interactions as
#'   well, with column names corresponding to the names of the main effects
#'   (e.g. \code{x1, x2, E}) and their interactions (e.g. \code{x1:E, x2:E}).
#'   All columns should be scaled to have mean 0 and variance 1; this is done
#'   internally by the \code{\link{sail}} function.
#' @param y response variable. For \code{family="gaussian"} should be a 1 column
#'   matrix or numeric vector. For \code{family="binomial"}, if the response is
#'   a vector it can be numeric with 0 for failure and 1 for success, or a
#'   factor with the first level representing "failure" and the second level
#'   representing "success". Alternatively, For binomial logistic regression,
#'   the response can be a matrix where the first column is the number of
#'   "successes" and the second column is the number of "failures".
#' @param family response type. see \code{y} for details. Currently only
#' \code{family = "gaussian"} is implemented.
#' @param main.effect.names character vector of main effects names. MUST be
#'   ordered in the same way as the column names of \code{x}. e.g. if the column
#'   names of \code{x} are \code{"x1","x2"} then \code{main.effect.names =
#'   c("x1","x2")}
#' @param interaction.names character vector of interaction names. MUST be
#'   separated by a colon (e.g. x1:x2), AND MUST be ordered in the same way as
#'   the column names of \code{x}
#' @param lambda.beta sequence of tuning parameters for the main effects. If
#'   \code{NULL} (default), this function will automatically calculate a
#'   sequence using the \code{\link{sail_once}} function which will be over a
#'   grid of tuning parameters for gamma as well. If the user specifies a
#'   sequence then this function will not automatically perform the serach over
#'   a grid. You will need to create the grid yourself e.g. repeat the
#'   lambda.gamma for each value of lambda.beta
#' @param lambda.gamma sequence of tuning parameters for the interaction
#'   effects. Default is \code{NULL} which means this function will
#'   automatically calculate a sequence of tuning paramters. See
#'   \code{\link{sail_once}} for details on how this sequence is calculated.
#' @param nlambda.gamma number of tuning parameters for gamma. This needs to be
#'   specified even for user defined inputs
#' @param nlambda.beta number of tuning parameters for beta. This needs to be
#'   specified even for user defined inputs
#' @param nlambda total number of tuning parameters. If \code{lambda.beta =
#'   NULL} and \code{lambda.gamma = NULL} then \code{nlambda} should be equal to
#'   \code{nlambda.beta x nlambda.gamma}. This is important to specify
#'   especially when a user defined sequence of tuning parameters is set.
#' @param thresh Convergence thresh for coordinate descent. Each
#'   coordinate-descent loop continues until the change in the objective
#'   function after all coefficient updates is less than thresh. Default
#'   value is \code{1e-4}.
#' @param maxit Maximum number of passes over the data for all tuning
#'   parameter values; default is 100.
#' @param cores The number of cores to use for certain calculations in the
#'   \code{\link{sail}} function, i.e. at most how many child processes will be
#'   run simultaneously using the \code{parallel} package. Must be at least one,
#'   and parallelization requires at least two cores. Default is 2.
#' @param initialization.type The procedure used to estimate the regression
#'   coefficients and used in the \code{\link{uni_fun}} function. If
#'   \code{"univariate"} then a series of univariate regressions is performed
#'   with the response variable \code{y}. If \code{"ridge"} then ridge
#'   regression is performed using the \code{\link[glmnet]{cv.glmnet}} function
#'   and the tuning parameter is chosen using 10 fold cross validation. The
#'   default is \code{"ridge"}.
#' @param center Should \code{x} and \code{y} be centered. Default is
#'   \code{TRUE}. Centering \code{y} applies to \code{family="gaussian"} only.
#' @param normalize Should \code{x} be scaled to have unit variance. Default is
#'   \code{TRUE}
#' @param verbose Should iteration number and vector of length \code{nlambda} be
#'   printed to console? Default is \code{TRUE}. 0 represents the algorithm has
#'   not converged for the pair of tuning parameters lambda.beta and
#'   lambda.gamma and 1 means it has converged
#' @param lambda.factor The factor for getting the minimal lambda in lambda
#'   sequence, where \code{min(lambda) = lambda.factor * max(lambda).
#'   max(lambda)} is the smallest value of lambda for which all coefficients are
#'   zero. The default depends on the relationship between \code{N} (the number
#'   of rows in the matrix of predictors) and \code{p} (the number of
#'   predictors). If \code{N > p}, the default is \code{1e-6}, close to zero. If
#'   \code{N < p}, the default is \code{0.01}. A very small value of
#'   lambda.factor will lead to a saturated fit.
#' @param weights observation weights. Can be total counts if responses are
#'   proportion matrices. Default is 1 for each observation. Currently NOT
#'   IMPLEMENTED
#' @details the index of the tuning parameters is as follows. If for example
#'   there are 10 lambda_gammas, and 20 lambda_betas, then the first
#'   lambda_gamma gets repeated 20 times. So the first twenty entries of tuning
#'   parameters correspond to 1 lambda_gamma and the 20 lambda_betas
#' @return An object with S3 class "sail" \describe{ \item{b0}{Intercept
#'   sequence of length \code{nlambda}} \item{beta}{A nvars x \code{nlambda}
#'   matrix of main effects (\eqn{\beta}) coefficients, stored in sparse column
#'   format \code{("CsparseMatrix")}} \item{alpha}{A nvars x \code{nlambda}
#'   matrix of interaction effects (\eqn{\alpha}) coefficients, stored in sparse
#'   column format \code{("CsparseMatrix")}} \item{gamma}{A nvars x
#'   \code{nlambda} matrix of (\eqn{\gamma}) coefficients, stored in sparse
#'   column format \code{("CsparseMatrix")}} \item{lambda.beta}{The sequence of
#'   tuning parameters used for the main effects} \item{lambda.gamma}{The
#'   sequence of tuning parameters used for the interaction effects}
#'   \item{tuning.parameters}{2 x nlambda matrix of tuning parameters. The first
#'   row corresponds to \code{lambda.beta} and the second row corresponds to
#'   \code{lambda.gamma}} \item{dfbeta}{list of length \code{nlambda} where each
#'   element gives the index of the nonzero \eqn{\beta} coefficients}
#'   \item{dfalpha}{list of length \code{nlambda} where each element gives the
#'   index of the nonzero \eqn{\alpha} coefficients} \item{x}{x matrix }
#'   \item{y}{response data} \item{bx}{column means of x matrix} \item{by}{mean
#'   of response} \item{sx}{column standard deviations of x matrix}
#'   \item{call}{the call to the function} \item{nlambda.gamma}{nlambda.gamma}
#'   \item{nlambda.beta}{nlambda.beta} \item{nlambda}{nlambda}
#'   \item{interaction.names}{interaction names} \item{main.effect.names}{main
#'   effect names} }
#' @note if the user specifies lambda.beta and lambda.gamma then they this will
#'   not take all possible combinations of lambda.beta and lambda.gamma. It will
#'   be the first element of each as a pair, and so on. This is done on purpose
#'   for use with the cv.sail function which uses the same lambda sequences for
#'   each fold.
#' @seealso \code{\link{sail_once}}
#' @examples
#' # number of observations
#' n <- 100
#'
#' # number of predictors
#' p <- 5
#'
#' # environment variable
#' e <- sample(c(0,1), n, replace = T)
#'
#' # main effects
#' x <- cbind(matrix(rnorm(n*p), ncol = p), e)
#'
#' # need to label columns
#' dimnames(x)[[2]] <- c("x1","x2","x3","x4","x5","e")
#'
#' # design matrix without intercept (can be user defined interactions)
#' X <- model.matrix(~(x1+x2+x3)*e+x1*x4+x3*x5-1, data = as.data.frame(x))
#'
#' # names must appear in the same order as X matrix
#' interaction_names <- grep(":", colnames(X), value = T)
#' main_effect_names <- setdiff(colnames(X), interaction_names)
#'
#' # response
#' Y <- X %*% rbinom(ncol(X), 1, 0.6) + 3*rnorm(n)
#'
#' # standardize data
#' data_std <- standardize(X,Y)
#'
#' result <- sail(x = data_std$x, y = data_std$y,
#'             main.effect.names = main_effect_names,
#'             interaction.names = interaction_names)
#' @author Sahir Bhatnagar
#'
#'   Maintainer: Sahir Bhatnagar \email{sahir.bhatnagar@@mail.mcgill.ca}
#'
#' @export

sail <- function(x, y, e, df = NULL, degree = 3, basis.intercept = FALSE,
                 center.x = TRUE,
                 group.penalty = c("gglasso", "MCP", "SCAD"),
                 family = c("gaussian", "binomial"),
                 weights, # observation weights
                 penalty.factor = rep(1, 1 + 2 * nvars), # predictor (adaptive lasso) weights, the last entry must be for the E variable
                 lambda.factor = ifelse(nobs < (1 + 2 * bscols * nvars), 0.01, 0.0001),
                 lambda = NULL,
                 alpha = 0.5,
                 nlambda = 100,
                 thresh = 1e-4,
                 maxit = 1000,
                 dfmax = 2 * nvars + 1,
                 exclude,
                 # initialization.type = c("ridge","univariate"),
                 # center = TRUE,
                 # normalize = FALSE,
                 verbose = TRUE) {

  # browser()

  ### Prepare all the generic arguments, then hand off to family functions

  family <- match.arg(family)
  if(alpha >= 1){
    warning("alpha >=1; set to 0.99")
    alpha <- 0.99
  }
  if(alpha <= 0){
    warning("alpha <= 0; set to 0.01")
    alpha <- 0.01
  }
  alpha <- as.double(alpha)

  # if (missing(df) & missing(degree)) stop("at least one of either df or degree must be supplied")

  # if(!missing(df)) df <- as.integer(df)
  # if(!missing(degree)) degree <- as.integer(degree)

  # if(df < degree){
  #   warning(sprintf("df too small, needs to be >= degree; set to %g",degree))
  #   df <- degree
  # }

  group.penalty <- match.arg(group.penalty)
  this.call <- match.call()
  nlam <- as.integer(nlambda)

  if (!is.matrix(x))
    stop("x has to be a matrix")
  if (any(is.na(x)))
    stop("Missing values in x not allowed")

  y <- drop(y)
  e <- drop(e)
  np <- dim(x)
  if (is.null(np) | (np[2] <= 1))
    stop("x should be a matrix with 2 or more columns")

  nobs <- as.integer(np[1])
  if (missing(weights))
    weights = rep(1, nobs)
  else if (length(weights) != nobs)
    stop(sprintf("number of elements in weights (%f) not equal to the number
                 of rows of x (%f)", length(weights), nobs))

  nvars <- as.integer(np[2])
  dimy <- dim(y)
  nrowy <- ifelse(is.null(dimy), length(y), dimy[1])
  dime <- dim(e)
  nrowe <- ifelse(is.null(dime), length(e), dime[1])

  if (nrowy != nobs)
    stop(paste("number of observations in y (", nrowy, ") not equal to the
               number of rows of x (", nobs, ")", sep = ""))

  if (nrowe != nobs)
    stop(paste("number of observations in e (", nrowe, ") not equal to the
               number of rows of x (", nobs, ")", sep = ""))

  if (length(y) != nobs)
    stop("x and y have different number of rows")

  if (length(e) != nobs)
    stop("x and e have different number of rows")

  if (!is.numeric(y))
    stop("The response y must be numeric. Factors must be converted to numeric")

  if (!is.numeric(e))
    stop("The environment variable e must be numeric. Factors must be converted to numeric")

  vnames <- colnames(x)
  if (is.null(vnames)) vnames <- paste("V",seq(nvars),sep="")

  ne <- as.integer(dfmax)
  if (missing(exclude)) exclude <- integer(0)
  if(any(penalty.factor == Inf)) {
    exclude <- c(exclude,seq(nvars)[penalty.factor==Inf])
    exclude <- sort(unique(exclude))
  }
  if(length(exclude)>0){
    jd <- match(exclude,seq(nvars),0)
    if(!all(jd>0)) stop("Some excluded variables out of range")
    penalty.factor[jd] <- 1 #ow can change lambda sequence
    jd <- as.integer(c(length(jd),jd))
  } else jd <- as.integer(0)
  vp <- as.double(penalty.factor)
  if (length(vp) < (1 + 2 * nvars)) stop("penalty.factor must be of length 1 + 2*ncol(x), and
                                     the first entry should correspond to the penalty.factor
                                     for X_E, the next ncol(x) correspond to main effects, and then
                                         interactions")
  we <- vp[1] # adaptive lasso weights for environment
  wj <- vp[(seq_len(nvars) + 1)] # adaptive lasso weights for main effects
  wje <- vp[(nvars + 2):length(vp)] # adaptive lasso weights for interactions

  thresh <- as.double(thresh)

  bscols <- ncol(splines::bs(x[,1], df = df, degree = degree, intercept = basis.intercept))

  if(is.null(lambda)){
    if(lambda.factor >= 1) stop("lambda.factor should be less than 1")
    flmin <- as.double(lambda.factor)
    ulam <- double(1)
  } else {
    flmin <- as.double(1)
    if(any(lambda<0)) stop("lambdas should be non-negative")
    ulam <- as.double(rev(sort(lambda)))
    nlam <- as.integer(length(lambda))
  }


  fit <- switch(family,
                gaussian = lspath(x = x,
                                  y = y,
                                  e = e,
                                  df = df,
                                  degree = degree,
                                  basis.intercept = basis.intercept,
                                  center.x = center.x,
                                  group.penalty = group.penalty,
                                  weights = weights,
                                  nlambda = nlam,
                                  thresh = thresh,
                                  maxit = maxit,
                                  verbose = verbose,
                                  alpha = alpha,
                                  nobs = nobs,
                                  nvars = nvars,
                                  jd = jd,
                                  vp = vp, # penalty.factor
                                  we = we,
                                  wj = wj,
                                  wje = wje,
                                  flmin = flmin, # lambda.factor
                                  vnames = vnames, #variable names
                                  ne = ne, # dfmax
                                  ulam = ulam)
  )

  fit$call <- this.call
  fit$nobs <- nobs
  class(fit) = c(class(fit), "sail")
  fit

  }


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

cv.sail <- function (x, y, e, df = NULL, degree = 3,
                     weights,
                     lambda = NULL,
                     type.measure = c("mse", "deviance", "class", "auc", "mae"),
                     nfolds = 10, foldid, grouped = TRUE, keep = FALSE, parallel = FALSE, ...) {

  if (missing(type.measure)) type.measure <- "default" else type.measure <- match.arg(type.measure)
  if (!is.null(lambda) && length(lambda) < 2)
    stop("Need more than one value of lambda for cv.sail")
  N <- nrow(x)
  if (missing(weights))
    weights <- rep(1, N)
  else weights <- as.double(weights)
  y <- drop(y)
  sail.call <- match.call(expand.dots = TRUE)
  which <- match(c("type.measure", "nfolds", "foldid", "grouped",
                  "keep"), names(sail.call), F)
  if (any(which))
    sail.call <- sail.call[-which]
  sail.call[[1]] <- as.name("sail")
  sail.object <- sail(x = x, y = y, e = e, df = df, degree = degree,
                     weights = weights, lambda = lambda, ...)

  sail.object$call <- sail.call

  ###Next line is commented out so each call generates its own lambda sequence
  # lambda <- sail.object$lambda

  # nz = sapply(predict(sail.object, type = "nonzero"), length)
  nz <- sapply(sail.object$active, length)
  if (missing(foldid)) foldid = sample(rep(seq(nfolds), length = N)) else nfolds = max(foldid)
  if (nfolds < 3)
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  outlist = as.list(seq(nfolds))
  if (parallel) {
    outlist = foreach(i = seq(nfolds), .packages = c("sail")) %dopar%
    {
      which = foldid == i

      if (length(dim(y)) > 1)
        y_sub = y[!which, ]
      else y_sub = y[!which]

      if (length(dim(e)) > 1)
        e_sub = e[!which, ]
      else e_sub = e[!which]

      sail(x = x[!which, , drop = FALSE], y = y_sub, e = e_sub, df = df, degree = degree,
           lambda = lambda, weights = weights[!which], ...)
    }
  } else {
    for (i in seq(nfolds)) {

      which = foldid == i

      if (is.matrix(y))
        y_sub = y[!which, ]
      else y_sub = y[!which]

      if (length(dim(e)) > 1)
        e_sub = e[!which, ]
      else e_sub = e[!which]

      outlist[[i]] = sail(x = x[!which, , drop = FALSE], y = y_sub, e = e_sub, df = df, degree = degree,
                          lambda = lambda, weights = weights[!which], ...)
    }
  }

  fun = paste("cv", class(sail.object)[[1]], sep = ".")
  lambda = sail.object$lambda

  cvstuff = do.call(fun, list(outlist, lambda, x, y, e, df, degree, weights,
                              foldid, type.measure, grouped, keep))
  cvm = cvstuff$cvm
  cvsd = cvstuff$cvsd
  nas = is.na(cvsd)
  if (any(nas)) {
    lambda = lambda[!nas]
    cvm = cvm[!nas]
    cvsd = cvsd[!nas]
    nz = nz[!nas]
  }
  cvname = cvstuff$name
  out = list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvup = cvm +
               cvsd, cvlo = cvm - cvsd, nzero = nz, name = cvname,
             sail.fit = sail.object)
  if (keep)
    out = c(out, list(fit.preval = cvstuff$fit.preval, foldid = foldid))
  lamin = if (cvname == "AUC")
    getmin(lambda, -cvm, cvsd)
  else getmin(lambda, cvm, cvsd)
  obj = c(out, as.list(lamin))
  class(obj) = "cv.sail"
  obj
}

