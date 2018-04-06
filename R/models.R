
#' @title Fit Sparse Additive Interaction Model with Strong Heredity
#' @description Function to fit the Sparse Additive Interaction Model with
#'   strong heredity for a sequence of tuning parameters. This is a penalized
#'   regression method that ensures the interaction term is non-zero only if its
#'   corresponding main-effects are non-zero. This model only considers the
#'   interactions between a single exposure (E) variable and a high-dimensional
#'   matrix (X). Additve (non-linear) main effects and interactions can be
#'   specified by the user. This can also be seen as a varying-coefficient
#'   model.
#' @param x input matrix of dimension \code{n x p}, where \code{n} is the number
#'   of subjects and p is number of X variables. Each row is an observation
#'   vector. Can be a high-dimensional (n < p) matrix. Can be a user defined
#'   design matrix of main effects only (without intercept) if
#'   \code{expand=FALSE}
#' @param y response variable. For \code{family="gaussian"} should be a 1 column
#'   matrix or numeric vector. For \code{family="binomial"}, should be a 1
#'   column matrix or numeric vector with -1 for failure and 1 for success.
#' @param e exposure or environment vector. Must be a numeric vector. Factors
#'   must be converted to numeric.
#' @param basis user defined basis expansion function. This function will be
#'   applied to every column in \code{x}. Specify \code{function(i) i} if no
#'   expansion is desired. Default: \code{function(i) splines::bs(i, df = 5)}.
#' @param group.penalty group lasso penalty. Can be one of \code{"gglasso"}
#'   (group lasso), \code{"grMCP"} (group MCP) or \code{"grSCAD"} (group SCAD).
#'   See references for details. Default: \code{"gglasso"}.
#' @param family response type. See \code{y} for details. Currently only
#'   \code{family = "gaussian"} is implemented. Default: \code{"gaussian"}.
#' @param center.x should the columns of \code{x} (after basis expansion) be
#'   centered (logical). Default: \code{TRUE}.
#' @param center.e should exposure variable \code{e} be centered. Default:
#'   \code{TRUE}.
#' @param expand should \code{basis} be applied to every column of \code{x}
#'   (logical). Set to \code{FALSE} if you want a user defined main effects
#'   design matrix. If \code{FALSE} the \code{group} membership argument must
#'   also be supplied. Default: \code{TRUE}.
#' @param group a vector of consecutive integers, starting from 1, describing
#'   the grouping of the coefficients. Only required when \code{expand=FALSE}.
#' @param weights observation weights. Default is 1 for each observation.
#'   Currently NOT IMPLEMENTED.
#' @param penalty.factor separate penalty factors can be applied to each
#'   coefficient. This is a number that multiplies lambda to allow differential
#'   shrinkage. Can be 0 for some variables, which implies no shrinkage, and
#'   that variable is always included in the model. Default is 1 for all
#'   variables. Must be of length \code{1 + 2*ncol(x)} where the first entry
#'   corresponds to the penalty.factor for \code{e}, the next \code{ncol(x)}
#'   entries correspond to main effects, and the following \code{ncol(x)}
#'   entries correspond to the interactions.
#' @param lambda.factor the factor for getting the minimal lambda in the lambda
#'   sequence, where \code{min(lambda) = lambda.factor * max(lambda)}.
#'   \code{max(lambda)} is the smallest value of lambda for which all
#'   coefficients are zero. The default depends on the relationship between
#'   \code{N} (the number of rows in the matrix of predictors) and \code{q} (the
#'   total number of predictors in the design matrix - including interactions).
#'   If \code{N > q}, the default is \code{1e-4}, close to zero. If \code{N <
#'   p}, the default is \code{0.01}. A very small value of lambda.factor will
#'   lead to a saturated fit.
#' @param lambda a user supplied lambda sequence. Typically, by leaving this
#'   option unspecified users can have the program compute its own lambda
#'   sequence based on \code{nlambda} and \code{lambda.factor}. Supplying a
#'   value of lambda overrides this. It is better to supply a decreasing
#'   sequence of lambda values than a single (small) value, if not, the program
#'   will sort user-defined lambda sequence in decreasing order automatically.
#'   Default: \code{NULL}.
#' @param alpha the mixing tuning parameter, with \eqn{0<\alpha<1}. It controls
#'   the penalization strength between the main effects and the interactions.
#'   The penalty is defined as \deqn{\lambda(1-\alpha)(w_e|\beta_e|+ \sum w_j
#'   ||\beta_j||_2) + \lambda\alpha(\sum w_{je} |\gamma_j|)}Larger values of
#'   \code{alpha} will favor selection of main effects over interactions.
#'   Smaller values of \code{alpha} will allow more interactions to enter the
#'   final model. Default: \code{0.5}
#' @param nlambda the number of lambda values. Default: 100
#' @param thresh convergence threshold for coordinate descent. Each
#'   coordinate-descent loop continues until the change in the objective
#'   function after all coefficient updates is less than \code{thresh}. Default:
#'   \code{1e-04}.
#' @param maxit maximum number of outer-loop iterations allowed at fixed lambda
#'   value. If models do not converge, consider increasing \code{maxit}.
#'   Default: 1000.
#' @param dfmax limit the maximum number of variables in the model. Useful for
#'   very large \code{q} (the total number of predictors in the design matrix -
#'   including interactions), if a partial path is desired. Default: \code{2 * p
#'   + 1} where p is the number of columns in \code{x}.
#' @param verbose display progress. Can be either 0,1 or 2. 0 will not display
#'   any progress, 2 will display very detailed progress and 1 is somewhere in
#'   between. Default: 1.
#' @return an object with S3 class "sail", \code{"*"}, where \code{"*"} is
#'   "lspath" or "logitreg". Results are provided for converged values of lambda
#'   only. \describe{\item{call}{the call that produced this object}
#'   \item{a0}{intercept sequence of length \code{nlambda}} \item{beta}{a (#
#'   main effects after basis expansion x \code{nlambda}) matrix of main effects
#'   coefficients, stored in sparse column format \code{("dgCMatrix")}}
#'   \item{alpha}{a (# interaction effects after basis expansion x
#'   \code{nlambda}) matrix of interaction effects coefficients, stored in
#'   sparse column format \code{("dgCMatrix")}} \item{gamma}{A \code{p x
#'   \code{nlambda}} matrix of (\eqn{\gamma}) coefficients, stored in sparse
#'   column format \code{("dgCMatrix")}} \item{bE}{exposure effect estimates of
#'   length \code{nlambda}} \item{active}{list of length \code{nlambda}
#'   containing character vector of selected variables} \item{lambda}{the actual
#'   sequence of lambda values used} \item{lambda2}{value for the mixing tuning
#'   parameter \eqn{0<\alpha<1}} \item{dfbeta}{the number of nonzero main effect
#'   coefficients for each value of lambda} \item{dfalpha}{the number of nonzero
#'   interaction coefficients for each value of lambda} \item{dfenviron}{the
#'   number of nonzero exposure (\code{e}) coefficients for each value of
#'   lambda} \item{dev.ratio}{the fraction of (null) deviance explained (for
#'   "lspath", this is the R-square). The deviance calculations incorporate
#'   weights if present in the model. The deviance is defined to be
#'   2*(loglike_sat - loglike), where loglike_sat is the log-likelihood for the
#'   saturated model (a model with a free parameter per observation). Hence
#'   dev.ratio=1-dev/nulldev.} \item{converged}{vector of logicals of length
#'   \code{nlambda} indicating if the algorithm converged} \item{nlambda}{number
#'   of converged lambdas} \item{design}{design matrix (X, E, X:E), of dimension
#'   \code{n x (2*ncols*p+1)} if \code{expand=TRUE}. This is used in the
#'   \code{predict} method.} \item{nobs}{number of observations}
#'   \item{nvars}{number of main effect variables} \item{vnames}{character of
#'   variable names for main effects (without expansion)} \item{ncols}{an
#'   integer of basis for each column of x if \code{expand=TRUE}, or an integer
#'   vector of basis for each variable if \code{expand=FALSE}}
#'   \item{center.x}{were the columns of x (after expansion) centered?}
#'   \item{center.e}{was \code{e} centered?} \item{basis}{user defined basis
#'   expansion function} \item{expand}{was the basis function applied to each
#'   column of x?} \item{group}{a vector of consecutive integers describing the
#'   grouping of the coefficients. Only if expand=FALSE}
#'   \item{interaction.names}{character vector of names of interaction
#'   variables} \item{main.effect.names}{character vector of names of main
#'   effect variables (with expansion)} }
#' @details The objective function for \code{family="gaussian"} is \deqn{RSS/2n
#'   + \lambda(1-\alpha)(w_e|\beta_e|+ \sum w_j ||\beta_j||_2) +
#'   \lambda\alpha(\sum w_{je} |\gamma_j|)} where \code{RSS} is the residual sum
#'   of squares and \code{n} is the number of observations. See Bhatnagar et al.
#'   (2018+) for details.
#'
#'   It is highly recommended to specify \code{center.x = TRUE} and
#'   \code{center.e = TRUE} for both convergence and interpretation reasons. If
#'   centered, the final estimates can be interpreted as the effect of the
#'   predictor on the response while holding all other predictors at their mean
#'   value. For computing speed reasons, if models are not converging or running
#'   slow, consider increasing \code{thresh}, decreasing \code{nlambda}, or
#'   increasing \code{lambda.factor} before increasing \code{maxit}. Then try
#'   increasing the value of \code{alpha} (which translates to more penalization
#'   on the interactions).
#'
#'   By default, \code{sail} uses the group lasso penalty on the basis
#'   expansions of \code{x}. To use the group MCP and group SCAD penalties (see
#'   Breheny and Huang 2015), the \code{grpreg} package must be installed.
#'
#' @examples
#' f.basis <- function(i) splines::bs(i, degree = 3)
#' data("sailsim")
#' fit <- sail(x = sailsim$x[,1:20], y = sailsim$y, e = sailsim$e,
#'             basis = f.basis, nlambda = 10)
#'
#' # estimated coefficients at each value of lambda
#' coef(fit)
#'
#' # predicted response at each value of lambda
#' predict(fit)
#'
#' #predicted response at a specific value of lambda
#' predict(fit, s = 1.5)
#' \dontrun{
#' if(interactive()){
#' # plot solution path for main effects and interactions
#' plot(fit)
#' # plot solution path only for main effects
#' plot(fit, type = "main")
#' # plot solution path only for interactions
#' plot(fit, type = "interaction")
#'  }
#' }
#'
#' @references Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010).
#'   Regularization Paths for Generalized Linear Models via Coordinate Descent.
#'   Journal of Statistical Software, 33(1), 1-22.
#'   \url{http://www.jstatsoft.org/v33/i01/}.
#' @references Breheny P and Huang J (2015). Group descent algorithms for
#'   nonconvex penalized linear and logistic regression models with grouped
#'   predictors. Statistics and Computing, 25: 173-187.
#' @references Yang Y, Zou H (2015). A fast unified algorithm for solving
#'   group-lasso penalize learning problems. Statistics and Computing. Nov
#'   1;25(6):1129-41.
#'   \url{http://www.math.mcgill.ca/yyang/resources/papers/STCO_gglasso.pdf}
#' @references Bhatnagar SR, Yang Y, Greenwood CMT. Sparse additive interaction
#'   models with the strong heredity property (2018+). Preprint.
#' @seealso \code{\link[splines]{bs}} \code{\link{cv.sail}}
#' @author Sahir Bhatnagar
#'
#'   Maintainer: Sahir Bhatnagar \email{sahir.bhatnagar@@gmail.com}
#' @rdname sail
#' @export
sail <- function(x, y, e,
                 basis = function(i) splines::bs(i, df = 5),
                 group.penalty = c("gglasso", "grMCP", "grSCAD"),
                 family = c("gaussian", "binomial"),
                 center.x = TRUE, # if true, this centers X
                 center.e = TRUE, # if true, this centers E
                 expand = TRUE, # if true, use basis to expand X's, else user should provide main effects design with group membership
                 group,
                 weights, # observation weights
                 penalty.factor = rep(1, 1 + 2 * nvars), # predictor (adaptive lasso) weights, the first entry must be for the E variable, then Xs, then X:E (gammas)
                 lambda.factor = ifelse(nobs < (1 + 2 * bscols * nvars), 0.01, 0.0001),
                 lambda = NULL,
                 alpha = 0.5,
                 nlambda = 100,
                 thresh = 1e-4,
                 fdev = 1e-5,
                 maxit = 1000,
                 dfmax = 2 * nvars + 1,
                 verbose = 1) {

  # browser()

  ### Prepare all the generic arguments, then hand off to family functions

  family <- match.arg(family)
  if (alpha >= 1) {
    warning("alpha >=1; set to 0.99")
    alpha <- 0.99
  }
  if (alpha <= 0) {
    warning("alpha <= 0; set to 0.01")
    alpha <- 0.01
  }
  alpha <- as.double(alpha)

  if (!expand & missing(group)) stop("group argument must be supplied when expand = FALSE")
  if (expand & !is.function(basis)) stop("basis needs to be a valid function when expand = TRUE")
  # tryCatch(is.function(function(i) splines::bs(i, df = 5)), error = function(e) FALSE)

  group.penalty <- match.arg(group.penalty)
  if (group.penalty %in% c("grMCP", "grSCAD")) {
    if (!requireNamespace("grpreg", quietly = TRUE)) {
      stop("Package \"grpreg\" needed for this function to work. Please install it.",
           call. = FALSE)
    }
  }

  this.call <- match.call()
  nlam <- as.integer(nlambda)

  if (!is.matrix(x)) {
    stop("x has to be a matrix")
  }
  if (any(is.na(x))) {
    stop("Missing values in x not allowed")
  }

  y <- drop(y)
  e <- drop(e)
  np <- dim(x)
  if (is.null(np) | (np[2] <= 1)) {
    stop("x should be a matrix with 2 or more columns")
  }

  nobs <- as.integer(np[1])
  if (missing(weights)) {
    weights <- rep(1, nobs)
  } else if (length(weights) != nobs) {
    stop(sprintf("number of elements in weights (%f) not equal to the number
                 of rows of x (%f)", length(weights), nobs))
  }

  if (!expand) {
    # nvars needs to be the number of original X variables.
    # if user provides their own design matrix, then we need to derive this from
    # the number of unique groups
    nvars <- length(unique(group))
  } else {
    nvars <- as.integer(np[2])
  }

  dimy <- dim(y)
  nrowy <- ifelse(is.null(dimy), length(y), dimy[1])
  dime <- dim(e)
  nrowe <- ifelse(is.null(dime), length(e), dime[1])

  if (nrowy != nobs) {
    stop(paste("number of observations in y (", nrowy, ") not equal to the
               number of rows of x (", nobs, ")", sep = ""))
  }

  if (nrowe != nobs) {
    stop(paste("number of observations in e (", nrowe, ") not equal to the
               number of rows of x (", nobs, ")", sep = ""))
  }

  if (length(y) != nobs) {
    stop("x and y have different number of rows")
  }

  if (length(e) != nobs) {
    stop("x and e have different number of rows")
  }

  if (!is.numeric(y)) {
    stop("The response y must be numeric. Factors must be converted to numeric")
  }

  if (!is.numeric(e)) {
    stop("The environment variable e must be numeric. Factors must be converted to numeric")
  }

  vnames <- colnames(x)
  if (is.null(vnames) & expand) {
    vnames <- paste("V", seq(nvars), sep = "")
  } else if (is.null(vnames) & !expand) {
    stop("x must have column names when expand = FALSE")
  }

  ne <- ifelse(expand, as.integer(dfmax), 2 * length(group) + 1)

  vp <- as.double(penalty.factor)
  if (length(vp) != (1 + 2 * nvars)) stop("penalty.factor must be of length 1 + 2*ncol(x)", call. = FALSE)
  we <- vp[1] # adaptive lasso weights for environment
  wj <- vp[(seq_len(nvars) + 1)] # adaptive lasso weights for main effects
  wje <- vp[(nvars + 2):length(vp)] # adaptive lasso weights for interactions

  thresh <- as.double(thresh)
  fdev <- as.double(fdev)


  if (!expand) {
    # this is for the user defined design matrix
    lambda.factor <- ifelse(nobs < (1 + 2 * length(group)), 0.01, 0.0001)
  } else {
    bscols <- ncol(basis(x[, 1, drop = FALSE])) # used for total number of variables for lambda.factor
  }

  if (is.null(lambda)) {
    if (lambda.factor >= 1) stop("lambda.factor should be less than 1")
    flmin <- as.double(lambda.factor)
    # ulam <- double(1)
    ulam <- NULL
  } else {
    flmin <- as.double(1)
    if (any(lambda < 0)) stop("lambdas should be non-negative")
    ulam <- as.double(rev(sort(lambda)))
    nlam <- as.integer(length(lambda))
  }

  fit <- switch(family,
    gaussian = lspath(
      x = x,
      y = y,
      e = e,
      basis = basis,
      center.x = center.x,
      center.e = center.e,
      expand = expand,
      group = group,
      group.penalty = group.penalty,
      weights = weights,
      nlambda = nlam,
      thresh = thresh,
      fdev = fdev,
      maxit = maxit,
      verbose = verbose,
      alpha = alpha,
      nobs = nobs,
      nvars = nvars,
      vp = vp, # penalty.factor
      we = we, # we, wj, wje are subsets of vp
      wj = wj,
      wje = wje,
      flmin = flmin, # lambda.factor
      vnames = vnames, # variable names
      ne = ne, # dfmax
      ulam = ulam
    )
  )

  fit$call <- this.call
  fit$nobs <- nobs
  class(fit) <- c(class(fit), "sail")
  fit
}



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
