
#' @title Fit Sparse Additive Interaction Model with Strong Heredity
#' @description Function to fit the Sparse Additive Interaction Model with
#'   strong heredity for a sequence of tuning parameters. This is a penalized
#'   regression method that ensures the interaction term is non-zero only if its
#'   corresponding main-effects are non-zero. This model only considers the
#'   interactions between a single exposure (E) variable and a high-dimensional
#'   matrix (X). Additive (non-linear) main effects and interactions can be
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
#' @param strong Flag specifying strong hierarchy (TRUE) or weak hierarchy
#'   (FALSE). Default FALSE.
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
#' @param weights observation weights, a vector of size n.
#'   Default is 1 for each observation.
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
#' @param fdev minimum fractional change in deviance for stopping path. Default:
#'   \code{1e-5}.
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
#'   nlambda} matrix of (\eqn{\gamma}) coefficients, stored in sparse column
#'   format \code{("dgCMatrix")}} \item{bE}{exposure effect estimates of length
#'   \code{nlambda}} \item{active}{list of length \code{nlambda} containing
#'   character vector of selected variables} \item{lambda}{the actual sequence
#'   of lambda values used} \item{lambda2}{value for the mixing tuning parameter
#'   \eqn{0<\alpha<1}} \item{dfbeta}{the number of nonzero main effect
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
#' @examples
#' f.basis <- function(i) splines::bs(i, degree = 3)
#' # we specify dfmax to early stop the solution path to
#' # limit the execution time of the example
#' fit <- sail(x = sailsim$x, y = sailsim$y, e = sailsim$e,
#'             basis = f.basis, nlambda = 50, dfmax = 5)
#'
#' # estimated coefficients at each value of lambda
#' coef(fit)
#'
#' # predicted response at each value of lambda
#' predict(fit)
#'
#' #predicted response at a specific value of lambda
#' predict(fit, s = 0.5)
#' if(interactive()){
#' # plot solution path for main effects and interactions
#' plot(fit)
#' # plot solution path only for main effects
#' plot(fit, type = "main")
#' # plot solution path only for interactions
#' plot(fit, type = "interaction")
#'  }
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
                 strong = TRUE,
                 group.penalty = c("gglasso", "grMCP", "grSCAD"),
                 family = c("gaussian", "binomial"),
                 center.x = TRUE, # if true, this centers X
                 center.e = TRUE, # if true, this centers E
                 expand = TRUE, # if true, use basis to expand X's, else user should provide main effects design with group membership
                 group,
                 weights,
                 penalty.factor = rep(1, 1 + 2 * nvars), # predictor (adaptive lasso) weights, the first entry must be for the E variable, then Xs, then X:E (gammas)
                 lambda.factor = ifelse(nobs < (1 + 2 * bscols * nvars), 0.00001, 0.00000001),
                 lambda = NULL,
                 alpha = 0.5,
                 nlambda = 100,
                 thresh = 1e-6,
                 fdev = 1e-6,
                 maxit = 200,
                 dfmax = 2 * nvars + 1,
                 verbose = 0) {

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
        call. = FALSE
      )
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

  weights=weights*length(weights) / sum(weights)


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
    flmin<- as.double(lambda.factor)
    # ulam <- double(1)
    ulam <- NULL
  } else {
    flmin <- as.double(1)
    if (any(lambda < 0)) stop("lambdas should be non-negative")
    ulam <- as.double(rev(sort(lambda)))
    nlam <- as.integer(length(lambda))
  }

  if (strong) {
      fit <- switch(family,
                    gaussian = lspathweights(
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
    } else {
      fit <- switch(family,
                    gaussian = lspathweak(
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
    }

  fit$call <- this.call
  fit$nobs <- nobs
  class(fit) <- c(class(fit), "sail")
  fit
}
