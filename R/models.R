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
#'   internally by the \code{\link{shim}} function.
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
#'   sequence using the \code{\link{shim_once}} function which will be over a
#'   grid of tuning parameters for gamma as well. If the user specifies a
#'   sequence then this function will not automatically perform the serach over
#'   a grid. You will need to create the grid yourself e.g. repeat the
#'   lambda.gamma for each value of lambda.beta
#' @param lambda.gamma sequence of tuning parameters for the interaction
#'   effects. Default is \code{NULL} which means this function will
#'   automatically calculate a sequence of tuning paramters. See
#'   \code{\link{shim_once}} for details on how this sequence is calculated.
#' @param nlambda.gamma number of tuning parameters for gamma. This needs to be
#'   specified even for user defined inputs
#' @param nlambda.beta number of tuning parameters for beta. This needs to be
#'   specified even for user defined inputs
#' @param nlambda total number of tuning parameters. If \code{lambda.beta =
#'   NULL} and \code{lambda.gamma = NULL} then \code{nlambda} should be equal to
#'   \code{nlambda.beta x nlambda.gamma}. This is important to specify
#'   especially when a user defined sequence of tuning parameters is set.
#' @param threshold Convergence threshold for coordinate descent. Each
#'   coordinate-descent loop continues until the change in the objective
#'   function after all coefficient updates is less than threshold. Default
#'   value is \code{1e-4}.
#' @param max.iter Maximum number of passes over the data for all tuning
#'   parameter values; default is 100.
#' @param cores The number of cores to use for certain calculations in the
#'   \code{\link{shim}} function, i.e. at most how many child processes will be
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
#' @return An object with S3 class "shim" \describe{ \item{b0}{Intercept
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
#'   for use with the cv.shim function which uses the same lambda sequences for
#'   each fold.
#' @seealso \code{\link{shim_once}}
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
#' result <- shim(x = data_std$x, y = data_std$y,
#'             main.effect.names = main_effect_names,
#'             interaction.names = interaction_names)
#' @author Sahir Bhatnagar
#'
#'   Maintainer: Sahir Bhatnagar \email{sahir.bhatnagar@@mail.mcgill.ca}
#'
#' @export

funshim <- function(x, y, main.effect.names, interaction.names,
                    family = c("gaussian", "binomial"),
                    weights,
                    lambda.factor = ifelse(nobs < nvars, 0.01, 0.001),
                    lambda.beta = NULL,
                    lambda.gamma = NULL,
                    nlambda.gamma = 10,
                    nlambda.beta = 10,
                    nlambda = 100,
                    threshold = 1e-4,
                    max.iter = 100,
                    initialization.type = c("ridge","univariate"),
                    center = TRUE,
                    normalize = TRUE,
                    verbose = TRUE,
                    cores = 2) {

  if(missing(main.effect.names)) stop("main.effect.names cannot be missing")
  if(missing(interaction.names)) stop("interaction.names cannot be missing")

  initialization.type <- match.arg(initialization.type)
  family <- match.arg(family)
  this.call <- match.call()

  if (!is.matrix(x))
    stop("x has to be a matrix")
  if (any(is.na(x)))
    stop("Missing values in x not allowed")

  y <- drop(y)
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
  if (nrowy != nobs)
    stop(paste("number of observations in y (", nrowy, ") not equal to the
               number of rows of x (", nobs, ")", sep = ""))

  if (length(y) != nobs)
    stop("x and y have different number of rows")
  if (!is.numeric(y))
    stop("The response y must be numeric. Factors must be converted to numeric")

  vnames <- colnames(x)
  if (!all(c(main.effect.names, interaction.names) %in% vnames))
    stop("Some variables specified in main.effect.names were not found in
         the columnames of x")

  if (any(c(is.null(lambda.beta) & !is.null(lambda.gamma),
            !is.null(lambda.beta) & is.null(lambda.gamma))))
    stop("lambda.beta or lambda.gamma is NULL while the other is not NULL. Both
         should be NULL, or both should be specified")

  if (any(c(is.null(lambda.beta),is.null(lambda.gamma)))) {
    if (lambda.factor >= 1)
      stop("lambda.factor should be less than 1")
  } else {
    flmin = as.double(1)

    nDimLambdaBeta <- dim(lambda.beta)
    nDimLambdaGamma <- dim(lambda.gamma)

    if (!is.null(nDimLambdaBeta)) {
      if (nDimLambdaBeta[2] > 1)
        stop("lambda.beta must be a vector or 1 column matrix") }

    if (!is.null(nDimLambdaGamma)) {
      if (nDimLambdaGamma[2] > 1)
        stop("lambda.gamma must be a vector or 1 column matrix")}

    if (any(c(lambda.beta < 0, lambda.gamma < 0)))
      stop("lambda.beta and lambda.gamma should be non-negative")

    if (length(lambda.beta) != length(lambda.gamma))
      stop("length of lambda.beta needs to be the same as length or lambda.gamma")
  }


  fit <- switch(family,
                gaussian = lspath(x = x, y = y, main.effect.names = main.effect.names,
                                  interaction.names = interaction.names,
                                  lambda.beta = lambda.beta, lambda.gamma = lambda.gamma,
                                  lambda.factor = lambda.factor,
                                  nlambda.gamma = nlambda.gamma,
                                  nlambda.beta = nlambda.beta,
                                  nlambda = nlambda,
                                  threshold = threshold, max.iter = max.iter,
                                  initialization.type = initialization.type,
                                  center = center, normalize = normalize,
                                  verbose = verbose,
                                  cores = cores)
  )

  fit$call <- this.call
  fit$nobs <- nobs
  class(fit) = c(class(fit), "funshim")
  fit

  }


#' Cross-validation for shim
#'
#' @description Does k-fold cross-validation for shim and determines the optimal
#'   pair of tuning parameters (\eqn{\lambda_\beta} and \eqn{\lambda_\gamma})
#'
#' @inheritParams shim
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
#' @details The function runs shim nfolds+1 times; the first to get the tuning
#'   parameter sequences, and then the remainder to compute the fit with each of
#'   the folds omitted. The error is accumulated, and the average error and
#'   standard deviation over the folds is computed. Note also that the results
#'   of cv.shim are random, since the folds are selected at random using the
#'   \code{\link{createfolds}} function. Users can reduce this randomness by
#'   running cv.shim many times, and averaging the error curves.
#' @author Sahir Bhatnagar
#'
#'   Maintainer: Sahir Bhatnagar \email{sahir.bhatnagar@@mail.mcgill.ca}
#' @export

cv.shim <- function(x, y, main.effect.names, interaction.names,
                    weights,
                    lambda.beta = NULL, lambda.gamma = NULL,
                    nlambda.gamma = 10,
                    nlambda.beta = 10,
                    nlambda = 100,
                    parallel = TRUE,
                    type.measure = c("mse"),
                    nfolds = 10, ...) {


  # x = X; y = Y; main.effect.names = main_effect_names;
  # interaction.names = interaction_names;
  # lambda.beta = NULL ; lambda.gamma = NULL
  # threshold = 1e-4 ; max.iter = 500 ; initialization.type = "ridge";
  # nlambda.gamma = 5; nlambda.beta = 20; cores = 1;
  # center=TRUE; normalize=TRUE
  #
  # nfolds = 5
  # grouped = TRUE; keep = FALSE; parallel = TRUE

  #initialization.type <- match.arg(initialization.type)
  type.measure <- match.arg(type.measure)

  if (!is.null(lambda.beta) && length(lambda.beta) < 2)
    stop("Need more than one value of lambda.beta for cv.shim")
  if (!is.null(lambda.gamma) && length(lambda.gamma) < 2)
    stop("Need more than one value of lambda.gamma for cv.shim")

  N <- nrow(x)
  if (missing(weights))
    weights = rep(1, N)
  else weights = as.double(weights)
  y <- drop(y)
  shim.call <- match.call(expand.dots = TRUE)
  which <- match(c("type.measure", "nfolds"), names(shim.call), F)
  if (any(which))
    shim.call = shim.call[-which]
  shim.call[[1]] = as.name("shim")

  shim.object <- shim(x = x, y = y,
                      main.effect.names = main.effect.names,
                      interaction.names = interaction.names,
                      weights = weights,
                      lambda.beta = lambda.beta,
                      lambda.gamma = lambda.gamma,
                      nlambda = nlambda,
                      nlambda.gamma = nlambda.gamma,
                      nlambda.beta = nlambda.beta, ...)

  shim.object$call = shim.call

  nz.main = sapply(predict(shim.object, type = "nonzero")[["main"]], length)
  nz.interaction = sapply(predict(shim.object, type = "nonzero")[["interaction"]], length)

  foldid <- createfolds(y, k = nfolds)
  if (nfolds < 3)
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  outlist = as.list(seq(nfolds))

  pb <- progress::progress_bar$new(
    format = "  performing cross validation [:bar] :percent eta: :eta",
    total = nfolds, clear = FALSE, width= 90)
  pb$tick(0)

  if (parallel) {
    outlist = foreach(i = seq(nfolds), .packages = c("glmnet")) %dopar%
    {
      which = foldid == i
      if (is.matrix(y)) y_sub = y[!which, ] else y_sub = y[!which]
      #print(paste("Foldid = ",i))
      pb$tick()
      shim(x = x[!which, , drop = FALSE],
           y = y_sub,
           main.effect.names = main.effect.names,
           interaction.names = interaction.names,
           weights = weights[!which],
           lambda.beta = shim.object$lambda.beta,
           lambda.gamma = shim.object$lambda.gamma,
           nlambda = shim.object$nlambda,
           nlambda.gamma = nlambda.gamma,
           nlambda.beta = nlambda.beta, ...)

    }
  } else {
    for (i in seq(nfolds)) {
      which = foldid == i

      if (is.matrix(y))
        y_sub = y[!which, ] else y_sub = y[!which]

        pb$tick()

        outlist[[i]] = shim(x[!which, , drop = FALSE],
                            y = y_sub,
                            main.effect.names = main.effect.names,
                            interaction.names = interaction.names,
                            weights = weights[!which],
                            lambda.beta = shim.object$lambda.beta,
                            lambda.gamma = shim.object$lambda.gamma,
                            nlambda = shim.object$nlambda,
                            nlambda.gamma = nlambda.gamma,
                            nlambda.beta = nlambda.beta, ...)
    }
  }

  lambda.beta <- shim.object$lambda.beta
  lambda.gamma <- shim.object$lambda.gamma
  fun <- paste("cv", class(shim.object)[[1]], sep = "_")
  cvstuff <- do.call(fun, list(outlist = outlist,
                               x = x, y = y, foldid = foldid,
                               nlambda = shim.object$nlambda,
                               nlambda.beta = shim.object$nlambda.beta,
                               nlambda.gamma = shim.object$nlambda.gamma))

  cvm <- cvstuff$cvm
  cvsd <- cvstuff$cvsd
  # the following checks is any of the tunining parameter pairs
  # have a cvsd==NA or did not converge. cvstuff$converged should be
  # equal to the number of folds because it is the sum of booleans
  # that converged over all folds.
  nas <- (is.na(cvsd) + (cvstuff$converged != nfolds)) != 0
  if (any(nas)) {
    lambda.beta = lambda.beta[!nas]
    lambda.gamma = lambda.gamma[!nas]
    cvm = cvm[!nas]
    cvsd = cvsd[!nas]
    # this is the total number of non-zero parameters (both betas and alphas)
    nz.main = nz.main[!nas]
    nz.interaction = nz.interaction[!nas]
  }
  cvname <- cvstuff$name

  df <- as.data.frame(cbind(lambda.beta = lambda.beta,
                            lambda.gamma = lambda.gamma,
                            mse = cvm,
                            upper = cvm + cvsd,
                            lower = cvm - cvsd,
                            nz.main = nz.main,
                            nz.interaction = nz.interaction,
                            log.gamma = round(log(lambda.gamma),2)))

  rownames(df) <- gsub("\\.(.*)", "",rownames(df))

  out <- list(lambda.beta = lambda.beta, lambda.gamma = lambda.gamma,
              cvm = cvm, cvsd = cvsd, cvup = cvm + cvsd,
              cvlo = cvm - cvsd, nz.main = nz.main, name = cvname,
              nz.interaction = nz.interaction,
              shim.fit = shim.object, converged = cvstuff$converged,
              cvm.mat.all = cvstuff$cvm.mat.all,
              df = df)

  lamin.beta = if (cvname == "AUC")
    getmin_type(lambda.beta, -cvm, cvsd, type = "beta") else getmin_type(lambda.beta, cvm, cvsd, type = "beta")

  obj <- c(out, as.list(lamin.beta))
  class(obj) <- "cv.shim"
  obj
}


