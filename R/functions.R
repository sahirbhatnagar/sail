#' Univariate regressions
#'
#' @description Function used to create initial estimates in fitting algorithm
#'   of the strong heredity interaction model implemented in the
#'   \code{\link{sail}} function
#' @param x Design matrix of dimension \code{n x q}, where \code{n} is the
#'   number of subjects and q is the total number of variables; each row is an
#'   observation vector. This must include all main effects and interactions as
#'   well, with column names corresponding to the names of the main effects
#'   (e.g. \code{x1, x2, E}) and their interactions (e.g. \code{x1:E, x2:E}).
#'   All columns should be scaled to have mean 0 and variance 1; this is done
#'   internally by the \code{\link{sail}} function.
#' @param y response variable (matrix form) of dimension \code{n x 1}
#' @param type The procedure used to estimate the regression coefficients. If
#'   \code{"univariate"} then a series of univariate regressions is performed
#'   with the response variable \code{y}. If \code{"ridge"} then ridge
#'   regression is performed using the \code{\link[glmnet]{cv.glmnet}} function
#'   and the tuning parameter is chosen using 10 fold cross validation. The
#'   default is \code{"ridge"}.
#' @param variables character vector of variable names for which you want the
#'   univariate regression estimate. Must be contained in the column names of
#'   \code{x}. This only applies if \code{type="univariate"}.
#' @param include.intercept Should intercept be fitted (default is
#'   \code{FALSE}). Should be set to \code{TRUE} if \code{y} is not centered
#' @return Regression coefficients as a \code{q x 1 data.frame}
#' @note \code{p} is defined as the number of main effects. I have introduced
#'   \code{q} as being the total number of variables (e.g. the number of columns
#'   in the design matrix).
#'
#' @author
#' Sahir Bhatnagar
#'
#' Maintainer: Sahir Bhatnagar \email{sahir.bhatnagar@@mail.mcgill.ca}
#'
#' @seealso \code{\link{sail}}, \code{\link[glmnet]{cv.glmnet}}
#'
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
#' dimnames(x)[[2]] <- c(paste0("x",1:p), "e")
#'
#' # design matrix without intercept
#' X <- model.matrix(~(x1+x2+x3+x4+x5)*e-1, data = as.data.frame(x))
#'
#' # response
#' Y <- X %*% rbinom(ncol(X), 1, 0.2) + 3*rnorm(n)
#'
#' uni_fun(X, Y)
#' @export

uni_fun <- function(x, y, type = c("ridge", "univariate"),
                    variables, include.intercept = FALSE) {
  type <- match.arg(type)
  res <- switch(type,
    univariate = {
      plyr::ldply(variables, function(i) {
        if (include.intercept) {
          fit <- stats::lm.fit(x = methods::cbind2(rep(1, nrow(x)), x[, i, drop = F]), y = y)
          fit$coefficients[2]
        } else {
          fit <- stats::lm.fit(x = x[, i, drop = F], y = y)
          fit$coefficients[1]
        }
      }) %>%
        magrittr::set_rownames(variables) %>%
        magrittr::set_colnames("univariate_beta") %>%
        as.matrix()
    },
    ridge = {
      # fit the ridge to get betas and alphas
      fit <- glmnet::cv.glmnet(
        x = x, y = y, alpha = 0,
        standardize = F,
        intercept = include.intercept
      )
      # remove intercept (even if include.intercept is FALSE,
      # coef.glmnet returns
      # an intercept set to 0)
      as.matrix(coef(fit, s = "lambda.min")[-1, ])
    }
  )

  return(res)
}




#' Calculate Sequence of Tuning Parameters
#'
#' @description Function to calculate the sequence of tuning parameters based on
#'   the design matrix \code{x} and the response variable {y}. This is used in
#'   the \code{\link{sail_once}} function to calculate the tuning parameters
#'   applied to the main effects
#'
#' @inheritParams uni_fun
#' @param weights Separate penalty factors can be applied to each coefficient.
#'   This is a number that multiplies lambda to allow differential shrinkage,
#'   and can be used to apply adaptive LASSO. Can be 0 for some variables, which
#'   implies no shrinkage, and that variable is always included in the model.
#'   Default is 1 for all variables (and implicitly infinity for variables
#'   listed in exclude). Note: the penalty factors are internally rescaled to
#'   sum to nvars, and the lambda sequence will reflect this change.
#' @param lambda.factor The factor for getting the minimal lambda in lambda
#'   sequence, where \code{min(lambda) = lambda.factor * max(lambda).
#'   max(lambda)} is the smallest value of lambda for which all coefficients are
#'   zero. The default depends on the relationship between \code{N} (the number
#'   of rows in the matrix of predictors) and \code{p} (the number of
#'   predictors). If \code{N > p}, the default is \code{1e-6}, close to zero. If
#'   \code{N < p}, the default is \code{0.01}. A very small value of
#'   lambda.factor will lead to a saturated fit.
#' @param nlambda the number of lambda values - default is 100.
#' @param scale_x should the columns of x be scaled - default is FALSE
#' @param center_y should y be mean centered - default is FALSE
#' @return numeric vector of length \code{q}
#' @details The maximum lambda is calculated using the following inequality:
#'   \deqn{(N*w_j)^-1 | \sum x_ij y_i | \le \lambda_max}
#'
#'   The minimum lambda is given by lambda.factor*lambda_max. The sequence of
#'   nlambda values are decreasing from lambda_max to lambda_min on the log
#'   scale.
#'
#'   The penalty factors are internally rescaled to sum to the number of
#'   predictor variables in glmnet. Therefore, to get the correct sequence of
#'   lambdas when there are weights, this function first rescales the weights
#'   and then calclated the sequence of lambdas.
#'
#'   This formula is taken from section 2.5 of the \code{glmnet} paper in the
#'   Journal of Statistical Software (see references for details)
#'
#' @author
#' Sahir Bhatnagar
#'
#' Maintainer: Sahir Bhatnagar \email{sahir.bhatnagar@@mail.mcgill.ca}
#'
#' @references Friedman, J., Hastie, T. and Tibshirani, R. (2008)
#'   \emph{Regularization Paths for Generalized Linear Models via Coordinate
#'   Descent}, \url{http://www.stanford.edu/~hastie/Papers/glmnet.pdf}
#'   \emph{Journal of Statistical Software, Vol. 33(1), 1-22 Feb 2010}
#'   \url{http://www.jstatsoft.org/v33/i01/}
#'
#'   Yang, Y., & Zou, H. (2015). A fast unified algorithm for solving
#'   group-lasso penalize learning problems. \emph{Statistics and Computing},
#'   25(6), 1129-1141.
#'   \url{http://www.math.mcgill.ca/yyang/resources/papers/gglasso.pdf}
#'
#'
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
#' dimnames(x)[[2]] <- c(paste0("x",1:p), "e")
#'
#' # design matrix without intercept
#' X <- model.matrix(~(x1+x2+x3+x4+x5)*e-1, data = as.data.frame(x))
#'
#' # response
#' Y <- X %*% rbinom(ncol(X), 1, 0.2) + 3*rnorm(n)
#'
#' lambda_sequence(X,Y)
#' @export

lambda_sequence <- function(x, y, weights = NULL,
                            # lambda.factor = ifelse(nobs < nvars, 0.01, 1e-06),
                            lambda.factor,
                            nlambda, scale_x = F, center_y = F) {

  # when scaling, first you center then you standardize
  if (any(as.vector(weights) < 0)) stop("Weights must be positive")
  np <- dim(x)
  nobs <- as.integer(np[1])
  nvars <- as.integer(np[2])

  if (!is.null(weights) & length(as.vector(weights)) < nvars) {
    stop("You must provide weights for every column of x")
  }

  # scale the weights to sum to nvars
  w <- if (is.null(weights)) rep(1, nvars) else as.vector(weights) / sum(as.vector(weights)) * nvars

  sx <- if (scale_x) apply(x, 2, function(i) scale(i, center = TRUE, scale = mysd(i))) else x
  sy <- if (center_y) as.vector(scale(y, center = T, scale = F)) else as.vector(y)
  lambda.max <- max(abs(colSums(sy * sx) / w)) / nrow(sx)

  rev(exp(seq(log(lambda.factor * lambda.max), log(lambda.max), length.out = nlambda)))
}

#' Calculate Adaptive Weights based on Ridge Regression
#'
#' @description uses ridge regression from \code{glmnet} package to calculate
#'   the adaptive weights used in the fitting algorithm implemented in the
#'   \code{sail} function.
#' @param x Design matrix of dimension \code{n x q}, where \code{n} is the
#'   number of subjects and q is the total number of variables; each row is an
#'   observation vector. This must include all main effects and interactions as
#'   well, with column names corresponding to the names of the main effects
#'   (e.g. \code{x1, x2, E}) and their interactions (e.g. \code{x1:E, x2:E}).
#'   All columns should be scaled to have mean 0 and variance 1; this is done
#'   internally by the \code{\link{sail}} function.
#' @param y response variable (matrix form) of dimension \code{n x 1}
#' @param main.effect.names character vector of main effects names
#' @param interaction.names character vector of interaction names. MUST be
#'   separated by a colon (e.g. \code{x1:e, x2:e})
#' @param include.intercept logical if intercept should be fitted. Default is
#'   \code{FALSE}. Should be set to \code{TRUE} if \code{y} is not centered
#' @return \code{q x 1} matrix of weights for the main effects and interaction
#'   terms
#' @details Ridge regression is performed using the
#'   \code{\link[glmnet]{cv.glmnet}} function and the tuning parameter is chosen
#'   using 10 fold cross validation
#'
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
#' main_effect_names <- c(paste0("x",1:p), "e")
#' interaction_names <- paste0("x",1:p, ":e")
#' dimnames(x)[[2]] <- main_effect_names
#'
#' # design matrix without intercept
#' X <- model.matrix(~(x1+x2+x3+x4+x5)*e-1, data = as.data.frame(x))
#'
#' # response
#' Y <- X %*% rbinom(ncol(X), 1, 0.2) + 3*rnorm(n)
#'
#' # standardize data
#' data_std <- standardize(X,Y)
#'
#' ridge_weights(x = data_std$x, y = data_std$y,
#'               main.effect.names = main_effect_names,
#'               interaction.names = interaction_names)
#'
#' @author Sahir Bhatnagar
#'
#'   Maintainer: Sahir Bhatnagar \email{sahir.bhatnagar@@mail.mcgill.ca}
#'
#' @references Friedman, J., Hastie, T. and Tibshirani, R. (2008)
#'   \emph{Regularization Paths for Generalized Linear Models via Coordinate
#'   Descent}, \url{http://www.stanford.edu/~hastie/Papers/glmnet.pdf}
#'   \emph{Journal of Statistical Software, Vol. 33(1), 1-22 Feb 2010}
#'   \url{http://www.jstatsoft.org/v33/i01/}
#'
#'
#' @export

ridge_weights <- function(x,
                          y,
                          group,
                          main.effect.names.list,
                          interaction.names.list,
                          include.intercept = F) {

  # include.intercept=F
  # main.effect.names = c(main_effect_names, "X_E")
  # interaction.names = interaction_names
  # ===========================

  # fit the ridge to get betas and alphas
  fit <- glmnet::cv.glmnet(
    x = x, y = y, alpha = 0,
    standardize = F,
    intercept = include.intercept
  )

  # remove intercept (even if include.intercept is FALSE, coef.glmnet returns
  # an intercept set to 0)
  betas.and.alphas <- as.matrix(coef(fit, s = "lambda.1se")[-1, ])

  # create output matrix
  weights <- matrix(nrow = 2 * length(unique(group)) + 1)
  dimnames(weights)[[1]] <- c(paste0("X", unique(group)), "X_E", paste0("X", unique(group), ":X_E"))


  # l2 norm of theta hat
  norm_theta_hat <- sapply(
    seq_along(main.effect.names.list),
    function(k) l2norm(betas.and.alphas[main.effect.names.list[[k]], ])
  )

  # l2 norm of alpha hat
  norm_alpha_hat <- sapply(
    seq_along(interaction.names.list),
    function(k) l2norm(betas.and.alphas[interaction.names.list[[k]], ])
  )

  # beta_E hat
  beta_E_hat <- betas.and.alphas["X_E", ]

  # main effects weights
  weights[paste0("X", unique(group)), ] <- 1 / norm_theta_hat

  weights["X_E", ] <- abs(1 / beta_E_hat)

  weights[paste0("X", unique(group), ":X_E"), ] <- abs(beta_E_hat * norm_theta_hat / norm_alpha_hat)

  return(weights)
}


#' Soft Thresholding Function
#'
#' @inheritParams uni_fun
#' @param beta vector of regression coefficients to be thresholded
#' @param lambda tuning parameters
#' @param weight vector of weights for each beta
#' @return matrix of thresholded regression coefficients
#' @note user must supply \code{x} AND \code{y}, or \code{beta}, but not both. I
#'   set it up this way because to get the sequence of lambdas, I use the
#'   \code{beta} argument so that I only compute this once. I use the \code{x,
#'   y} argument for the CV folds. \code{lambda} can be a vector and this
#'   functions will return each thresholded beta for each lambda
#' @examples
#'  soft(beta = 0.5, lambda = seq(0.001,0.65,length.out = 10), weight = 1)
#' @author
#' Sahir Bhatnagar
#'
#' Maintainer: Sahir Bhatnagar \email{sahir.bhatnagar@@mail.mcgill.ca}
#' @export

soft <- function(x, y, beta, lambda, weight) {
  if (missing(x) & missing(y) & missing(beta)) stop("user must supply x AND y, or beta but not both")
  if (missing(x) & missing(y)) return(list("beta" = sign(beta) * pmax(0, abs(beta) - lambda * weight)))
  if (missing(beta)) {
    beta <- coef(stats::lm.fit(x = x[, 1, drop = F], y = y))[1]

    b_lasso <- sign(beta) * pmax(0, abs(beta) - lambda * weight)

    # need to return a matrix, because this is used in the step to
    # calculate y_tilde in the sail function
    return(matrix(b_lasso, ncol = 1))
  }
}

# SoftThreshold <- function(x, lambda) {
#   sign(x)* pmax(0, abs(x) - lambda)
# }

SoftThreshold <- function(x, lambda) {
  # note: this works also if lam is a matrix of the same size as x.
  sign(x) * (abs(x) - lambda) * (abs(x) > lambda)
}


#' Fit Strong Heredity model with one iteration
#'
#' @inheritParams ridge_weights
#' @inheritParams lambda_sequence
#' @param nlambda.gamma number of tuning parameters for gamma
#' @param nlambda.beta number of tuning parameters for beta
#' @param initialization.type The procedure used to estimate the regression
#'   coefficients and used in the \code{\link{uni_fun}} function. If
#'   \code{"univariate"} then a series of univariate regressions is performed
#'   with the response variable \code{y}. If \code{"ridge"} then ridge
#'   regression is performed using the \code{\link[glmnet]{cv.glmnet}} function
#'   and the tuning parameter is chosen using 10 fold cross validation. The
#'   default is \code{"ridge"}.
#' @description This function runs the first iteration of the fitting algorithm
#'   just to get the sequence of \code{lambda_gamma} and \code{lambda_beta}
#' @seealso \code{\link{sail}}
#' @return list of length 2, first element is \code{lambda_gamma} and second
#'   element is \code{lambda_beta}
#' @details A unique sequence of tuning parameters for the main effects
#'   (lambda.beta) is calculated for each tuning parameter for the interaction
#'   terms (lambda.gamma)
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
#' sail_once(x = data_std$x, y = data_std$y,
#'           main.effect.names = main_effect_names,
#'           interaction.names = interaction_names,
#'           nlambda.gamma = 5, nlambda.beta = 5)
#' @author
#' Sahir Bhatnagar
#'
#' Maintainer: Sahir Bhatnagar \email{sahir.bhatnagar@@mail.mcgill.ca}
#' @export

sail_once <- function(x, y, main.effect.names, interaction.names,
                      initialization.type = c("ridge", "univariate"),
                      nlambda.gamma = 20,
                      nlambda.beta = 20,
                      lambda.factor = ifelse(nobs < nvars, 0.01, 1e-6)) {
  initialization.type <- match.arg(initialization.type)
  np <- dim(x)
  nobs <- np[1]
  nvars <- np[2]

  # total number of tuning parameters
  nlambda <- nlambda.gamma * nlambda.beta

  adaptive.weights <- ridge_weights(
    x = x, y = y,
    main.effect.names = main.effect.names,
    interaction.names = interaction.names
  )

  # initialization
  betas_and_alphas <- uni_fun(
    variables = colnames(x), x = x, y = y,
    include.intercept = F,
    type = initialization.type
  )

  # this converts the alphas to gammas
  uni_start <- convert(betas_and_alphas,
    main.effect.names = main.effect.names,
    interaction.names = interaction.names
  )

  beta_hat_previous <- uni_start[main.effect.names, , drop = F]
  gamma_hat_previous <- uni_start[interaction.names, , drop = F]

  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # get tuning parameters for gamma
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  # this is a nsubjects x lambda matrix for each tuning parameter stored in a list
  # each element of the list corresponds to a tuning parameter
  y_tilde <- y - x[, main.effect.names, drop = F] %*% beta_hat_previous

  # calculate x_tilde for each beta vector corresponding to a diffent tuning parameter
  x_tilde <- xtilde(
    interaction.names = interaction.names,
    data.main.effects = x[, main.effect.names, drop = F],
    beta.main.effects = beta_hat_previous
  )

  # get the sequence of lambda_gammas using the first iteration of
  # x_tilde and y_tilde
  # x_tilde only has the interaction columns, therefore, penalty.factor
  # must also only include the weights for the interaction terms

  lambda_gamma <- lambda_sequence(
    x = x_tilde,
    y = y_tilde,
    weights = adaptive.weights[colnames(x_tilde), ],
    nlambda = nlambda.gamma,
    scale_x = F, center_y = F
  )

  # dont use GLMNET to get lambda sequence because it truncates the sequence,
  # so you may end up with a sequence that is less than nlambda.gamma
  # pass lambda_gamma to glmnet to get coefficients
  fit_gamma_hat_glmnet <- glmnet::glmnet(
    x = x_tilde,
    y = y_tilde,
    lambda = lambda_gamma,
    nlambda = nlambda.gamma,
    penalty.factor = adaptive.weights[colnames(x_tilde), ],
    standardize = F, intercept = F,
    lambda.min.ratio = lambda.factor
  )

  # use lambda gammas produced by lambda_sequence function
  lambda_gamma_glmnet <- lambda_gamma

  # get gamma coefficients and remove intercept
  # this results in a matrix of size p*(p-1)/2 x nlambda_gamma i.e.
  # the number of interaction variables by the number of lambda_gammas
  gamma_hat_next <- as.matrix(coef(fit_gamma_hat_glmnet))[-1, , drop = F]


  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # get tuning parameters for beta
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  beta_hat_next <- beta_hat_previous

  # for the lambda_beta sequence, calculate the sequences for each
  # beta, and then take the sequence that contains the maximum
  # this will ensure that all betas will be 0 under the maximum lambda

  # this is to store the lambda_beta sequences for each main.effect.names
  # then we will determine the maximum and minimum lambda_beta across
  # all main effects, for each of the nlambda.beta*nlambda.gamma combinations
  lambda_beta_temp <- vector("list", length(main.effect.names))
  names(lambda_beta_temp) <- main.effect.names


  # for each main effect, and for each nlambda_gamma,
  # we get a sequence of nlambda_beta tuning parameters
  # only store the max and min for each main effect for each lambda_gamma
  lambda_beta_seq_for_every_lambda_gamma <-
    replicate(nlambda.gamma,
      matrix(
        nrow = 2,
        ncol = length(main.effect.names),
        dimnames = list(
          paste0("lambda_", c("min", "max")),
          main.effect.names
        )
      ),
      simplify = "array"
    )

  # index data.frame to figure out which j < j'
  index <- data.frame(main.effect.names, seq_along(main.effect.names),
    stringsAsFactors = F
  )
  colnames(index) <- c("main.effect.names", "index")

  for (k in seq_len(ncol(gamma_hat_next))) {
    for (j in main.effect.names) {

      # determine the main effects not in j
      j_prime_not_in_j <- setdiff(main.effect.names, j)

      y_tilde_2 <- y -
        x[, j_prime_not_in_j, drop = F] %*% beta_hat_next[j_prime_not_in_j, , drop = F] -
        as.matrix(rowSums(xtilde_mod(
          beta.main.effects = beta_hat_next[j_prime_not_in_j, , drop = F],
          gamma.interaction.effects = gamma_hat_next[, k, drop = F],
          interaction.names = interaction.names[-grep(paste0("\\b", j, "\\b"), interaction.names)],
          data.main.effects = x[, j_prime_not_in_j, drop = F]
        )), ncol = 1)

      # j' less than j (main effects)
      j.prime.less <- index[
        which(index[, "index"] < index[which(index$main.effect.names == j), 2]),
        "main.effect.names"
      ]

      # need to make sure paste(j.prime.less,j,sep=":") are variables in x matrix
      # this is to get around situations where there is only interactions with the E variable
      j.prime.less.interaction <- intersect(paste(j.prime.less, j, sep = ":"), colnames(x))

      # need to get the main effects that are in j.prime.greater.interaction
      j.prime.less <- gsub("\\:(.*)", "", j.prime.less.interaction)

      # the if conditions in term1 and term2 are to check if there are
      # any variables greater or less than j
      # lapply is faster than mclapply here
      term_1 <- if (length(j.prime.less.interaction) != 0) {
        x[, j.prime.less.interaction] %*%
          (gamma_hat_next[j.prime.less.interaction, k, drop = F] *
            beta_hat_next[j.prime.less, , drop = F])
      } else {
        matrix(rep(0, nrow(x)), ncol = 1)
      }

      # j' greater than j
      j.prime.greater <- index[
        which(index[, "index"] > index[which(index$main.effect.names == j), 2]),
        "main.effect.names"
      ]

      # need to make sure j.prime.greater is a variable in x matrix
      # this is to get around situations where there is only interactions with the E variable
      j.prime.greater.interaction <- intersect(paste(j, j.prime.greater, sep = ":"), colnames(x))

      # need to get the main effects that are in j.prime.greater.interaction
      j.prime.greater <- if (all(gsub("\\:(.*)", "", j.prime.greater.interaction) == j)) {
        gsub("(.*)\\:", "", j.prime.greater.interaction)
      } else {
        gsub("\\:(.*)", "", j.prime.greater.interaction)
      }


      term_2 <- if (length(j.prime.greater.interaction) != 0) {
        x[, j.prime.greater.interaction] %*%
          (gamma_hat_next[j.prime.greater.interaction, k, drop = F] *
            beta_hat_next[j.prime.greater, , drop = F])
      } else {
        matrix(rep(0, nrow(x)), ncol = 1)
      }

      x_tilde_2 <- x[, j, drop = F] + term_1 + term_2

      lambda_beta_seq <- lambda_sequence(x_tilde_2,
        y_tilde_2,
        weights = adaptive.weights[colnames(x_tilde_2), ],
        nlambda = nlambda.beta
      )

      # the seqeunce of lambda_betas for variable j
      lambda_beta_seq_for_every_lambda_gamma[, j, k] <- c(min(lambda_beta_seq), max(lambda_beta_seq))
    }
  }


  # it is possible that the sequence of lambdas
  return(list(
    lambda_gamma = lambda_gamma_glmnet,
    lambda_beta = lapply(seq_len(length(lambda_gamma_glmnet)), function(i)
      rev(exp(seq(log(min(lambda_beta_seq_for_every_lambda_gamma[, , i])),
        log(max(lambda_beta_seq_for_every_lambda_gamma[, , i])),
        length.out = nlambda.beta
      ))))
  ))
}



#' Convert alphas to gammas
#'
#' @description function that takes a vector of betas (which are the main
#'   effects) and alphas (which are the interaction effects) and converts the
#'   alphas to gammas.
#' @param betas.and.alphas q x 1 data.frame or matrix of main effects and
#'   interaction estimates. For example the output from the \code{uni_fun}
#'   function. The rownames must be appropriately labelled because these labels
#'   will be used in other functions
#' @param main.effect.names character vector of main effects names. MUST be
#'   ordered in the same way as the column names of \code{x}. e.g. if the column
#'   names of \code{x} are \code{\"x1\",\"x2\"} then \code{main.effect.names =
#'   c("x1","x2")}
#' @param interaction.names character vector of interaction names. MUST be
#'   separated by a colon (e.g. x1:x2), AND MUST be
#'   ordered in the same way as the column names of \code{x}
#' @param epsilon threshold to avoid division by a very small beta e.g. if any
#'   of the main effects are less than epsilon, set gamma to zero. This should
#'   not really be an important parameter because this function is only used in
#'   the initialization step, where the intial estimates are from OLS or ridge
#'   regression and therefor should not be very close to 0
#' @details note that \deqn{y = \beta_0 + \beta_1 x_1 + ... + \beta_p x_p +
#'   \alpha_{12} x_1 x_2 + ... + \alpha_{p-1,p} x_p x_{p-1} } and
#'   \deqn{\alpha_{ij} = \gamma_{ij} * \beta_i*\beta_j , i < j}
#'
#'   This function is used because the fitting algorithm estimates the gammas,
#'   and furthermore, the L1 penalty is placed on the gammas. It is used only in
#'   the initialization step in the \code{\link{sail}} function
#'
#' @seealso \code{\link{sail}}, \code{\link{Q_theta}}
#' @return a labelled q x 1 data.frame of betas and gammas

# I thought I needed this.. but maybe not.. Its to convert
# alphas to gammas for initialization, but I never end up using the initialization
# gammas except for in the Q function of Choi et al. So perhpas I don't need it
# its causing issues because its no longer a simple transformation. I need to use l2norms
convert <- function(betas.and.alphas,
                    # main.effect.names, interaction.names,
                    epsilon = 1e-5,
                    group = group,
                    main.effect.names.list,
                    interaction.names.list) {
  browser()
  betas_and_gammas <- matrix(nrow = nrow(betas.and.alphas)) %>%
    magrittr::set_rownames(rownames(betas.and.alphas))

  sapply(unique(group), function(gr)
    l2norm(betas.and.alphas[interaction.names.list[[gr]], ]) / (l2norm(betas.and.alphas[main.effect.names.list[[gr]], ]) * as.vector(betas.and.alphas["X_E", ])))



  # add back the main effects which dont need to be transformed
  for (j in main.effect.names) {
    betas_and_gammas[j, ] <- betas.and.alphas[j, ]
  }

  return(betas_and_gammas)
}

#' Convert gammas to alphas
#'
#' @description function that takes a vector of betas (which are the main
#'   effects) and gammas and converts the alphas to gammas. This function is
#'   used to calculate the linear predictor of the likelihood function (the Q
#'   function in the fitting algorithm)
#' @param betas.and.gammas q x 1 data.frame or matrix of betas and gamma
#'   estimates. For example the output from the \code{convert} function. The
#'   rownames must be appropriately labelled because these labels will be used
#'   in other functions
#' @inheritParams convert
#' @return a labelled q x 1 data.frame of betas and alphas

convert2 <- function(beta,
                     gamma,
                     main.effect.names,
                     interaction.names,
                     group = group,
                     intercept = NULL) {
  betas.and.gammas <- methods::rbind2(beta, gamma)

  # create output matrix
  betas.and.alphas <- matrix(nrow = length(beta) * 2 - 1)
  dimnames(betas.and.alphas)[[1]] <- c(
    as.vector(do.call(c, main.effect.names)), "X_E",
    as.vector(do.call(c, interaction.names))
  )

  alphas <- do.call(rbind, lapply(unique(group), function(ind) {
    as.matrix(gamma[ind, ] * beta["X_E", ] * beta[main.effect.names[[ind]], , drop = FALSE])
  }))

  rownames(alphas) <- paste(rownames(alphas), "X_E", sep = ":")

  betas.and.alphas[rownames(beta), 1] <- beta
  betas.and.alphas[rownames(alphas), 1] <- alphas

  # add back intercept if it is non-NULL
  if (!is.null(intercept)) betas.and.alphas["(Intercept)", ] <- intercept

  return(betas.and.alphas)
}



#' Calculate working X's to update Gammas.
#'
#' @description function used to calculate working X's (xtilde) in step 3 of
#'   algorithm
#' @param interaction.names.list list of interaction names where each element of
#'   the list corresponds to the original X. Each list element should be of
#'   length df
#' @param x data frame or matrix containing all the mains effects, the
#'   environment, and the interaction of the basis expansions
#' @param beta.main.effects  matrix containing the coefficients of main effects
#' @param nlambda number of tuning parameters
#' @return matrix of working X's (xtilde)

xtilde <- function(x,
                   group,
                   main.effect.names.list,
                   interaction.names.list,
                   beta.main.effects) {

  # note that x[,interaction.names.list[[jj]]] contains the product term X_E:Phi_j, so
  # you dont need to multiply by X_E
  xtildas <- lapply(
    seq_along(unique(group)),
    function(jj)
      as.vector(beta.main.effects["X_E", ]) *
        (x[, interaction.names.list[[jj]]] %*% beta.main.effects[main.effect.names.list[[jj]], , drop = F])
  )


  xtildas <- do.call(cbind, xtildas)

  dimnames(xtildas)[[2]] <- paste0("X", unique(group))

  return(xtildas)
}



#' Calculate working X's to update Betas
#'
#' @description function used to calculate working X's (xtilde) in step 4 of
#'   algorithm
#' @param interaction.names character vector of interaction names. must be
#'   separated by a ':' (e.g. x1:x2)
#' @param data.main.effects data frame or matrix containing the main effects
#'   data
#' @param beta.main.effects data frame or matrix containing the coefficients of
#'   main effects
#' @param gamma.interaction.effects data frame or matrix containing the gamma
#'   parameters
#' @return matrix of working X's (xtilde) of dimension n x (p*(p-1)/2)
#' @note this function is a modified x_tilde for step 4 because we thought maybe
#'   there was a typo. Math and results suggests that there is a typo in the
#'   original paper.

xtilde_mod <- function( # interaction.names, data.main.effects, beta.main.effects,
                       # gamma.interaction.effects,
                       x,
                       group,
                       main.effect.names.list,
                       interaction.names.list,
                       beta.main.effects,
                       gamma.interaction.effects) {

  # create output matrix. no pipe is faster
  xtildas <- matrix(
    ncol = length(interaction.names),
    nrow = nrow(data.main.effects)
  )
  colnames(xtildas) <- interaction.names

  for (k in interaction.names) {

    # get names of main effects corresponding to interaction
    main <- strsplit(k, ":")[[1]]

    # step 4 to calculate x tilda
    xtildas[, k] <- prod(beta.main.effects[main, ]) * gamma.interaction.effects[k, ] *
      data.main.effects[, main[1], drop = F] *
      data.main.effects[, main[2], drop = F]
  }



  xtildas <- lapply(
    seq_along(unique(group)),
    function(j)
      as.vector(beta.main.effects["X_E", ]) *
        x[, "X_E"] *
        (x[, interaction.names.list[[j]]] %*% beta.main.effects[main.effect.names.list[[j]], , drop = F])
  )


  xtildas <- do.call(cbind, xtildas)



  return(xtildas)
}


# j = "x10"
# j_prime_not_in_j <- setdiff(main.effect.names,j)
# interaction.names = interaction_names
# Rprof(tmp <- tempfile())
# (beta.main.effects = beta_hat_next[j_prime_not_in_j, , drop = F]);
# (gamma.interaction.effects = gamma_hat_next);
# (interaction.names = interaction_names[-grep(paste0("\\b",j,"\\b"), interaction_names)]);
# (data.main.effects = x[,j_prime_not_in_j, drop = F])
#
# # create output matrix. no pipe is faster
# xtildas <- matrix(ncol = length(interaction.names),
#                   nrow = nrow(data.main.effects))
# colnames(xtildas) <- interaction.names
# dim(xtildas)
#
# for (k in interaction.names) {
#   k = "x1:e"
#   # get names of main effects corresponding to interaction
#   main <- unlist(stringr::str_split(k, ":"))
# main <- c("x1","x4")
#   strsplit(k, ":")[[1]]
#
#   # step 4 to calculate x tilda
#   xtildas[,k] <- prod(beta.main.effects[main,]) * gamma.interaction.effects[k,] *
#     data.main.effects[,main[1],drop = F] *
#     data.main.effects[,main[2],drop = F]
# }
#
# gamma.interaction.effects[1,1] <- 1
# v1 <- prod(beta.main.effects[main,],gamma.interaction.effects[1,1]) *
#   data.main.effects[,main[1],drop = F] *
#   data.main.effects[,main[2],drop = F]
# v2 <- prod(beta.main.effects[main,],gamma.interaction.effects[1,1])*apply(data.main.effects[,main], 1, prod)
#
# xtildas[,k] <- prod(beta.main.effects[main,],gamma.interaction.effects[1,1])*apply(data.main.effects[,main], 1, prod)
# xtildas[,k] <- prod(beta.main.effects[main,]) * gamma.interaction.effects[1,1] *
#   data.main.effects[,main[1],drop = F] *
#   data.main.effects[,main[2],drop = F]
# str(v1)
# all.equal(as.numeric(v1),as.numeric(v2))
#
# compare <- microbenchmark("v1" = {
#   xtildas[,3] <- prod(beta.main.effects[main,],gamma.interaction.effects[1,1]) *
#     data.main.effects[,main[1],drop = F] *
#     data.main.effects[,main[2],drop = F]
# },
# "v3" = {
#   xtildas[,1] <- prod(beta.main.effects[main,])*gamma.interaction.effects[1,1]*
#     data.main.effects[,main[1],drop = F] *
#     data.main.effects[,main[2],drop = F]
# } ,times = 2000)
# autoplot(compare)
#
#
#
#
# library(microbenchmark)
# library(ggplot2)
# compare <- microbenchmark("v1" = { y_tilde_2 <- y -
#   x[,j_prime_not_in_j, drop = F] %*%
#   beta_hat_next[j_prime_not_in_j, , drop = F] -
#   rowSums(
#     xtilde_mod(beta.main.effects = beta_hat_next[j_prime_not_in_j, , drop = F],
#                gamma.interaction.effects = gamma_hat_next,
#                interaction.names = interaction.names[-grep(paste0("\\b",j,"\\b"), interaction.names)],
#                data.main.effects = x[,j_prime_not_in_j, drop = F]))},
#   "v2" = { y_tilde_v2 <- y -
#     x[,j_prime_not_in_j, drop = F] %*%
#     beta_hat_next[j_prime_not_in_j, , drop = F] -
#     myrowSums(
#       xtilde_mod(beta.main.effects = beta_hat_next[j_prime_not_in_j, , drop = F],
#                  gamma.interaction.effects = gamma_hat_next,
#                  interaction.names = interaction.names[-grep(paste0("\\b",j,"\\b"), interaction.names)],
#                  data.main.effects = x[,j_prime_not_in_j, drop = F]))}, times = 1000L)
# autoplot(compare)
#
#
# all.equal(y_tilde_2,y_tilde_v2)
#
# Rprof()
# summaryRprof(tmp)



# myrowSums <- function (x, na.rm = FALSE, dims = 1L) {
#   dn <- dim(x)
#   p <- prod(dn[-(id <- seq_len(dims))])
#   dn <- dn[id]
#   z <- .Internal(rowSums(x, prod(dn), p, na.rm))
#   if (length(dn) > 1L) {
#     dim(z) <- dn
#     dimnames(z) <- dimnames(x)[id]
#   }
#   else names(z) <- dimnames(x)[[1L]]
#   z
# }




#' Likelihood function
#'
#' @description calculates likelihood function. Used to assess convergence of
#'   fitting algorithm. This corresponds to the Q(theta) function in the paper
#'
#' @inheritParams uni_fun
#' @param beta p x 1 matrix of main effect estimates
#' @param gamma p*(p-1)/2 x 1 matrix of gamma estimates
#' @param weights adaptive weights calculated by \code{ridge_weights} function
#'   with rownames corresponding to column names of x
#' @param lambda.beta a single tuning parameter for main effects
#' @param lambda.gamma a single tuning parameter for gammas
#' @param main.effect.names character vector of main effects names
#' @param interaction.names character vector of interaction names. must be
#'   separated by a colon (e.g. \code{x1:E})
#' @return value of likelihood function
#' @note you dont use the intercept in the calculation of the Q function
#' because its not being penalized

Q_theta <- function(R, nobs, lambda, alpha,
                    we, wj, wje,
                    betaE, theta_list, gamma) {

  # browser()
  (1 / (2 * nobs)) * crossprod(R) +
    lambda * (1 - alpha) * (
      we * abs(betaE) +
        sum(sapply(seq_along(theta_list), function(i) l2norm(theta_list[[i]]) * wj[i]))
    ) +
    lambda * alpha * sum(wje * abs(gamma))
}



#' Internal eclust functions
#'
#' @description Internal eclust helper functions
#'
#' @details These functions are not intended for use by users.
#'
#' @author Sahir Bhatnagar
#'
#'   Maintainer: Sahir Bhatnagar \email{sahir.bhatnagar@@mail.mcgill.ca}
#'
#' @name eclust-internal
NULL


#' @description \code{repcol} is to return a matrix with \code{x} being repeated
#'   into \code{n} columns
#' @param x is numeric vector to be repeated
#' @param n how many times \code{x} needs to be repeated
#' @rdname eclust-internal

repcol <- function(x, n) {
  s <- NCOL(x)
  matrix(x[, rep(1:s, each = n)], nrow = NROW(x), ncol = NCOL(x) * n)
}


#' @description \code{mysd} is to calculate standard deviation but with divisor of n
#'   and not n-1
#' @param i a vector of numerics
#' @rdname eclust-internal

mysd <- function(i) sqrt(crossprod(i - mean(i)) / length(i))

#' @description \code{check_col_0} is to check how many columns are 0
#' @param M is a matrix
#' @rdname eclust-internal

check_col_0 <- function(M) {
  M[, colSums(abs(M)) != 0, drop = F]
}







#' @description \code{nonzero} is to determine which coefficients are non-zero
#' @param beta vector or 1 column matrix of regression coefficients
#' @rdname eclust-internal
#' @export
nonzero <- function(beta, bystep = FALSE) {
  ### bystep = FALSE means which variables were ever nonzero
  ### bystep = TRUE means which variables are nonzero for each step
  beta <- as.matrix(beta)
  nr <- nrow(beta)
  if (nr == 1) {
    if (bystep) {
      apply(beta, 2, function(x) if (abs(x) > 0) {
          1
        } else {
          NULL
        })
    } else {
      if (any(abs(beta) > 0)) {
        1
      } else {
        NULL
      }
    }
  }
  else {
    beta <- abs(beta) > 0
    which <- seq(nr)
    ones <- rep(1, ncol(beta))
    nz <- as.vector((beta %*% ones) > 0)
    which <- which[nz]
    if (bystep) {
      if (length(which) > 0) {
        beta <- as.matrix(beta[which, , drop = FALSE])
        nzel <- function(x, which) if (any(x)) {
            which[x]
          } else {
            NULL
          }
        which <- apply(beta, 2, nzel, which)
        if (!is.list(which)) {
          which <- data.frame(which)
        }
        which
      }
      else {
        dn <- dimnames(beta)[[2]]
        which <- vector("list", length(dn))
        names(which) <- dn
        which
      }
    }
    else {
      which
    }
  }
}



#' Standardize Data
#'
#' @description Function that standardizes the data before running the fitting
#'   algorithm. This is necessary in all penalization methods so that the effect
#'   of a given penalty is the same for each predictor. This is used in the
#'   \code{\link{sail}} function
#' @inheritParams uni_fun
#' @param intercept Should \code{x} and \code{y} be centered. Default is
#'   \code{TRUE}
#' @param normalize Should \code{x} be scaled to have unit variance. Default is
#'   \code{TRUE}
#' @return list of length 5:
#' \describe{
#'   \item{x}{centered and normalized \code{x} matrix}
#'   \item{y}{centered \code{y} numeric vector}
#'   \item{bx}{numeric vector of column means of \code{x} matrix}
#'   \item{by}{mean of \code{y}}
#'   \item{sx}{standard deviations (using a divisor of \code{n}
#'   observations) of columns of \code{x} matrix}
#' }
#' @author
#' Sahir Bhatnagar
#'
#' Maintainer: Sahir Bhatnagar \email{sahir.bhatnagar@@mail.mcgill.ca}
#' @export

standardize <- function(x, y, center = TRUE, normalize = FALSE) {
  x <- as.matrix(x)
  # y <- as.numeric(y)
  n <- nrow(x)
  p <- ncol(x)

  if (center) {
    bx <- colMeans(x)
    # by <- mean(y)
    x <- scale(x, bx, FALSE)
    # y <- y - mean(y)
  } else {
    bx <- rep(0, p)
    by <- 0
  }
  if (normalize) {
    sx <- sqrt(colSums(x^2) / n)
    x <- scale(x, FALSE, sx)
  } else {
    sx <- rep(1, p)
  }

  return(list(
    x = x,
    # y = y,
    bx = bx, by = by, sx = sx
  ))
}

#' @description \code{\%ni\%} is the opposite of \code{\%in\%}
#' @rdname eclust-internal
"%ni%" <- Negate("%in%")

#' Compute cross validation error
#'
#' @description functions used to calculate cross validation error and used by
#'   the \code{\link{cv.sail}} function
#'
#' @param outlist list of cross validated fitted models. List is of length equal
#'   to \code{nfolds} argument in \code{\link{cv.sail}} function
#' @param foldid numeric vector indicating which fold each observation belongs
#'   to
#' @param nlambda.beta number of tuning paramters for the main effect
#' @param nlambda.gamma number of tuning parameters for the interaction effects
#' @param nlambda total number of tuning parameter pairs. This includes those
#'   pairs of tuning parameters that didn't converge in the CV folds.
#' @param outlist list containing results from cv.sail function. each element of
#'   the list is a run of sail_multiple_faster for each CV fold
#' @inheritParams uni_fun
#' @seealso \code{\link{cv.sail}}
#' @details The output of the \code{cv_lspath} function only returns values for those tuning
#'   paramters that DID converge

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







#' @describeIn cv_lspath Computations for crossvalidation error
#' @export
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

#' @describeIn cv_lspath Function that returns the tuning paramter corresponding
#'   to the minimum cross validated error and cross validated error within 1
#'   standard error of the minimum
#' @param type should be one of "beta" or "gamma"
getmin_type <- function(lambda, cvm, cvsd, type) {
  cvmin <- min(cvm, na.rm = TRUE)
  idmin <- cvm <= cvmin
  lambda.min <- lambda[idmin][which.max(lambda[idmin])]
  idmin <- match(lambda.min, lambda)
  semin <- (cvm + cvsd)[idmin]
  idmin <- cvm <= semin
  lambda.1se <- lambda[idmin][which.max(lambda[idmin])]

  # this is to get the index of the tuning parameter pair which is labelled by "s"
  # e.g. "s25" corresponds to the 25th pair of tuning parameters
  lambda.min.name <- gsub("\\.(.*)", "", names(lambda.min))
  lambda.1se.name <- gsub("\\.(.*)", "", names(lambda.1se))
  res <- list(
    lambda.min = lambda.min, lambda.min.name,
    lambda.1se = lambda.1se, lambda.1se.name
  )
  names(res) <- c(
    paste0("lambda.min.", type), "lambda.min.name",
    paste0("lambda.1se.", type), "lambda.1se.name"
  )
  res
}


#' @describeIn cv_lspath Interpolation function.
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


#' @describeIn cv_lspath Interpolation function.
lambda.interp <- function(lambda, s) {
  if (length(lambda) == 1) {
    nums <- length(s)
    left <- rep(1, nums)
    right <- left
    sfrac <- rep(1, nums)
  } else {
    s[s > max(lambda)] <- max(lambda)
    s[s < min(lambda)] <- min(lambda)
    k <- length(lambda)
    sfrac <- (lambda[1] - s) / (lambda[1] - lambda[k])
    lambda <- (lambda[1] - lambda) / (lambda[1] - lambda[k])
    coord <- stats::approx(lambda, seq(lambda), sfrac)$y
    left <- floor(coord)
    right <- ceiling(coord)
    sfrac <- (sfrac - lambda[right]) / (lambda[left] - lambda[right])
    sfrac[left == right] <- 1
  }
  list(left = left, right = right, frac = sfrac)
}

#' @describeIn cv_lspath Interpolation function.
lamfix <- function(lam) {
  llam <- log(lam)
  lam[1] <- exp(2 * llam[2] - llam[3])
  lam
}

#' Update Weights based on betas and gammas.
#'
#' @description uses betas and gammas to update weights. this is used to update
#'   the weights at each iteration of the fitting algorithm in the \code{sail}
#'   function
#' @param betas.and.gammas q x 1 data.frame or matrix of betas and gamma
#'   estimates. The rownames must be appropriately labelled because these labels
#'   are be used in this function and must match those in the arguments
#'   \code{main.effect.names} and \code{interaction.names}
#' @param main.effect.names character vector of main effects names
#' @param interaction.names character vector of interaction names. must be
#'   separated by a colon (e.g. x1:x2)
#' @return q x 1 matrix of weights

update_weights <- function(betas,
                           # gammas,
                           alphas,
                           main.effect.names,
                           interaction.names,
                           group,
                           epsilon = 1e-7) {

  # browser()

  betas.and.alphas <- rbind(betas, alphas)

  # create output matrix
  weights <- matrix(nrow = 2 * length(unique(group)) + 1)
  dimnames(weights)[[1]] <- c(paste0("X", unique(group)), "X_E", paste0("X", unique(group), ":X_E"))


  # l2 norm of theta hat
  norm_theta_hat <- sapply(
    seq_along(main.effect.names),
    function(k) l2norm(betas.and.alphas[main.effect.names[[k]], ])
  )

  # l2 norm of alpha hat
  norm_alpha_hat <- sapply(
    seq_along(interaction.names),
    function(k) l2norm(betas.and.alphas[interaction.names[[k]], ])
  )

  # beta_E hat
  beta_E_hat <- betas.and.alphas["X_E", ]

  # main effects weights
  weights[paste0("X", unique(group)), ] <- ifelse(norm_theta_hat < epsilon, 1e6, 1 / norm_theta_hat)

  weights["X_E", ] <- if (abs(as.vector(beta_E_hat)) < epsilon) 1e6 else abs(1 / as.vector(beta_E_hat))

  weights[paste0("X", unique(group), ":X_E"), ] <- ifelse(norm_alpha_hat < epsilon, 1e6, abs(beta_E_hat * norm_theta_hat / norm_alpha_hat))

  return(weights)
}

#' @rdname eclust-internal
isEmpty <- function(x) {
  return(length(x) == 0)
}


#' @description \code{checkargs.xy} function to check inputs of sail function
#'
#' @rdname eclust-internal
checkargs.xy <- function(x, y) {
  if (missing(x)) stop("x is missing")
  if (is.null(x) || !is.matrix(x)) stop("x must be a matrix")
  if (missing(y)) stop("y is missing")
  if (is.null(y) || !is.numeric(y)) stop("y must be numeric")
  if (ncol(x) == 0) stop("There must be at least one predictor
                         [must have ncol(x) > 0]")
  if (checkcols(x)) stop("x cannot have duplicate columns")
  if (length(y) == 0) stop("There must be at least one data point
                           [must have length(y) > 0]")
  if (length(y) != nrow(x)) stop("Dimensions don't match
                                 [length(y) != nrow(x)]")
}


#' @rdname eclust-internal
checkargs.misc <- function(sigma = NULL, alpha = NULL, k = NULL,
                           gridrange = NULL, gridpts = NULL, griddepth = NULL,
                           mult = NULL, ntimes = NULL,
                           beta = NULL, lambda = NULL, tol.beta = NULL, tol.kkt = NULL,
                           bh.q = NULL) {
  if (!is.null(sigma) && sigma <= 0) stop("sigma must be > 0")
  if (!is.null(lambda) && lambda < 0) stop("lambda must be >= 0")
  if (!is.null(alpha) && (alpha <= 0 || alpha >= 1)) stop("alpha must be
                                                          between 0 and 1")
  if (!is.null(k) && length(k) != 1) stop("k must be a single number")
  if (!is.null(k) && (k < 1 || k != floor(k))) stop("k must be an integer >= 1")
  if (!is.null(gridrange) && (length(gridrange) != 2 ||
    gridrange[1] > gridrange[2])) {
    stop("gridrange must be an interval of the form c(a,b) with a <= b")
  }
  if (!is.null(gridpts) && (gridpts < 20 || gridpts != round(gridpts))) {
    stop("gridpts must be an integer >= 20")
  }
  if (!is.null(griddepth) && (griddepth > 10 || griddepth != round(griddepth))) {
    stop("griddepth must be an integer <= 10")
  }
  if (!is.null(mult) && mult < 0) stop("mult must be >= 0")
  if (!is.null(ntimes) && (ntimes <= 0 || ntimes != round(ntimes))) {
    stop("ntimes must be an integer > 0")
  }
  if (!is.null(beta) && sum(beta != 0) == 0) stop("Value of lambda too large,
                                                  beta is zero")
  if (!is.null(lambda) && length(lambda) != 1) stop("lambda must be a single
                                                    number")
  if (!is.null(lambda) && lambda < 0) stop("lambda must be >=0")
  if (!is.null(tol.beta) && tol.beta <= 0) stop("tol.beta must be > 0")
  if (!is.null(tol.kkt) && tol.kkt <= 0) stop("tol.kkt must be > 0")
}

#' @rdname eclust-internal
tabular <- function(df, ...) {
  stopifnot(is.data.frame(df))

  align <- function(x) if (is.numeric(x)) "r" else "l"
  col_align <- vapply(df, align, character(1))

  cols <- lapply(df, format, ...)
  contents <- do.call(
    "paste",
    c(cols, list(sep = " \\tab ", collapse = "\\cr\n  "))
  )

  paste("\\tabular{", paste(col_align, collapse = ""), "}{\n  ",
    contents, "\n}\n",
    sep = ""
  )
}


#' @rdname eclust-internal
error.bars <- function(x, upper, lower, width = 0.02, ...) {
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x, upper, x, lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
  range(upper, lower)
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


l2norm <- function(x) sqrt(sum(x^2))



#' Simulate Data
#'
#' @description function to simulate data
#'
gendata <- function(n, p, df, degree, intercept = FALSE,
                    # E = stats::rnorm(n = n, sd = 0.5),
                    E = rbinom(n = n, size = 1, prob = 0.5),
                    # E = runif(n=n),
                    betaE = 2, SNR = 1) {

  # covariates
  X <- replicate(n = p, stats::runif(n))

  ncols <- degree + intercept
  # coefficients: each is a vector of length df and corresponds to the expansion of X_j
  b1 <- stats::rnorm(n = ncols)
  b2 <- stats::rnorm(n = ncols)
  b3 <- stats::rnorm(n = ncols)
  b4 <- stats::rnorm(n = ncols)
  b5 <- stats::rnorm(n = ncols)
  bE1 <- stats::rnorm(n = ncols)
  bE2 <- stats::rnorm(n = ncols)

  # b1 <- stats::runif(n = df, 0.5, 1.5)
  # b2 <- stats::runif(n = df, 1, 1.5)
  # b3 <- stats::runif(n = df, 0.1, 0.5)
  # b4 <- stats::runif(n = df, -1.5, -0.5)
  # b5 <- stats::runif(n = df, -2.0, -1.0)
  # bE1 <- stats::runif(n = df, 0.3, 1.5)
  # bE2 <- stats::runif(n = df, 0.8, 1.5)

  # error
  error <- stats::rnorm(n)

  Y.star <- splines::bs(X[, 1], df = df, degree = degree) %*% b1 +
    splines::bs(X[, 2], df = df, degree = degree) %*% b2 +
    # splines::bs(X[, 3], df = df, degree = degree) %*% b3 +
    # splines::bs(X[, 4], df = df, degree = degree) %*% b4 +
    # bs(X[,5], df = df, degree = degree) %*% b5 +
    betaE * E +
    E * splines::bs(X[, 1], df = df, degree = degree) %*% bE1 +
    E * splines::bs(X[, 2], df = df, degree = degree) %*% bE2

  k <- sqrt(stats::var(Y.star) / (SNR * stats::var(error)))

  Y <- Y.star + as.vector(k) * error

  return(list(
    x = X, y = Y, e = E, df = df, b1 = b1, b2 = b2, b3 = b3, b4 = b4,
    b5 = b5,
    bE1 = bE1, bE2 = bE2
  ))
}

gendata2 <- function(n, p, corr = 0,
                     E = truncnorm::rtruncnorm(n, a = -2, b = 2),
                     betaE = 2, SNR = 1, hierarchy = TRUE) {
  # this is modified from "VARIABLE SELECTION IN NONPARAMETRIC ADDITIVE MODEL" huang et al, Ann Stat.
  # n = 200
  # p = 10
  # corr = 1


  # covariates
  W <- replicate(n = p, truncnorm::rtruncnorm(n, a = 0, b = 1))
  U <- truncnorm::rtruncnorm(n, a = 0, b = 1)
  V <- truncnorm::rtruncnorm(n, a = 0, b = 1)

  X1 <- (W[, 1] + corr * U) / (1 + corr)
  X2 <- (W[, 2] + corr * U) / (1 + corr)
  X3 <- (W[, 3] + corr * U) / (1 + corr)
  X4 <- (W[, 4] + corr * U) / (1 + corr)

  X <- (W[, 5:p] + corr * V) / (1 + corr)

  Xall <- cbind(X1, X2, X3, X4, X)

  colnames(Xall) <- paste0("X", seq_len(p))

  # see "Variable Selection in NonParametric Addditive Model" Huang Horowitz and Wei
  f1 <- function(t) 5 * t
  f2 <- function(t) 4.5 * (2 * t - 1)^2
  f3 <- function(t) 4 * sin(2 * pi * t) / (2 - sin(2 * pi * t))
  f4 <- function(t) 6 * (0.1 * sin(2 * pi * t) + 0.2 * cos(2 * pi * t) +
      0.3 * sin(2 * pi * t)^2 + 0.4 * cos(2 * pi * t)^3 +
      0.5 * sin(2 * pi * t)^3)

  # error
  error <- stats::rnorm(n)

  if (hierarchy) {
    Y.star <- f1(X1) +
      f2(X2) +
      f3(X3) +
      f4(X4) +
      betaE * E +
      E * f3(X3) +
      E * f4(X4)
  } else {
    Y.star <- f1(X1) +
      f2(X2) +
      # f3(X3) +
      # f4(X4) +
      betaE * E +
      E * f3(X3) +
      E * f4(X4)
  }

  k <- sqrt(stats::var(Y.star) / (SNR * stats::var(error)))

  Y <- Y.star + as.vector(k) * error

  return(list(
    x = Xall, y = Y, e = E, f1 = f1(X1),
    f2 = f2(X2), f3 = f3(X3), f4 = f4(X4), betaE = betaE,
    f1.f = f1, f2.f = f2, f3.f = f3, f4.f = f4
  ))
}


gendataPaper <- function(n, p, corr = 0,
                         E = truncnorm::rtruncnorm(n, a = -1, b = 1),
                         # E = rbinom(n,1,0.5),
                         betaE = 2, SNR = 2, hierarchy = c("strong", "weak", "none"),
                         nonlinear = TRUE, interactions = TRUE) {
  # this is modified from "VARIABLE SELECTION IN NONPARAMETRIC ADDITIVE MODEL" huang et al, Ann Stat.
  # n = 200
  # p = 10
  # corr = 1

  hierarchy <- match.arg(hierarchy)

  # covariates
  W <- replicate(n = p, truncnorm::rtruncnorm(n, a = 0, b = 1))
  U <- truncnorm::rtruncnorm(n, a = 0, b = 1)
  V <- truncnorm::rtruncnorm(n, a = 0, b = 1)

  # W <- replicate(n = p, rnorm(n))
  # U <- rnorm(n)
  # V <- rnorm(n)

  X1 <- (W[, 1] + corr * U) / (1 + corr)
  X2 <- (W[, 2] + corr * U) / (1 + corr)
  X3 <- (W[, 3] + corr * U) / (1 + corr)
  X4 <- (W[, 4] + corr * U) / (1 + corr)

  X <- (W[, 5:p] + corr * V) / (1 + corr)

  Xall <- cbind(X1, X2, X3, X4, X)

  colnames(Xall) <- paste0("X", seq_len(p))

  # see "Variable Selection in NonParametric Addditive Model" Huang Horowitz and Wei
  f1 <- function(t) 5 * t
  f2 <- function(t) 4.5 * (2 * t - 1)^2
  f3 <- function(t) 4 * sin(2 * pi * t) / (2 - sin(2 * pi * t))
  f4 <- function(t) 6 * (0.1 * sin(2 * pi * t) + 0.2 * cos(2 * pi * t) +
      0.3 * sin(2 * pi * t)^2 + 0.4 * cos(2 * pi * t)^3 +
      0.5 * sin(2 * pi * t)^3)

  # error
  error <- stats::rnorm(n)

  if (!nonlinear) {
    # linear scenario; obeys hierachy. Scenario 2
    # Y.star <- 2 * (X1 - 1)  +
    #   2.5 * (X2 + 2) +
    #   2.7 * (X3) +
    #   3 * X4 +
    #   betaE * E +
    #   2 * E * X3 +
    #   2.5 * E * X4

    Y.star <- -1.5 * (X1 - 2) +
      1 * (X2 + 1) +
      1.5 * (X3) -
      2 * X4 +
      betaE * E +
      E * X3 -
      1.5 * E * X4


    scenario <- "2"
  } else {
    if (!interactions) {
      # main effects only; non-linear Scenario 3
      Y.star <- f1(X1) +
        f2(X2) +
        f3(X3) +
        f4(X4) +
        betaE * E
      scenario <- "3"
    } else {
      if (hierarchy == "none" & interactions) {
        # interactions only; non-linear
        Y.star <- E * f3(X3) +
          E * f4(X4)
        scenario <- "1c"
      } else if (hierarchy == "strong" & interactions) {
        # strong hierarchy; non-linear
        Y.star <- f1(X1) +
          f2(X2) +
          f3(X3) +
          f4(X4) +
          betaE * E +
          E * f3(X3) +
          E * f4(X4)
        scenario <- "1a"
      } else if (hierarchy == "weak" & interactions) {
        # weak hierarchy; linear
        Y.star <- f1(X1) +
          f2(X2) +
          # f3(X3) +
          # f4(X4) +
          betaE * E +
          E * f3(X3) +
          E * f4(X4)
        scenario <- "1b"
      }
    }
  }

  k <- sqrt(stats::var(Y.star) / (SNR * stats::var(error)))

  Y <- Y.star + as.vector(k) * error

  return(list(
    x = Xall, y = Y, e = E, Y.star = Y.star, f1 = f1(X1),
    f2 = f2(X2), f3 = f3(X3), f4 = f4(X4), betaE = betaE,
    f1.f = f1, f2.f = f2, f3.f = f3, f4.f = f4,
    X1 = X1, X2 = X2, X3 = X3, X4 = X4, scenario = scenario
  ))
}




gendata4 <- function(n = 100, p = 100, E = stats::rnorm(n), betaE = 1, SNR = 1) {
  # this is modified from SPAM Ravikumar et all JRSSB
  # n = 200
  # p = 10
  # corr = 1

  # covariates
  X <- replicate(n = p, stats::runif(n, -2.5, 2.5))
  colnames(X) <- paste0("X", seq_len(p))

  f1 <- function(t) -sin(1.5 * t)
  f2 <- function(t) t^3 + 1.5 * (t - 0.5)^2
  f3 <- function(t) -stats::dnorm(t, 0.5, 0.8)
  f4 <- function(t) sin(exp(-0.5 * t))

  X1 <- X[, 1]
  X2 <- X[, 2]
  X3 <- X[, 3]
  X4 <- X[, 4]

  # error
  error <- stats::rnorm(n)

  Y.star <- f1(X1) +
    f2(X2) +
    f3(X3) +
    f4(X4) +
    betaE * E +
    E * f3(X3) +
    E * f4(X4)

  k <- sqrt(stats::var(Y.star) / (SNR * stats::var(error)))

  Y <- Y.star + as.vector(k) * error

  return(list(
    x = X, y = Y, e = E, f1 = f1(X1),
    f2 = f2(X2), f3 = f3(X3), f4 = f4(X4), betaE = betaE,
    f1.f = f1, f2.f = f2, f3.f = f3, f4.f = f4
  ))
}


# from radchenko
gendata3 <- function(n = 300, p = 50, betaE = 1, SNR = 1) {

  # n = 200
  # p = 10
  # corr = 1
  # ===================

  # covariates
  # browser()

  W <- replicate(n = p, stats::runif(n, min = 0, max = 1))

  X1 <- W[, 1]
  X2 <- W[, 2]
  X3 <- W[, 3]
  X4 <- W[, 4]
  X5 <- W[, 5]
  E <- stats::runif(n, min = 0, max = 1)

  Xall <- W

  colnames(Xall) <- paste0("X", seq_len(p))

  # see "Variable Selection in NonParametric Addditive Model" Huang Horowitz and Wei
  # f1 <- function(t) t
  f2 <- function(t) 3 / (1 + t)^2
  f3 <- function(t) 2 * sin(t)
  f4 <- function(t) 4 * exp(t)
  f5 <- function(t) 6 * t^3

  # error
  error <- stats::rnorm(n)

  # Y.star <- sqrt(0.5) * (f1(X1)  +
  #                          f2(X2) +
  #                          f3(X3) +
  #                          f4(X4) +
  #                          f5(X5) +
  #                          betaE * E +
  #                          E * f1(X1) +
  #                          E * f2(X2))

  Y.star <-
    # f1(X1)  +
    f2(X2) +
      f3(X3) +
      f4(X4) +
      f5(X5) +
      betaE * E +
      # E * f1(X1) +
      E * f5(X5)

  k <- sqrt(stats::var(Y.star) / (SNR * stats::var(error)))

  Y <- Y.star + as.vector(k) * error

  return(list(
    x = Xall, y = Y, e = E,
    # f1 = f1(X1),
    f2 = f2(X2), f3 = f3(X3),
    f4 = f4(X4), f5 = f5(X5)
  ))
}



bic <- function(eta, sigma2, beta, eigenvalues, x, y, nt, c, df_lambda) {
  -2 * crossprod(y - x %*% betas.and.alphas) + c * df_lambda
}



# expan2 <- function(expr, xname = "x") {
#   sexpr <- substitute(expr)
#   if (is.name(sexpr)) {
#     expr <- call(as.character(sexpr), as.name(xname))
#   }
#   else {
#     if (!((is.call(sexpr) || is.expression(sexpr)) && xname %in%
#           all.vars(sexpr)))
#       stop(gettextf("'expr' must be a function, or a call or an expression containing '%s'",
#                     xname), domain = NA)
#     expr <- sexpr
#   }
#
#   Xmat = replicate(200, rnorm(100))
#   ll <- list(x = Xmat)
#   names(ll) <- xname
#   y <- eval(expr, envir = ll, enclos = parent.frame())
#
#   return(y)
#
# }
#
# set.seed(123)
# (mat <- replicate(4, rnorm(10)))
#
# fit <- function(x, expr = splines::bs(i, df = 5)) {
#   sexpr <- substitute(expr)
#   sexpr[[2]] <- substitute(x[,j])
#
#   lapply(seq_len(ncol(x)), function(j) eval(sexpr))
#
# }
#
# result <- fit(x = mat)
# lapply(result, head)
#
#
# set.seed(123)
# (mat <- replicate(4, rnorm(10)))
#
# fit <- function(x, expr = function(i) splines::bs(i, df = 5)) {
#
#   nvars <- ncol(x)
#   x <- scale(x, center = TRUE, scale = FALSE)
#
#   design <- design_mat(x = x, expr = expr, nvars = nvars)
#   design
#   # then fit some model on design
#
# }
#
# design_mat <- function(x, expr, nvars) {
#
#   lapply(seq_len(nvars), function(j) expr(x[, j]))
#
# }
#
# fit(x = mat)
#
#
# expan(x = replicate(5, rnorm(10)))
# splines::bs(rnorm(20), df = 5)
# expan(expr = splines::bs(y, intercept = TRUE), xname = "y")

design_sail <- function(x, e, expand, group, basis, nvars, vnames, center.x, center.e) {
  if (center.e) {
    e <- drop(standardize(e, center = TRUE, normalize = FALSE)$x)
  }

  if (!expand) {
    # Dont Expand X's if expand=FALSE. use user supplied design matrix
    if (center.x) {
      Phi_j_list <- lapply(split(seq(group), group), function(j) standardize(x[, j, drop = FALSE],
          center = TRUE
        )$x)
    } else {
      Phi_j_list <- lapply(split(seq(group), group), function(j) x[, j, drop = FALSE])
    }

    Phi_j <- do.call(cbind, Phi_j_list)
    main_effect_names <- vnames
    dimnames(Phi_j)[[2]] <- main_effect_names

    # X_E x Phi_j
    XE_Phi_j_list <- lapply(Phi_j_list, function(i) e * i)
    XE_Phi_j <- do.call(cbind, XE_Phi_j_list)
    interaction_names <- paste(main_effect_names, "E", sep = ":")
    dimnames(XE_Phi_j)[[2]] <- interaction_names
  } else {

    # Expand X's
    if (center.x) {
      Phi_j_list <- lapply(
        seq_len(nvars),
        function(j) standardize(basis(x[, j, drop = FALSE]),
            center = TRUE
          )$x
      )
    } else {
      Phi_j_list <- lapply(
        seq_len(nvars),
        function(j) basis(x[, j, drop = FALSE])
      )
    }

    ncols <- ncol(Phi_j_list[[1]]) # this is to get the number of columns for each expansion
    Phi_j <- do.call(cbind, Phi_j_list)
    main_effect_names <- paste(rep(vnames, each = ncols), rep(seq_len(ncols), times = nvars), sep = "_")
    dimnames(Phi_j)[[2]] <- main_effect_names

    # E x Phi_j
    XE_Phi_j_list <- lapply(Phi_j_list, function(i) e * i)
    XE_Phi_j <- do.call(cbind, XE_Phi_j_list)
    interaction_names <- paste(main_effect_names, "E", sep = ":")
    dimnames(XE_Phi_j)[[2]] <- interaction_names
  }

  # this is used for the predict function
  design <- cbind(Phi_j, "E" = e, XE_Phi_j)

  return(list(
    Phi_j_list = Phi_j_list, Phi_j = Phi_j,
    XE_Phi_j_list = XE_Phi_j_list, XE_Phi_j = XE_Phi_j,
    main_effect_names = main_effect_names, interaction_names = interaction_names,
    design = design, ncols = if (expand) ncols else sapply(Phi_j_list, ncol)
  ))
}


#' Not in convenience function
"%ni%" <- Negate("%in%")


#' Color Blind Palette
cbbPalette <- function() {
  c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
}

#' An alternative to \code{summaryRprof()}
#'
#' \code{proftools} parses a profiling file and prints an easy-to-understand
#' table showing the most time-intensive function calls.
#'
#' Line numbers are included if \code{Rprof()} was run with
#' \code{line.numbering=TRUE}. If it was run with \code{memory.profiling=TRUE},
#' this function will probably break.
#'
#' Below the table are printed any files identified if line numbering is true,
#' the total time recorded by \code{Rprof()}, and the "parent call".  The
#' parent call consists of the parent call stack of all the call stacks in the\
#' table. Note that this is the parent call stack of only the printed lines,
#' not of all stacks recorded by \code{Rprof()}. This makes the table easier to read and fit into the console.
#'
#' @param file A profiling file generated by \code{Rprof()}
#' @param lines The number of lines (call stacks) you want returned. Lines are
#' printed from most time-intensive to least.
proftable <- function(file, lines = 10) {
  profdata <- readLines(file)
  interval <- as.numeric(strsplit(profdata[1L], "=")[[1L]][2L]) / 1e+06
  filelines <- grep("#File", profdata)
  files <- profdata[filelines]
  profdata <- profdata[-c(1, filelines)]
  total.time <- interval * length(profdata)
  ncalls <- length(profdata)
  profdata <- gsub("\\\"| $", "", profdata)
  calls <- lapply(profdata, function(x) rev(unlist(strsplit(x, " "))))
  stacktable <- as.data.frame(table(sapply(calls, function(x) paste(x, collapse = " > "))) / ncalls * 100, stringsAsFactors = FALSE)
  stacktable <- stacktable[order(stacktable$Freq[], decreasing = TRUE), 2:1]
  colnames(stacktable) <- c("PctTime", "Call")
  stacktable <- head(stacktable, lines)
  shortcalls <- strsplit(stacktable$Call, " > ")
  shortcalls.len <- range(sapply(shortcalls, length))
  parent.call <- unlist(lapply(seq(shortcalls.len[1]), function(i) Reduce(intersect, lapply(shortcalls, "[[", i))))
  shortcalls <- lapply(shortcalls, function(x) setdiff(x, parent.call))
  stacktable$Call <- sapply(shortcalls, function(x) paste(x, collapse = " > "))
  if (length(parent.call) > 0) {
    parent.call <- paste(paste(parent.call, collapse = " > "), "> ...")
  } else {
    parent.call <- "None"
  }
  frac <- sum(stacktable$PctTime)
  attr(stacktable, "total.time") <- total.time
  attr(stacktable, "parent.call") <- parent.call
  attr(stacktable, "files") <- files
  attr(stacktable, "total.pct.time") <- frac
  print(stacktable, row.names = FALSE, right = FALSE, digits = 3)
  if (length(files) > 0) {
    cat("\n")
    cat(paste(files, collapse = "\n"))
    cat("\n")
  }
  cat(paste("\nParent Call:", parent.call))
  cat(paste("\n\nTotal Time:", total.time, "seconds\n"))
  cat(paste0("Percent of run time represented: ", format(frac, digits = 3)), "%")

  invisible(stacktable)
}
