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


SoftThreshold <- function(x, lambda) {
  # note: this works also if lam is a matrix of the same size as x.
  sign(x) * (abs(x) - lambda) * (abs(x) > lambda)
}


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

standardize <- function(x, center = TRUE, normalize = FALSE) {
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

#' Dopar function so that foreach stays in Imports and not Depends
`%dopar%` <- foreach::`%dopar%`


#' Color Blind Palette
cbbPalette <- function() {
  c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
}



#' @describeIn cv.lspath Computations for crossvalidation error
#' @note taken verbatim from glmnet
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



#' @describeIn cv.lspath Interpolation function.
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
