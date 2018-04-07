#' Internal sail functions
#'
#' @description Internal sail helper functions
#'
#' @details These functions are not intended for use by users.
#'
#' @name sail-internal
NULL


#' @rdname sail-internal
SoftThreshold <- function(x, lambda) {
  # note: this works also if lam is a matrix of the same size as x.
  sign(x) * (abs(x) - lambda) * (abs(x) > lambda)
}


#' @description \code{\%ni\%} is the opposite of \code{\%in\%}
#' @rdname sail-internal
"%ni%" <- Negate("%in%")


#' @rdname sail-internal
l2norm <- function(x) sqrt(sum(x^2))



#' @rdname sail-internal
`%dopar%` <- foreach::`%dopar%`



#' @description \code{cbbPalette} gives a Color Blind Palette
#' @rdname sail-internal
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


#' @rdname sail-internal
error.bars <- function(x, upper, lower, width = 0.02, ...) {
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x, upper, x, lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
  range(upper, lower)
}

#' @description \code{nonzero} is to determine which coefficients are non-zero
#' @param beta vector or 1 column matrix of regression coefficients
#' @rdname sail-internal
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






#' Likelihood function
#'
#' @description calculates likelihood function. Used to assess convergence of
#'   fitting algorithm. This corresponds to the Q(theta) function in the paper
#'
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
