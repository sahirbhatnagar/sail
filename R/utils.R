#' Internal sail functions
#'
#' @description Internal sail helper functions
#'
#' @details These functions are not intended for use by users.
#'
#' @name sail-internal
NULL


#' @rdname sail-internal
#' @param x numeric value of a coefficient
#' @param lambda tuning parameter value
SoftThreshold <- function(x, lambda) {
  # note: this works also if lam is a matrix of the same size as x.
  sign(x) * (abs(x) - lambda) * (abs(x) > lambda)
}



"%ni%" <- Negate("%in%")


#' @rdname sail-internal
l2norm <- function(x) sqrt(sum(x^2))



`%dopar%` <- foreach::`%dopar%`



#' @description \code{cbbPalette} gives a Color Blind Palette
#' @rdname sail-internal
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


fix.lam <- function(lam){
  if(length(lam)>2){
    llam=log(lam)
    lam[1]=exp(2*llam[2]-llam[3])
  }
  lam
}


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
#' @param bystep \code{bystep = FALSE} means which variables were ever nonzero.
#'   \code{bystep = TRUE} means which variables are nonzero for each step
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



#' @description \code{check_col_0} is to check how many columns are 0 and is
#'   used in the fitting functions \code{lspath}
#' @param M is a matrix
#' @rdname sail-internal
check_col_0 <- function(M) {
  M[, colSums(abs(M)) != 0, drop = F]
}


#' Objective function
#'
#' @description calculates likelihood function. Used to assess convergence of
#'   fitting algorithm. This corresponds to the Q(theta) function in the paper
#'
#' @param R residual
#' @param nobs number of observations
#' @inheritParams sail
#' @param we penalty factor for exposure variable
#' @param wj penalty factor for main effects
#' @param wje penalty factor for interactions
#' @param weights observational weights, default is 1 for all observations
#' @param betaE estimate of exposure effect
#' @param theta_list estimates of main effects
#' @param gamma estimates of gamma parameter
#' @return value of the objective function
Q_theta <- function(R, nobs, lambda, alpha,weights,
                    we, wj, wje,
                    betaE, theta_list, gamma) {
  if (missing(weights)) {
    weights <- rep(1, nobs)
  } else if (length(weights) != nobs) {
    stop(sprintf("number of elements in weights (%f) not equal to the number
                 of rows of x (%f)", length(weights), nobs))
  }
  (1 / (2 * nobs)) * crossprod(sqrt(weights)*R) +
    lambda * (1 - alpha) * (
      we * abs(betaE) +
        sum(sapply(seq_along(theta_list), function(i) l2norm(theta_list[[i]]) * wj[i]))
    ) +
    lambda * alpha * sum(wje * abs(gamma))
}


#' Standardize Data
#'
#' @description Function that standardizes the data before running the fitting
#'   algorithm. This is used in the \code{\link{sail}} function
#' @param x data to be standardized
#' @param center Should \code{x} be centered. Default is \code{TRUE}
#' @param normalize Should \code{x} be scaled to have unit variance. Default is
#'   \code{FALSE}
#' @return list of length 3: \describe{ \item{x}{centered and possibly
#'   normalized \code{x} matrix} \item{bx}{numeric vector of column means of
#'   \code{x} matrix} \item{sx}{standard deviations (using a divisor of
#'   \code{n} observations) of columns of \code{x} matrix} }
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
    # y = y, by = by,
    bx = bx, sx = sx
  ))
}


#' @title Sail design matrix
#' @description Create design matrix used in \code{\link{sail}} function
#' @inheritParams sail
#' @param nvars number of variables
#' @param vnames variable names
design_sail <- function(x, e, expand, group, basis, nvars, vnames, center.x, center.e) {
    # e <- drop(standardize(e, center = TRUE, normalize = FALSE)$x)
    e <- drop(scale(e, center = center.e, scale = FALSE))
    me <- attr(e, "scaled:center") # mean of X_E



  if (!expand) {
    # Dont Expand X's if expand=FALSE. use user supplied design matrix
    Phi_j_list <- lapply(split(seq(group), group),
                         function(j) scale(x[, j, drop = FALSE],
                                           center = center.x, scale = FALSE))

    Phi_j <- do.call(cbind, Phi_j_list)
    main_effect_names <- vnames
    dimnames(Phi_j)[[2]] <- main_effect_names


    # mPhi_j is the column means of basis expansion of X
    mPhi_j <- sapply(Phi_j_list, attr, which = "scaled:center")



    # X_E x Phi_j
    XE_Phi_j_list <- lapply(Phi_j_list, function(i) scale(e * i, center = center.x, scale = FALSE))
    XE_Phi_j <- do.call(cbind, XE_Phi_j_list)
    interaction_names <- paste(main_effect_names, "E", sep = ":")
    dimnames(XE_Phi_j)[[2]] <- interaction_names

    # mXEPhi_j is the column means of the product of X and E
    mXE_Phi_j <- sapply(XE_Phi_j_list, attr, which = "scaled:center")



  } else {

    # Expand X's
    Phi_j_list <- lapply(
      seq_len(nvars),
      function(j) scale(basis(x[, j, drop = FALSE]),
                        center = center.x, scale = FALSE
          )
      )

    # mPhi_j is the column means of basis expansion of X
    mPhi_j <- sapply(Phi_j_list, attr, which = "scaled:center")

    ncols <- ncol(Phi_j_list[[1]]) # this is to get the number of columns for each expansion
    Phi_j <- do.call(cbind, Phi_j_list)
    main_effect_names <- paste(rep(vnames, each = ncols), rep(seq_len(ncols), times = nvars), sep = "_")
    dimnames(Phi_j)[[2]] <- main_effect_names

    # E x Phi_j
    # at this point E and Phi_j have already been centered
    # in the following line we center the product of centered variables
    # this is what we've seen being done in hierNet
    # https://github.com/cran/hierNet/blob/5ea5ef716d134879d16ff497f620174f1cc5d813/R/funcs.R#L35-L40
    XE_Phi_j_list <- lapply(Phi_j_list, function(i) scale(e * i, center = center.x, scale = FALSE))
    # XE_Phi_j_list <- lapply(Phi_j_list, function(i) e * i)
    XE_Phi_j <- do.call(cbind, XE_Phi_j_list)
    interaction_names <- paste(main_effect_names, "E", sep = ":")
    dimnames(XE_Phi_j)[[2]] <- interaction_names

    # mXEPhi_j is the column means of the product of X and E
    mXE_Phi_j <- sapply(XE_Phi_j_list, attr, which = "scaled:center")

  }

  # this is used for the predict function
  design <- cbind(Phi_j, "E" = e, XE_Phi_j)

  return(list(
    Phi_j_list = Phi_j_list,
    Phi_j = Phi_j,
    mPhi_j = if (center.x) mPhi_j else 0,
    XE_Phi_j_list = XE_Phi_j_list,
    XE_Phi_j = XE_Phi_j,
    mXE_Phi_j = if (center.x) mXE_Phi_j else 0,
    main_effect_names = main_effect_names,
    interaction_names = interaction_names,
    E = e, # this is the possibly centered e
    mE = me,
    design = design,
    ncols = if (expand) ncols else sapply(Phi_j_list, ncol)
  ))
}



