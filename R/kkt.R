# nocov start
# rm(list = ls())
# source("R/sim-data.R")
# l2norm <- function(x) sqrt(sum(x^2))
#
# alpha <- 0.5
# w_e <- 1 ; w_j <- 1
#
# R <- Y - b0
#
# # this is not the right formula with weights.. but ignore for now. need to figure out the weights after
# lambda_max <- (1 / (n * (1 - alpha))) * max( (1 / w_e) * (crossprod(e, R) ),
#                                              max((1 / w_j) * sapply(Phi_j_list, function(i) l2norm(crossprod(i, R)))))
#
# R_without_j <- Y - b0
#
#
#
# design_array <- array(1:8, dim = c(2,2,2))
#
# apply(design_array, c(1,2), FUN = function(i) i * 2)
#
# phi <- design_array[,,1]
#
# all.equal(phi %*% matrix(b1),
#      phi %*% b1)
#
#
# Wj <- E * phi
#
# all.equal(phi %*% b1 + 2 * bE * Wj %*% b1,
#           (phi + 2 * bE * Wj) %*% b1)
#


dlogit <- function(r, delta) {
  dl <- -1 / (1 + exp(r))
  dl
}

dls <- function(r, delta) {
  dl <- -r
  dl
}

dwls <- function(r, delta,weights) {
  dl <- -r
  dl
}




# m1 <- gglasso(x = bardet$x, y = bardet$y, group = group1, loss = "ls")
# violations <- gglasso:::KKT(b0 = m1$b0, beta = m1$beta, y = bardet$y, x = bardet$x, lambda = m1$lambda,
#                       pf = rep(sqrt(5),20), group = group1, thr = 1e-3, loss = "ls")
# names(fit)
#
# tt <- margin(b0 = fit$a0, betaE = fit$bE, beta = fit$beta, gamma = fit$gamma,
#              alpha = fit$alpha,
#              y = DT$y, phij = fit$Phi_j, xe_phij = fit$XE_Phi_j,
#              e = DT$e, df = fit$df, loss = "ls")

# this gives -R = -(Y - hat(Y))

margin <- function(b0, betaE, beta, gamma, alpha, y, phij, xe_phij, e, df, loss = c("ls", "logit",'wls')) {
  loss <- match.arg(loss)
  nobs <- length(as.vector(y))
  beta <- as.matrix(beta)
  alpha <- as.matrix(alpha)
  gamma <- as.matrix(gamma)
  b0MAT <- matrix(rep(b0, nobs), nrow = nobs, byrow = TRUE)

  # dim(beta)
  # dim(b0MAT)
  # dim(xe_phij)
  # dim(phij)
  # dim(alpha)
  # dim(beta)
  # b0[1:5]
  # b0MAT[1:5,1:5]
  # browser()

  link <- b0MAT + phij %*% beta + matrix(e) %*% matrix(betaE, nrow = 1) + xe_phij %*% alpha

  if (loss %in% c("logit")) {
    r <- y * link
  }
  else {
    r <- as.vector(y) - link
  }

  fun <- paste("d", loss, sep = "")

  # this is a matrix of -R = -(Y - hat(Y)) of dimension nobs x nlambda
  dMat <- apply(r, c(1, 2), eval(fun))

  #   e + dim(xe_phij %*% beta)
  # dim(gamma)
  #
  #
  #   if (loss %in% c("logit", "sqsvm", "hsvm")) {
  #     yxdMat <- t(x) %*% (dMat * y)/nobs
  #   } else yxdMat <- t(x) %*% dMat/nobs
  #   yxdMat

  return(dMat)
}




KKT <- function(b0, betaE, beta, gamma, alpha, y, phij, xe_phij, e, df,weights,
                lambda, lambda2, group, we, wj, wje, thr, loss = c("ls", "logit",'wls')) {
  loss <- match.arg(loss)
  bn <- as.integer(max(group))
  nobs <- length(as.vector(y))

  # this gives -R = -(Y - hat(Y)) for loss = "ls"
  dl <- margin(
    b0 = b0, betaE = betaE, beta = beta, gamma = gamma, alpha = alpha, y = y,
    phij = phij, xe_phij = xe_phij, e = e, df = df, loss = loss
  )

  # KKT for beta0 -----------------------------------------------------------
  # this is the gradient for beta0, this should be of length nlambda
  B0 <- t(dl) %*% ( matrix(1, nrow = nobs) )/ nobs

    ctr <- 0

  for (l in 1:length(lambda)) {
    if (abs(B0[l, ]) > thr) {
      message("violate at b0 ", B0[l, ], " lambda=", lambda[l], "\n")
      ctr <- ctr + 1
    }
  }

  warning("% of violations for beta0", ctr / length(lambda), "\n")


  # KKT for betaE -----------------------------------------------------------


  # results for betaE
  # BE <- matrix(NA, ncol = length(lambda))

  ctr <- 0
  for (l in 1:length(lambda)) {
    xdMat_betaE <- e + rowSums(do.call(cbind, lapply(seq_along(unique(group)), function(j) {
      index <- (group == j)
      as.matrix(gamma[j, l] * (xe_phij[, index, drop = F] %*% beta[index, l, drop = F]))
    })))

    dl_norm_betaE <- t(xdMat_betaE) %*% (dl[, l, drop = FALSE]) / nobs

    if (betaE[l] == 0) {
      BE <- dl_norm_betaE / (-lambda[l] * (1 - lambda2) * we)
      if (abs(BE) >  thr) {
        warning("violate at bE = 0", abs(BE), " lambda=", lambda[l], "\n")
        ctr <- ctr + 1
      }
    } else {
      BE <- as.vector(dl_norm_betaE + lambda[l] * (1 - lambda2) * we * sign(betaE[l]))
      if (abs(BE) > thr) {
        warning("violate at bE != 0", abs(BE), " lambda=", lambda[l], "\n")
        ctr <- ctr + 1
      } # else {warning("no violation at bE != 0", BE, " lambda=",lambda[l], "\n")}
    }
  }

  warning("% of violations for betaE", ctr / length(lambda), "\n")
  # return(list(v = ctr/length(lambda)))

  # browser()


  # KKT for gamma -----------------------------------------------------------

  ctr <- 0
  for (l in 1:length(lambda)) {
    for (g in 1:bn) {

      # browser()
      ind <- (group == g)

      xdMat_gammaj <- as.matrix(betaE[l] * (xe_phij[, ind, drop = F] %*% beta[ind, l, drop = F]))

      dl_norm_gammaj <- t(xdMat_gammaj) %*% (dl[, l, drop = FALSE] )/ nobs

      if (gamma[g, l] == 0) {
        BE <- dl_norm_gammaj / (-lambda[l] * (1 - lambda2) * we)
        if (abs(BE) > thr) {
          warning("violate at gamma_j = 0", BE, " lambda=", lambda[l], "\n")
          ctr <- ctr + 1
        }
      } else {
        BE <- as.vector(dl_norm_gammaj + lambda[l] * (1 - lambda2) * we * sign(gamma[g, l]))
        if (abs(BE) > thr) {
          warning("violate at gamma_j != 0", BE, " lambda=", lambda[l], "\n")
          ctr <- ctr + 1
        } # else {warning("no violation at bE != 0", BE, " lambda=",lambda[l], "\n")}
      }
    }
  }

  warning("% of violations for gamma", ctr / length(lambda), "\n")



  # KKT for theta -----------------------------------------------------------

  # results for beta (aka theta)
  # B <- matrix(NA, ncol = length(lambda))

  ctr <- 0
  for (l in 1:length(lambda)) {
    for (g in 1:bn) {

      # browser()
      ind <- (group == g)

      # this is the pre-multiplier of R in eq (17)
      xdMat <- phij[, ind, drop = F] + gamma[g, l] * betaE[l] * xe_phij[, ind, drop = F]

      # this is the first part of eq (17)
      dl_prod <- t(xdMat) %*% (dl[, l, drop = FALSE]) / nobs

      dl_norm <- l2norm(dl_prod)

      # l2 norm of beta for the expansion of covariate corresponding to group g
      b_norm <- l2norm(beta[ind, l])

      if (b_norm != 0) {
        AA <- dl_prod + beta[ind, l] * lambda[l] * (1 - lambda2) * wj[g] / b_norm
        # warning(AA,"\n")
        if (sum(abs(AA)) >= thr) {
          warning("violate at bX != 0", sum(abs(AA)), " lambda=", lambda[l], "\n")
          ctr <- ctr + 1
        }
      } else {
        BB <- dl_norm - lambda[l] * (1 - lambda2) * wj[g]
        # warning(BB,"\n")
        if (BB > thr) {
          warning("violate at bX = 0", BB, " lambda=", lambda[l], "\n")
          ctr <- ctr + 1
        }
      }
    }
  }
  warning("% of violations for bTheta", ctr / length(lambda), "\n")
  return(ctr / length(lambda))
}
# nocov end
