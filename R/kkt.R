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
  dl <- -1/(1 + exp(r))
  dl
}

dls <- function(r, delta) {
  dl <- -r
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
margin <- function(b0, betaE, beta, gamma, alpha, y, phij, xe_phij, e, df, loss = c("ls", "logit")) {

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
  } else { r <- as.vector(y) - link }

  fun <- paste("d", loss, sep = "")

  # this is a matrix of -R = -(Y - hat(Y)) of dimension nobs x nlambda
  dMat <- apply(r, c(1, 2), eval(fun))

  # this is the gradient for beta0, this should be of length nlambda
  # t(dMat) %*% matrix(1, nrow = nobs) / nobs

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


tt <- KKT(b0 = fit$a0, betaE = fit$bE, beta = fit$beta, gamma = fit$gamma,
          alpha = fit$alpha, y = DT$y, phij = fit$Phi_j, xe_phij = fit$XE_Phi_j,
          e = DT$e, df = fit$df,
          lambda = fit$lambda, lambda2 = fit$lambda2, group = fit$group,
          we = fit$we, wj = fit$wj, wje = fit$wje, thr = 1e-5, loss = "ls")

KKT <- function(b0, betaE, beta, gamma, alpha, y, phij, xe_phij, e, df,
                lambda, lambda2, group, we, wj, wje, thr, loss = c("ls","logit")) {

  loss <- match.arg(loss)
  bn <- as.integer(max(group))
  dl <- margin(b0 = b0, betaE = betaE, beta = beta, gamma = gamma, alpha = alpha,
               y = y, phij = phij, xe_phij = xe_phij, e = e, df = df, loss = loss)

  browser()
  B <- matrix(NA, ncol = length(lambda))
  ctr <- 0
  for (l in 1:length(lambda)) {
    for (g in 1:bn) {
      ind <- (group == g)
      dl_norm <- sqrt(crossprod(dl[ind, l], dl[ind, l]))
      b_norm <- sqrt(crossprod(beta[ind, l], beta[ind, l]))
      if (b_norm != 0) {
        AA <- dl[ind, l] + beta[ind, l] * lambda[l] * as.vector(pf[g]/b_norm)
        if (sum(abs(AA)) >= thr) {
          cat("violate at b != 0", sum(abs(AA)), "\n")
          ctr <- ctr + 1
        }
      } else {
        BB <- dl_norm - pf[g] * lambda[l]
        if (BB > thr) {
          cat("violate at b = 0", BB, "\n")
          ctr <- ctr + 1
        }
      }
    }
  }
  # cat("# of violations", ctr/length(lambda), "\n")
  return(ctr/length(lambda))
}


