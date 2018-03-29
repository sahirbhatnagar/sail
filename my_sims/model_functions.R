## @knitr models

# devtools::load_all("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/")

make_easy_sail_model <- function(n, p, df, SNR, betaE) {

  # n = 400;p=10;df=5;SNR=4;betaE=2

  #==============

  # group membership
  group <- rep(seq_len(p), each = df)

  # covariates
  X <- replicate(n = p, runif(n))
  E <- rnorm(n = n, sd = 0.5)

  # Expand X's
  Phi_j <- do.call(cbind, lapply(seq_len(p), function(j) splines::bs(X[,j], df = df)))
  XE_Phi_j <- E * Phi_j

  main_effect_names <- paste(paste0("X", group), rep(seq_len(df), times = p), sep = "_")
  interaction_names <- paste(main_effect_names, "X_E", sep = ":")

  dimnames(Phi_j)[[2]] <- main_effect_names
  dimnames(XE_Phi_j)[[2]] <- interaction_names

  design <- cbind(Phi_j, "X_E" = E, XE_Phi_j)

  true_beta <- matrix(rep(0, ncol(design), ncol = 1))
  dimnames(true_beta)[[1]] <- colnames(design)
  # the first 5 main effects and the first 2 interactions are active
  true_beta[c(1:(5*df),(p * df + 2):(p * df + 1 + 2 * df) ),1] <- rnorm(n = 7 * df)
  true_beta["X_E",] <- betaE

  causal <- rownames(true_beta[which(true_beta[,1]!=0),,drop=F])
  not_causal <- setdiff(rownames(true_beta), causal)


  new_model(name = "sail_easy",
            label = sprintf("n = %s, p = %s, df = %s, SNR = %s, betaE = %s", n, p, df, SNR, betaE),
            params = list(n = n, p = p, df = df, SNR = SNR, betaE = betaE,
                          true_beta = true_beta, causal = causal, design = design, X = X, E = E,
                          not_causal = not_causal),
            simulate = function(n, design, true_beta, SNR, nsim) {
              # error <- stats::rnorm(n)
              error <- MASS::mvrnorm(nsim, mu = rep(0, n), Sigma = diag(n))
              Y.star <- as.numeric(design %*% true_beta)
              k <- sqrt(stats::var(Y.star) / (SNR * apply(error, 1, var)))
              error2 <- sweep(t(error), 2, k, FUN = "*")
              y <- Y.star + error2
              return(split(y, col(y))) # make each col its own list element
            })
}



make_gendata_Paper <- function(n, p, corr, betaE, SNR, lambda.type, parameterIndex) {

  main <- paste0("X", seq_len(p))
  vnames <- c(main, "E", paste0(main,":E"))

  if (parameterIndex == 1) { # 1a
    hierarchy = "strong" ; nonlinear = TRUE ; interactions = TRUE
    causal <- c("X1","X2","X3","X4","E","X3:E","X4:E")
  } else if (parameterIndex == 2) { # 1b
    hierarchy = "weak" ; nonlinear = TRUE ; interactions = TRUE
    causal <- c("X1","X2","E","X3:E","X4:E")
  } else if (parameterIndex == 3) { # 1c
    hierarchy = "none" ; nonlinear = TRUE ; interactions = TRUE
    causal <- c("X3:E","X4:E")
  } else if (parameterIndex == 4) { # 2
    hierarchy = "strong"; nonlinear = FALSE; interactions = TRUE
    causal <- c("X1","X2","X3","X4","E","X3:E","X4:E")
  } else if (parameterIndex == 5) { # 3
    hierarchy = "strong" ; nonlinear = TRUE ; interactions = FALSE
    causal <- c("X1","X2","X3","X4","E")
  }

  not_causal <- setdiff(vnames, causal)

  DT <- gendataPaper(n = n, p = p, SNR = SNR, betaE = betaE,
                     hierarchy = hierarchy, nonlinear = nonlinear, interactions = interactions,
                     corr = corr, E = truncnorm::rtruncnorm(n, a = -1, b = 1))

  # used for glmnet and lasso backtracking
  X_linear_design <- design_sail(x = DT$x, e = DT$e, nvars = p,
                                 vnames = paste0("X",1:p), degree = 1,
                                 center.x = FALSE, basis.intercept = FALSE)$design

  new_model(name = "gendata_Paper",
            label = sprintf("n = %s, p = %s, corr = %s, betaE = %s, SNR = %s, hierarchy = %s,
                            nonlinear = %s, interactions = %s, scenario = %s",
                            n, p, corr, betaE, SNR, hierarchy,
                            nonlinear, interactions, parameterIndex),
            params = list(n = n, p = p, corr = corr, betaE = betaE, SNR = SNR, lambda.type = lambda.type,
                          hierarchy = hierarchy, nonlinear = nonlinear, vnames = vnames,
                          interactions = interactions, causal = causal, X_linear_design = X_linear_design,
                          not_causal = not_causal, x = DT$x, e = DT$e, Y.star = DT$Y.star, EX = cbind(E=DT$e, DT$x)),
            simulate = function(n, Y.star, nsim) {
              error <- MASS::mvrnorm(nsim, mu = rep(0, n), Sigma = diag(n))
              # Y.star <- as.numeric(design %*% true_beta)
              k <- sqrt(stats::var(Y.star) / (SNR * apply(error, 1, var)))
              error2 <- sweep(t(error), 2, k, FUN = "*")
              y <- Y.star + error2
              return(split(y, col(y))) # make each col its own list element
            })

}



make_gendata_Paper_not_simulator <- function(n, p, corr, betaE, SNR, lambda.type, parameterIndex) {

  main <- paste0("X", seq_len(p))
  vnames <- c(main, "E", paste0(main,":E"))

  if (parameterIndex == 1) { # 1a
    hierarchy = "strong" ; nonlinear = TRUE ; interactions = TRUE
    causal <- c("X1","X2","X3","X4","E","X3:E","X4:E")
  } else if (parameterIndex == 2) { # 1b
    hierarchy = "weak" ; nonlinear = TRUE ; interactions = TRUE
    causal <- c("X1","X2","E","X3:E","X4:E")
  } else if (parameterIndex == 3) { # 1c
    hierarchy = "none" ; nonlinear = TRUE ; interactions = TRUE
    causal <- c("X3:E","X4:E")
  } else if (parameterIndex == 4) { # 2
    hierarchy = "strong"; nonlinear = FALSE; interactions = TRUE
    causal <- c("X1","X2","X3","X4","E","X3:E","X4:E")
  } else if (parameterIndex == 5) { # 3
    hierarchy = "strong" ; nonlinear = TRUE ; interactions = FALSE
    causal <- c("X1","X2","X3","X4","E")
  }

  not_causal <- setdiff(vnames, causal)

  DT <- gendataPaper(n = n, p = p, SNR = SNR, betaE = betaE,
                     hierarchy = hierarchy, nonlinear = nonlinear, interactions = interactions,
                     corr = corr, E = truncnorm::rtruncnorm(n, a = -1, b = 1))

  return(DT)
  # # used for glmnet and lasso backtracking
  # X_linear_design <- design_sail(x = DT$x, e = DT$e, nvars = p,
  #                                vnames = paste0("X",1:p), degree = 1,
  #                                center.x = FALSE, basis.intercept = FALSE)$design
  #
  # new_model(name = "gendata_Paper",
  #           label = sprintf("n = %s, p = %s, corr = %s, betaE = %s, SNR = %s, hierarchy = %s,
  #                           nonlinear = %s, interactions = %s, scenario = %s",
  #                           n, p, corr, betaE, SNR, hierarchy,
  #                           nonlinear, interactions, parameterIndex),
  #           params = list(n = n, p = p, corr = corr, betaE = betaE, SNR = SNR, lambda.type = lambda.type,
  #                         hierarchy = hierarchy, nonlinear = nonlinear, vnames = vnames,
  #                         interactions = interactions, causal = causal, X_linear_design = X_linear_design,
  #                         not_causal = not_causal, x = DT$x, e = DT$e, Y.star = DT$Y.star, EX = cbind(E=DT$e, DT$x)),
  #           simulate = function(n, Y.star, nsim) {
  #             error <- MASS::mvrnorm(nsim, mu = rep(0, n), Sigma = diag(n))
  #             # Y.star <- as.numeric(design %*% true_beta)
  #             k <- sqrt(stats::var(Y.star) / (SNR * apply(error, 1, var)))
  #             error2 <- sweep(t(error), 2, k, FUN = "*")
  #             y <- Y.star + error2
  #             return(split(y, col(y))) # make each col its own list element
  #           })

}



# nsim = 10;n=100;SNR=3
# error <- MASS::mvrnorm(nsim, mu = rep(0, n), Sigma = diag(n))
# Y.star <- rnorm(n)
# k <- sqrt(stats::var(Y.star) / (SNR * apply(error, 1, var)))
# error2 <- sweep(t(error), 2, k, FUN = "*")
# y <- Y.star + error2
# split(y, col(y))
#
# all.equal(tr[,1], t(error)[,1]*k[1])
