## @knitr models

# devtools::load_all()

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


# nsim = 10;n=100;SNR=3
# error <- MASS::mvrnorm(nsim, mu = rep(0, n), Sigma = diag(n))
# Y.star <- rnorm(n)
# k <- sqrt(stats::var(Y.star) / (SNR * apply(error, 1, var)))
# error2 <- sweep(t(error), 2, k, FUN = "*")
# y <- Y.star + error2
# split(y, col(y))
#
# all.equal(tr[,1], t(error)[,1]*k[1])
