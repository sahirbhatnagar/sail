# pacman::p_load(splines)
# pacman::p_load(microbenchmark)
# pacman::p_load(ggplot2)
#
# set.seed(12345)
# p <- 50
# n <- 200
# df <- 5
#
# # covariates
# X <- replicate(n = p, runif(n))
#
# # environment
# E <- rnorm(n = n, sd = 0.5)
#
# # coefficients: each is a vector of length df and corresponds to the expansion of X_j
# b0 <- 1.5
# b1 <- rnorm(n = df)
# b2 <- rnorm(n = df)
# b3 <- rnorm(n = df)
# b4 <- rnorm(n = df)
# b5 <- rnorm(n = df)
# bE1 <- rnorm(n = df)
# bE2 <- rnorm(n = df)
#
# # beta for environment
# bE <- 2
#
#
# # error
# error <- rnorm(n = n)
#
# Y <- b0 +
#   bs(X[,1], df = df) %*% b1  +
#   bs(X[,2], df = df) %*% b2 +
#   bs(X[,3], df = df) %*% b3 +
#   bs(X[,4], df = df) %*% b4 +
#   bs(X[,5], df = df) %*% b5 +
#   bE * E +
#   E * bs(X[,1], df = df) %*% bE1 +
#   E * bs(X[,2], df = df) %*% bE2 +
#   error
#
#
#
# y <- drop(Y)
# e <- drop(E)
# np <- dim(X)
# nobs <- as.integer(np[1])
# nvars <- as.integer(np[2])
#
# # group membership
# group <- rep(seq_len(nvars), each = df)
#
# # Expand X's
# Phi_j_list <- lapply(seq_len(nvars), function(j) splines::bs(X[,j], df = df))
# design_array <- array(NA, dim = c(n, df, p))
#
# for (i in 1:p) {
#   design_array[,,i] <- splines::bs(X[,i], df = df)
# }
#
# design_array[,,1][1:5,1:5]
# Phi_j_list[[1]][1:5,1:5]
#
# Phi_j <- do.call(cbind, Phi_j_list)
#
# tt <- microbenchmark(docall = {
#   Phi_j_list <- lapply(seq_len(nvars), function(j) splines::bs(X[,j], df = df))
#   Phi_j <- do.call(cbind, Phi_j_list)
#
# },
# array = {
#   design_array <- array(NA, dim = c(n, df, p))
#
#   for (i in 1:p) {
#     design_array[,,i] <- splines::bs(X[,i], df = df)
#   }
# }, times = 500)
#
# autoplot(tt)
#
# # Phi_j <- do.call(cbind, lapply(seq_len(nvars), function(j) splines::ns(x[,j], df = df)))
# head(Phi_j)
# main_effect_names <- paste(paste0("X", group), rep(seq_len(df), times = nvars), sep = "_")
# dimnames(Phi_j)[[2]] <- main_effect_names
#
# # X_E x Phi_j
# XE_Phi_j <- e * Phi_j
# interaction_names <- paste(main_effect_names, "X_E", sep = ":")
# dimnames(XE_Phi_j)[[2]] <- interaction_names
#
# design <- cbind(Phi_j, "X_E" = e, XE_Phi_j)
