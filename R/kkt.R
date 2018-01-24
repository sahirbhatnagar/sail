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
