######################################
#' R Source code file for creating simulated dataset to be included in the sail package
#' Simulation scenario modified from
#' "Variable Selection in NonParametric Addditive Model" Huang, Horowitz and Wei, Ann. Stat
#' Author: Sahir Bhatnagar
#' Created: April 5, 2018
#' Updated:
#####################################

if (!require("pacman")) install.packages("pacman")
pacman::p_load(truncnorm)

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
                     corr = corr)
  # , E = truncnorm::rtruncnorm(n, a = -1, b = 1))

  return(DT)
}

gendataPaper <- function(n, p, corr = 0,
                         E = truncnorm::rtruncnorm(n, a = -1, b = 1),
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

set.seed(123456)
DT <- make_gendata_Paper_not_simulator(n = 100, p = 200, corr = 0.0,
                                       betaE = 2, SNR = 2,
                                       parameterIndex = 1)

sailsim <- list(x = DT$x, y = DT$y, e = DT$e)
devtools::use_data(sailsim, overwrite = TRUE)
