# pacman::p_load_gh('sahirbhatnagar/sail')
devtools::load_all()
gendataIntro <- function (n, p, corr = 0, E = truncnorm::rtruncnorm(n, a = -1, b = 1),
                          betaE = 2, SNR = 2, hierarchy = c("strong", "weak", "none"),
                          nonlinear = TRUE, interactions = TRUE, causal,
                          not_causal) {
  if (!requireNamespace("truncnorm", quietly = TRUE)) {
    stop("Package \"truncnorm\" needed for this function to simulate data. Please install it.",
         call. = FALSE)
  }

  hierarchy <- match.arg(hierarchy)
  W <- replicate(n = p, truncnorm::rtruncnorm(n, a = 0, b = 1))
  U <- truncnorm::rtruncnorm(n, a = 0, b = 1)
  V <- truncnorm::rtruncnorm(n, a = 0, b = 1)
  X1 <- (W[, 1] + corr * U)/(1 + corr)
  X2 <- (W[, 2] + corr * U)/(1 + corr)
  X3 <- (W[, 3] + corr * U)/(1 + corr)
  X4 <- (W[, 4] + corr * U)/(1 + corr)
  X <- (W[, 5:p] + corr * V)/(1 + corr)
  Xall <- cbind(X1, X2, X3, X4, X)
  colnames(Xall) <- paste0("X", seq_len(p))

  f1 <- function(x) 2 * x
  f2 <- function(x) -2 * (2 * x - 1)^2
  # f3 <- function(x) 2 * sin(2 * pi * x)/(2 - sin(2 * pi *
  #                                                  x))
  # f4 <- function(x) -2 * (0.1 * sin(2 * pi * x) + 0.2 *
  #                          cos(2 * pi * x) + 0.3 * sin(2 * pi * x)^2 + 0.4 *
  #                          cos(2 * pi * x)^3 + 0.5 * sin(2 * pi * x)^3)

  f3 <- function(x) 2 * sin(x)
  f4 <- function(x) -2 * exp(x)

  f3.inter <- function(x, e) e * f3(x)
  f4.inter <- function(x, e) e * f4(x)

  error <- stats::rnorm(n)

  Y.star <- f3(X3) + f4(X4) + betaE * E + f4.inter(X4, E)

  k <- sqrt(stats::var(Y.star)/(SNR * stats::var(error)))
  Y <- Y.star + as.vector(k) * error
  return(list(x = Xall, y = Y, e = E, Y.star = Y.star, f1 = f1(X1),
              f2 = f2(X2), f3 = f3(X3), f4 = f4(X4), betaE = betaE,
              f1.f = f1, f2.f = f2, f3.f = f3, f4.f = f4, X1 = X1,
              X2 = X2, X3 = X3, X4 = X4))
}

DT <- gendataIntro(n = 100, p = 20, corr = 0, SNR = 2, betaE = 2)
f.basis <- function(i) splines::bs(i, degree = 3)
fit <- sail(x = DT$x, y = DT$y, e = DT$e,
            basis = f.basis,
            alpha = 0.2,
            verbose = 1)

plot(fit)

library(doMC)
registerDoMC(cores = 8)
# data("sailsim")
f.basis <- function(i) splines::bs(i, degree = 3)
# f.basis <- function(i) i
cvfit <- cv.sail(x = DT$x, y = DT$y, e = DT$e,
                 alpha = 0.2,
                 verbose = 1,
                 dfmax = 20,
                 basis = f.basis, nfolds = 10, parallel = TRUE)
plot(cvfit)
predict(cvfit, s="lambda.min", type = "nonzero")

predict(cvfit, s="lambda.1se", type = "nonzero")

plot(cvfit$sail.fit)
abline(v = log(cvfit[["lambda.1se"]]))

cvfit$lambda



