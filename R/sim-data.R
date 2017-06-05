pacman::p_load(splines)

set.seed(12345)
p <- 50
n <- 300
df <- 5

# covariates
X <- replicate(n = p, runif(n))

# environment
E <- rnorm(n = n, sd = 0.5)

# coefficients: each is a vector of length df and corresponds to the expansion of X_j
b0 <- 1
b1 <- rnorm(n = df)
b2 <- rnorm(n = df)
b3 <- rnorm(n = df)
b4 <- rnorm(n = df)
b5 <- rnorm(n = df)
bE1 <- rnorm(n = df)
bE2 <- rnorm(n = df)

# beta for environment
bE <- 2


# error
error <- rnorm(n = n)

Y <- b0 +
  bs(X[,1], df = df) %*% b1  +
  bs(X[,2], df = df) %*% b2 +
  bs(X[,3], df = df) %*% b3 +
  bs(X[,4], df = df) %*% b4 +
  bs(X[,5], df = df) %*% b5 +
  bE * E +
  E * bs(X[,1], df = df) %*% bE1 +
  E * bs(X[,2], df = df) %*% bE2 +
  error


