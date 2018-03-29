# x: the matrix containing the predictors y: the response N: the number of slicing
# schemes. If missing, the slicing schemes with G=3, ..., [log(n)]+1 is used

# return: a vector containing K^G.

k.filter <- function(x, y, N = NULL, nslices = NULL, slicing.scheme = NULL, response.type = "continuous",
                     method = "fused") {
  method <- match.arg(arg = method, choices = c("fused", "single"))
  response.type <- match.arg(arg = response.type, choices = c(
    "continuous", "discrete",
    "categorical"
  ))

  if (!is.null(slicing.scheme)) {
    if (method == "fused") {
      warning("A slicing scheme is given. Using a single Kolmogorov filter.")
    }
    obj <- k.filter.single(x = x, y = y, slicing.scheme = slicing.scheme)
  } else {
    if (response.type == "categorical") {
      y <- factor(y)
      obj <- k.filter.single(x = x, y = y)
    }
    if (response.type == "continuous" | response.type == "discrete") {
      n <- nrow(x)
      if (is.null(N)) {
        N <- ceiling(log(n)) - 2
      }
      if (is.null(nslices)) {
        nslices <- 3:(N + 2)
      }
      obj <- fused.k.filter(x = x, y = y, N = N, nslices = nslices)
    }
  }
  obj
}


k.filter.single <- function(x, y, slicing.scheme = NULL, nslices = NULL) {
  if (is.factor(y)) {
    y.dm <- y
  }
  if (!is.factor(y)) {
    if (is.null(slicing.scheme)) {
      if (is.null(nslices)) {
        stop("If y is not a factor, either slicing.scheme or nslices should be specified")
      }
      slicing.scheme <- quantile(quantile(y, seq(1 / nslices, 1 - 1 / nslices, 1 / nslices)))
    }
    K <- length(slicing.scheme) + 1
    slicing.scheme <- c(min(y) - 1, slicing.scheme, max(y) + 1)
    y.dm <- cut(y, slicing.scheme, label = 1:K, right = F)
  }
  y.dm <- as.numeric(y.dm)
  K <- length(unique(y.dm))
  p <- ncol(x)
  ks.stat <- matrix(0, p, K * (K - 1) / 2)

  nclass <- 0
  for (j in 1:(K - 1)) {
    for (l in (j + 1):K) {
      nclass <- nclass + 1
      for (i in 1:p) {
        ks.stat[i, nclass] <- ks.test(x[y.dm == j, i], x[y.dm == l, i])$statistic
      }
    }
  }
  ks.stat.max0 <- apply(ks.stat, 1, max)
  k.rank <- rank(-ks.stat.max0, ties.method = "max")
  list(k.stat = ks.stat.max0, k.rank = k.rank)
}


fused.k.filter <- function(x, y, N = NULL, nslices = NULL) {

  # browser()
  n <- nrow(x)
  p <- ncol(x)

  if (is.null(N)) {
    N <- ceiling(log(n)) - 2
  }
  if (is.null(nslices)) {
    nslices <- 3:(N + 2)
  }

  ks.stat.single <- matrix(0, N, p)
  ks.stat.max <- rep(0, p)

  for (K in nslices) {
    slicing.scheme <- quantile(y, seq(0, 1, 1 / K))
    slicing.scheme[1] <- slicing.scheme[1] - 1
    slicing.scheme[K + 1] <- slicing.scheme[K + 1] + 1
    y.dm <- cut(y, slicing.scheme, labels = c(1:K), right = F)
    ks.stat <- matrix(0, p, K * (K - 1) / 2)

    nclass <- 0
    for (j in 1:(K - 1)) {
      for (l in (j + 1):K) {
        nclass <- nclass + 1
        for (i in 1:p) {
          ks.stat[i, nclass] <- ks.test(x[y.dm == j, i], x[y.dm == l, i])$statistic
        }
      }
    }
    ks.stat.max0 <- apply(ks.stat, 1, max)
    ks.stat.single[K - 2, ] <- ks.stat.max0
    ks.stat.max <- ks.stat.max + ks.stat.max0
  }

  k.rank <- rank(-ks.stat.max, ties.method = "max")

  list(
    k.stat = ks.stat.max, k.stat.single = ks.stat.single, N = N, nslices = nslices,
    k.rank = k.rank
  )
}
