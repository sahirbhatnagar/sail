expect_matrix <- function(x) inherits(x,"matrix")

context("generate simulated data")

DT1 <- gendata(n = 75, p = 100, corr = 0, betaE = 2, SNR = 1, parameterIndex = 1)
DT2 <- gendata(n = 75, p = 100, corr = 0, betaE = 2, SNR = 1, parameterIndex = 2)
DT3 <- gendata(n = 75, p = 100, corr = 0, betaE = 2, SNR = 1, parameterIndex = 3)
DT4 <- gendata(n = 75, p = 100, corr = 0, betaE = 2, SNR = 1, parameterIndex = 4)
DT5 <- gendata(n = 75, p = 100, corr = 0, betaE = 2, SNR = 1, parameterIndex = 5)

test_that("testing gendata", {

  expect_true(expect_matrix(DT1$x))
  expect_true(expect_matrix(DT2$x))
  expect_true(expect_matrix(DT3$x))
  expect_true(expect_matrix(DT4$x))
  expect_true(expect_matrix(DT5$x))

})
