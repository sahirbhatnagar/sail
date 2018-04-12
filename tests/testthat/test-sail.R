library(sail)
context("sail model fit")

data("sailsim")
data("oasis")
f.basis <- function(i) splines::bs(i, degree = 3)

test_that("no error in fitting sail for both simulated and real data", {

  fit_sim <- try(sail(x = sailsim$x, y = sailsim$y, e = sailsim$e, basis = f.basis),
                silent = TRUE)

  fit_oasis <- try(sail(x = oasis$x, y = oasis$y, e = oasis$e, basis = f.basis),
                   silent = TRUE)

  expect_false(inherits(fit_sim, "try-error"))
  expect_false(inherits(fit_oasis, "try-error"))

})
