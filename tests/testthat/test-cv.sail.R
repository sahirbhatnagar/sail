library(doParallel)
doParallel::registerDoParallel(cores = 2)
expect_matrix <- function(x) inherits(x,"matrix")
set.seed(1234) # we set the seed so that the cv error curves remain identical when testing (randomness is introduced due to CV folds)

context("cv.sail model fit, parallel, predict, plot with both packaged datasets")

testthat::skip_on_cran()
testthat::skip_on_appveyor()

data("sailsim")
data("oasis")
f.basis <- function(i) splines::bs(i, degree = 3)

cvfit_sim <- try(cv.sail(x = sailsim$x, y = sailsim$y, e = sailsim$e, basis = f.basis,
                         dfmax = 5, nfolds = 3),
                 silent = TRUE)

cvfit_oasis <- try(cv.sail(x = oasis$x, y = oasis$y, e = oasis$e, basis = f.basis,
                           dfmax = 5, nfolds = 6, alpha = 0.8),
                   silent = TRUE)

cvfit_parallel <- try(cv.sail(x = sailsim$x, y = sailsim$y, e = sailsim$e,
                              dfmax = 5,
                              basis = f.basis, nfolds = 3, parallel = TRUE),
                      silent = FALSE)

cvfit_parallel_weak <- try(cv.sail(x = sailsim$x, y = sailsim$y, e = sailsim$e, strong = TRUE,
                              dfmax = 5,
                              basis = f.basis, nfolds = 3, parallel = TRUE),
                      silent = FALSE)

new_x <- replicate(20, rnorm(50))
new_e <- rnorm(50, sd = 0.5)

test_that("no error in fitting cv.sail and parallel version for both simulated and real data", {

  expect_false(inherits(cvfit_sim, "try-error"))
  expect_false(inherits(cvfit_oasis, "try-error"))
  expect_false(inherits(cvfit_parallel, "try-error"))
  expect_false(inherits(cvfit_parallel_weak, "try-error"))
  expect_is(cvfit_sim, "cv.sail")
  expect_is(cvfit_oasis, "cv.sail")
  expect_is(cvfit_parallel, "cv.sail")

})


test_that("no error in predict for cv.sail", {

  expect_true(expect_matrix(predict(cvfit_sim, type = "nonzero", s = "lambda.min")))
  expect_true(expect_matrix(predict(cvfit_oasis, type = "nonzero", s = "lambda.min")))
  expect_true(expect_matrix(predict(cvfit_parallel, type = "nonzero", s = "lambda.min")))

  expect_true(expect_matrix(predict(cvfit_sim)))
  expect_true(expect_matrix(predict(cvfit_oasis)))
  expect_true(expect_matrix(predict(cvfit_parallel)))

  expect_equivalent(predict(cvfit_sim, s = Inf),
                    matrix(rep(cvfit_sim$sail.fit$a0[1], cvfit_sim$sail.fit$nobs)))

  expect_equivalent(predict(cvfit_oasis, s = Inf),
                    matrix(rep(cvfit_oasis$sail.fit$a0[1], cvfit_oasis$sail.fit$nobs)))

  expect_equivalent(predict(cvfit_parallel, s = Inf),
                    matrix(rep(cvfit_parallel$sail.fit$a0[1], cvfit_parallel$sail.fit$nobs)))

  expect_matrix(predict(cvfit_sim, newx = new_x, newe = new_e))
  expect_equal(dim(predict(cvfit_sim, newx = new_x, newe = new_e, s = c(5,3,1,0.5))), c(nrow(new_x),4))

  expect_equivalent(class(coef(cvfit_sim)), "dgCMatrix")
  expect_equal(dim(coef(cvfit_sim))[[1]], dim(cvfit_sim$sail.fit$design)[[2]] + 1)

})



context("cv plot")

test_that("check plots for cv curve", {

  disp_cv_curve <- function() plot(cvfit_sim)
  vdiffr::expect_doppelganger("cv.sail cv curve", disp_cv_curve)

})
