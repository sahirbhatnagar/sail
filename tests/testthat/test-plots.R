set.seed(1234) # we set the seed so that the cv error curves remain identical when testing (randomness is introduced due to CV folds)
data("sailsim")

f.basis <- function(i) splines::bs(i, degree = 5)

context("effect plot")

test_that("plot for main effects and interaction effects", {

  testthat::skip_on_cran()
  testthat::skip_on_appveyor()

  library(doParallel)
  doParallel::registerDoParallel(cores = 2)
  cvfit_sim <- cv.sail(x = sailsim$x, y = sailsim$y, e = sailsim$e, basis = f.basis,
                           nfolds = 10)

  disp_main_curve <- function() plotMain(cvfit_sim$sail.fit, x = sailsim$x, xvar = "X4",
                                         s = cvfit_sim$lambda.min, f.truth = sailsim$f4,
                                         legend.position = "topright")
  vdiffr::expect_doppelganger("main effect curve", disp_main_curve)


  disp_inter_curve <- function() plotInter(cvfit_sim$sail.fit, x = sailsim$x, xvar = "X4",
                                          f.truth = sailsim$f4.inter,
                                          s = cvfit_sim$lambda.min,
                                          title_z = "Estimated")
  vdiffr::expect_doppelganger("interaction persp plot", disp_inter_curve)


})
