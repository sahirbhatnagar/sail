library(splines)

context("sail model fit strong and weak hierarchy with both packaged datasets")

data("sailsim")
data("oasis")

# data for weak hierarchy
DT <- gendata(n = 200, p = 10, corr = 0, SNR = 2, betaE = 2, parameterIndex = 2)

f.basis <- function(i) splines::bs(i, degree = 3)


fit_sim <- try(sail(x = sailsim$x, y = sailsim$y, e = sailsim$e, basis = f.basis),
               silent = TRUE)

fit_oasis <- try(sail(x = oasis$x, y = oasis$y, e = oasis$e, basis = f.basis),
                 silent = TRUE)

fit_weak <- try(sail(x = DT$x, y = DT$y, e = DT$e, basis = f.basis, strong = FALSE),
                silent = TRUE)

test_that("no error in fitting sail for both simulated and real data", {

  expect_false(inherits(fit_sim, "try-error"))
  expect_false(inherits(fit_oasis, "try-error"))
  expect_false(inherits(fit_weak, "try-error"))
  expect_is(fit_sim, "sail")
  expect_is(fit_oasis, "sail")
  expect_is(fit_weak, "sail")

  expect_equivalent(class(coef(fit_sim)), "dgCMatrix")
  expect_equal(dim(coef(fit_sim))[[1]], dim(fit_sim$design)[[2]] + 1)
  expect_equal(dim(coef(fit_sim))[[2]], sum(fit_sim$converged))

})


context("sail model fit with penalty factor")

test_that("no error in fitting sail with different penalty.factor", {

  fit_pf <- try(sail(x = sailsim$x, y = sailsim$y, e = sailsim$e, basis = f.basis,
                      penalty.factor = c(0, 1, 0.4, 0.6,0.7, rep(1, 2*ncol(sailsim$x) - 4))),
                 silent = TRUE)

  fit_pf_weak <- try(sail(x = sailsim$x, y = sailsim$y, e = sailsim$e, basis = f.basis, strong = FALSE,
                     penalty.factor = c(0, 1, 0.4, 0.6,0.7, rep(1, 2*ncol(sailsim$x) - 4))),
                silent = TRUE)

  expect_false(inherits(fit_pf, "try-error"))
  expect_false(inherits(fit_pf_weak, "try-error"))

})




context("sail model fit with user defined design")

test_that("no error in fitting sail with user defined design", {


  x_df <- as.data.frame(sailsim$x)
  x_df$race <- factor(sample(1:5, nrow(x_df), replace = TRUE))
  x <- model.matrix(~ 0 +  bs(X1) + bs(X2) + ns(X3, 5) + poly(X4, 6) +
                      X5 + poly(X6,2) + race, data = x_df)

  fit_design <- try(sail(x = x, y = sailsim$y, e = sailsim$e, expand = FALSE,
                         group = attr(x, "assign")),
                    silent = TRUE)

  fit_design_weak <- try(sail(x = x, y = sailsim$y, e = sailsim$e, expand = FALSE,
                         group = attr(x, "assign"), strong = FALSE),
                    silent = TRUE)

  fit_design_pf <- try(sail(x = x, y = sailsim$y, e = sailsim$e, expand = FALSE,
                            group = attr(x, "assign"),
                            penalty.factor = c(1, 0.0, 0.0, rep(1, 2*length(unique(attr(x, "assign"))) - 3),.8)),
                       silent = TRUE)

  fit_design_pf_weak <- try(sail(x = x, y = sailsim$y, e = sailsim$e, expand = FALSE, strong = FALSE,
                            group = attr(x, "assign"),
                            penalty.factor = c(1, 0.0, 0.0, rep(1, 2*length(unique(attr(x, "assign"))) - 3),.8)),
                       silent = TRUE)

  expect_false(inherits(fit_design, "try-error"))
  expect_false(inherits(fit_design_weak, "try-error"))
  expect_false(inherits(fit_design_pf, "try-error"))
  expect_false(inherits(fit_design_pf_weak, "try-error"))


})



context("solution path plot")

test_that("plot solution path", {
  disp_solution_path <- function() plot(fit_sim)
  vdiffr::expect_doppelganger("sail solution path", disp_solution_path)
})



test_that("plot solution path weak", {
  disp_solution_path <- function() plot(fit_weak)
  vdiffr::expect_doppelganger("sail solution path weak", disp_solution_path)
})




