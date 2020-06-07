#' OASIS Brain Data
#'
#' A dataset containing a small subset of the OASIS Brain Data project. Noise
#' variables have been added to increase the number of predictors and show the
#' utility of the sail package.
#'
#' @format A list with 3 elements: \describe{ \item{x}{a matrix of dimension
#'   \code{136 x 30} with the following colums:\describe{\item{Age}{Patient's
#'   age} \item{EDUC}{Education} \item{MMSE}{Mini-Mental State Exam}
#'   \item{eTIV}{Estimated Total Intracranial Volume} \item{nWBV}{Normalized
#'   Whole Brain Volume} \item{ASF}{Atlas scaling factor}
#'   \item{noise[1-24]}{24 independent standard normal noise variables}}}\item{y}{a
#'   numeric vector of length 136 representing the right Hippocampal volume for
#'   each patient}\item{e}{a binary 0/1 vector of length 136, representing
#'   Dementia status. 0: Non-demented, 1: Demented} }
#' @source \url{https://github.com/stnava/RMI/tree/master/tomfletcher}
#' @source \url{http://www.oasis-brains.org/}
#' @examples
#' oasis
"oasis"


#' Simulated Data Used in Bhatnagar et al. (2018+) Paper
#'
#' A dataset containing simulated data used in the accompanying paper to this
#' package
#'
#' @details The code used to simulate the data is available at
#'   \url{https://github.com/sahirbhatnagar/sail/blob/master/data-raw/SIMULATED_data.R}.
#'   See \code{\link{gendata}} for more details. The true model is given by
#'   \deqn{Y = f1(X1) + f2(X2) +  f3(X3) + f4(X4) +  E * (2 +  f3(X3) +
#'   f4(X4))} where \describe{\item{}{f1(t)=5t}\item{}{f2(t)=3(2t -
#'   1)^2}\item{}{f3(t)= 4sin(2pi*t) /
#'   (2-sin(2pi*t)}\item{}{f4(t)=6(0.1sin(2pi*t) + 0.2cos(2pi*t) +
#'   0.3sin(2pi*t)^2 + 0.4cos(2pi*t)^3 + 0.5sin(2pi*t)^3)}}
#' @format A list with 7 elements: \describe{ \item{x}{a matrix of dimension
#'   \code{100 x 20} where rows are observations and columns are
#'   predictors}\item{y}{a numeric response vector of length 100 }\item{e}{a
#'   numeric exposure vector of length 100}\item{f1,f2,f3,f4}{the true
#'   functions} }
#' @references Lin, Y., & Zhang, H. H. (2006). Component selection and smoothing
#'   in multivariate nonparametric regression. The Annals of Statistics, 34(5),
#'   2272-2297.
#' @references Huang J, Horowitz JL, Wei F. Variable selection in nonparametric
#'   additive models (2010). Annals of statistics. Aug 1;38(4):2282.
#' @references Bhatnagar SR, Yang Y, Greenwood CMT. Sparse additive interaction
#'   models with the strong heredity property (2018+). Preprint.
#' @examples
#' sailsim
"sailsim"
