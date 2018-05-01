rm(list=ls())
dev.off()
devtools::load_all("/home/sahir/git_repositories/sail/")
# source("my_sims/model_functions.R")

# set.seed(123)
# DT <- make_gendata_Paper_not_simulator(n = 100, p = 20, corr = 0.0,
#                                        betaE = 2, SNR = 2,
#                                        parameterIndex = 1)
# DT$e
#
# DT <- gendata4(n = 100, p = 20, E = stats::rbinom(100, 2,0.5), betaE = 2, SNR = 4)
DT <- sail::gendata(n = 200, p = 20, corr = 0, SNR = 2, betaE = 1, parameterIndex = 2)
DT$f1.f
DT$x
DT$not_causal
DT$causal
cbind(DT$x, E = DT$e)

trainIndex <- drop(caret::createDataPartition(DT$y, p = 0.5, list = FALSE, times = 1))



library(doMC)
registerDoMC(cores = 8)
f.basis <- function(i) splines::bs(i, degree = 5)
fit <- sail(x = DT$x, y = DT$y, e = DT$e,
            # alpha = 0.5,
            strong = FALSE,
            verbose = 2,
            basis = f.basis)

plot(fit)

# f.basis <- function(i) i
cvfit <- cv.sail(x = DT$x, y = DT$y, e = DT$e,
                 strong = FALSE,
                 # alpha = 0.5,
                 verbose = 1,
                 basis = f.basis, nfolds = 10, parallel = TRUE)

cvfit2 <- cv.sail(x = DT$x, y = DT$y, e = DT$e,
                 strong = TRUE,
                 # alpha = 0.5,
                 verbose = 2,
                 basis = f.basis, nfolds = 5, parallel = TRUE)

cvfit$sail.fit
# plot cv-error curve
plot(cvfit)
plot(cvfit2)

cvfit$cvm[which(cvfit$lambda.min==cvfit$lambda)]
cvfit2$cvm[which(cvfit2$lambda.min==cvfit2$lambda)]

# non-zero estimated coefficients at lambda.1se
predict(cvfit, type = "non")
predict(cvfit, type = "nonzero", s="lambda.min")



new_x <- replicate(20, rnorm(50))
new_e <- rnorm(50, sd = 0.5)
predict(cvfit, newx = new_x, newe = new_e, s = c(5,3,1,0.5))


expect_equal(apply(linpred, 2, median), object$linear.predictors,
             tolerance = tol,
             check.attributes = FALSE)


# plot main effect for X1
plotMain(cvfit$sail.fit, x = sailsim$x, xvar = "X3",legend.position = "topright",
         s = cvfit$lambda.min, f.truth = sailsim$f3)

plotInter(cvfit$sail.fit, x = sailsim$x, xvar = "X4",
          f.truth = sailsim$f4.inter,
          s = cvfit$lambda.min,
          title_z = "Estimated")


# with basis function -----------------------------------------------------

f.basis <- function(i) splines::bs(i, degree = 3)
# f.basis <- function(i) i
# DT$y <- scale(DT$y, center = TRUE, scale = FALSE) centering response does change the solution path
data("sailsim")
system.time(fit <- sail(x = sailsim$x[,1:10,drop=F], y = sailsim$y, e = sailsim$e,
                        basis = f.basis,
                        # lambda.factor = 0.001,
                        # thresh = 1e-3,
                        # dfmax = 2,
                        nlambda = 100, verbose = 1))

plot(fit)
predict(fit, s=c(2.15, 0.32, 0.40), type="nonzero")
pacman::p_load(doMC)
registerDoMC(cores = 8)
system.time(cvfit <- cv.sail(x = sailsim$x[,1:20,drop=F], y = sailsim$y, e = sailsim$e,
                             basis = f.basis,type.measure = "mse",
                             nlambda = 100, verbose = 1, nfolds = 10,
                             parallel = TRUE, grouped = TRUE))
plot(fit)
plot(cvfit)

str(fit)
coef(fit)
predict(fit)
coef(cvfit, s = "lambda.1se")
coef(cvfit, s = "lambda.min")
coef(cvfit, s = "lambda.min")[nonzero(coef(cvfit, s = "lambda.min")),,drop=F]

devtools::load_all("/home/sahir/git_repositories/sail/")
system.time(fit <- sail(x = DT$x, y = DT$y, e = DT$e, basis = f.basis, alpha = 0.5,
                        # lambda = c(2.10166227446906, 1.82791866705889, 1.66553141063233, 1.51757019050867,
                        # 1.31990515959739, 1.04600226643003, 0.953078307983438)#,
                        penalty.factor = c(0, 1, 1, rep(1, 2*ncol(DT$x) - 3),1)
))
plot(fit)
fit
coef(fit)
fit$lambda[c(3,6,8,10,13,18,20)] %>% dput

registerDoMC(cores = 8)
system.time(cvfit <- cv.sail(x = DT$x, y = DT$y, e = DT$e, basis = f.basis, alpha = 0.5,
                             parallel = TRUE, nfolds = 10, verbose = T,
                             nlambda = 100,
                             lambda = c(2.10166227446906, 1.82791866705889, 1.66553141063233, 1.51757019050867,
                                        1.31990515959739, 1.04600226643003, 0.953078307983438)))#,
# penalty.factor = c(1, 0.0, 0.0, rep(1, 2*ncol(DT$x) - 3),.8)))
cvfit$lambda[c(3,6,8,10,13,18,20)] %>% dput
plot(cvfit)
plot(cvfit$sail.fit)
coef(cvfit, s = "lambda.min")
predict(cvfit, s = 0.3)
# with pre-specified design -----------------------------------------------

x_df <- as.data.frame(DT$x)
head(x_df)
x_df$race <- factor(sample(1:5, nrow(x_df), replace = TRUE))
str(x_df)
x <- model.matrix(~ 0 +  bs(X1) + bs(X2) + ns(X3, 5) + poly(X4, 6) +
                    X5 + poly(X6,2) + race, data = x_df)

head(x)
devtools::load_all("/home/sahir/git_repositories/sail/")
fit <- sail(x = x, y = DT$y, e = DT$e, expand = FALSE, group = attr(x, "assign"),
            alpha = 0.5, penalty.factor = c(1, 0.0, 0.0, rep(1, 2*length(unique(attr(x, "assign"))) - 3),.8))
coef(fit)
plot(fit)
fit
fit$active

split(attr(x, "assign"), attr(x, "assign"))
registerDoMC(cores = 5)
cvfit <- cv.sail(x = x, y = DT$y, e = DT$e, expand = FALSE, alpha = 0.5, group = attr(x, "assign"),
                 parallel = TRUE, nfolds = 5, verbose = T, nlambda = 25,
                 penalty.factor = c(1, 0.0, 0.0, rep(1, 2*length(unique(attr(x, "assign"))) - 3),.8))

plot(cvfit)
coef(cvfit, s = "lambda.1se")

# DT <- gendata4(n = 200, p = 50, E = rnorm(200), betaE = 2, SNR = 3)
# Rprof(tmp <- tempfile())

# if df = 1 and degree = 1, then dont do any bsplines expansion, use original data

# design_sail(x = DT$x, e = DT$e, nvars = 25, vnames = paste0("X",1:25), df = 5, degree = 3)

registerDoMC(cores = 5)

DT$y <- scale(DT$y, center = TRUE, scale = FALSE)
foldid <- sample(1:10,size=length(DT$y),replace=TRUE)

cvfit.8 <- cv.sail(x = DT$x, y = DT$y, e = DT$e, df = 5, degree = 3, thresh = 1e-4, maxit = 1000,alpha = .8,
                   parallel = TRUE, foldid = foldid, nfolds = 10, verbose = T, nlambda = 100)


x <- runif(100)
B <- splines::bs(x[order(x)], degree = 5, intercept = T, Boundary.knots = c(0,1))
matplot(x, B, type="l", lwd=2)

xx <- replicate(10, rnorm(20))
colnames(xx) <- paste0("V",1:ncol(xx))
xx_df <- as.data.frame(xx)
colnames(xx_df)
str(xx_df)
x_des <- model.matrix(~0+ bs(V1) + V2 + ns(V3, 5), data = xx_df)
head(x_des)
split(x_des, split(attr(x_des, "assign"), attr(x_des, "assign")))
group <- attr(x_des, "assign")
lapply(split(seq(group), group), function(i) x_des[,i,drop=FALSE])

# x <- seq(0, 1, length=n)
# x2 <- truncnorm::rtruncnorm(n, a = -2, b = 2)
# x <- rnorm(n)
# B <- bs(x, df = NULL, degree=5, intercept = FALSE)


# ncol(B)
# attr(B,"knots")
# attr(B,"Boundary.knots")
#
# B1 <- bs(x, df = NULL, degree=5, intercept = FALSE)
# B1 <- apply(B1, 2, function(i) i[order(x)])
# matplot(x[order(x)], B1, type="n", lwd=2, main = "df=NULL, degree=5")
# for (i in 1:5) lines(x[order(x)], B1[,i], col = cbbPalette()[i],lwd=2)

n <- 1000
x <- truncnorm::rtruncnorm(n, a = 0, b = 1)
B2 <- bs(x, df = 5, degree=3, intercept = FALSE)
knots <- attr(B2,"knots")
bound.knots <- attr(B2,"Boundary.knots")
B2 <- apply(B2, 2, function(i) i[order(x)])
matplot(x[order(x)], B2, type="n", lwd=2, main = sprintf("df=5, degree=3, inner.knots at c(33.33%%, 66.66%%) percentile"),
        ylab = "bs(x)", cex.lab=1.5, xlab = "x")
for (i in 1:5) lines(x[order(x)], B2[,i], col = cbbPalette()[c(2,3,4,6,7)][i],lwd=3)
abline(v = knots, lty = 2, col = "gray", lwd = 4)


# dgp <- sin(2*pi*x[order(x)])
dgp <- bs(x[order(x)], df = 5, degree = 3, intercept = FALSE)
dgp <- bs(x[order(x)], df = 5, degree = 3, intercept = FALSE)

dgp <-  -dnorm(x[order(x)], 0.5, 0.8)
b1 <- rnorm(ncol(dgp))
b2 <- rnorm(ncol(dgp))
# y <- dgp + rnorm(n, sd=.1)
y <- dgp %*% betas  + rnorm(n, sd=.5)

model <- lm(y~B1)
summary(model)
plot(x[order(x)], y, cex=.25, col="grey")
lines(x[order(x)], fitted(model), lwd=2)
# lines(x[order(x)], dgp, col="red", lty=2, lwd = 2)
lines(x[order(x)], dgp %*% betas, col="red", lty=2, lwd = 2)

model <- lm(y~B2)
summary(model)
plot(x[order(x)], y, cex=.25, col="grey")
lines(x[order(x)], fitted(model), lwd=2)
# lines(x[order(x)], dgp, col="red", lty=2, lwd = 2)
lines(x[order(x)], dgp %*% betas, col="red", lty=2, lwd = 2)

# Gendata -----------------------------------------------------------------

# set.seed(1234)
# DT <- gendata(n = 200, p = 25, SNR = 3, betaE = 2, df = NULL, degree = 5)
DT <- gendata2(n = 200, p = 20, SNR = 2,
               # E = rbinom(n = 200, size = 1, prob = 0.5),
               # E = rnorm(200),
               betaE = 2, corr = 0.2)#, df = 5)
# DT <- gendata4(n = 200, p = 50, E = rnorm(200), betaE = 2, SNR = 3)
# Rprof(tmp <- tempfile())

# if df = 1 and degree = 1, then dont do any bsplines expansion, use original data

# design_sail(x = DT$x, e = DT$e, nvars = 25, vnames = paste0("X",1:25), df = 5, degree = 3)

registerDoMC(cores = 5)

DT$y <- scale(DT$y, center = TRUE, scale = FALSE)
foldid <- sample(1:10,size=length(DT$y),replace=TRUE)

system.time(
  cvfit.8 <- cv.sail(x = DT$x, y = DT$y, e = DT$e, df = 5, degree = 3, thresh = 1e-4, maxit = 1000,alpha = .8,
                     parallel = TRUE, foldid = foldid, nfolds = 10, verbose = T, nlambda = 100)
)

system.time(
  cvfit.5 <- cv.sail(x = DT$x, y = DT$y, e = DT$e, df = 5, degree = 3, thresh = 1e-4, maxit = 1000,alpha = .5,
                     parallel = TRUE, foldid = foldid, nfolds = 10, verbose = T, nlambda = 100)
)

system.time(
  cvfit.2 <- cv.sail(x = DT$x, y = DT$y, e = DT$e, df = 5, degree = 3, thresh = 1e-4, maxit = 1000, alpha = .2,
                     parallel = TRUE, basis.intercept = FALSE,
                     foldid = foldid, nfolds = 10, verbose = T, nlambda = 100)
)

dev.off()

c(bottom, left, top, right)
ylim <- range(cvfit.8$cvm, cvfit.5$cvm,cvfit.2$cvm)
xlim <- range(log(cvfit.8$lambda), log(cvfit.5$lambda),log(cvfit.2$lambda), 3)

pdf("figure/alpha.pdf", width = 11, height = 8)
par(mfrow=c(2,2))
par(mar=c(5.1, 4.1, 2.6, 2.1))
plot(cvfit.8);text(-2,35,expression(paste(alpha,"=",0.8)), cex = 2);
plot(cvfit.5);text(-2.8,35,expression(paste(alpha,"=",0.5)), cex = 2);
plot(cvfit.2);text(-2.7,38,expression(paste(alpha,"=",0.2)), cex = 2);
plot(log(cvfit.8$lambda),cvfit.8$cvm,pch=19,col=cbbPalette()[c(2)],xlab="log(Lambda)",ylab=cvfit.8$name,
     ylim = ylim, xlim = xlim)
points(log(cvfit.5$lambda),cvfit.5$cvm,pch=19,col=cbbPalette()[c(3)])
points(log(cvfit.2$lambda),cvfit.2$cvm,pch=19,col=cbbPalette()[c(4)])
legend("bottomright",legend=c(expression(paste(alpha,"=",0.8)),
                              expression(paste(alpha,"=",0.5)),
                              expression(paste(alpha,"=",0.2))),
       pch=19,col=cbbPalette()[c(2:4)],
       cex = 1.5, bty = "n")
dev.off()




c(cvfit.8$cvm, cvfit.5$cvm,cvfit.2$cvm)[which.min(c(cvfit.8$cvm, cvfit.5$cvm,cvfit.2$cvm))];which.min(c(cvfit.8$cvm, cvfit.5$cvm,cvfit.2$cvm))

cvfit <- cvfit.5
coef(cvfit, s = "lambda.min")[nonzero(coef(cvfit, s = "lambda.min")),,drop=F]
coef(cvfit, s = "lambda.1se")[nonzero(coef(cvfit, s = "lambda.1se")),,drop=F]

dev.off()

pdf("figure/gendata2_cvfit.pdf", width = 11, height = 8)
plot(cvfit)
dev.off()

pdf("figure/gendata2_fit.pdf", width = 11, height = 8)
plot(cvfit$sail.fit)
dev.off()

dev.off()
par(mfrow=c(2,2))
for (i in 1:4){
  xv <- paste0("X",i)
  ind <- cvfit$sail.fit$group == which(cvfit$sail.fit$vnames == xv)
  design.mat <- cvfit$sail.fit$design[,cvfit$sail.fit$main.effect.names[ind],drop = FALSE]
  # f.truth <- design.mat %*% DT$b1
  f.truth <- DT[[paste0("f",i)]]
  plotMain(object = cvfit$sail.fit, xvar = xv, s = cvfit$lambda.min, f.truth = f.truth, legend.position = "topleft")
}
# f.truth <- design.mat %*% DT$b2
# coef(cvfit$sail.fit, s = 0.08)

cvfit
cvfit$lambda.min
cvfit$lambda.1se
coef(cvfit)
plot(cvfit)

f3.persp = function(X, E) {
  # E * as.vector(DT$betaE) + DT$f3.f(X) +
  E * DT$f3.f(X)
}

f4.persp = function(X, E) {
  # E * as.vector(DT$betaE) + DT$f4.f(X) +
  E * DT$f4.f(X)
}

i = 3
xv <- paste0("X",i)
plotInter(object = cvfit$sail.fit, xvar = xv, s = cvfit$lambda.min,
          f.truth = f3.persp,
          simulation = T,
          npoints = 40)

i = 4
xv <- paste0("X",4)
plotInter(object = cvfit$sail.fit, xvar = xv, s = cvfit$lambda.min,
          f.truth = f4.persp,
          simulation = T,
          npoints = 40)




plot(DT$y, predict(cvfit, s = "lambda.min"))
abline(a=0,b=1)
crossprod(DT$y - predict(cvfit, s = "lambda.1se"))
cor(DT$y, predict(cvfit, s = "lambda.min"))^2

class(fit)
fit <- cvfit$sail.fit
plot(fit)
fit
tt <- KKT(b0 = fit$a0, betaE = fit$bE, beta = fit$beta, gamma = fit$gamma,
          alpha = fit$alpha, y = DT$y, phij = fit$Phi_j, xe_phij = fit$XE_Phi_j,
          e = DT$e, df = fit$df,
          lambda = fit$lambda, lambda2 = fit$lambda2, group = fit$group,
          we = fit$we, wj = fit$wj, wje = fit$wje, thr = 1e-1, loss = "ls")





system.time(
  fit <- sail(x = DT$x, y = DT$y, e = DT$e, df = 1, degree = 1, thresh = 1e-3,
              maxit = 1000,
              alpha = .1,
              # dfmax = 15,
              verbose = TRUE, nlambda = 100)
)

fit
plot(fit)
fit$active

tt <- KKT(b0 = fit$a0, betaE = fit$bE, beta = fit$beta, gamma = fit$gamma,
          alpha = fit$alpha, y = DT$y, phij = fit$Phi_j, xe_phij = fit$XE_Phi_j,
          e = DT$e, df = fit$df,
          lambda = fit$lambda, lambda2 = fit$lambda2, group = fit$group,
          we = fit$we, wj = fit$wj, wje = fit$wje, thr = 1e-1, loss = "ls")


# pacman::p_load(glinternet)
# dx <- cbind(DT$e,DT$x)
# glfit <- glinternet.cv(dx, DT$y, interactionCandidates=1, numLevels = rep(1, ncol(dx)))
# plot(glfit)
# coef(glfit)


# Hierbasis example -------------------------------------------------------


set.seed(1)

# Generate the points x.
n <- 100
p <- 30

x <- matrix(rnorm(n*p), ncol = p)

# A simple model with 3 non-zero functions.
e <- rnorm(n = n, sd=0.5)
y <- rnorm(n, sd = 0.1) + e * sin(x[, 1]) + x[, 2] + (x[, 3])^3

devtools::load_all()

system.time(
  cvfit <- cv.sail(x = x, y = y, e = e, df = 3, degree = 3, thresh = 1e-4,
                   maxit = 1000,
                   alpha = .1,
                   parallel = TRUE,
                   nfolds = 10,
                   # dfmax = 15,
                   verbose = FALSE, nlambda = 100)
)

plot(cvfit)
coef(cvfit, s = "lambda.min")[nonzero(coef(cvfit, s = "lambda.min")),,drop=F]
coef(cvfit, s = "lambda.1se")[nonzero(coef(cvfit, s = "lambda.1se")),,drop=F]
fit <- cvfit$sail.fit

tt <- KKT(b0 = fit$a0, betaE = fit$bE, beta = fit$beta, gamma = fit$gamma,
          alpha = fit$alpha, y = y, phij = fit$Phi_j, xe_phij = fit$XE_Phi_j,
          e = e, df = fit$df,
          lambda = fit$lambda, lambda2 = fit$lambda2, group = fit$group,
          we = fit$we, wj = fit$wj, wje = fit$wje, thr = 1e-1, loss = "ls")
fit
plot(fit)
fit$active
fit$a0


# Gendata3 ----------------------------------------------------------------


devtools::load_all()

DT <- gendata3(n = 100, p = 50, betaE = 1.5, SNR = 3)
# DT <- gendata2(n = 200, p = 20, SNR = 2, betaE = 2)#, df = 5)
# Rprof(tmp <- tempfile())

# if df = 1 and degree = 1, then dont do any bsplines expansion, use original data

system.time(
  fit <- sail(x = DT$x, y = DT$y, e = DT$e, df = 5, degree = 3, thresh = 1e-4,
              maxit = 1000,
              alpha = .1,
              # dfmax = 15,
              verbose = TRUE, nlambda = 100)
)

system.time(
  cvfit <- cv.sail(x = DT$x, y = DT$y, e = DT$e, df = 5, degree = 3, thresh = 1e-4,
                   maxit = 1000,
                   alpha = .01,
                   parallel = TRUE,
                   nfolds = 10,
                   # dfmax = 15,
                   verbose = FALSE, nlambda = 100)
)


tt <- KKT(y = DT$y, e = DT$e, b0 = fit$a0, betaE = fit$bE, beta = fit$beta, gamma = fit$gamma,
          alpha = fit$alpha, phij = fit$Phi_j, xe_phij = fit$XE_Phi_j, df = fit$df,
          lambda = fit$lambda, lambda2 = fit$lambda2, group = fit$group,
          we = fit$we, wj = fit$wj, wje = fit$wje, thr = 1e-1, loss = "ls")

fit
plot(fit)
fit$active



# Hierbasis ---------------------------------------------------------------

pacman::p_load_gh("asadharis/HierBasis")

require(Matrix)

set.seed(1)

# Generate the points x.
n <- 100
p <- 30

x <- matrix(rnorm(n*p), ncol = p)

# A simple model with 3 non-zero functions.
y <- rnorm(n, sd = 0.1) + sin(x[, 1]) + x[, 2] + (x[, 3])^3

mod <- AdditiveHierBasis(x, y, nbasis = 50, max.lambda = 30,
                         beta.mat = NULL,
                         nlam = 50, alpha = 0.5,
                         lam.min.ratio = 1e-4, m.const = 3,
                         max.iter = 300, tol = 1e-4)

# Plot the individual functions.
xs <- seq(-3,3,length = 300)
plot(mod,1,30, type  ="l",col = "red", lwd = 2, xlab = "x", ylab = "f_1(x)",
     main = "Estimating the Sine function")
lines(xs, sin(xs), type = "l", lwd = 2)
legend("topleft", c("Estimated Function", "True Function"),
       col = c("red", "black"), lwd = 2, lty = 1)

plot(mod,2,30, type  ="l",col = "red", lwd = 2, xlab = "x", ylab = "f_2(x)",
     main = "Estimating the Linear function")
lines(xs, xs, type = "l", lwd = 2)
legend("topleft", c("Estimated Function", "True Function"),
       col = c("red", "black"), lwd = 2, lty = 1)

plot(mod,3,30, type  ="l",col = "red", lwd = 2, xlab = "x", ylab = "f_3(x)",
     main = "Estimating the cubic polynomial")
lines(xs, xs^3, type = "l", lwd = 2)
legend("topleft", c("Estimated Function", "True Function"),
       col = c("red", "black"), lwd = 2, lty = 1)



# other datasets ----------------------------------------------------------

pacman::p_load(TCGA2STAT)

rnaseq_os.ov <- getTCGA(disease="OV", data.type="RNASeq", type="RPKM", clinical=TRUE)

rnaseq_os.ov$clinical %>% colnames

Y <- rnaseq_os.ov$clinical[,"daystodeath"] %>% as.character() %>% as.numeric()
Y <- Y[!is.na(Y)]

rnaseq_os.ov$dat %>% str

top <- names(sort(rowSds(rnaseq_os.ov$dat), decreasing = T)[1:100])

X <- t(rnaseq_os.ov$dat[top,])
rownames(X)
rnaseq_os.ov$merged.dat %>% str
rnaseq_os.ov$clinical

library(MASS)
data("Boston")
help(Boston)
skimr::skim(as.data.frame(Boston))
Y <- Boston$medv
E <- Boston$ptratio
hist(E)
hist(Y)
X <- as.matrix(Boston[, -which(colnames(Boston) %in% c("medv","ptratio", "black","chas","zn",
                                                       "crim","rad"))])
skimr::skim(as.data.frame(X))

f.basis <- function(i) splines::bs(i, degree = 2)
fit_boston <- sail(x = X, y = Y, e = E, basis = f.basis, alpha = 0.1, center.e = FALSE)
plot(fit_boston)
coef(fit_boston)
fit_boston

registerDoMC(cores = 8)
f.basis <- function(i) splines::bs(i, degree = 2)
fit_boston <- cv.sail(x = X, y = Y, e = E, basis = f.basis, alpha = 0.5,
                      parallel = TRUE, nfolds = 10)
plot(fit_boston)
coef(fit_boston, s = "lambda.min")




# pacman::p_load(splines)
# pacman::p_load(microbenchmark)
# pacman::p_load(ggplot2)
#
# set.seed(12345)
# p <- 50
# n <- 200
# df <- 5
#
# # covariates
# X <- replicate(n = p, runif(n))
#
# # environment
# E <- rnorm(n = n, sd = 0.5)
#
# # coefficients: each is a vector of length df and corresponds to the expansion of X_j
# b0 <- 1.5
# b1 <- rnorm(n = df)
# b2 <- rnorm(n = df)
# b3 <- rnorm(n = df)
# b4 <- rnorm(n = df)
# b5 <- rnorm(n = df)
# bE1 <- rnorm(n = df)
# bE2 <- rnorm(n = df)
#
# # beta for environment
# bE <- 2
#
#
# # error
# error <- rnorm(n = n)
#
# Y <- b0 +
#   bs(X[,1], df = df) %*% b1  +
#   bs(X[,2], df = df) %*% b2 +
#   bs(X[,3], df = df) %*% b3 +
#   bs(X[,4], df = df) %*% b4 +
#   bs(X[,5], df = df) %*% b5 +
#   bE * E +
#   E * bs(X[,1], df = df) %*% bE1 +
#   E * bs(X[,2], df = df) %*% bE2 +
#   error
#
#
#
# y <- drop(Y)
# e <- drop(E)
# np <- dim(X)
# nobs <- as.integer(np[1])
# nvars <- as.integer(np[2])
#
# # group membership
# group <- rep(seq_len(nvars), each = df)
#
# # Expand X's
# Phi_j_list <- lapply(seq_len(nvars), function(j) splines::bs(X[,j], df = df))
# design_array <- array(NA, dim = c(n, df, p))
#
# for (i in 1:p) {
#   design_array[,,i] <- splines::bs(X[,i], df = df)
# }
#
# design_array[,,1][1:5,1:5]
# Phi_j_list[[1]][1:5,1:5]
#
# Phi_j <- do.call(cbind, Phi_j_list)
#
# tt <- microbenchmark(docall = {
#   Phi_j_list <- lapply(seq_len(nvars), function(j) splines::bs(X[,j], df = df))
#   Phi_j <- do.call(cbind, Phi_j_list)
#
# },
# array = {
#   design_array <- array(NA, dim = c(n, df, p))
#
#   for (i in 1:p) {
#     design_array[,,i] <- splines::bs(X[,i], df = df)
#   }
# }, times = 500)
#
# autoplot(tt)
#
# # Phi_j <- do.call(cbind, lapply(seq_len(nvars), function(j) splines::ns(x[,j], df = df)))
# head(Phi_j)
# main_effect_names <- paste(paste0("X", group), rep(seq_len(df), times = nvars), sep = "_")
# dimnames(Phi_j)[[2]] <- main_effect_names
#
# # X_E x Phi_j
# XE_Phi_j <- e * Phi_j
# interaction_names <- paste(main_effect_names, "X_E", sep = ":")
# dimnames(XE_Phi_j)[[2]] <- interaction_names
#
# design <- cbind(Phi_j, "X_E" = e, XE_Phi_j)



