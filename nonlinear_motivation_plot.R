sim_nonlinear <- function(n = 400, p = 100, E = rnorm(n), betaE = 1, SNR = 1) {
  # this is modified from SPAM Ravikumar et all JRSSB
  # n = 200
  # p = 10
  # corr = 1
  E = rbinom(n, 1, 0.5)
  # covariates
  X <- replicate(n = p, runif(n, 0, 1))
  colnames(X) <- paste0("X", seq_len(p))

  f1 <- function(x) sin((7/8)* pi * x)
  f2 <- function(x) x^3 + 1.5 * (x - 0.5)^2
  f3 <- function(x) -dnorm(x, 0.5, 0.8)
  f4 <- function(x) sin(exp(-0.5 * x))

  X1 <- X[,1]
  X2 <- X[,2]
  X3 <- X[,3]
  X4 <- X[,4]

  # error
  error <- stats::rnorm(n)

  Y.star <- f1(X1)  +
    f2(X2) +
    f3(X3) +
    f4(X4) +
    betaE * E +
    E * f3(X3) +
    E * f4(X4)

  Y.star <- 18.5 + 10*(X1) + 1.4*E + 10 * E * f1(X1)

  k <- sqrt(stats::var(Y.star) / (SNR * stats::var(error)))

  Y <- Y.star + as.vector(k) * error

  return(list(x = X, y = Y, e = E, f1 = f1(X1),Y.star = Y.star,
              f2 = f2(X2), f3 = f3(X3), f4 = f4(X4), betaE = betaE,
              f1.f = f1, f2.f = f2, f3.f = f3, f4.f = f4))

}

# f1 <- function(x) sin((7/8)* pi * x)
# f2 <- function(x) x^3 + 1.5 * (x - 0.5)^2
# f3 <- function(x) -dnorm(x, 0.5, 0.8)
# f4 <- function(x) sin(exp(-0.5 * x))

DT <- sim_nonlinear()
x1 <- DT$x[,1]
# f1 <- DT$f1
e <- DT$e
# f1 <- function(x) 22*sin((7/8)* pi * x)
# curve(f1)

lin_pred <- DT$Y.star
unexposed_index <- which(e==0)
exposed_index <- which(e==1)

e0 <- lin_pred[unexposed_index]
e1 <- lin_pred[exposed_index]

x1e0 <- x1[unexposed_index]
x1e1 <- x1[exposed_index]

min.length.top <- range(lin_pred)[1] ; max.length.top <- range(lin_pred)[2]
# dev.off()
#
pdf(file="mcgillsims/figures/nonlinear_motivation.pdf",width=11,height=8)
par(mfrow=c(1,1), tcl=-0.5, family="serif",omi=c(0.2,0.2,0,0))
par(mai=c(1,1,1,0.2))
plot(x1, lin_pred,
     pch = 19,
     ylab = "Obesity",
     xlab = "DNA Methylation",
     col = cbbPalette()[c(4,7)],
     bty="n",
     xaxt="n",
     type = "n",
     cex.lab = 2,
     cex.axis = 2,
     cex = 2,
     main = "",
     cex.main = 2.5,
     ylim = c(min.length.top-3, max.length.top+5))
axis(1, labels = T, cex.axis = 2)
lines(x1e0[order(x1e0)],e0[order(x1e0)], col = cbbPalette()[c(4)], lwd = 4)
lines(x1e1[order(x1e1)],e1[order(x1e1)], col = cbbPalette()[c(7)], lwd = 4)
legend("topleft", c("Controls", "GD Affected Pregnancy"),
                   col = cbbPalette()[c(4,7)], pch = 19, cex = 2, bty = "n")
rug(x1)
dev.off()

# mai = c(1,1,2,0))
# oma = c(1,1,2,1))
plot_apoe_inter(cv_obj = cvfit, X = X, original_name = "X60", sail_name = "X19", xlab =  "supramarginal gyrus right", legend = T,
                apoe = FALSE)
plot_apoe_inter(cv_obj = cvfit, X = X, original_name = "X60", sail_name = "X19", xlab =  "supramarginal gyrus right", legend = F,
                apoe = TRUE, ylab = "", main = "APOE e4 = 1")
dev.off()
