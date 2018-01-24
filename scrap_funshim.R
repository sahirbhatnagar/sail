.datatable.aware=TRUE
rm(list=ls())
devtools::load_all()
data.table:::cedta()
# devtools::document()
set.seed(12345)
DT <- gendata(n = 400, p = 10, df = 5, SNR = 2)

pacman::p_load(mgcv)
dat <- data.frame(Y = DT$y, DT$x, E = DT$e)
colnames(dat)
b <- gam(Y ~ E + s(X1) + s(X2) + s(X3)+ s(X4)+ s(X5)+ s(X6)+ s(X7) +
           s(X8) + s(X9) + s(X10) +
           # s(X11) + s(X12) + s(X13)+ s(X14)+ s(X15)+ s(X16)+ s(X17)+
           # s(X18) + s(X19) + s(X20) +
           ti(X1, E) + ti(X2, E) + ti(X3, E) + ti(X4, E) + ti(X5, E) + ti(X6, E) +
           ti(X7, E) + ti(X8, E) + ti(X9, E) + ti(X10, E),  #+
           # ti(X11, E) + ti(X12, E) + ti(X13, E) + ti(X14, E) + ti(X15, E) + ti(X16, E) +
           # ti(X17, E) + ti(X18, E) + ti(X19, E) + ti(X20, E)
           data = dat, select=TRUE, method="REML")

summary(b)
coef(b)

b$model
# # DT <- gendata2(n = 400, p = 10, corr = 0, SNR = 3.11)
# cvfit <- cv.sail(x = DT$x, y = DT$y, e = DT$e, df = DT$df, maxit = 500, cores = 5,
#                  nfolds = 5,
#                  nlambda.gamma = 10, nlambda.beta = 10, nlambda = 100,
#                  thresh = 1e-5, center=TRUE, normalize=FALSE, verbose = T)
#
# plot(cvfit)
# coef(cvfit)
# coef(cvfit, s = "lambda.1se")
# coef(cvfit, s = "lambda.min")
# saveRDS(cvfit, file = "data/cvfit_gendata_df5_p10_SNR2_n400.rds")
# saveRDS(DT, file = "data/DT_gendata_df5_p10_SNR2_n400.rds")

cvfit <- readRDS("data/cvfit_gendata_df5_p10_SNR2_n400.rds")
DT <- readRDS("data/DT_gendata_df5_p10_SNR2_n400.rds")

# fit <- sail(x = DT$x, y = DT$y, e = DT$e, df = DT$df,
#                maxit = 500, nlambda.gamma = 5, nlambda.beta = 5,
#                # group.penalty = "SCAD",
#                nlambda = 25,
#                thresh = 1e-5, center=TRUE, normalize=FALSE, verbose = T)
#
# fit$lambda.beta
# coef(fit)
# plot.grpreg(fit)
# rm(fit)
#
#
#
# fit$
#
# matplot(t(as.matrix(fit$beta[1:5,])))
# library(magrittr)
# as.matrix(fit$beta[1:5,])
#
# args(cv.funshim)
# pacman::p_load(doParallel)
# doParallel::registerDoParallel(cores = 5)
#
#
# library(magrittr)
# cvfit %>% class
# dev.off()
#
# plot(cvfit)

# source("R/plot.R")

# fit = cvfit$funshim.fit
#
# BIC_fit <- vector(mode = "numeric", length = fit$nlambda)
# EBIC_fit <- vector(mode = "numeric", length = fit$nlambda)
# n = nrow(DT$x)
# p = ncol(DT$x) * DT$df
#
# for (lam in seq_len(fit$nlambda)) {
#
#   # lam = 1
#   #================
#   betas <- fit$beta[,lam, drop = F]
#   alphas <- fit$alpha[,lam, drop = F]
#   a0 <- fit$b0[lam]
#   betas.and.alphas <- rbind(a0, betas, alphas)
#
#   # fit$design %>% colnames()
#
#   BIC_fit[lam] <- log(crossprod(as.matrix(DT$y - cbind(1,fit$design) %*% betas.and.alphas))) +
#     log(n) / n * (fit$dfbeta[lam,] + fit$dfalpha[lam,])
#   EBIC_fit[lam] <- log(crossprod(as.matrix(DT$y - cbind(1,fit$design) %*% betas.and.alphas))) +
#     (log(n) / n) * (fit$dfbeta[lam,] + fit$dfalpha[lam,]) + 0.5 * (log(p) / n) * (fit$dfbeta[lam,] + fit$dfalpha[lam,])
# }
#
# dev.off()
# plot(BIC_fit, pch = 19, col = "red")
# plot(EBIC_fit, pch = 19, col = "red")
# cbind(rbind2(fit$beta[,which.min(BIC_fit), drop = F], fit$alpha[,which.min(BIC_fit), drop = F])[nonzero(rbind2(fit$beta[,which.min(BIC_fit), drop = F], fit$alpha[,which.min(BIC_fit), drop = F])),,drop=F],
# rbind2(fit$beta[,which.min(EBIC_fit), drop = F], fit$alpha[,which.min(EBIC_fit), drop = F])[nonzero(rbind2(fit$beta[,which.min(EBIC_fit), drop = F], fit$alpha[,which.min(EBIC_fit), drop = F])),,drop=F])

# pacman::p_load(mgcv)
# dat <- data.frame(Y = DT$y, DT$x, E = DT$e)
# colnames(dat)
#
# b <- gam(Y ~ E + s(X1) + s(X2) + s(X3)+ s(X4)+ s(X5)+ s(X6)+ s(X7) +
#            s(X8) + s(X9) + s(X10) +
#            s(X11) + s(X12) + s(X13)+ s(X14)+ s(X15)+ s(X16)+ s(X17)+
#            s(X18) + s(X19) + s(X20) +
#            ti(X1, E) + ti(X2, E) + ti(X3, E) + ti(X4, E) + ti(X5, E) + ti(X6, E) +
#            ti(X7, E) + ti(X8, E) + ti(X9, E) + ti(X10, E) +
#            ti(X11, E) + ti(X12, E) + ti(X13, E) + ti(X14, E) + ti(X15, E) + ti(X16, E) +
#            ti(X17, E) + ti(X18, E) + ti(X19, E) + ti(X20, E) , data = dat, select=TRUE, method="REML")
#
#
#
# summary(b)
# gam()
# plot(b)




# plot for gendata mai eff--------------------------------------------------------
# rm(list=ls())
# cvfit <- readRDS("data/cvfit_gendata_df3_p10_SNR2.rds")
# fit <- cvfit$funshim.fit
# set.seed(12345)
# DT <- gendata(n = 400, p = 10, df = 5, SNR = 4)
# pacman::p_load(doParallel)
# doParallel::registerDoParallel(cores = 5)
# cvfit <- cv.funshim(x = DT$x, y = DT$y, e = DT$e, df = DT$df, maxit = 1000, cores = 5,
#                     nfolds = 5,
#                     nlambda.gamma = 12, nlambda.beta = 12, nlambda = 144,
#                     thresh = 1e-5, center=TRUE, normalize=FALSE, verbose = T)
# saveRDS(cvfit, file = "data/cvfit_gendata_df5_p10_SNR4.rds")
# cvfit <- readRDS("data/cvfit_gendata_df5_p20_SNR4_n400.rds")
# DT <- readRDS("data/DT_gendata_df5_p20_SNR4_n400.rds")
cvfit <- readRDS("data/cvfit_gendata_df5_p10_SNR2_n400.rds")
DT <- readRDS("data/DT_gendata_df5_p10_SNR2_n400.rds")

fit <- cvfit$sail.fit
# source("R/plot.R")
pacman::p_load(cowplot)


# png(filename="gendata_cvfit.png",width=11,height=8,units="in",res=150)

pdf(file="~/Dropbox/jobs/hec/talk/gendata_cvfit.pdf",width=11,height=8)#,units="in",res=150)
plot(cvfit) + cowplot::panel_border()
dev.off()

colnames(DT$x) <- paste0("X",1:ncol(DT$x))
(true.coefs <- cbind(DT$b1,DT$b2,DT$b3,DT$b4,DT$b5, rep(0,5),rep(0,5),rep(0,5),rep(0,5),rep(0,5),
                     DT$bE1,DT$bE2,rep(0,5),rep(0,5),rep(0,5),rep(0,5),rep(0,5),rep(0,5),rep(0,5),rep(0,5)))
colnames(true.coefs) <- c(paste0("X",1:10), paste0("X",1:10,":X_E"))

dev.off()

lambda_type <- "lambda.1se"
trop <- RSkittleBrewer::RSkittleBrewer("trop")

i = "X1"

min.length.top = min(sapply(paste0("X",1:5), function(i) {
  range(as.vector(fit$design[,paste(i,seq_len(DT$df), sep = "_")] %*% true.coefs[,i,drop=F]),
        as.vector(fit$design[,paste(i,seq_len(DT$df), sep = "_")] %*%
                    coef(cvfit, s = lambda_type)[paste(i,seq_len(DT$df), sep = "_"),,drop=F]))[1]}))

max.length.top = max(sapply(paste0("X",1:5), function(i) {
  range(as.vector(fit$design[,paste(i,seq_len(DT$df), sep = "_")] %*% true.coefs[,i,drop=F]),
        as.vector(fit$design[,paste(i,seq_len(DT$df), sep = "_")] %*%
                    coef(cvfit, s = lambda_type)[paste(i,seq_len(DT$df), sep = "_"),,drop=F]))[2]}))

min.length.bottom = min(sapply(paste0("X",6:10), function(i) {
  range(as.vector(fit$design[,paste(i,seq_len(DT$df), sep = "_")] %*% true.coefs[,i,drop=F]),
        as.vector(fit$design[,paste(i,seq_len(DT$df), sep = "_")] %*%
                    coef(cvfit, s = lambda_type)[paste(i,seq_len(DT$df), sep = "_"),,drop=F]))[1]}))

max.length.bottom = max(sapply(paste0("X",6:10), function(i) {
  range(as.vector(fit$design[,paste(i,seq_len(DT$df), sep = "_")] %*% true.coefs[,i,drop=F]),
        as.vector(fit$design[,paste(i,seq_len(DT$df), sep = "_")] %*%
                    coef(cvfit, s = lambda_type)[paste(i,seq_len(DT$df), sep = "_"),,drop=F]))[2]}))



# main effects, red is estimated, black is truth, lambda.min
dev.off()

# png(filename="gendata_main_eff.png",width=11,height=8,units="in",res=150)
# mai = c(bottom, left, top, right)
pdf(file="~/Dropbox/jobs/hec/talk/gendata_main_eff.pdf",width=11,height=8)#,units="in",res=150)
par(mfrow=c(2,5), tcl=-0.5, family="serif", omi=c(0.2,0.2,0,0))
i = "X1"
par(mai=c(0.55,0.61,0.1,0.2))
plot(DT$x[,i], fit$design[,paste(i,seq_len(DT$df), sep = "_")] %*%
       coef(cvfit, s = lambda_type)[paste(i,seq_len(DT$df), sep = "_"),,drop=F],
     pch = 19,
     ylab = TeX(sprintf("$f(x_%s)$",as.numeric(gsub("X","",i)))),
     xlab = TeX(sprintf("$x_%s$",as.numeric(gsub("X","",i)))),
     col = trop[1],
     bty="n",
     xaxt="n",
     cex.lab = 2,
     cex.axis = 2,
     ylim = c(min.length.top, max.length.top))
axis(1, labels = F)
points(DT$x[,i], fit$design[,paste(i,seq_len(DT$df), sep = "_")] %*% true.coefs[,i,drop=F], pch = 19, col = trop[2])
legend(0.0,-1.3, c("Truth","Estimated"), col = trop[2:1], pch = 19, cex = 2, bty = "n")

for(i in paste0("X",2:5)) {
plot(DT$x[,i], fit$design[,paste(i,seq_len(DT$df), sep = "_")] %*%
       coef(cvfit, s = lambda_type)[paste(i,seq_len(DT$df), sep = "_"),,drop=F],
     pch = 19,
     ylab = TeX(sprintf("$f(x_%s)$",as.numeric(gsub("X","",i)))),
     xlab = TeX(sprintf("$x_%s$",as.numeric(gsub("X","",i)))),
     col = trop[1],
     bty="n",
     yaxt="n",
     xaxt="n",
     cex.lab = 2,
     cex.axis = 2,
     ylim = c(min.length.top, max.length.top))
points(DT$x[,i], fit$design[,paste(i,seq_len(DT$df), sep = "_")] %*% true.coefs[,i,drop=F], pch = 19, col = trop[2])
axis(1, labels = F)
axis(2, labels = F)
}


i = "X6"
par(mai=c(0.55,0.6,0.1,0.2))
plot(DT$x[,i], fit$design[,paste(i,seq_len(DT$df), sep = "_")] %*%
       coef(cvfit, s = lambda_type)[paste(i,seq_len(DT$df), sep = "_"),,drop=F],
     pch = 19,
     ylab = TeX(sprintf("$f(x_%s)$",as.numeric(gsub("X","",i)))),
     xlab = TeX(sprintf("$x_%s$",as.numeric(gsub("X","",i)))),
     col = trop[1],
     bty="n",
     cex.lab = 2,
     cex.axis = 2,
     ylim = c(min.length.top, max.length.top))
points(DT$x[,i], fit$design[,paste(i,seq_len(DT$df), sep = "_")] %*% true.coefs[,i,drop=F], pch = 19, col = trop[2])
# legend("bottomright", c("Truth","Estimated"), col = trop[2:1], pch = 19, cex = 1.5, bty = "n")

for(i in paste0("X",7:10)) {
  plot(DT$x[,i], fit$design[,paste(i,seq_len(DT$df), sep = "_")] %*%
         coef(cvfit, s = lambda_type)[paste(i,seq_len(DT$df), sep = "_"),,drop=F],
       pch = 19,
       ylab = TeX(sprintf("$f(x_{%s})$",as.numeric(gsub("X","",i)))),
       xlab = TeX(sprintf("$x_{%s}$",as.numeric(gsub("X","",i)))),
       col = trop[1],
       bty="n",
       yaxt="n",
       cex.lab = 2,
       cex.axis = 2,
       # xaxt="n",
       ylim = c(min.length.top, max.length.top))
  points(DT$x[,i], fit$design[,paste(i,seq_len(DT$df), sep = "_")] %*% true.coefs[,i,drop=F], pch = 19, col = trop[2])
  axis(2, labels = F)
}

# coef(cvfit, s = "lambda.1se")[nonzero(coef(cvfit, s = "lambda.1se")),,drop=F]
dev.off()


# gendata plot interactions effects ---------------------------------------------
i = "X1"

x <- seq(range(DT$x[,i])[1], range(DT$x[,i])[2], length.out = 30)
e <- seq(range(DT$e)[1], range(DT$e)[2], length.out = 30)
f.est <- function(x, e) { e * splines::bs(x, DT$df) %*%
    as.matrix(coef(cvfit, s = "lambda.1se")[paste0(i,"_", seq_len(DT$df),":X_E"),]) }
f.truth <- function(x, e) { e * splines::bs(x, DT$df) %*%
    true.coefs[,paste0(i,":X_E"), drop = F] }
z.est <- outer(x, e, f.est)
z.truth <- outer(x, e, f.truth)

op <- par(bg = "white")
z_range <- c(min(z.est, z.truth), max(z.est, z.truth))
# z_range <- c(-1.868998, 2.118003)

dev.off()
# png(filename="gendata_inter_X1.png",width=11,height=8,units="in",res=150)

pdf(file="~/Dropbox/jobs/hec/talk/gendata_inter_X1.pdf",width=11,height=8)#,units="in",res=150)
par(mfrow=c(1,2), tcl=-0.5, family="serif", omi=c(0.2,0.2,0,0))
par(mai=c(0.,0.,0.3,0.))
persp(x, e, z.truth,
      zlim = z_range,
      theta=30, phi=30,
      ltheta = 120, expand = 0.5,
      r=2, shade=0.3, axes=TRUE,scale=TRUE, box=T,
      nticks=5,
      # ticktype="detailed",
      col=trop[2],
      xlab=sprintf("f(x_%s)",as.numeric(gsub("X","",i))),
      ylab="X_E",
      zlab="Y", main="Truth")
persp(x, e, z.est,
      theta=30, phi=30,
      ltheta = 120, expand = 0.5,
      r=2, shade=0.3, axes=TRUE,scale=TRUE, box=T,
      nticks=5,
      zlim = z_range,
      # ticktype="detailed",
      col=trop[1],
      xlab=sprintf("f(x_%s)",as.numeric(gsub("X","",i))),
      ylab="X_E",
      zlab="Y", main="Estimated X_E*f(X_1)")
dev.off()




pdf(file="~/Dropbox/jobs/hec/talk/non_linear_example.pdf",width=11,height=8)#,units="in",res=150)
par(tcl=-0.5, family="serif", omi=c(0.2,0.2,0,0))
par(mai=c(0.,0.,0.3,0.))
persp(x, e, z.truth,
      zlim = z_range,
      theta=30, phi=30,
      ltheta = 120, expand = 0.5,
      r=2, shade=0.3, axes=TRUE,scale=TRUE, box=T,
      nticks=5,
      # ticktype="detailed",
      col=trop[2],
      cex.lab = 2.5,
      cex.main = 2.4,
      xlab="Variable explicative",
      ylab="Environnement",
      zlab="Variable réponse", main="Interaction non-linéaire")
dev.off()





coef(cvfit, s = lambda_type)[nonzero(coef(cvfit, s = lambda_type)),,drop=T]
true.coefs
true.coefs




# plot  for gendata2 main eff--------------------------------------------------------
cvfit <- readRDS("data/cvfit_gendata2_df5_p10_SNR3.rds")
rm(list=ls())
devtools::load_all()
pacman::p_load(truncnorm)
source("R/plot.R")
class(cvfit) <- "cv.sail"

coef(cvfit, s = "lambda.1se")
coef(cvfit, s = "lambda.min")[nonzero(coef(cvfit, s = "lambda.min")),,drop=F]
coef(cvfit, s = "lambda.1se")[nonzero(coef(cvfit, s = "lambda.1se")),,drop=F]

plot(cvfit)
set.seed(12345)
DT <- gendata2(n = 400, p = 10, corr = 0, SNR = 3.11)
fit <- cvfit$funshim.fit

lambda_type = "lambda.1se"

min.length.top = min(range(DT$f1(DT$x[,"X1"]),
                           as.vector(fit$design[,paste("X1",1:5, sep = "_")] %*%
                                       coef(cvfit, s = lambda_type)[paste("X1",1:5, sep = "_"),,drop=F])),
                     range(DT$f2(DT$x[,"X2"]),
                           as.vector(fit$design[,paste("X2",1:5, sep = "_")] %*%
                                       coef(cvfit, s = lambda_type)[paste("X2",1:5, sep = "_"),,drop=F])))
max.length.top = max(range(DT$f1(DT$x[,"X1"]),
                           as.vector(fit$design[,paste("X1",1:5, sep = "_")] %*%
                                       coef(cvfit, s = lambda_type)[paste("X1",1:5, sep = "_"),,drop=F])),
                     range(DT$f2(DT$x[,"X2"]),
                           as.vector(fit$design[,paste("X2",1:5, sep = "_")] %*%
                                       coef(cvfit, s = lambda_type)[paste("X2",1:5, sep = "_"),,drop=F])))

min.length.bot = min(range(DT$f3(DT$x[,"X3"]),
                           as.vector(fit$design[,paste("X3",1:5, sep = "_")] %*%
                                       coef(cvfit, s = lambda_type)[paste("X3",1:5, sep = "_"),,drop=F])),
                     range(DT$f4(DT$x[,"X4"]),
                           as.vector(fit$design[,paste("X4",1:5, sep = "_")] %*%
                                       coef(cvfit, s = lambda_type)[paste("X4",1:5, sep = "_"),,drop=F])))
max.length.bot = max(range(DT$f3(DT$x[,"X3"]),
                           as.vector(fit$design[,paste("X3",1:5, sep = "_")] %*%
                                       coef(cvfit, s = lambda_type)[paste("X3",1:5, sep = "_"),,drop=F])),
                     range(DT$f4(DT$x[,"X4"]),
                           as.vector(fit$design[,paste("X4",1:5, sep = "_")] %*%
                                       coef(cvfit, s = lambda_type)[paste("X4",1:5, sep = "_"),,drop=F])))


dev.off()
# png(filename="gendata2_main_eff.png",width=11,height=8,units="in",res=150)
pdf(file="~/Dropbox/jobs/hec/talk/gendata2_main_eff.pdf",width=17,height=9)#,units="in",res=150)
par(mfrow=c(2,2), tcl=-0.5, family="serif", omi=c(0.2,0.2,0,0))
trop <- RSkittleBrewer::RSkittleBrewer("tropical")
fit = cvfit$funshim.fit
i = "X1"
# mai = c(bottom, left, top, right)
par(mai=c(0.7,0.8,0.4,0))
plot(DT$x[,i], fit$design[,paste(i,1:5, sep = "_")] %*%
       coef(cvfit, s = lambda_type)[paste(i,1:5, sep = "_"),,drop=F],
     pch = 19,
     ylab = latex2exp::TeX("$f(x_1)$"),
     xlab = latex2exp::TeX("$x_1$"),
     xaxt="n", bty="n",
     col = trop[1],
     cex.axis=2,
     cex.lab=2,cex.main=2,
     main = latex2exp::TeX("$f(x_1) = 5x_1$"),
     ylim = c(min.length.top,max.length.top+0.15)
          # ylim = range(DT$f1(DT$x[,i]),
                  # as.vector(fit$design[,paste(i,1:5, sep = "_")] %*%
                              # coef(cvfit, s = lambda_type)[paste(i,1:5, sep = "_"),,drop=F]))
)
points(DT$x[,i], DT$f1(DT$x[,i]), pch = 19, col=trop[2])
axis(1, labels=F)
legend(0.5,.8, c("Truth","Estimated"), col = trop[2:1], pch = 19, cex = 2, bty = "n")


i = "X2"
plot(DT$x[,i], fit$design[,paste(i,1:5, sep = "_")] %*%
       coef(cvfit, s = lambda_type)[paste(i,1:5, sep = "_"),,drop=F],
     pch = 19,
     ylab = latex2exp::TeX("$f(x_2)$"),
     xlab = latex2exp::TeX("$x_2$"),
     xaxt="n", bty="n",
     yaxt="n",
     col = trop[1],
     cex.axis=2,
     cex.lab=2,cex.main=2,
     main = latex2exp::TeX("$f(x_2) = 3 (2x_2 - 1)^2$"),
     ylim = c(min.length.top,max.length.top+0.15)
     # ylim = range(DT$f2(DT$x[,i]),
     #              as.vector(fit$design[,paste(i,1:5, sep = "_")] %*%
     #                          coef(cvfit, s = lambda_type)[paste(i,1:5, sep = "_"),,drop=F]))
)
points(DT$x[,i], DT$f2(DT$x[,i]), pch = 19, col = trop[2])
axis(1, labels=F)
axis(2, labels=F)
# legend("bottomright", c("Truth","Estimated"), col = trop[2:1], pch = 19)


i = "X3"
par(mai=c(0.7,0.8,0.8,0))
plot(DT$x[,i], fit$design[,paste(i,1:5, sep = "_")] %*%
       coef(cvfit, s = lambda_type)[paste(i,1:5, sep = "_"),,drop=F],
     pch = 19,
     ylab = latex2exp::TeX("$f(x_3)$"),
     xlab = latex2exp::TeX("$x_3$"),
     bty="n",
     col = trop[1],
     cex.axis=2,
     cex.lab=2,cex.main=2,
     main = latex2exp::TeX("$f(x_3) = \\frac{4\\sin(2\\pi x_3)}{2 - \\sin(2\\pi x_3)}$"),
     ylim = c(min.length.bot,max.length.bot+0.15)
     # ylim = range(DT$f3(DT$x[,i]),
     #              as.vector(fit$design[,paste(i,1:5, sep = "_")] %*%
     #                          coef(cvfit, s = lambda_type)[paste(i,1:5, sep = "_"),,drop=F]))
)
points(DT$x[,i], DT$f3(DT$x[,i]), pch = 19, col = trop[2])
# legend("topright", c("Truth","Estimated"), col = trop[2:1], pch = 19)


i = "X4"
par(mai=c(0.7,0.8,0.8,0))
title <- as.list(expression(paste("f(x", phantom()[{
  paste("4")
}], "",") = 6(0.1sin(2",pi,"x", phantom()[{
  paste("4")
}], "",") + 0.2cos(2",pi,"x", phantom()[{
  paste("4")
}], "",") + 0.3sin(2",pi,"x", phantom()[{
  paste("4")
}], "",")","",
phantom()^{paste("2")},"+ "),
paste("     ","0.4cos(2",pi,"x", phantom()[{
  paste("4")
}], "",")","",
phantom()^{paste("3")} ,"+ 0.5sin(2",pi,"x", phantom()[{
  paste("4")
}], "",")","",
phantom()^{paste("3")},")")))

plot(DT$x[,i], fit$design[,paste(i,1:5, sep = "_")] %*%
       coef(cvfit, s = lambda_type)[paste(i,1:5, sep = "_"),,drop=F],
     pch = 19,
     ylab = latex2exp::TeX("$f(x_4)$"),
     xlab = latex2exp::TeX("$x_4$"),
     yaxt="n",
     bty="n",
     col = trop[1],
     main = NULL,
     cex.axis=2,
     cex.lab=2,cex.main=2,
     # main = latex2exp::TeX("$f(x) = 6(0.1\\sin(2\\pi x) + 0.2\\cos(2\\pi x) + 0.3\\sin(2\\pi x)^2 + 0.4\\cos(2\\pi x)^3 + 0.5\\sin(2\\pi x)^3)$"),
     ylim = c(min.length.bot,max.length.bot+0.15)
     # ylim = range(DT$f4(DT$x[,i]),
     #              as.vector(fit$design[,paste(i,1:5, sep = "_")] %*%
     #                          coef(cvfit, s = lambda_type)[paste(i,1:5, sep = "_"),,drop=F]))
)
mtext(do.call(expression, title ),side=3, line = c(1,-1) , cex = 2)
points(DT$x[,i], DT$f4(DT$x[,i]), pch = 19, col=trop[2])
axis(2, labels=F)
dev.off()
# legend("bottomright", c("Truth","Estimated"), col = trop[2:1], pch = 19)


# plot  for gendata2 interaction X4--------------------------------------------------------

i = "X4"
x <- seq(range(DT$x[,i])[1], range(DT$x[,i])[2], length.out = 30)
e <- seq(range(DT$e)[1], range(DT$e)[2], length.out = 30)
f.est <- function(x, e) { e * splines::bs(x, df = 5) %*%
    as.matrix(coef(cvfit, s = "lambda.1se")[paste0(i,"_", seq_len(5),":X_E"),]) }
f.truth <- function(x, e) { e * DT$f4(x)  }
z.est <- outer(x, e, f.est)
z.truth <- outer(x, e, f.truth)

op <- par(bg = "white")
z_range <- c(min(z.est, z.truth), max(z.est, z.truth))
# z_range <- c(-1.868998, 2.118003)

dev.off()
# png(filename="gendata2_inter_X4.png",width=11,height=8,units="in",res=150)
pdf(file="~/Dropbox/jobs/hec/talk/gendata2_inter_X4.pdf",width=11,height=8)#,units="in",res=150)
par(mfrow=c(1,2), tcl=-0.5, family="serif", omi=c(0.2,0.2,0,0))
par(mai=c(0.,0.,0.3,0.))
persp(x, e, z.truth,
      zlim = z_range,
      theta=30, phi=30,
      ltheta = 120, expand = 0.5,
      r=2, shade=0.3, axes=TRUE,scale=TRUE, box=T,
      nticks=5,
      # ticktype="detailed",
      col=trop[2],
      xlab=sprintf("f(x_%s)",as.numeric(gsub("X","",i))),
      ylab="X_E",
      zlab="Y", main="Truth")
persp(x, e, z.est,
      theta=30, phi=30,
      ltheta = 120, expand = 0.5,
      r=2, shade=0.3, axes=TRUE,scale=TRUE, box=T,
      nticks=5,
      zlim = z_range,
      # ticktype="detailed",
      col=trop[1],
      xlab=sprintf("f(x_%s)",as.numeric(gsub("X","",i))),
      ylab="X_E",
      zlab="Y", main="Estimated X_E*f(X_4)")
dev.off()


# plot  for gendata2 interaction X3--------------------------------------------------------


i = "X3"
x <- seq(range(DT$x[,i])[1], range(DT$x[,i])[2], length.out = 30)
e <- seq(range(DT$e)[1], range(DT$e)[2], length.out = 30)
f.est <- function(x, e) { e * splines::bs(x, df = 5) %*%
    as.matrix(coef(cvfit, s = "lambda.1se")[paste0(i,"_", seq_len(5),":X_E"),]) }
f.truth <- function(x, e) { e * DT$f3(x)  }
z.est <- outer(x, e, f.est)
z.truth <- outer(x, e, f.truth)

op <- par(bg = "white")
z_range <- c(min(z.est, z.truth), max(z.est, z.truth))
# z_range <- c(-1.868998, 2.118003)

dev.off()
# png(filename="gendata2_inter_X4.png",width=11,height=8,units="in",res=150)
pdf(file="~/Dropbox/jobs/hec/talk/gendata2_inter_X3.pdf",width=11,height=8)#,units="in",res=150)
par(mfrow=c(1,2), tcl=-0.5, family="serif", omi=c(0.2,0.2,0,0))
par(mai=c(0.,0.,0.3,0.))
persp(x, e, z.truth,
      zlim = z_range,
      theta=30, phi=30,
      ltheta = 120, expand = 0.5,
      r=2, shade=0.3, axes=TRUE,scale=TRUE, box=T,
      nticks=5,
      # ticktype="detailed",
      col=trop[2],
      xlab=sprintf("f(x_%s)",as.numeric(gsub("X","",i))),
      ylab="X_E",
      zlab="Y", main="Truth")
persp(x, e, z.est,
      theta=30, phi=30,
      ltheta = 120, expand = 0.5,
      r=2, shade=0.3, axes=TRUE,scale=TRUE, box=T,
      nticks=5,
      zlim = z_range,
      # ticktype="detailed",
      col=trop[1],
      xlab=sprintf("f(x_%s)",as.numeric(gsub("X","",i))),
      ylab="X_E",
      zlab="Y", main="Estimated X_E*f(X_3)")
dev.off()




