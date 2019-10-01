## ---- load-support ----
library(psych)
load("/Users/tianyuan.lu/Desktop/sail/real_data/methods_comparison_error_summary.RData")
par(family="serif")

## ---- support-error-cross ----
error.crosses(errorSummary[c(6:10),],
              errorSummary[c(1:5),],
              labels=unique(errorSummary$group1),
              xlab="Number of Active Variables",
              main = "SUPPORT Data: Means (+/- 1 SD) from 200 Train/Validate/Test Splits",
              sd = TRUE,
              cex.lab = 1.4,
              cex.axis = 1.4,
              cex.main = 1.5,
              offset = -0.8,
              xlim = c(0, 120),
              ylab="Test Set AUROC",
              colors = sail:::cbbPalette[c(1,3,4,7,2)],
              pch=16,cex=1)

## ---- dzclass-interaction ----
load("/Users/tianyuan.lu/Desktop/sail/real_data/sail200Bootstrap.RData")
par(mfrow=c(1,2), tcl=-0.5, family="serif",
    cex.lab = 1.6, font.lab = 1.6, cex.axis = 1.6,
    cex.main = .1,
    omi=c(0.2,0.2,0,0),
    mar = c(4, 4, 1, 1.1) + 0.1)
object = sailfittrain
x = usedatcomplete
varname = "num.co"
xvar = paste0("bs(",varname,", degree = 3)", 1:3)
design = designmat2
s = optLambda
e = usedatcomplete$dzclass

ind <- which(object$vnames %in% xvar)
ctPred <- list()
ARFPred <- list()
for (ite in 1:200) {
  allCoefs <- COEFlist[[ite]]
  a0 <- allCoefs[1, ]
  betas <- as.matrix(allCoefs[object$main.effect.names[ind], , drop = FALSE])
  alphas <- as.matrix(allCoefs[object$interaction.names[ind], , drop = FALSE])
  betaE <- as.matrix(allCoefs["E", , drop = FALSE])
  design.mat.main <- design[, object$main.effect.names[ind], drop = FALSE]
  design.mat.int <- design[, object$main.effect.names[ind], drop = FALSE] * e
  originalE <- e
  originalX <- x[,varname]
  f.hat <- drop(originalE * as.vector(betaE) + design.mat.main %*% betas + design.mat.int %*% alphas)
  ylims <- range(f.hat)
  control = drop(unique(originalE))[1]
  ARF = drop(unique(originalE))[2]
  cont_index <- which(originalE==control)
  ARF_index <- which(originalE==ARF)
  cont_pred <- f.hat[cont_index]
  ARF_pred <- f.hat[ARF_index]
  cont_pred[order(originalX[cont_index])] -> ctPred[[ite]]
  ARF_pred[order(originalX[ARF_index])] -> ARFPred[[ite]]
}
matrix(unlist(ctPred),byrow=T,nrow=200) -> ctPredmat
matrix(unlist(ARFPred),byrow=T,nrow=200) -> ARFPredmat
apply(ctPredmat,2,function(x) quantile(x,0.95)) -> ctPredmax
apply(ctPredmat,2,function(x) quantile(x,0.05)) -> ctPredmin
apply(ARFPredmat,2,function(x) quantile(x,0.95)) -> ARFPredmax
apply(ARFPredmat,2,function(x) quantile(x,0.05)) -> ARFPredmin
originalX[cont_index][order(originalX[cont_index])] -> polygonXct
originalX[ARF_index][order(originalX[ARF_index])] -> polygonXarf

allCoefs <- coef(object, s = s)
a0 <- allCoefs[1, ]
betas <- as.matrix(allCoefs[object$main.effect.names[ind], , drop = FALSE])
alphas <- as.matrix(allCoefs[object$interaction.names[ind], , drop = FALSE])
betaE <- as.matrix(allCoefs["E", , drop = FALSE])
design.mat.main <- design[, object$main.effect.names[ind], drop = FALSE]
design.mat.int <- design[, object$main.effect.names[ind], drop = FALSE] * e
originalE <- e # this is the centered E

originalX <- x[,varname]

f.hat <- drop(originalE * as.vector(betaE) + design.mat.main %*% betas + design.mat.int %*% alphas)
ylims <- range(f.hat)
control = drop(unique(originalE))[1]
ARF = drop(unique(originalE))[2]
cont_index <- which(originalE==control)
ARF_index <- which(originalE==ARF)
cont_pred <- f.hat[cont_index]
ARF_pred <- f.hat[ARF_index]
min.length.top <- range(f.hat)[1] ; max.length.top <- range(f.hat)[2]
cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
color = cbbPalette[c(6,4,7)]
plot(originalX, f.hat,
     pch = 19,
     col = color[1],
     bty="n",
     xlim=c(0,10),
     ylim=c(-0.4,3),
     xaxt="n",
     type = "n",
     main = "",
     xlab = "Number of comorbidities",
     ylab = "Marginal risk")
axis(1, labels = T)
ctPredmat <- ctPredmat[apply(ctPredmat,1,function(x) sum(x > ctPredmax)==0),]
ctPredmat <- ctPredmat[apply(ctPredmat,1,function(x) sum(x < ctPredmin)==0),]

ARFPredmat <- ARFPredmat[apply(ARFPredmat,1,function(x) sum(x > ARFPredmax)==0),]
ARFPredmat <- ARFPredmat[apply(ARFPredmat,1,function(x) sum(x < ARFPredmin)==0),]

rgb(red=0,green=114,blue=178,max=255,alpha=25) -> Col1
rgb(red=0,green=158,blue=115,max=255,alpha=25) -> Col2

apply(ctPredmat,1,function(x) lines(polygonXct,-x,col=Col1,lty=1))
apply(ARFPredmat,1,function(x) lines(polygonXarf,-x,col=Col2,lty=1))

lines(originalX[cont_index][order(originalX[cont_index])], -cont_pred[order(originalX[cont_index])], col = color[1], lwd = 3)
lines(originalX[ARF_index][order(originalX[ARF_index])], -ARF_pred[order(originalX[ARF_index])], col = color[2], lwd = 3)
graphics::rug(originalX, side = 1)

object = sailfittrain
x = usedatcomplete
varname = "age"
xvar = paste0("bs(",varname,", degree = 3)", 1:3)
design = designmat2
s = optLambda
e = usedatcomplete$dzclass
ind <- which(object$vnames %in% xvar)
ctPred <- list()
ARFPred <- list()
for (ite in 1:200) {
  allCoefs <- COEFlist[[ite]]
  a0 <- allCoefs[1, ]
  betas <- as.matrix(allCoefs[object$main.effect.names[ind], , drop = FALSE])
  alphas <- as.matrix(allCoefs[object$interaction.names[ind], , drop = FALSE])
  betaE <- as.matrix(allCoefs["E", , drop = FALSE])
  design.mat.main <- design[, object$main.effect.names[ind], drop = FALSE]
  design.mat.int <- design[, object$main.effect.names[ind], drop = FALSE] * e
  originalE <- e
  originalX <- x[,varname]
  f.hat <- drop(originalE * as.vector(betaE) + design.mat.main %*% betas + design.mat.int %*% alphas)
  ylims <- range(f.hat)
  control = drop(unique(originalE))[1]
  ARF = drop(unique(originalE))[2]
  cont_index <- which(originalE==control)
  ARF_index <- which(originalE==ARF)
  cont_pred <- f.hat[cont_index]
  ARF_pred <- f.hat[ARF_index]
  cont_pred[order(originalX[cont_index])] -> ctPred[[ite]]
  ARF_pred[order(originalX[ARF_index])] -> ARFPred[[ite]]
}
matrix(unlist(ctPred),byrow=T,nrow=200) -> ctPredmat
matrix(unlist(ARFPred),byrow=T,nrow=200) -> ARFPredmat
apply(ctPredmat,2,function(x) quantile(x,0.95)) -> ctPredmax
apply(ctPredmat,2,function(x) quantile(x,0.05)) -> ctPredmin
apply(ARFPredmat,2,function(x) quantile(x,0.95)) -> ARFPredmax
apply(ARFPredmat,2,function(x) quantile(x,0.05)) -> ARFPredmin
originalX[cont_index][order(originalX[cont_index])] -> polygonXct
originalX[ARF_index][order(originalX[ARF_index])] -> polygonXarf

allCoefs <- coef(object, s = s)
a0 <- allCoefs[1, ]
betas <- as.matrix(allCoefs[object$main.effect.names[ind], , drop = FALSE])
alphas <- as.matrix(allCoefs[object$interaction.names[ind], , drop = FALSE])
betaE <- as.matrix(allCoefs["E", , drop = FALSE])
design.mat.main <- design[, object$main.effect.names[ind], drop = FALSE]
design.mat.int <- design[, object$main.effect.names[ind], drop = FALSE] * e
originalE <- e # this is the centered E

originalX <- x[,varname]

f.hat <- drop(originalE * as.vector(betaE) + design.mat.main %*% betas + design.mat.int %*% alphas)
ylims <- range(f.hat)
control = drop(unique(originalE))[1]
ARF = drop(unique(originalE))[2]
cont_index <- which(originalE==control)
ARF_index <- which(originalE==ARF)
cont_pred <- f.hat[cont_index]
ARF_pred <- f.hat[ARF_index]
min.length.top <- range(f.hat)[1] ; max.length.top <- range(f.hat)[2]
cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
color = cbbPalette[c(6,4,7)]
plot(originalX, f.hat,
     pch = 19,
     col = color[1],
     bty="n",
     ylim=c(-0.1,1.3),
     xaxt="n",
     type = "n",
     main = "",
     xlab = "Age",
     ylab = "")
axis(1, labels = T)
ctPredmat <- ctPredmat[apply(ctPredmat,1,function(x) sum(x > ctPredmax)==0),]
ctPredmat <- ctPredmat[apply(ctPredmat,1,function(x) sum(x < ctPredmin)==0),]

ARFPredmat <- ARFPredmat[apply(ARFPredmat,1,function(x) sum(x > ARFPredmax)==0),]
ARFPredmat <- ARFPredmat[apply(ARFPredmat,1,function(x) sum(x < ARFPredmin)==0),]

rgb(red=0,green=114,blue=178,max=255,alpha=25) -> Col1
rgb(red=0,green=158,blue=115,max=255,alpha=25) -> Col2

apply(ctPredmat,1,function(x) lines(polygonXct,-x,col=Col1,lty=1))
apply(ARFPredmat,1,function(x) lines(polygonXarf,-x,col=Col2,lty=1))

lines(originalX[cont_index][order(originalX[cont_index])], -cont_pred[order(originalX[cont_index])], col = color[1], lwd = 3)
lines(originalX[ARF_index][order(originalX[ARF_index])], -ARF_pred[order(originalX[ARF_index])], col = color[2], lwd = 3)
legend("topleft", c("non-ARF/MOSF", "ARF/MOSF"),
       col = color[1:3], pch = 19, cex = 1.6, bty = "n")
graphics::rug(originalX, side = 1)

