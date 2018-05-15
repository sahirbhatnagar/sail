rm(list=ls())
# devtools::document()
devtools::load_all()
# source("R/ksfilter.R")
pacman::p_load(data.table)
pacman::p_load(magrittr)
pacman::p_load(tidyverse)
pacman::p_load(doParallel)
pacman::p_load(splines)
pacman::p_load(foreach)
pacman::p_load(doMC)
pacman::p_load(glmnet)
pacman::p_load_current_gh('sahirbhatnagar/glinternet')
pacman::p_load(DataExplorer)
pacman::p_load(LassoBacktracking)
# pacman::p_load(mice)


# DT <- xlsx::read.xlsx("rda/For Sahir.xlsx", sheetIndex = 1)
DT <- readRDS("rda/scott/scott.rds") %>% as.data.table()
rm_vars <- paste0("std2_",c("com_100m","com_200m","com_300m",
             "gov_100m",
             "ind_100m",
             "highway_100m", "highway_200m", "highway_300m",
             "water_100m", "water_200m", "water_300m",
             "rail_100m","rail_200m"))
DT[, hist(ros_summer)]
DT[, hist(ros_annual)]
x_vars <- DT %>% colnames() %>% grep("std2", . , value = TRUE)
x_vars <- setdiff(x_vars, rm_vars)
X <- DT[, ..x_vars]
# str(X)
Xc <- X[complete.cases(X),] %>% as.matrix()

Y <- DT[complete.cases(X)]$ros_annual
# Y <- DT[complete.cases(X)]$ros_summer

help(glmnet)
set.seed(123456)
cvfitlasso <- cv.glmnet(Xc, Y, alpha = 1, nfolds = 10)
plot(cvfitlasso)
# coef(cvfitlasso, s = "lambda.min")
sprintf("%.2f (%.2f, %.2f)",
        cvfitlasso$cvm[which(cvfitlasso$lambda.min==cvfitlasso$lambda)],
        cvfitlasso$cvlo[which(cvfitlasso$lambda.min==cvfitlasso$lambda)],
        cvfitlasso$cvup[which(cvfitlasso$lambda.min==cvfitlasso$lambda)])

cor(predict(cvfitlasso, s="lambda.min", newx = Xc), Y)^2

set.seed(123456)
which(colnames(Xc)=="std2_pop_100m")
cvglinternet <- glinternet.cv(X = Xc, Y = Y,nFolds = 10,
                              numLevels = rep(1, ncol(Xc)),
                              family = "gaussian",
                              interactionCandidates = c(71),
                              verbose = T)
tc <- coef(cvglinternet, lambdaType = "lambdaHat")

(mains <- colnames(Xc)[tc$mainEffects$cont])
(inters <- paste0(colnames(Xc)[tc$interactions$contcont[,2]],":E"))
(cvmse <- cvglinternet$cvErr[which(cvglinternet$lambdaHat==cvglinternet$lambda)])
plot(cvglinternet)

sprintf("%.2f (%.2f, %.2f)",
cvglinternet$cvErr[which(cvglinternet$lambdaHat==cvglinternet$lambda)],
cvglinternet$cvErr[which(cvglinternet$lambdaHat==cvglinternet$lambda)]-cvglinternet$cvErrStd[which(cvglinternet$lambdaHat==cvglinternet$lambda)],
cvglinternet$cvErr[which(cvglinternet$lambdaHat==cvglinternet$lambda)]+cvglinternet$cvErrStd[which(cvglinternet$lambdaHat==cvglinternet$lambda)])
abline(h=310.96)


plot((predict(cvglinternet, X = Xc) - Y)^2)
abline(h=0)
plot(cvglinternet)
abline(h=78.255)
plot(predict(cvglinternet, X = Xc), Y)
abline(a=0,b=1)
cor(predict(cvglinternet, X = Xc), Y)^2





set.seed(123456)
fitBT <- cvLassoBT(x = Xc,
                   y = Y, iter_max=10, nperms=1, nfolds = 5, nlambda = 100)
coefBT <- as.matrix(predict(fitBT$BT_fit, type = "coef",
                            s = fitBT$lambda[fitBT$cv_opt[1]], iter = fitBT$cv_opt[2]))
nzcoef <- coefBT[sail:::nonzero(coefBT),,drop=F]
fitBT$cv_opt_err

cor(predict(fitBT$BT_fit, type = "resp",
            s = fitBT$lambda[fitBT$cv_opt[1]], iter = fitBT$cv_opt[2]), Y)^2
plot(predict(fitBT$BT_fit, type = "resp",
            s = fitBT$lambda[fitBT$cv_opt[1]], iter = fitBT$cv_opt[2]), Y)
abline(a=0,b=1)

# f.basis <- function(i) splines::bs(i, degree = 3)
f.basis <- function(i) i
E <- DT[complete.cases(X)]$std2_pop_100m
Xcc <- Xc[, -which(colnames(Xc)=="std2_pop_100m")]
registerDoMC(cores = 10)
set.seed(123456)
cvfit03 <- cv.sail(x = Xcc, y = Y, e = E, basis = f.basis, alpha = 0.5,
                   maxit = 1000, thresh = 5e-03,
                   # fdev = 1e-10,
                   grouped = FALSE,
                   strong = FALSE, verbose = 2, nfolds = 10, parallel = TRUE)
# saveRDS(cvfit03, file = "rda/scott/cvfit_alpha01_groupedFalse_strongFalse_nfolds5_bs3_Estd2_pop_100m.rds")
plot(cvfit03)
plot(cvfit03$sail.fit)
predict(cvfit03, type="non", s = "lambda.min")

sprintf("%.2f (%.2f, %.2f)",
        cvfit03$cvm[which(cvfit03$lambda.min==cvfit03$lambda)],
        cvfit03$cvlo[which(cvfit03$lambda.min==cvfit03$lambda)],
        cvfit03$cvup[which(cvfit03$lambda.min==cvfit03$lambda)])

cor(predict(cvfit03, s="lambda.min"), Y)^2
plot(predict(cvfit03, s="lambda.min"), Y)
mean((predict(cvfit03, s="lambda.min")- Y)^2)

abline(a=0,b=1)

tf <- readRDS("rda/scott/cvfit_alpha03_groupedFalse_strongFalse_nfolds5_bs3_Estd2_pop_100m.rds")
plot(tf)
tf$cvm[which(tf$lambda.min==tf$lambda)]
cor(predict(tf, s="lambda.min"), Y)^2
# E=population density ----------------------------------------------------

# E <- DT[complete.cases(X)]$std2_d_airport_yyz
# Xcc <- Xc[, -which(colnames(Xc)=="std2_d_airport_yyz")]
# E <- DT[complete.cases(X)]$std2_highway_300m
# Xcc <- Xc[, -which(colnames(Xc)=="std2_highway_300m")]
# E <- DT[complete.cases(X)]$std2_pop_100m
# Xcc <- Xc[, -which(colnames(Xc)=="std2_pop_100m")]
E <- DT[complete.cases(X)]$std2_majrd_100m
Xcc <- Xc[, -which(colnames(Xc)=="std2_majrd_100m")]
# E <- DT[complete.cases(X)]$std2_build_height_1000m
# Xcc <- Xc[, -which(colnames(Xc)=="std2_build_height_1000m")]
# pacman::p_load(DataExplorer)
# DataExplorer::HistogramContinuous(DT[, grep("build_heigh", colnames(DT), value = T), with=F])
# apply(as.matrix(DT[, grep("build_heigh", colnames(DT), value = T), with=F]), 2, sd, na.rm=TRUE)

f.basis <- function(i) splines::bs(i, degree = 3)
# f.basis <- function(i) i
fit <- sail(x = Xcc, y = Y, e = E, basis = f.basis, alpha = 0.5,
            strong = FALSE, verbose = 2, nlambda = 100)
plot(fit)
coef(fit)[nonzero(coef(fit)),]

cvfit02 <- cv.sail(x = Xcc, y = Y, e = E, basis = f.basis, alpha = 0.2,
                   strong = FALSE, verbose = 2, nfolds = length(Y), parallel = TRUE)
cvfit03 <- cv.sail(x = Xcc, y = Y, e = E, basis = f.basis, alpha = 0.3,
                   strong = FALSE, verbose = 2, nfolds = length(Y), parallel = TRUE)
cvfit04 <- cv.sail(x = Xcc, y = Y, e = E, basis = f.basis, alpha = 0.4,
                   strong = TRUE, verbose = 2, nfolds = length(Y), parallel = TRUE)

cvfit03 <- cv.sail(x = Xcc, y = Y, e = E, basis = f.basis, alpha = 0.3,
                   grouped = FALSE,
                   strong = FALSE, verbose = 2, nfolds = 5, parallel = TRUE)
plot(cvfit03)

cvfit
plot(cvfit)
predict(cvfit, type="non", s = "lambda.min")
cvfit$cvm[which(cvfit$lambda.min==cvfit$lambda)]
cvfit03$cvm[which(cvfit03$lambda.min==cvfit03$lambda)]

cvfit03 <- readRDS("cvfit_alpha03_groupedFalse_strongFalse_nfolds5_bs3_Estd2_pop_100m.rds")
plot(cvfit03)
predict(cvfit03, type="non", s = "lambda.min")
cvfit03$cvm[which(cvfit03$lambda.min==cvfit03$lambda)]
cor(predict(cvfit03, s="lambda.min"), Y)^2

pdf("maineff_distance_to_aiport.pdf", width = 11, height = 8)
plotMain(cvfit03$sail.fit, x = Xcc, xvar = "std2_d_airport_yyz", s = cvfit03$lambda.min,
         ylab = "f(Distance to Airport (YYZ) (m))", xlab = "Distance to Airport (YYZ) (m)")
dev.off()

pdf("maineff_traffic_counts.pdf", width = 11, height = 8)
plotMain(cvfit03$sail.fit, x = Xcc, xvar = "std2_traffic_500", s = cvfit03$lambda.min,
         ylab = "f(Traffic Counts (500 m))", xlab = "Traffic Counts (500 m)")
dev.off()

help(plotMain)

pdf("interaction_distance_density.pdf", width = 11, height = 11)
plotInter(cvfit03$sail.fit, x = Xcc, xvar = "std2_d_airport_yyz", s = cvfit03$lambda.min,
          xlab = "Distance to Airport (YYZ) (m)", title_z = "Interaction Effect",
          ylab = "Population Density (100 m)", zlab = "Annual ROS")
dev.off()
help(plotMain)
help(plotInter)




# E=building height -------------------------------------------------------

# DT <- readRDS("rda/scott.rds") %>% as.data.table()
DT <- readRDS("rda/scott_raw.rds") %>% as.data.table()
str(DT)
DT[, hist(ros_summer)]
DT[, hist(ros_annual)]
colnames(DT)
head(DT)
rm_vars <- c("com_100m","com_200m","com_300m",
             "gov_100m",
             "ind_100m",
             "highway_100m", "highway_200m", "highway_300m",
             "water_100m", "water_200m", "water_300m",
             "rail_100m","rail_200m")

x_vars <- setdiff(colnames(DT), c("sample_id", "labelid", "siteid", "ros_summer", "ros_annual",
                                  "bc_ngm3_w", "bc_ngm3_s", "bc_annual_mean", rm_vars))


X <- DT[, ..x_vars]
str(X)
Xc <- X[complete.cases(X),] %>% as.matrix()
Y <- DT[complete.cases(X)]$ros_annual

E <- standardize(DT[complete.cases(X)]$build_height_500m, center = TRUE, normalize = TRUE)$x
Xcc <- standardize(Xc[, -which(colnames(Xc)=="build_height_500m")], center = TRUE, normalize = TRUE)$x

f.basis <- function(i) splines::bs(i, degree = 3)
# f.basis <- function(i) i
# f.basis <- function(i) gamsel::basis.gen(i)
# pacman::p_load(genefilter)
# which(genefilter::rowSds(t(Xcc))==0)
# genefilter::rowSds(t(Xcc)) %>% plot
#
#
#
# standardize(Xcc, center = TRUE, normalize = TRUE)
fit <- sail(x = Xcc, y = Y, e = E, basis = f.basis, alpha = 0.1, center.x = FALSE, center.e = FALSE,
            strong = FALSE, verbose = 2, nlambda = 100, dfmax = 5)
plot(fit)
registerDoMC(cores = 5)

cvfit03 <- cv.sail(x = Xcc, y = Y, e = E, basis = f.basis, alpha = 0.3,
                   grouped = FALSE,
                   strong = FALSE, verbose = 2, nfolds = 5, parallel = TRUE)

plot(cvfit03)
predict(cvfit03, type="non", s = "lambda.min")
cvfit03$cvm[which(cvfit03$lambda.min==cvfit03$lambda)]
cor(predict(cvfit03, s="lambda.min"), Y)^2


DT$bc_ngm3_s
DT$std2_d_airport_ytz #(marker for city center), closer
DT$std2_d_airport_yyz
DT$std2_majrd_1000m
DT$std2_highway_300m
DT$std2_pop_100m
DT$std2_traffic_300

E = DT$std2_build_height_1000m #(average building height)
#(max height for max building height)
DT$std2_build_300m





dev.off()
DT$bc_ngm3_s %>% hist
