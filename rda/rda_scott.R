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
pacman::p_load(DataExplorer)
# pacman::p_load(mice)


# DT <- xlsx::read.xlsx("rda/For Sahir.xlsx", sheetIndex = 1)
DT <- readRDS("rda/scott.rds") %>% as.data.table()
DT[, hist(ros_summer)]
DT[, hist(ros_annual)]
x_vars <- DT %>% colnames() %>% grep("std2", . , value = TRUE)
X <- DT[, ..x_vars]
# str(X)
Xc <- X[complete.cases(X),] %>% as.matrix()
Y <- DT[complete.cases(X)]$bc_ngm3_s

cvfitlasso <- cv.glmnet(Xc, Y, alpha = 0.5, nfolds = 5)
plot(cvfitlasso)
coef(cvfitlasso, s = "lambda.min")
cvfitlasso$cvm[which(cvfitlasso$lambda.min==cvfitlasso$lambda)]


# E=population density ----------------------------------------------------



# E <- DT[complete.cases(X)]$std2_d_airport_yyz
# Xcc <- Xc[, -which(colnames(Xc)=="std2_d_airport_yyz")]
# E <- DT[complete.cases(X)]$std2_highway_300m
# Xcc <- Xc[, -which(colnames(Xc)=="std2_highway_300m")]
# E <- DT[complete.cases(X)]$std2_pop_100m
# Xcc <- Xc[, -which(colnames(Xc)=="std2_pop_100m")]
E <- DT[complete.cases(X)]$std2_highway_100m
Xcc <- Xc[, -which(colnames(Xc)=="std2_highway_100m")]
# pacman::p_load(DataExplorer)
# DataExplorer::HistogramContinuous(DT[, grep("build_heigh", colnames(DT), value = T), with=F])
# apply(as.matrix(DT[, grep("build_heigh", colnames(DT), value = T), with=F]), 2, sd, na.rm=TRUE)

f.basis <- function(i) splines::bs(i, degree = 3)
# f.basis <- function(i) i
fit <- sail(x = Xcc, y = Y, e = E, basis = f.basis, alpha = 0.5,
            strong = TRUE, verbose = 2, nlambda = 100)
plot(fit)
coef(fit)[nonzero(coef(fit)),]
registerDoMC(cores = 5)
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
