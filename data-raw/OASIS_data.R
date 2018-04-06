######################################
# R Source code file for creating dataset to be included in the sail package
# data taken from:
# https://github.com/stnava/RMI/tree/master/tomfletcher
# which is from the OASIS brain project
# http://www.oasis-brains.org/
# Code modified from https://github.com/stnava/RMI/blob/master/tomfletcher/ModelSelection.Rnw
# Author: Sahir Bhatnagar
# Created: April 5, 2018
# Updated:
#####################################

## Read in the longitudinal OASIS data

if (!require("pacman")) install.packages("pacman")
pacman::p_load(magrittr)
pacman::p_load(dplyr)

clinical <- read.csv("https://raw.githubusercontent.com/stnava/RMI/master/tomfletcher/oasis_longitudinal.csv")
hippo <- read.csv("https://raw.githubusercontent.com/stnava/RMI/master/tomfletcher/oasis_longitudinal_hippocampus.csv")
ldat <- merge(hippo, clinical, by.x = c("ID", "Visit"), by.y = c("Subject.ID", "Visit"))

## To simplify things, we'll remove subjects that converted to dementia
converts <- unique(ldat[ldat$Group == "Converted",]$ID)
ldat <- ldat[!(ldat$ID %in% converts),]
ldat$Group <- factor(ldat$Group)

visit1 <- ldat[ldat$Visit==1,,drop=F]
visit1 <- visit1[ , !(names(visit1) %in% c("Visit","MR.Delay","SES", "Hand"))]
visit1 <- visit1[complete.cases(visit1),,drop=FALSE]

visit2 <- ldat[ldat$Visit==2,,drop=F]
visit2 <- visit2[ , !(names(visit2) %in% c("Visit","MR.Delay", "SES","Hand"))]
visit2 <- visit2[complete.cases(visit2),,drop=FALSE]

X <- visit1 %>% dplyr::select(Age, EDUC,MMSE,eTIV,nWBV, ASF) %>% as.matrix()
noise <- replicate(24, rnorm(nrow(X)))
colnames(noise) <- paste0("noise", seq(ncol(noise)))
X <- cbind(X, noise)

X2 <- visit2 %>% dplyr::select(Age, EDUC,MMSE,eTIV,nWBV, ASF) %>% as.matrix()
noise <- replicate(24, rnorm(nrow(X2)))
colnames(noise) <- paste0("noise", seq(ncol(noise)))
X2 <- cbind(X2, noise)

Y <- visit1 %>% dplyr::pull(RightHippoVol)
Y2 <- visit2 %>% dplyr::pull(RightHippoVol)

E <- (visit1 %>% dplyr::pull(Group) %>% as.numeric())  - 1
E2 <- (visit2 %>% dplyr::pull(Group) %>% as.numeric())  - 1

oasis <- list(x = X, y = Y, e = E)
devtools::use_data(oasis, overwrite = TRUE)
# Function to plot raw longitudinal data (assumes two groups)
# long.plot <- function(data, yname, idname, agename, groupname, ylab = yname, main = "",
#                       pch = 19, cex.main = 1.5, cex.lab = 1.25, cex.axis = 1.25, alpha = 0.5){
#   cols = c(rgb(0,0,1,alpha), rgb(1,0,0,alpha))
#   y = data[,yname]
#   age = data[,agename]
#   yrange = c(min(y), max(y))
#   plot(y ~ age, ylab = ylab, xlab = "Age", main = main, ylim = yrange,
#        pch = pch, cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis,
#        col = cols[data[,groupname]])
#
#   ids = unique(data[,idname])
#
#   for(id in ids)
#   {
#     x = data[data[,idname] == id,]
#     lines(x[,yname] ~ x[,agename], col = cols[x[,groupname]])
#   }
# }

# long.plot(ldat, "RightHippoVol", "ID", "Age", "Group",
#           main = "OASIS Longitudinal Hippocampus Data")
# legend("topright", c("Nondemented", "Demented"),
#        col = c("red", "blue"), pch = 19)
# ldat$ID %>% unique() %>% length


# ## Read in the cross-sectional OASIS data, including hippocampal volumes
# x = read.csv("~/Downloads/RMI/tomfletcher/oasis_cross-sectional.csv")
# y = read.csv("~/Downloads/RMI/tomfletcher/oasis_hippocampus.csv")
# cdat = merge(x, y, by = "ID")
#
# ## Let's look at only elderly subjects
# # cdat = cdat[cdat$Age >= 60,]
# rownames(cdat) = NULL
#
# ## Remove columns we won't need
# cdat = cdat[ , !(names(cdat) %in% c("Delay", "Hand", "SES", "Educ", "ASF"))]
#
# cdat <- cdat[complete.cases(cdat),,drop=FALSE]

# plot_missing(visit1)
# plot_bar(visit1)
# plot_histogram(visit1)
# plot_histogram(visit2)
# skimr::skim(as.data.frame(visit1))
# skimr::skim(as.data.frame(visit2))
# dput(colnames(visit1))
#
# visit1
# Xm <- model.matrix(~0+bs(Age,3)+bs(EDUC,3)+bs(MMSE,3)+bs(eTIV,3)+bs(nWBV,3)+bs(ASF,3)+CDR+M.F,
#                    data = visit1)
# group <- attr(Xm, "assign")
# Xm <- Xm[,-which(colnames(Xm)=="M.FF")]
# group <- group[-length(group)]
# E <- (visit1 %>% dplyr::pull(Group) %>% as.numeric())  - 1
# head(Xm)
# skimr::skim(as.data.frame(X))
# rm(list=ls())
# ## Read in the cross-sectional OASIS data, including hippocampal volumes
# x = read.csv("data-raw/oasis_cross-sectional.csv")
# y = read.csv("data-raw/oasis_hippocampus.csv")
# cdat = merge(x, y, by = "ID")
#
# ## Let's look at only elderly subjects
# cdat = cdat[cdat$Age >= 60,]
# rownames(cdat) = NULL
# cdat$CDR[cdat$CDR>0] <- 1
# cdat$CDR %>% table
# ## Remove columns we won't need
# cdat <- cdat[ , !(names(cdat) %in% c("Delay", "Hand", "SES","ID"))]
# # head(cdat)
#
# cdat <- cdat %>%
#   dplyr::select(Age, Educ,MMSE,eTIV,nWBV, ASF, M.F, LeftHippoVol,RightHippoVol, CDR) %>%
#   mutate(M.F = as.numeric(M.F)-1)
# cdat <- cdat[complete.cases(cdat),,drop=FALSE]
#
# # head(cdat)
#
# X <- cdat %>%
#   dplyr::select(Age, Educ,MMSE,eTIV,nWBV, ASF, M.F) %>%
#   # dplyr::select(CDR, Educ,MMSE,eTIV,nWBV, ASF, M.F) %>%
#   # dplyr::select(Age, Educ,MMSE,eTIV,nWBV, ASF, CDR) %>%
#   as.matrix()
# # head(X)
# noise <- replicate(43, rnorm(nrow(X)))
# colnames(noise) <- paste0("noise", seq(ncol(noise)))
# X <- cbind(X, noise)
#
# E <- cdat %>% dplyr::pull(CDR)
#
# Y <- cdat %>% dplyr::pull(RightHippoVol)
# # Y <- cdat %>% dplyr::pull(LeftHippoVol)
#
# oasis2 <- list(x = X, y = Y, e = E)
# devtools::use_data(oasis2, overwrite = TRUE)



