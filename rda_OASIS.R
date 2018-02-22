rm(list=ls())
# devtools::document()
devtools::load_all()
# source("R/ksfilter.R")
pacman::p_load(data.table)
pacman::p_load(magrittr)
pacman::p_load(genefilter)
pacman::p_load(tidyverse)
pacman::p_load(doParallel)
pacman::p_load(splines)
pacman::p_load(foreach)
pacman::p_load(doMC)
pacman::p_load(DataExplorer)


## Read in the cross-sectional OASIS data, including hippocampal volumes
x = read.csv("~/Downloads/RMI/tomfletcher/oasis_cross-sectional.csv")
y = read.csv("~/Downloads/RMI/tomfletcher/oasis_hippocampus.csv")
cdat = merge(x, y, by = "ID")

## Let's look at only elderly subjects
# cdat = cdat[cdat$Age >= 60,]
rownames(cdat) = NULL

## Remove columns we won't need
cdat = cdat[ , !(names(cdat) %in% c("Delay", "Hand", "SES", "Educ", "ASF"))]

cdat <- cdat[complete.cases(cdat),,drop=FALSE]

## Read in the longitudinal OASIS data
clinical = read.csv("~/Downloads/RMI/tomfletcher/oasis_longitudinal.csv")
hippo = read.csv("~/Downloads/RMI/tomfletcher/oasis_longitudinal_hippocampus.csv")
ldat = merge(hippo, clinical, by.x = c("ID", "Visit"), by.y = c("Subject.ID", "Visit"))

## To simplify things, we'll remove subjects that converted to dementia
converts = unique(ldat[ldat$Group == "Converted",]$ID)
ldat = ldat[!(ldat$ID %in% converts),]
ldat$Group = factor(ldat$Group)

# Function to plot raw longitudinal data (assumes two groups)
long.plot <- function(data, yname, idname, agename, groupname, ylab = yname, main = "",
           pch = 19, cex.main = 1.5, cex.lab = 1.25, cex.axis = 1.25, alpha = 0.5){
    cols = c(rgb(0,0,1,alpha), rgb(1,0,0,alpha))
    y = data[,yname]
    age = data[,agename]
    yrange = c(min(y), max(y))
    plot(y ~ age, ylab = ylab, xlab = "Age", main = main, ylim = yrange,
         pch = pch, cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis,
         col = cols[data[,groupname]])

    ids = unique(data[,idname])

    for(id in ids)
    {
      x = data[data[,idname] == id,]
      lines(x[,yname] ~ x[,agename], col = cols[x[,groupname]])
    }
}


long.plot(ldat[ldat$Visit==2,,drop=F], "RightHippoVol", "ID", "Age", "Group",
          main = "OASIS Longitudinal Hippocampus Data")
legend("topright", c("Nondemented", "Demented"),
       col = c("red", "blue"), pch = 19)
ldat$ID %>% unique() %>% length

visit1 <- ldat[ldat$Visit==1,,drop=F]

visit1 <- visit1[ , !(names(visit1) %in% c("Visit","MR.Delay", "SES","Hand"))]

plot_missing(visit1)
plot_bar(visit1)
plot_histogram(visit1)
visit1 <- visit1[complete.cases(visit1),,drop=FALSE]

X <- visit1 %>% dplyr::select(Age, EDUC,MMSE,CDR,eTIV,nWBV, ASF) %>% as.matrix()
Y <- visit1 %>% dplyr::pull(RightHippoVol)
Y <- drop(scale(Y, center = TRUE, scale = FALSE))
hist(Y)
E <- (visit1 %>% dplyr::pull(Group) %>% as.numeric())  - 1

system.time(
  fit <- sail(x = X, y = Y, e = E, df = 3, degree = 3, thresh = 1e-4,
              maxit = 1000,
              alpha = .3,
              # dfmax = 15,
              verbose = TRUE, nlambda = 100)
)


plot(fit)
fit

registerDoMC(cores = 8)
foldid=sample(1:5,size=length(Y),replace=TRUE)
foldid %>% table
sample(rep(seq(10), length = length(Y)))

# cv1=cv.glmnet(x,y,foldid=foldid,alpha=1)
# cv.5=cv.glmnet(x,y,foldid=foldid,alpha=.5)
# cv0=cv.glmnet(x,y,foldid=foldid,alpha=0)

system.time(
  cvfit.9 <- cv.sail(x = X, y = Y, e = E, df = 3, degree = 3, thresh = 1e-3,
                   maxit = 1000,
                   # foldid = foldid,
                   alpha = .1,
                   parallel = TRUE,
                   nfolds = 3,
                   # dfmax = 15,
                   verbose = TRUE, nlambda = 100)
)

plot(cvfit.9)
coef(cvfit, s = "lambda.min")[nonzero(coef(cvfit, s = "lambda.min")),,drop=F]
coef(cvfit, s = "lambda.1se")[nonzero(coef(cvfit, s = "lambda.1se")),,drop=F]
xv <- "MMSE"
ind <- cvfit$sail.fit$group == which(cvfit$sail.fit$vnames == xv)
design.mat <- cvfit$sail.fit$design[,cvfit$sail.fit$main.effect.names[ind],drop = FALSE]
# f.truth <- design.mat %*% DT$b4
# f.truth <- DT$f2
plotMain(object = cvfit$sail.fit, xvar = xv, s = cvfit$lambda.min, legend.position = "topleft")


