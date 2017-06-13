rm(list = ls())
devtools::load_code()
source("R/ksfilter.R")
pacman::p_load(data.table)
pacman::p_load(magrittr)
pacman::p_load(genefilter)
pacman::p_load(tidyverse)
load("data/adni/ADNI_pheno.RData") # old one 340 subjects




X <- as.matrix(amy_mat)
E <- covr_endo$Age_bl
Y <- stage

fmla <- as.formula(paste0("~0+",paste(colnames(X), collapse = "+"), "+ E +",paste(colnames(X),":E", collapse = "+") ))

str(X)

Xmat <- model.matrix(fmla, data = data.frame(X,E))

Ftest <- colFtests(X, factor(Y))
Ftest$name <- rownames(Ftest)

Ftest %>%
  dplyr::arrange(p.value) %>%
  filter(p.value < 1e-4) %>%
  select(name) %>% unlist

coef(cvfit, s="lambda.min")

nonz_coef <- coef(cvfit, s="lambda.min")[nonzero(coef(cvfit, s="lambda.min")),,drop=F]

plot()


rbind(Ftest$name, paste0("X",seq_len(ncol(X))))


top <- Ftest %>%
  dplyr::arrange(p.value) %>%
  filter(p.value < 1e-4) %>%
  select(name) %>% unlist(use.names = FALSE)

dev.off()
plot(density(amy_mat[,1]), type = "n", ylim = c(0, 10), xlim = c(0.3, 2.5))
n = dim(amy_mat)[2]-1
for(i in 1:n){
  lines(density(c(amy_mat[,i])))
}



fit <- funshim(x = X, y = Y, e = E, df = 3, maxit = 200, nlambda.gamma = 5, nlambda.beta = 5,
               nlambda = 100,lambda.factor = 0.0001,
               thresh = 1e-3, center=TRUE, normalize=F, verbose = T)

coef(fit)[nonzero(coef(fit)),]
pacman::p_load(doParallel)
doParallel::registerDoParallel(cores = 5)
getDoParWorkers()
devtools::load_all()
source("R/plot.R")
source("R/fitting.R")

cvfit <- cv.funshim(x = X,
                    y = Y,
                    e = E,
                    df = 2,
                    maxit = 200,
                    cores = 5,
                    nfolds = 5,
                    nlambda.gamma = 7,
                    nlambda.beta = 7,
                    nlambda = 49,
                    # lambda.factor = 0.0001,
                    thresh = 1e-3, center=TRUE, normalize=TRUE, verbose = T)
plot(cvfit)
coef(cvfit, s = "lambda.min")[nonzero(coef(cvfit, s = "lambda.min")),]





pacman::p_load(gamsel)
bases=pseudo.bases(X,degree=10,df=6)
# Gaussian gam
gamsel.out=gamsel(X,Y,bases=bases)
par(mfrow=c(1,2),mar=c(5,4,3,1))
summary(gamsel.out)
gamsel.cv=cv.gamsel(X,Y,bases=bases)
plot(gamsel.cv)
coef(gamsel.cv)
gamsel::getActive(gamsel.cv)

par(mfrow=c(2,2))
plot(gamsel.out, newx=X[,getActive(gamsel.out, index = 24, type = "nonlinear")$l24],index=24)
preds<-predict(gamsel.out, X, index = 24, type = "terms")

preds[1:5,1,1:5]


ds[order(ds$obj.k.rank)[1:20],]


trop <- RSkittleBrewer::RSkittleBrewer("tropical")
par(mai=c(0.2,0.4,0.2,0))
plot(X[,69], splines::bs(X[,69],3) %*%
       as.matrix(coef(cvfit, s = "lambda.min")[paste("X1",1:3, sep = "_"),,drop=F]),
     pch = 19,
     ylab = sprintf("f(%s)","X69"),
     xlab = "X69",
     xaxt="n", bty="n",
     col = trop[1]#,
     # main = latex2exp::TeX("$f(x_1) = 5x_1$") #,
     # ylim = c(min.length.top,max.length.top+0.15)
     # ylim = range(DT$f1(DT$x[,i]),
     # as.vector(fit$design[,paste(i,1:5, sep = "_")] %*%
     # coef(cvfit, s = "lambda.min")[paste(i,1:5, sep = "_"),,drop=F]))
)
points(DT$x[,i], DT$f1(DT$x[,i]), pch = 19, col=trop[2])
axis(1, labels=F)
legend("bottomright", c("Truth","Estimated"), col = trop[2:1], pch = 19, cex = 2, bty = "n")



# ks filter adni amy_beta


n<-200
p<-5000

x<-matrix(runif(n*p),n,p)
#x<-x%*%sigma.sqrt

y <- 4*x[,1]+2*tan(x[,2]*pi/2)+5*x[,3]^2+rnorm(n)


obj <- fused.k.filter(x=X,y=Y, nslices = 50)
obj <- k.filter.single(x=X*E,y=factor(Y))

max(obj$k.rank)

obj <- k.filter(x,y)

# the minimum number of predictors needed to keep all useful predictors
# in this case, all useful predictor are x1,x2 and x3
max(obj$k.rank[c(1:3)])
max(obj$k.rank[c(1:3)])

obj$k.rank

X[,obj$k.rank]
ds <- data.frame(colnames(X), obj$k.rank, stringsAsFactors = F)

ds[order(ds$obj.k.rank)[1:20],]$colnames.X.

# each row probably represents the number of slices used to categorize the
# response. each column of k.stat.single is a predictor
dim(obj$k.stat.single)
obj$k.stat.single[,1:10]


obj$k.stat



obj2 <- k.filter(x=X,y=Y)

obj2$k.rank





# fdg ---------------------------------------------------------------------

# global suvr is taken at baseline
# pheno_mat has Tau and Ptau, global Abeta, gloabl SUVR
# DT_info has APOE and diag at baseline
rm(list=ls())
load("~/Downloads/pheno_ADNI.RData")
DT_info <- read.csv("data/adni/fdg_info.csv", stringsAsFactors = FALSE)

str(DT_info)
# rownames(DT_info) <- DT_info$PTID
dim(DT_info)
colnames(DT_info)

# this is the fdg data
temp <- read.csv("data/adni/fdg_baseline.csv", stringsAsFactors = FALSE)
str(temp)
temp[1:5,1:6]
# rownames(temp) <- temp$PTID
dim(temp)
colnames(temp)

pheno_data <- data.frame(PTID = tpid.intr, pheno_mat, Diag_mat, cov_mat, stringsAsFactors = F)
str(pheno_data)
dim(pheno_data)
colnames(pheno_data)

library(tidyverse)
DT <- pheno_data %>% left_join(temp, by = "PTID") %>% left_join(DT_info, by = "PTID")
DT <- DT[complete.cases(DT),]



table(DT$Diag2)
table(DT$Diag1) # this has a 7 in it
table(DT$diag_3bl.y)

Y <- DT$diag_3bl.y
# E <- DT$APOE_bin
E <- DT$Ptau
X <- as.matrix(DT[,c(grep("X", colnames(DT), value = TRUE))])
colnames(X) <- paste0("X",1:ncol(X))
# X <- as.matrix(DT[,top])
# X <- as.matrix(DT[,SUVR])
devtools::load_all()

rm(list=ls())

X <- t(readRDS("~/Downloads/top.rds"))
cov <- readRDS("~/Downloads/cov.rds")
Y <- cov[,"X.12_asthma"]
E <- cov[,"current_smoker"]


fit <- funshim(x = X,
               y = Y,
               e = E,
               df = 3,
               maxit = 200,
               nlambda.gamma = 7,
               nlambda.beta = 7,
               nlambda = 49,
               lambda.factor = 0.00001,
               thresh = 1e-3, center=TRUE, normalize=TRUE, verbose = T)

nonzero(coef(fit))

cvfit <- cv.funshim(x = X,
                    y = Y,
                    e = E,
                    df = 3,
                    maxit = 200,
                    cores = 5,
                    nfolds = 5,
                    nlambda.gamma = 7,
                    nlambda.beta = 7,
                    nlambda = 49,
                    lambda.factor = 0.00001,
                    thresh = 1e-3, center=TRUE, normalize=TRUE, verbose = T)

dev.off()
coef(cvfit, s = "lambda.min")
source("R/plot.R")
plot(cvfit)
coef(fit)[nonzero(coef(fit)),]
fit$dfalpha

pacman::p_load(genefilter)

str(X)
colttests(X,Y)
colFtests(X,factor(Y))

Ftest <- colFtests(X, factor(Y))
Ftest$name <- rownames(Ftest)

top <- Ftest %>%
dplyr::arrange(p.value) %>%
  filter(p.value < 1e-4) %>%
  select(name) %>% unlist(use.names = FALSE)








# catherine laprise
asthma ~ methylation*smoking



