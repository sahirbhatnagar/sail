pacman::p_load(data.table)
pacman::p_load(magrittr)
rm(list=ls())
load("data/adni/ADNI_pheno.RData")
DT_info <- fread("data/adni/fdg_info.csv")
DT <- fread("data/adni/fdg_baseline.csv")

amy_mat %>% as.matrix()


dev.off()
plot(density(amy_mat[,1]), type = "n", ylim = c(0, 10), xlim = c(0.3, 2.5))
n = dim(amy_mat)[2]-1
for(i in 1:n){
  lines(density(c(amy_mat[,i])))
}

X <- as.matrix(amy_mat)
e <- apoe
Y <- stage

fit <- funshim(x = X, y = Y, e = e, df = 3, maxit = 200, nlambda.gamma = 1, nlambda.beta = 50,
               nlambda = 50,
               thresh = 1e-5, center=TRUE, normalize=FALSE, verbose = T)
coef(fit)
