pacman::p_load(splines)
pacman::p_load(magrittr)
pacman::p_load(foreach)
pacman::p_load(methods)
pacman::p_load(doMC)
# pacman::p_load(gamsel)

# rm(list=ls())
# dev.off()
devtools::load_all()


# Gendata -----------------------------------------------------------------

DT <- gendata2(n = 200, p = 1000, SNR = 2, betaE = 1, hierarchy = FALSE, corr = 0,
               E = truncnorm::rtruncnorm(200, a = -1, b = 1))
# DT <- gendata(n = 200, p = 25, SNR = 1, betaE = 2, df = 5, degree = 3)
# DT <- gendata3(n = 200, p = 50, betaE = 2, SNR = 3)
# DT <- gendata4(n = 100, p = 100, E = truncnorm::rtruncnorm(100, a = -1, b = 1), betaE = 2, SNR = 2)

DT$y <- scale(DT$y, center = TRUE, scale = FALSE)
# foldid <- sample(1:10,size=length(DT$y),replace=TRUE)

registerDoMC(cores = 10)

cvfit <- cv.sail(x = DT$x, y = DT$y, e = DT$e, df = 5, degree = 3, basis.intercept = FALSE,
                 thresh = 1e-4,
                 maxit = 1000,
                 alpha = .2,
                 parallel = TRUE,
                 # foldid = foldid,
                 nfolds = 10, verbose = T, nlambda = 100)
# plot(cvfit)
# plot(cvfit2)
# plot(cvfit$sail.fit)
# cvfit$sail.fit
# coef(cvfit, s = "lambda.min")[nonzero(coef(cvfit, s = "lambda.min")),,drop=F]
# coef(cvfit, s = "lambda.1se")[nonzero(coef(cvfit, s = "lambda.1se")),,drop=F]
# coef(cvfit2, s = "lambda.min")[nonzero(coef(cvfit2, s = "lambda.min")),,drop=F]
# coef(cvfit2, s = "lambda.1se")[nonzero(coef(cvfit2, s = "lambda.1se")),,drop=F]
# cvfit <- cvfit2
# par(mfrow=c(2,2))
# for (i in 1:4){
#   xv <- paste0("X",i)
#   ind <- cvfit$sail.fit$group == which(cvfit$sail.fit$vnames == xv)
#   design.mat <- cvfit$sail.fit$design[,cvfit$sail.fit$main.effect.names[ind],drop = FALSE]
#   # f.truth <- design.mat %*% DT$b1
#   f.truth <- DT[[paste0("f",i)]]
#   plotMain(object = cvfit$sail.fit, xvar = xv, s = cvfit$lambda.min, f.truth = f.truth, legend.position = "topleft")
# }


saveRDS(object = cvfit,
        file = tempfile(pattern = "cvfit_gendata2_n200_p1000_SNR2_betaE1_df5_degree3_alpha2_weak_hier_",
                        tmpdir = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/mcgillsims/gendata2_p1000_weak_hier",
                        fileext = ".rds")
)


# files = list.files(path = '/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/mcgillsims',
#                    pattern = '*.rds', full.names = TRUE)
# dat_list = lapply(files, function (x) readRDS(x))
# plot(dat_list[[2]])
