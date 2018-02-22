pacman::p_load(splines)
pacman::p_load(magrittr)
pacman::p_load(foreach)
pacman::p_load(methods)
pacman::p_load(doMC)
pacman::p_load(gamsel)

rm(list=ls())
# dev.off()
devtools::load_all()


# Gendata -----------------------------------------------------------------

DT <- gendata2(n = 200, p = 50, SNR = 3, betaE = 2)
# DT <- gendata(n = 200, p = 25, SNR = 1, betaE = 2, df = 5, degree = 3)
# DT <- gendata3(n = 200, p = 50, betaE = 2, SNR = 3)
# DT <- gendata4(n = 100, p = 100, E = truncnorm::rtruncnorm(100, a = -1, b = 1), betaE = 2, SNR = 2)

# DT$y <- scale(DT$y, center = TRUE, scale = FALSE)
# foldid <- sample(1:10,size=length(DT$y),replace=TRUE)

registerDoMC(cores = 10)

cvfit <- cv.sail(x = DT$x, y = DT$y, e = DT$e, df = 5, degree = 3, thresh = 1e-4, maxit = 1000,
                 alpha = .2,
                 parallel = TRUE,
                 # foldid = foldid,
                 nfolds = 10, verbose = T, nlambda = 100)
plot(cvfit)
plot(cvfit$sail.fit)
cvfit$sail.fit
coef(cvfit, s = "lambda.min")[nonzero(coef(cvfit, s = "lambda.min")),,drop=F]
coef(cvfit, s = "lambda.1se")[nonzero(coef(cvfit, s = "lambda.1se")),,drop=F]
saveRDS(object = cvfit,
        file = tempfile(pattern = "cvfit_gendata_n200_p50_SNR3_betaE2_df5_degree3_alpha2_",
                        tmpdir = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/mcgillsims/gendata",
                        fileext = ".rds")
)


# files = list.files(path = '/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/mcgillsims',
#                    pattern = '*.rds', full.names = TRUE)
# dat_list = lapply(files, function (x) readRDS(x))
# plot(dat_list[[2]])
