# for i in 1 ; do qsub -v index=$i mcgill-simulation.sh ; done

pacman::p_load(splines)
pacman::p_load(magrittr)
pacman::p_load(foreach)
pacman::p_load(methods)
pacman::p_load(doMC)
pacman::p_load(profvis)
# pacman::p_load(gamsel)

# rm(list=ls())
# dev.off()
# devtools::load_all("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/")
devtools::load_all("/home/sahir/git_repositories/sail/")


parameterIndex <- as.numeric(as.character(commandArgs(trailingOnly = T)[1]))
# parameterIndex = 1

if (parameterIndex == 1) { # 1a
  hierarchy = "strong" ; nonlinear = TRUE ; interactions = TRUE
} else if (parameterIndex == 2) { # 1b
  hierarchy = "weak" ; nonlinear = TRUE ; interactions = TRUE
} else if (parameterIndex == 3) { # 1c
  hierarchy = "none" ; nonlinear = TRUE ; interactions = TRUE
} else if (parameterIndex == 4) { # 2
  hierarchy = "strong"; nonlinear = FALSE; interactions = TRUE
} else if (parameterIndex == 5) { # 3
  hierarchy = "strong" ; nonlinear = TRUE ; interactions = FALSE
} else if (parameterIndex == 6) { # this is scenario 2 but we fit sail with degree=1
  hierarchy = "strong"; nonlinear = FALSE; interactions = TRUE
}

lambda.type <- "lambda.min"
# hierarchy = "strong", nonlinear = TRUE, interactions = TRUE, # scenario 1a
# hierarchy = "weak", nonlinear = TRUE, interactions = TRUE, # scenario 1b
# hierarchy = "none", nonlinear = TRUE, interactions = TRUE, # scenario 1c
# hierarchy = "strong", nonlinear = FALSE, interactions = TRUE, # scenario 2
# hierarchy = "strong", nonlinear = TRUE, interactions = FALSE, # scenario 3

# Simulate Data -----------------------------------------------------------

n = 200
p = 50

DT <- gendataPaper(n = n, p = p, SNR = 2, betaE = 1,
                   hierarchy = hierarchy, nonlinear = nonlinear, interactions = interactions,
                   corr = 0,
                   E = truncnorm::rtruncnorm(n, a = -1, b = 1))

DT$y <- scale(DT$y, center = TRUE, scale = FALSE)

registerDoMC(cores = 5)
# profvis(
  # system.time(
Rprof(tmp <- tempfile())
fit <- sail(x = DT$x, y = DT$y, e = DT$e, df = 5, degree = 3, basis.intercept = FALSE,
            thresh = 1e-3,
            maxit = 1000,
            alpha = .2,
            center.x = TRUE,
            verbose = T, nlambda = 100)
Rprof()
proftable(tmp)
plot(fit)
if (parameterIndex != 6) {
cvfit <- cv.sail(x = DT$x, y = DT$y, e = DT$e, df = 5, degree = 3, basis.intercept = FALSE,
                 thresh = 1e-4,
                 maxit = 1000,
                 alpha = .2,
                 parallel = TRUE,
                 center.x = TRUE,
                 # foldid = foldid,
                 nfolds = 10, verbose = T, nlambda = 100)
} else if (parameterIndex == 6){
  cvfit <- cv.sail(x = DT$x, y = DT$y, e = DT$e, degree = 1, basis.intercept = FALSE,
                   thresh = 1e-4,
                   maxit = 1000,
                   alpha = .05,
                   parallel = TRUE,
                   # foldid = foldid,
                   nfolds = 10, verbose = T, nlambda = 100)
  DT$scenario <- "6"
}
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

if (parameterIndex != 6) {
  saveRDS(object = cvfit,
          file = tempfile(pattern = sprintf("cvfit_gendata2_n200_p1000_SNR2_betaE1_df5_degree3_alpha05_%s_",DT$scenario),
                          tmpdir = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/mcgillsims/gendata2_p1000_1c_2_3_6",
                          fileext = ".rds")
  )
} else if (parameterIndex == 6) {
  saveRDS(object = cvfit,
          file = tempfile(pattern = sprintf("cvfit_gendata2_n200_p1000_SNR2_betaE1_df1_degree1_alpha05_%s_",DT$scenario),
                          tmpdir = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/mcgillsims/gendata2_p1000_1c_2_3_6",
                          fileext = ".rds")
  )
}
# files = list.files(path = '/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_lambda_branch/mcgillsims',
#                    pattern = '*.rds', full.names = TRUE)
# dat_list = lapply(files, function (x) readRDS(x))
# plot(dat_list[[2]])
