## ---- load-ADNI-data ----

source("/home/sahir/git_repositories/sail/my_sims/plot_functions_rda_ADNI.R")

amy_mat <- read.csv("/home/sahir/git_repositories/sail/data-nogit/adni_new/csf_amyloid_final.csv",
                    stringsAsFactors = FALSE)
covr <- read.csv("/home/sahir/git_repositories/sail/data-nogit/adni_new/covariates.csv",
                 stringsAsFactors = FALSE, sep = ";")

# these are used for plotting. we use the entire data set for the plots
DT <- dplyr::inner_join(amy_mat, covr, by = c("PTID" = "IID")) %>%
  select(-AV45_path_bl)
X <- DT %>% select(starts_with("X"), diag_3bl.x, APOE_bin) %>%
  mutate(diag_3bl.x = diag_3bl.x - 1) %>%
  as.matrix()
Xnorm <- sail:::standardize(X, center = TRUE, normalize = TRUE)$x
E <- Xnorm[, "diag_3bl.x"]


# load-results-adni
df <- readRDS("/home/sahir/git_repositories/sail/my_sims/rda_results/rda_ADNI_may_17_2018v2.rds")
DT_res <- df %>% as.data.table()


## ---- get-best-sail-model ----

draw_ind <- DT_res[Method=="sail"][which.min(mse)]$Draw
dat <- draws(sim)@draws[[draw_ind]]
fit <- sail(x = dat$xtrain, y = dat$ytrain, e = dat$etrain,
            expand = FALSE,
            center.x = F,
            center.e = T,
            group = dat$group,
            alpha = 0.1,
            maxit = 250,
            strong = TRUE,
            verbose = 1)
ytest_hat <- predict(fit, newx = dat$xtest, newe = dat$etest)
msetest <- colMeans((dat$ytest - ytest_hat)^2)
lambda.min.index <- as.numeric(which.min(msetest))
lambda.min <- fit$lambda[which.min(msetest)]
yvalid_hat <- predict(fit, newx = dat$xvalid, newe = dat$evalid, s = lambda.min)
msevalid <- mean((dat$yvalid - drop(yvalid_hat))^2)
nzcoef <- predict(fit, s = lambda.min, type = "nonzero")

design <- do.call(rbind, list(dat$xtrain, dat$xtest, dat$xvalid))



## ---- plot-main-adni ----

# names taken from email from lai (search for ADNI in gmail)
plotMainADNI(fit, x = DT[as.numeric(as.character(rownames(design))), "X175"],
             xvar = paste0("bs(X175, 3)",1:3),design = design, s = lambda.min,
             ylab = "f(X175)", xlab = "Cuneus right")
plotMainADNI(fit, x = DT[as.numeric(as.character(rownames(design))), "X196"],
             xvar = paste0("bs(X196, 3)",1:3),design = design, s = lambda.min,
             ylab = "f(X196)", xlab = "Lateral occipitotemporal gyrus left")
plotMainADNI(fit, x = DT[as.numeric(as.character(rownames(design))), "X154"],
             xvar = paste0("bs(X154, 3)",1:3),design = design, s = lambda.min,
             ylab = "f(X154)", xlab =  "Middle occipital gyrus left")
plotMainADNI(fit, x = DT[as.numeric(as.character(rownames(design))), "X24"],
             xvar = paste0("bs(X24, 3)",1:3),design = design, s = lambda.min,
             ylab = "f(X24)", xlab =  "Middle occipital gyrus left")

## ---- plot-inter-adni-175 ----

par(mfrow=c(1,2), tcl=-0.5, family="serif",
    omi=c(0.2,0.2,0,0))
plotInterADNI(fit, x = DT[as.numeric(as.character(rownames(design))), "X175"],
              xvar = paste0("bs(X175, 3)",1:3), design = design, s = lambda.min,
              e = E, apoe = FALSE, legend = T, legend.position = "bottomleft",
              ylab = "f(X175)", xlab = "Cuneus right", main = "APOE e4 = 0")
plotInterADNI(fit, x = DT[as.numeric(as.character(rownames(design))), "X175"],
              xvar = paste0("bs(X175, 3)",1:3), design = design, s = lambda.min,
              e = E, apoe = TRUE, legend = F,
              ylab = "f(X175)", xlab = "Cuneus right", main = "APOE e4 = 1")


## ---- plot-inter-adni-154 ----

par(mfrow=c(1,2), tcl=-0.5, family="serif",
    omi=c(0.2,0.2,0,0))
plotInterADNI(fit, x = DT[as.numeric(as.character(rownames(design))), "X154"],
              xvar = paste0("bs(X154, 3)",1:3), design = design, s = lambda.min,
              e = E, apoe = FALSE, legend = F,
              ylab = "f(X154)", xlab = "Middle occipital gyrus left", main = "APOE e4 = 0")
plotInterADNI(fit, x = DT[as.numeric(as.character(rownames(design))), "X154"],
              xvar = paste0("bs(X154, 3)",1:3), design = design, s = lambda.min,
              e = E, apoe = TRUE, legend = T,
              ylab = "f(X154)", xlab = "Middle occipital gyrus left", main = "APOE e4 = 1")


## ----error-crosses ----

affect.mat2 <- describeBy(DT_res[, c("mse","nactive")], DT_res$Method, mat = TRUE)

par(family="serif")
error.crosses(affect.mat2[c(6:10),],
              affect.mat2[c(1:5),],
              labels=unique(affect.mat2$group1),
              xlab="Number of Active Variables",
              main = "ADNI Data: Means (+/- 1 SD) from 200 Train/Validate/Test Splits",
              sd = TRUE,
              cex.lab = 1.4,
              cex.axis = 1.4,
              cex.main = 1.5,
              xlim = c(0, 34),
              ylab="Test Set MSE",
              colors = sail:::cbbPalette[c(1,3,4,7,2)],
              pch=16,cex=2)

# error.crosses(affect.mat2[c(1:5),],
#               affect.mat2[c(6:10),],
#               labels=unique(affect.mat2$group1),
#               ylab="Number of active variables",
#               sd = TRUE,
#               # xlim = c(0, 35),
#               xlab="Test set MSE",
#               colors = sail:::cbbPalette[c(1,3,4,7,2)],
#               pch=16,cex=2)
