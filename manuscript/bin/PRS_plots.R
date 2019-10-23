#' This is a script which is called by the manuscript.Rnw file It runs cv.sail
#' on the actual PRS dataset (imputed once to fill in missing IQs), and then
#' plots the selected model's interaction curve

rm(list=ls())
## ---- PRS-packages --------------------------------------------------------
pacman::p_load(psych)
pacman::p_load(magrittr)
pacman::p_load(data.table)
pacman::p_load(sail)
pacman::p_load(mice)
pacman::p_load(doMC)
pacman::p_load(here)
source(here::here("manuscript/bin/PRS_plot_functions.R"))


## ---- load-PRS-results ------------------------------------------------

# see the file 'run_sail_200_bootstrap_PRS.R' for script used to run the models on
# bootstrap samples
load(file = here::here("manuscript/results/PRS_results.RData"))
affect.mat2 <- describeBy(simulation_results_PRS[, c("mse","nactive")], simulation_results_PRS$Method, mat = TRUE)
affect.mat2 <- affect.mat2[which(affect.mat2$group1 %in% c("Adaptivesailweak", "lasso", "lassoBT", "sail","sailweak")),]

## ---- load-PRS-data ------------------------------------------------------

gen3pc <- read.table(here::here("manuscript/raw_data/Gen_3PC_scores.txt"), header = TRUE)
gen3pc <- cbind(gen3pc,rownames(gen3pc))
colnames(gen3pc)[4] <- "SentrixID"
iq_md <- read.table(here::here("manuscript/raw_data/IQ_and_mental_development_variables_for_Sahir_with_study_ID.txt"),
                    header=TRUE)
snp_prs_na <- read.table(here::here("manuscript/raw_data/NFP_170614_INFO08_nodup_hard09_noambi_GWAS_EduYears_Pooled_beta_withaf_5000pruned_noambi_16Jan2018.score"),
                         header = TRUE, sep = ",")

# Merge the iq_md (SentrixID), snp_prs_na (Label.1) and gen3pc (SentrixID) datasets together
m1 <- merge(iq_md,snp_prs_na, by.x = "SentrixID",by.y = "Label.1")
IQdat <- merge(m1,gen3pc, by.x = "SentrixID", by.y = "SentrixID")

# Keep relavant columns
IQ4y <- subset(IQdat,
               select=c("IQ_4yrs", #Y variable
                        "Tx_group_bin", #E variable
                        "PRS_0.0001","PRS_0.001","PRS_0.01",
                        "PRS_0.05","PRS_0.1","PRS_0.2","PRS_0.3",
                        "PRS_0.4","PRS_0.5"))   #X variables

pt <- mice::mice(IQ4y, m = 1, seed = 1234)
IQ4y <- mice::complete(pt)
IQ4y <- as.matrix(IQ4y)
rownames(IQ4y) <- IQdat[,"SentrixID"]
IQ4y <- IQ4y[!is.na(IQ4y[,"IQ_4yrs"]),]

# cleanup
rm(list = c("IQdat","m1","snp_prs_na","iq_md","pt","gen3pc"))

## ---- PRS-cvsail-imputation-1----
doMC::registerDoMC(cores = 10)
fitcv <- cv.sail(x = IQ4y[,-c(1,2)], y = IQ4y[,"IQ_4yrs"], e = IQ4y[, "Tx_group_bin"],
                 basis = function(i) splines::bs(i, degree = 3),
                 expand = TRUE,
                 center.x = T,
                 center.e = T,
                 #group = dat$group,
                 alpha = 0.1,
                 #maxit = 250,
                 strong = FALSE,
                 verbose = 1,
                 parallel = TRUE
)


## ---- PRS-intervention-interaction ----


plotInterPRS(object = fitcv$sail.fit,
             originalDataNotCentered = IQ4y, # complete original data which contains both "x" and "e"
             xvar = "PRS_0.0001",
             s = fitcv$lambda.min,
             design = fitcv$sail.fit$design, # this contains the design matrix from sail
             evar = "Tx_group_bin", # this is the name of the E column in 'originalDataNotCentered'
             xlab = "Polygenic risk score at 0.0001 level of significance",
             degree = 3,
             ylab = "Marginal Risk",
             legend.position = "bottomright",
             #         main = "Effect of Intervention and PRS
             # on IQ at 4 years of age",
             rug = FALSE,
             labels = c("-0.0015", "-0.001", "-0.0005", "0", "0.0005", "0.001"),
             at = c(-0.0015, -0.001, -5e-04, 0, 5e-04, 0.001),
             color = sail:::cbbPalette[c(2,4)],
             legend = TRUE)




## ---- PRS-cvsail-imputation-others----

gen3pc <- read.table(here::here("manuscript/raw_data/Gen_3PC_scores.txt"), header = TRUE)
gen3pc <- cbind(gen3pc,rownames(gen3pc))
colnames(gen3pc)[4] <- "SentrixID"
iq_md <- read.table(here::here("manuscript/raw_data/IQ_and_mental_development_variables_for_Sahir_with_study_ID.txt"),
                    header=TRUE)
snp_prs_na <- read.table(here::here("manuscript/raw_data/NFP_170614_INFO08_nodup_hard09_noambi_GWAS_EduYears_Pooled_beta_withaf_5000pruned_noambi_16Jan2018.score"),
                         header = TRUE, sep = ",")

# Merge the iq_md (SentrixID), snp_prs_na (Label.1) and gen3pc (SentrixID) datasets together
m1 <- merge(iq_md,snp_prs_na, by.x = "SentrixID",by.y = "Label.1")
IQdat <- merge(m1,gen3pc, by.x = "SentrixID", by.y = "SentrixID")

# Keep relavant columns
IQ4y <- subset(IQdat,
               select=c("IQ_4yrs", #Y variable
                        "Tx_group_bin", #E variable
                        "PRS_0.0001","PRS_0.001","PRS_0.01",
                        "PRS_0.05","PRS_0.1","PRS_0.2","PRS_0.3",
                        "PRS_0.4","PRS_0.5"))   #X variables

pt <- mice::mice(IQ4y, m = 5, seed = 1234)

IQ4y <- mice::complete(pt, action = "all")
IQ4y <- lapply(IQ4y, as.matrix)
IQ4y <- lapply(IQ4y, function(i){
  rownames(i) <- IQdat[,"SentrixID"]
  i
})

# cleanup
rm(list = c("IQdat","m1","snp_prs_na","iq_md","pt","gen3pc"))

doMC::registerDoMC(cores = 10)

fitcvALL <- lapply(IQ4y, function(DT) {
  cv.sail(x = DT[,-c(1,2)], y = DT[,"IQ_4yrs"], e = DT[, "Tx_group_bin"],
          basis = function(i) splines::bs(i, degree = 3),
          expand = TRUE,
          center.x = T,
          center.e = T,
          #group = dat$group,
          alpha = 0.1,
          #maxit = 250,
          strong = FALSE,
          verbose = 1,
          parallel = TRUE
  )
}
)

## ---- cvsail-bootstrap-first-imputation----

# here we impute once and then run bootstraps on that first imputed data
# run on tmux
pacman::p_load(psych)
pacman::p_load(magrittr)
pacman::p_load(data.table)
pacman::p_load(sail)
pacman::p_load(mice)
pacman::p_load(doMC)
source("../raw_data/plot_functions_rda.R")

gen3pc <- read.table("../raw_data/Gen_3PC_scores.txt", header = TRUE)
gen3pc <- cbind(gen3pc,rownames(gen3pc))
colnames(gen3pc)[4] <- "SentrixID"
iq_md <- read.table("../raw_data/IQ_and_mental_development_variables_for_Sahir_with_study_ID.txt",
                    header=TRUE)
snp_prs_na <- read.table("../raw_data/NFP_170614_INFO08_nodup_hard09_noambi_GWAS_EduYears_Pooled_beta_withaf_5000pruned_noambi_16Jan2018.score",
                         header = TRUE, sep = ",")

# Merge the iq_md (SentrixID), snp_prs_na (Label.1) and gen3pc (SentrixID) datasets together
m1 <- merge(iq_md,snp_prs_na, by.x = "SentrixID",by.y = "Label.1")
IQdat <- merge(m1,gen3pc, by.x = "SentrixID", by.y = "SentrixID")

# Keep relavant columns
IQ4y <- subset(IQdat,
               select=c("IQ_4yrs", #Y variable
                        "Tx_group_bin", #E variable
                        "PRS_0.0001","PRS_0.001","PRS_0.01",
                        "PRS_0.05","PRS_0.1","PRS_0.2","PRS_0.3",
                        "PRS_0.4","PRS_0.5"))   #X variables

pt <- mice::mice(IQ4y, m = 1, seed = 1234)
IQ4y <- mice::complete(pt, action = "all")
IQ4y <- lapply(IQ4y, as.matrix)
IQ4y <- lapply(IQ4y, function(i){
  rownames(i) <- IQdat[,"SentrixID"]
  i
})

# cleanup
rm(list = c("IQdat","m1","snp_prs_na","iq_md","pt","gen3pc"))

doMC::registerDoMC(cores = 10)
fitcvALL <- lapply(1:199, function(i) {
  DT <- IQ4y$`1`
  DT <- DT[sample(1:nrow(DT),replace = TRUE),]
  message(paste("Bootstrap sample ",i))
  tryCatch({
    fitcv <- cv.sail(x = DT[,-c(1,2)], y = DT[,"IQ_4yrs"], e = DT[, "Tx_group_bin"],
                     basis = function(i) splines::bs(i, degree = 3),
                     expand = TRUE,
                     center.x = T,
                     center.e = T,
                     #group = dat$group,
                     alpha = 0.1,
                     #maxit = 250,
                     strong = FALSE,
                     verbose = 1,
                     parallel = TRUE
    )

    return(list(DT = DT, fitcv = fitcv))

  },
  error = function(err) {
    return(NA)
  }
  )
}
)

saveRDS(object = fitcvALL, file = "rda/PRS/sail_bootstrap_impute.rds")
ls()

fitcvALL2 <- readRDS(file = "rda/PRS/sail_bootstrap_impute.rds")
which(is.na(fitcvALL2))
fitcvALL2 <- fitcvALL2[!is.na(fitcvALL2)]
fitcvALL2[[99]]$fitcv %>% coef(s="lambda.min")


jj=15
plotInterPRS(object = fitcvALL2[[jj]]$fitcv$sail.fit,
             originalDataNotCentered = fitcvALL2[[jj]]$DT, # complete original data which contains both "x" and "e"
             xvar = "PRS_0.0001",
             s = fitcvALL2[[jj]]$fitcv$lambda.min,
             design = fitcvALL2[[jj]]$fitcv$sail.fit$design, # this contains the design matrix from sail
             evar = "Tx_group_bin", # this is the name of the E column in 'originalDataNotCentered'
             xlab = "PRS_0.0001",
             degree = 3,
             ylab = "Marginal Risk",
             legend.position = "bottomright",
             main = "Effect of Intervention and PRS
     on IQ at 4 years of age",
             rug = TRUE,
             color = sail:::cbbPalette[c(2,4)],
             legend = TRUE)


selected_vars <- as.matrix(do.call(cbind, lapply(fitcvALL2, function(i) coef(i$fitcv, s = "lambda.min"))))[-1,]
colnames(selected_vars) <- paste("Bootstrap",seq_along(fitcvALL2),sep = " ")
# create heatmap using pheatmap
pacman::p_load(pheatmap)


paletteLength <- 31 # should be odd numbered
myColor <- c(
  colorRampPalette(RColorBrewer::brewer.pal(9, "Reds")[5:9])(paletteLength/2),
  "white","white",
  colorRampPalette(RColorBrewer::brewer.pal(9, "Blues")[5:9])(paletteLength/2)
)
myBreaks <- c(
  seq(min(selected_vars[selected_vars<0]),max(selected_vars[selected_vars<0]), length.out = ceiling(paletteLength/2) + 1),
  max(selected_vars[selected_vars<0])*0.8, 0, min(selected_vars[selected_vars>0])*0.8, # this part is white
  seq(min(selected_vars[selected_vars>0]), max(selected_vars), length.out = floor(paletteLength/2))
)

# dev.off()
# RColorBrewer::brewer.pal(9, "Reds")[5:9]
# viridis::viridis()

pheatmap(selected_vars, cluster_rows = FALSE, cluster_cols = FALSE,
         color=myColor, breaks=myBreaks,show_colnames = FALSE)
# color = colorRampPalette(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))(100))

pacman::p_load(UpSetR)



## ---- PRS-intervention-interaction-others ----

par(mfrow = c(2,3))
lapply(seq_along(fitcvALL), function(i){
  plotInterPRS(object = fitcvALL[[i]]$sail.fit,
               originalDataNotCentered = IQ4y[[i]], # complete original data which contains both "x" and "e"
               xvar = "PRS_0.0001",
               s = fitcvALL[[i]]$lambda.min,
               design = fitcvALL[[i]]$sail.fit$design, # this contains the design matrix from sail
               evar = "Tx_group_bin", # this is the name of the E column in 'originalDataNotCentered'
               xlab = "PRS_0.0001",
               degree = 3,
               ylab = "Marginal Risk",
               legend.position = "bottomright",
               main = "Effect of Intervention and PRS
     on IQ at 4 years of age",
               rug = TRUE,
               color = sail:::cbbPalette[c(2,4)],
               legend = FALSE)
}
)
color = sail:::cbbPalette[c(2,4)]
graphics::plot.new()
legend("center",c("Intervention","Control"),col=c(color[2],color[1]),lwd=c(3,3))
# lapply(fitcvALL, plot)
# lapply(fitcvALL, coef, s = "lambda.min")

## ---- PRS-model-selection ----

selected_vars <- as.matrix(do.call(cbind, lapply(fitcvALL, coef, s = "lambda.min")))[-1,]
colnames(selected_vars) <- paste("Imputation",1:5,sep = " ")
# create heatmap using pheatmap
pacman::p_load(pheatmap)


paletteLength <- 31 # should be odd numbered
myColor <- c(
  colorRampPalette(RColorBrewer::brewer.pal(9, "Reds")[5:9])(paletteLength/2),
  "white","white",
  colorRampPalette(RColorBrewer::brewer.pal(9, "Blues")[5:9])(paletteLength/2)
)
myBreaks <- c(
  seq(min(selected_vars[selected_vars<0]),max(selected_vars[selected_vars<0]), length.out = ceiling(paletteLength/2) + 1),
  max(selected_vars[selected_vars<0])*0.8, 0, min(selected_vars[selected_vars>0])*0.8, # this part is white
  seq(min(selected_vars[selected_vars>0]), max(selected_vars), length.out = floor(paletteLength/2))
  )

# dev.off()
# RColorBrewer::brewer.pal(9, "Reds")[5:9]
# viridis::viridis()

pheatmap(selected_vars, cluster_rows = FALSE, cluster_cols = FALSE,
         color=myColor, breaks=myBreaks)
# color = colorRampPalette(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))(100))








## ---- PRS-error-crosses-plots ----

par(family="serif")
error.crosses(affect.mat2[grep("nactive", rownames(affect.mat2)),],
              affect.mat2[grep("mse", rownames(affect.mat2)),],
              labels=unique(affect.mat2$group1),
              xlab="Number of Active Variables",
              main = "IQ Data: Means (+/- 1 SE) from 200 Train/Validate/Test Splits",
              sd = FALSE,
              cex.lab = 1.4,
              cex.axis = 1.4,
              cex.main = 1.5,
              xlim = c(0, 10),
              # ylim=c(100,250),
              #ylim=c(-1000,1000),
              ylab="Test Set MSE",
              colors = sail:::cbbPalette[c(3,4,7,2)],
              pch=16, cex=2)







