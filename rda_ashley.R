# MMSCORE  = APOE * amyloid beta + Age
#
# MAVAN,
#
# Results after birth ~ Prenatal Depression * snps

#' Here are the variables included:
#'
#' "PSCID" is the id of the child
#'
#' "true_id" is the unique id for PSCID x time
#'
#'
#' "time_24" if 1, then its 24 months, otherwise its 18 months (important to
#' include as covariate to account for time differences in the intercept).
#'
#'
#' "ITSEA_Att Child" Attention of the child, mother reported, using the ITSEA
#' (The outcome)
#'
#'
#' "Pren_CESD HWB6_CESD HWB12_CESD" are prenatal, 6 months and 12 months CESD
#' maternal depression scores respectively (From 0 to 60).
#'
#'
#' "looatp overf overp insmilef insmilep reachf reachp kissf kissp phywotyf
#' phywotyp matvocalf pokef pokep mouthf mouthp laughf laughp imitatef imitatep
#' awayf awayp" are microanalytic measures of mother-child interaction (how
#' often mother look at child, how much mother play with child with toys, ect;
#' see BEST Documentation.ppt for more details). These are either frequencies
#' (variable ending with f) or percentage of time spent doing it (variable
#' ending with p). I recommend z-scores all of these as they are extremely
#' skewed.
#'
#'
#' "B_DRD2 B_DRD4_78 B_DAT B_BDNF_r B_COMT_alt_r" Child genes with coding we
#' used for our papers:
#'
#' B_DRD2=1 when having at least one A (=A1=T in literature),
#'
#' B_DRD4_78=1 when 6, 7 or 8 repeat (6 extremely rare so mostly irrelevant),
#'
#' DAT=1 when 10/10 (10 associated with low dopamine),
#'
#' BDNF_r = 1 when GG (=Val),
#'
#' B_COMT_alt_r=1 when AA (Met/Met, worrier type, associated with higher
#' baseline dopamine but less stress resilience).
#'
#'
#' "B_DRD2_cat_ B_DAT_cat_ B_BDNF_cat_ B_COMT_cat_" Raw child genes coding, in
#' case you want to code them differently.
#'
#'
#' "Factor3_6m mom_age_birth gender_male" Important covariates we adjusted for,
#' but you could interact those too. Factor3_6m is the child emotional
#' self-regulation at 6 months (mother reported), mom_age_birth_x is mother age
#' at birth (centered) and gender_male = 1 is child is a boy.
#'

rm(list=ls())
devtools::document()
devtools::load_all()
# source("R/ksfilter.R")
pacman::p_load(data.table)
pacman::p_load(magrittr)
pacman::p_load(genefilter)
pacman::p_load(tidyverse)
pacman::p_load(doParallel)
pacman::p_load(mice)

DT <- read.csv("~/git_repositories/sail/data/ashley/Celia_Sahir_data.csv", stringsAsFactors = FALSE)
# dput(colnames(DT))
col_names <- c("PSCID", "true_id", "time_24",
               "ITSEA_Att", # outcome: Attention of the child
               "Pren_CESD", # exposure: are prenatal, 6 months and 12 months CESD maternal depression scores
               "HWB6_CESD", "HWB12_CESD",
               # "B_DRD2_cat_", "B_DAT_cat_", "B_BDNF_cat_", "B_COMT_cat_",
               "B_DRD2", "B_DRD4_78", "B_DAT", "B_BDNF_r", "B_COMT_alt_r", # child genes
               "looatp",
               "overf",
               # "overp",
               "insmilef",
               # "insmilep",
               "reachf",
               # "reachp",
               "kissf",
               # "kissp",
               "phywotyf",
               # "phywotyp",
               "matvocalf",
               "pokef",
               # "pokep",
               "mouthf",
               # "mouthp",
               "laughf",
               # "laughp",
               "imitatef",
               # "imitatep",
               "awayf",
               # "awayp",
               "Factor3_6m", # the child emotional self-regulation at 6 months (mother reported)
               "mom_age_birth_x", # mother age at birth (centered)
               "gender_male" )# gender of child)

DT$Pren_CESD %>% table(useNA = "always")
DT$HWB6_CESD %>% table(useNA = "always")
DT$HWB12_CESD %>% table(useNA = "always")


Hmisc::describe(DT)
DT <- DT %>% select(col_names)
# DT$B_DRD2 %>% table
# DT <- DT %>% mutate()
set.seed(12345)
DT.m <- mice::mice(DT)
mice:::densityplot.mids(DT.m)
DT.c <- complete(DT.m, 2)

# DT.c <- DT[complete.cases(DT),]

DT.18 <- DT.c %>% filter(time_24==0)
DT.24 <- DT.c %>% filter(time_24==1)

complete.cases(DT.18) %>% sum
complete.cases(DT.24) %>% sum

Y <- as.numeric(DT.18$ITSEA_Att)
Y <- as.numeric(DT.24$ITSEA_Att)
table(Y, useNA = "always")

E <- as.numeric(DT.18$Pren_CESD)
E <- as.numeric(DT.18$HWB12_CESD)
E <- as.numeric(DT.18$Factor3_6m)

E <- as.numeric(DT.24$Pren_CESD)
E <- as.numeric(DT.24$HWB6_CESD)
E <- as.numeric(DT.24$Factor3_6m)

colnames(DT.18)
X <- DT.18 %>% select(-ITSEA_Att, -PSCID, -ends_with("CESD"),-true_id,-time_24, -Factor3_6m) %>% as.matrix()
X <- DT.24 %>% select(-ITSEA_Att, -PSCID, -ends_with("CESD"),-true_id,-time_24, -Factor3_6m) %>% as.matrix()
X <- DT.24 %>% select(-ITSEA_Att, -PSCID,-true_id,-time_24, -Factor3_6m, -HWB6_CESD, -HWB12_CESD) %>% as.matrix()
dim(X)

summary(lm(Y ~ X*E))


doParallel::registerDoParallel(cores = 5)
getDoParWorkers()
cvfit <- cv.sail(x = X, y = Y, e = E, df = 3, maxit = 500,
                    nlambda.gamma = 5, nlambda.beta = 5,
                    nlambda = 25,
                 group.penalty = "SCAD",
                    # lambda.beta = exp(seq(log(0.01 * 1000), log(1000), length.out = 25)),
                    # lambda.gamma =rep(1000,25),
                    lambda.factor = 0.0001,
                    nfolds = 5,
                    thresh = 1e-5, center=T, normalize=F, verbose = T)
plot(cvfit)
coef(cvfit, s="lambda.min")
cvfit$sail.fit$beta

fmla <- as.formula(paste0("~0+",paste(colnames(X), collapse = "+"), "+ E +",paste(colnames(X),":E", collapse = "+") ))
Xmat <- model.matrix(fmla, data = data.frame(X,E))
glfit <- cv.glmnet(x = Xmat, y = Y, nfolds = 10)
plot(glfit)
coef(glfit)



