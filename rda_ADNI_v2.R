table(apoe)
table(stage, )
library(magrittr)

amy_pheno <- xlsx::read.xlsx("~/Downloads/DrCelia_data.xlsx", sheetIndex = 1, stringAsFactors = FALSE)
amy_mat <- read.csv("~/git_repositories/sail/data/adni_new/csf_amyloid_final.csv", stringsAsFactors = FALSE)
covr <- read.csv("~/git_repositories/sail/data/adni_new/covariates.csv", stringsAsFactors = FALSE, sep = ";")
surv <- read.csv("~/git_repositories/sail/data/adni_new/fdg_info.csv", stringsAsFactors = FALSE, sep = ",")



sum(as.character(amy_pheno$PTID) %in% amy_mat$PTID)

sum(as.character(covr$IID) %in% amy_mat$PTID)



amy_mat$PTID %>% length()
rownames(amy_mat)

amy_pheno$MMSCORE
