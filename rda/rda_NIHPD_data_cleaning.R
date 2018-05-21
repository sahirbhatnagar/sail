# rm(list=ls())

# load packages -----------------------------------------------------------

# pacman::p_load(data.table)
# pacman::p_load(dplyr)
# pacman::p_load(zoo)
# `%ni%` = Negate(`%in%`)


nihdata <- function(nprobes, phenoVariable, exposure, filter, data) {

  # load raw-data -----------------------------------------------------------

  # dat <- R.matlab::readMat("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/nihpd_data/NIHP_Cortical_Thickness.mat")
  # saveRDS(dat$Matrix.All.81924, file = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/nihpd_data/NIHP_Cortical_Thickness81924.rds")
  # saveRDS(dat$Matrix.All.AAL, file = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/nihpd_data/NIHP_Cortical_Thickness.rds")

  if (data == "AAL") {
    dat <- readRDS(file = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/nihpd_data/NIHP_Cortical_Thickness.rds")
  } else {
    dat <- readRDS(file = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/nihpd_data/NIHP_Cortical_Thickness81924.rds")
  }

  # head(dat$Matrix.All.AAL)
  # head(dat$Matrix.All.81924)
  # Clean cortical thickness values -----------------------------------------

  DT <- dat %>%
    as.data.table() %>%
    setnames(c("V1","V2","V3","V4"), c("ID", "age","F","M")) %>%
    setkey(ID)
  DT[,`:=`(age=age/365.25)]
  DT[, newid := paste(ID, 1:.N, sep="_"), by = ID]
  #unique id for each row
  DT[, unique_id := 1:.N]
  # categorize the age as per Cereb.Cortex paper by Budha
  # DT[, fivenum(age)]
  DT[, age_binary := cut(age, c(4.8, 11.3, Inf))]
  DT[, age_binary2 := cut(age, c(4.8, 8.4, 11.3, 14.7 ,Inf))]
  DT[, table(age_binary)]
  DT[, table(age_binary2)]
  DT[,`:=`(freq=.N), by = ID]
  # DT[, table(freq, useNA = "always")]
  # DT[,1:4, with=F]

  brain_probes <- grep("V\\d*", colnames(DT),  value = TRUE)
  # brain_probes
  # length(brain_probes)
  # setdiff(colnames(DT), brain_probes)

  setcolorder(DT, c("unique_id","ID","age","age_binary","age_binary2","F","M","freq", brain_probes))

  # get phenotypes ----------------------------------------------------------

  pheno <- as.data.table(readxl::read_excel("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/nihpd_data/NIHPD_Demographics_IQ.xls"))
  setkey(pheno, Subject_ID)
  pheno[, newid := paste(Subject_ID, 1:.N, sep="_"), by = Subject_ID]
  data.table::setnames(pheno, "QC_050514: Scores (0 = Failed; 1 = Check; 2 = Passed)", "QC")
  data.table::setnames(pheno, "Age_Days_Date_of_Visit_to_DOB", "age")
  pheno[, age:=NULL]
  pheno <- pheno[!is.na(Timepoint_Label)]
  # pheno[, time:=NULL]
  pheno[, time:=vector("integer", nrow(pheno))]
  pheno[, time := if(Timepoint_Label %in% c("v1","V1")) 1L
        else if (Timepoint_Label %in% c("v2","V2")) 2L
        else if (Timepoint_Label %in% c("v3","V3","v4")) 3L, by=list(Subject_ID, Timepoint_Label)]
  pheno[, Scanner_Serial_T1:=as.numeric(Scanner_Serial_T1)]
  pheno[, `:=`(WASI_Full_Scale_IQ=as.numeric(WASI_Full_Scale_IQ),
               WASI_Performance_IQ=as.numeric(WASI_Performance_IQ),
               WASI_Verbal_IQ=as.numeric(WASI_Verbal_IQ),
               Household_Income_Level=as.numeric(Household_Income_Level))]

  # income is in the dataset at time 1, but sometimes at time 2 or 3,
  # so need to use the na.locf to carry obseervations forwards or backward
  # pheno[, income_binary:=NULL]
  pheno[, table(HUD_Adjusted_Family_Income, useNA = "always")]
  pheno[HUD_Adjusted_Family_Income=="."]
  pheno[HUD_Adjusted_Family_Income==".", HUD_Adjusted_Family_Income := NA]

  pheno[, table(HUD_Adjusted_Family_Income, useNA = "always")]
  pheno[, HUD_Adjusted_Family_Income := na.locf(HUD_Adjusted_Family_Income, na.rm = FALSE, fromLast = FALSE), by = Subject_ID]
  pheno[, HUD_Adjusted_Family_Income := na.locf(HUD_Adjusted_Family_Income, na.rm = FALSE, fromLast = TRUE), by = Subject_ID]
  pheno[, table(HUD_Adjusted_Family_Income, useNA = "always")]
  pheno[, table(Subject_ID, HUD_Adjusted_Family_Income, useNA="always")]

  pheno[, income := factor(HUD_Adjusted_Family_Income, levels = c("Low", "Medium", "High")) ]
  pheno[income %in% c("Low","Medium"), income_binary := "Low" ]
  pheno[income %in% "High", income_binary := "High"]
  pheno[, income_binary := factor(income_binary, levels = c("Low", "High"))]
  pheno[, table(income_binary, useNA = "always")]

  # also try to create binary income variable based on Household_Income_Level variable
  # low: score from 1-7, high: score from 8-10
  # this leads to perfectly balanced groups.. 164 in the low, 165 in the high

  pheno[, table(Household_Income_Level, useNA="always")]
  pheno[, Household_Income_Level := na.locf(Household_Income_Level, na.rm = FALSE, fromLast = FALSE), by = Subject_ID]
  pheno[, Household_Income_Level := na.locf(Household_Income_Level, na.rm = FALSE, fromLast = TRUE), by = Subject_ID]

  with(pheno,xtabs(~Household_Income_Level+income_binary))
  pheno[Household_Income_Level %in% 1:7, income_binary2 := "Low"]
  pheno[Household_Income_Level %in% 8:10, income_binary2 := "High"]
  pheno[, income_binary2 := factor(income_binary2, levels = c("Low","High"))]
  pheno[, xtabs(~Household_Income_Level + income_binary2)]


  # Merge pheno with data ---------------------------------------------------

  DT_with_pheno <- merge(DT, pheno, by = "newid")
  DT_with_pheno %>% dim

  setcolorder(DT_with_pheno, c(setdiff(colnames(DT_with_pheno), brain_probes), brain_probes))

  DT_with_pheno[, Subject_Gender := factor(Subject_Gender, levels = c("Male","Female"))]
  DT_with_pheno[, Site_Location := factor(Site_Location)]
  # DT_with_pheno[, "Site_Location", with = F] %>% str
  # str(DT_with_pheno$Subject_Gender)
  DT_with_pheno[, gender_num:=ifelse(Subject_Gender=="Male", -1, 1)]
  DT_with_pheno[, table(gender_num)]

  # convert age_binary into a numeric so that the fitting functions work
  DT_with_pheno[, E := as.numeric(age_binary) - 1]
  DT_with_pheno[, E2 := as.numeric(age_binary2) - 1]

  DT_with_pheno[, income_binary := as.numeric(income_binary) - 1]
  DT_with_pheno[, income_binary2 := as.numeric(income_binary2) - 1]

  DT_with_pheno[, table(E, age_binary)]
  DT_with_pheno[, table(E2, age_binary2)]
  DT_with_pheno[E2 %in% 1:2, E2:=NA]
  DT_with_pheno[E2 %in% 3, E2:=1]
  DT_with_pheno[, table(E2, age_binary2, useNA = "always")]


  # Filter out probes -------------------------------------------------------

  # Response is WASI_Full_Scale_IQ,  p-values per probe based on
  # lm controlling for Subject_Gender + Site_Location
  load("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/nihpd_git/uni_results.RData")
  # load("/home/bhatnaga/coexpression/rda/nihpd/uni_results.RData")

  # 10,000 most variable
  # filtered_probes <- uni_results[order(sd, decreasing = TRUE)[1:topn],]$probe
  if (data != "AAL") {
    if (filter == "pvalue") {
      filtered_probes <- uni_results[order(get(filter), decreasing = FALSE)[seq(nprobes)],]$probe
    } else if (filter == "sd") {
      # this is for filtering by sd
      filtered_probes <- uni_results[order(get(filter), decreasing = TRUE)[seq(nprobes)],]$probe
    }
  } else {
    # dont filter probes for AAL data
    filtered_probes <- brain_probes
  }

  # intersect(uni_results[order(get(filter_var), decreasing = FALSE)[1:topn],]$probe,
  #           uni_results[order(pvalue, decreasing = FALSE)[1:topn],]$probe)
  # length(filtered_probes)
  # filtered_probes <- brain_probes
  # uni_results[probe %in% filtered_probes, plot(sd, -log10(pvalue))]

  null_probes <- setdiff(brain_probes, filtered_probes)

  for (jj in null_probes){
    set(DT_with_pheno, i = NULL, j = jj, value = NULL)
  }

  dim(DT_with_pheno)
  # all(filtered_probes %in% colnames(DT_with_pheno))


  # Split data by visit -----------------------------------------------------


  # data for the first timepoint (this means that its the first obseration for that person)
  first_visit <- DT_with_pheno[unique(DT_with_pheno, by="ID"), unique_id , mult = "first"]
  DT1 <- DT_with_pheno[unique_id %in% first_visit]
  DT1 %>% dim
  DT1[, .N, by = "ID"]$N %>% table

  not_first_visit <- setdiff(seq(nrow(DT_with_pheno)), first_visit)

  second_visit <- DT_with_pheno[unique_id %ni% first_visit][unique(DT_with_pheno[unique_id %ni% first_visit], by = "ID"), unique_id, mult="first"]
  second_visit %>% length()
  DT2 <- DT_with_pheno[unique_id %in% second_visit]
  dim(DT2)

  third_visit <- DT_with_pheno[unique_id %ni% c(DT1$unique_id, DT2$unique_id), unique_id]
  length(third_visit)
  DT3 <- DT_with_pheno[unique_id %in% third_visit]


  # split data into X, Y, E -------------------------------------------------



  # this is the cortical thickness data for the 1st timepoint only
  # all the data for all timepoints is in dat$Thick.81924
  xtrain <- DT1[!is.na(get(phenoVariable))][!is.na(get(exposure))][, filtered_probes, with = F] %>% as.matrix()
  ytrain <- DT1[!is.na(get(phenoVariable))][!is.na(get(exposure))][, phenoVariable, with = F] %>% as.matrix() %>% drop
  etrain <- DT1[!is.na(get(phenoVariable))][!is.na(get(exposure))][, exposure, with = F] %>% as.matrix() %>% drop
  xtrain_lasso <- cbind(E=etrain, xtrain)


  # for tuning parameter. i call this validate in the paper.. but used train/test/valid in simulator
  xtest <- DT2[!is.na(get(phenoVariable))][!is.na(get(exposure))][, filtered_probes, with = F] %>% as.matrix()
  ytest <- DT2[!is.na(get(phenoVariable))][!is.na(get(exposure))][, phenoVariable, with = F] %>% as.matrix() %>% drop
  etest <- DT2[!is.na(get(phenoVariable))][!is.na(get(exposure))][, exposure, with = F] %>% as.matrix() %>% drop
  xtest_lasso <- cbind(E=etest, xtest)


  xvalid <- DT3[!is.na(get(phenoVariable))][!is.na(get(exposure))][, filtered_probes, with = F] %>% as.matrix()
  yvalid <- DT3[!is.na(get(phenoVariable))][!is.na(get(exposure))][, phenoVariable, with = F] %>% as.matrix() %>% drop
  evalid <- DT3[!is.na(get(phenoVariable))][!is.na(get(exposure))][, exposure, with = F] %>% as.matrix() %>% drop
  xvalid_lasso <- cbind(E=evalid, xvalid)

  return(list(xtrain = xtrain, etrain = etrain, ytrain = ytrain, xtrain_lasso = xtrain_lasso,
              xtest = xtest, etest = etest, ytest = ytest, xtest_lasso = xtest_lasso,
              xvalid = xvalid, evalid = evalid, yvalid = yvalid, xvalid_lasso = xvalid_lasso))
}
