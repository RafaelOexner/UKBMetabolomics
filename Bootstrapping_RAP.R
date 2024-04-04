#Revision bootstrapping RAP

#EN models revision

#setup
install.packages("caret")
install.packages("doParallel")
install.packages("hdnom")
install.packages("readxl")
install.packages("dplyr")
install.packages("survival")
install.packages("boot")
install.packages("pROC")

library(caret)
library(doParallel)
library(hdnom)
library(readxl)
library(dplyr)
library(survival)
library(glmnet)
library(boot)
library(survival)
library(pROC)

registerDoParallel(detectCores())

#fetch df
system("dx download 'RevisionMetabolomics/Before_exclusion.tsv'")
system("dx download 'RevisionMetabolomics/UKB metabolite field ID.xlsx'")

Complete <- read.delim("Before_exclusion.tsv", sep = '\t')
head(Complete)

#exclude missing values in PCP categories
PCP <- Complete[complete.cases(Complete[c(169,170,172,173,195,337,336,212,335,208,209,210)]), ]
#change colnames
annotation <- read_excel("UKB metabolite field ID.xlsx")
for (i in 1:ncol(PCP)) {
  for (j in 1:nrow(annotation)) {
    if (grepl(annotation$field_id[j], colnames(PCP)[i])) {
      colnames(PCP)[i] <- annotation$title[j]
      break #stop searching after first match
    }
  }
}
colnames(PCP)[c(169,170,172,173,195,208,209,210,212,335,336,337)] <- c("Age", "Sex", "Ethnicity", "Smoking status", "BMI", "Cholesterol", "HDL", "LDL", "Glucose (non-NMR)", "DM Medication", "HT Medication", "SBP")

#clinical features
#variables to include: 169 (Age), 170 (sex), 172 (Ethnicity), 173 (Smoking currently), 195 (BMI), 337 (SBPmin), 336 (BP meds), 212 (Glucose), 335 (DM meds), 208 (Cholesterol), 209 (HDL), 210 (LDL)
#numeric variables: 169, 195, 337, 212, 208, 209, 210
#PCP-HF numerical
PCP <- PCP %>% mutate_at(vars(c(169, 195, 337, 212, 208, 209, 210)), scale)
#categorical: 170 (M/F), 172 (White/Black), 173 (T/F), 336 (T/F), 335 (T/F), put categorical values as 0/1
PCP$Sex <- ifelse(PCP$Sex == "M", 1, 0)
PCP <- PCP[(PCP$Ethnicity=="B" | PCP$Ethnicity=="W"),] #exclude heterogenous non B/W ethnic backgrounds
PCP$Ethnicity <- ifelse(PCP$Ethnicity == "W", 1, 0)
PCP$`Smoking status` <- ifelse(PCP$`Smoking status` == TRUE, 1, 0)
PCP$`DM Medication` <- ifelse(PCP$`DM Medication` == TRUE, 1, 0)
PCP$`HT Medication` <- ifelse(PCP$`HT Medication` == TRUE, 1, 0)
#exclude patients at baseline risk
PCP <- PCP[PCP$HF_at_base==FALSE,]
PCP <- PCP[PCP$CAD_at_base==FALSE,]
PCP <- PCP[PCP$LLM_Meds==FALSE,]

#metabolomics
head(PCP)
PCP <- PCP %>% mutate_at(vars(1:168), log)
PCP <- PCP %>% mutate_at(vars(1:168), scale)
#exclude metabolite outliers
sd_values <- apply(PCP[,1:168], 2, sd)
mean_values <- apply(PCP[,1:168], 2, mean)
PCP[,1:168][PCP[,1:168] < (mean_values - 5*sd_values) | PCP[,1:168] > (mean_values + 5*sd_values)] <- NA
PCP <- PCP[complete.cases(PCP[,1:168]),]

#split dataset
set.seed(123)
trainIndex <- createDataPartition(PCP$HF, p = 0.8, list = FALSE)
trainData <- PCP[trainIndex, ]
testData <- PCP[-trainIndex, ]

#fit models on training dataset (hdnom(glmnet), 10fold CV, alpha optimisation, defaults to optimisation of deviance using partial likelihood)
fit_agesex <- fit_enet(as.matrix(trainData[,c(169:170)]), Surv(trainData$followup_HF, trainData$HF), nfolds = 10, alphas = seq(0.1, 0.9, 0.1), rule = "lambda.1se", seed = 123, parallel = TRUE)
fit_agesex_met <- fit_enet(as.matrix(trainData[,c(1:170)]), Surv(trainData$followup_HF, trainData$HF), nfolds = 10, alphas = seq(0.1, 0.9, 0.1), rule = "lambda.1se", seed = 123, parallel = TRUE)
fit_PCP <- fit_enet(as.matrix(trainData[,c(169:170,172,173,195,208,209,210,212,335,336,337)]), Surv(trainData$followup_HF, trainData$HF), nfolds = 10, alphas = seq(0.1, 0.9, 0.1), rule = "lambda.1se", seed = 123, parallel = TRUE) #includes age and sex
fit_PCP_met <- fit_enet(as.matrix(trainData[,c(1:170,172,173,195,208,209,210,212,335,336,337)]), Surv(trainData$followup_HF, trainData$HF), nfolds = 10, alphas = seq(0.1, 0.9, 0.1), rule = "lambda.1se", seed = 123, parallel = TRUE) #includes age and sex

#calculate LPs
testData$LP_agesex <- predict(fit_agesex$model, newx = as.matrix(testData[,c(169:170)]), type = "link")
testData$LP_agesex_met <- predict(fit_agesex_met$model, newx = as.matrix(testData[,c(1:170)]), type = "link")
testData$LP_PCP <- predict(fit_PCP$model, newx = as.matrix(testData[,c(169:170,172,173,195,208,209,210,212,335,336,337)]), type = "link")
testData$LP_PCP_met <- predict(fit_PCP_met$model, newx = as.matrix(testData[,c(1:170,172,173,195,208,209,210,212,335,336,337)]), type = "link")

trainData$LP_agesex <- predict(fit_agesex$model, newx = as.matrix(trainData[,c(169:170)]), type = "link")
trainData$LP_agesex_met <- predict(fit_agesex_met$model, newx = as.matrix(trainData[,c(1:170)]), type = "link")
trainData$LP_PCP <- predict(fit_PCP$model, newx = as.matrix(trainData[,c(169:170,172,173,195,208,209,210,212,335,336,337)]), type = "link")
trainData$LP_PCP_met <- predict(fit_PCP_met$model, newx = as.matrix(trainData[,c(1:170,172,173,195,208,209,210,212,335,336,337)]), type = "link")

#Calibration dfs
Calibration <- testData[!(testData$followup_HF < 10 & testData$HF == FALSE),] 
Calibration$HF10y <- ifelse(Calibration$followup_HF <= 10 & Calibration$HF == TRUE, TRUE, FALSE)
Calibration$risk_agesex <- predict(fit_agesex, as.matrix(trainData[,c(169:170)]), Surv(trainData$followup_HF, trainData$HF), newx = as.matrix(Calibration[,c(169:170)]), pred.at = 10)
Calibration$risk_agesex_met <- predict(fit_agesex_met, as.matrix(trainData[,c(1:170)]), Surv(trainData$followup_HF, trainData$HF), newx = as.matrix(Calibration[,c(1:170)]), pred.at = 10)
Calibration$risk_PCP <- predict(fit_PCP, as.matrix(trainData[,c(169:170,172,173,195,208,209,210,212,335,336,337)]), Surv(trainData$followup_HF, trainData$HF), newx = as.matrix(Calibration[,c(169:170,172,173,195,208,209,210,212,335,336,337)]), pred.at = 10)
Calibration$risk_PCP_met <- predict(fit_PCP_met, as.matrix(trainData[,c(1:170,172,173,195,208,209,210,212,335,336,337)]), Surv(trainData$followup_HF, trainData$HF), newx = as.matrix(Calibration[,c(1:170,172,173,195,208,209,210,212,335,336,337)]), pred.at = 10)

Calibration_train <- trainData[!(trainData$followup_HF < 10 & trainData$HF == FALSE),] 
Calibration_train$HF10y <- ifelse(Calibration_train$followup_HF <= 10 & Calibration_train$HF == TRUE, TRUE, FALSE)
Calibration_train$risk_agesex <- predict(fit_agesex, as.matrix(trainData[,c(169:170)]), Surv(trainData$followup_HF, trainData$HF), newx = as.matrix(Calibration_train[,c(169:170)]), pred.at = 10)
Calibration_train$risk_agesex_met <- predict(fit_agesex_met, as.matrix(trainData[,c(1:170)]), Surv(trainData$followup_HF, trainData$HF), newx = as.matrix(Calibration_train[,c(1:170)]), pred.at = 10)
Calibration_train$risk_PCP <- predict(fit_PCP, as.matrix(trainData[,c(169:170,172,173,195,208,209,210,212,335,336,337)]), Surv(trainData$followup_HF, trainData$HF), newx = as.matrix(Calibration_train[,c(169:170,172,173,195,208,209,210,212,335,336,337)]), pred.at = 10)
Calibration_train$risk_PCP_met <- predict(fit_PCP_met, as.matrix(trainData[,c(1:170,172,173,195,208,209,210,212,335,336,337)]), Surv(trainData$followup_HF, trainData$HF), newx = as.matrix(Calibration_train[,c(1:170,172,173,195,208,209,210,212,335,336,337)]), pred.at = 10)

#bootstraping performance metrics

#define NRI and IDI functions
calculate_nri <- function(df, pred_col1, pred_col2, outcome_col) {
  delta <- df[[pred_col1]] - df[[pred_col2]]
  appropriate <- ifelse(df[[outcome_col]] == TRUE, delta < 0, delta > 0)
  inappropriate <- ifelse(df[[outcome_col]] == TRUE, delta > 0, delta < 0)
  events_nri <- sum(appropriate[df[[outcome_col]] == TRUE]) / sum(df[[outcome_col]] == TRUE) - 
    sum(inappropriate[df[[outcome_col]] == TRUE]) / sum(df[[outcome_col]] == TRUE)
  non_events_nri <- sum(appropriate[df[[outcome_col]] == FALSE]) / sum(df[[outcome_col]] == FALSE) - 
    sum(inappropriate[df[[outcome_col]] == FALSE]) / sum(df[[outcome_col]] == FALSE)
  nri <- events_nri + non_events_nri
  return(list(nri = nri, events_nri = events_nri, non_events_nri = non_events_nri))
}
calculate_idi <- function(df, pred_col1, pred_col2, outcome_col) {
  mean_risk_HF <- 1-mean(df[[pred_col1]][df[[outcome_col]] == TRUE])
  mean_risk_noHF <- 1-mean(df[[pred_col1]][df[[outcome_col]] == FALSE])
  mean_risk_HF_Metabolomics <- 1-mean(df[[pred_col2]][df[[outcome_col]] == TRUE])
  mean_risk_noHF_Metabolomics <- 1-mean(df[[pred_col2]][df[[outcome_col]] == FALSE])
  slope_old <- mean_risk_HF - mean_risk_noHF
  slope_new <- mean_risk_HF_Metabolomics - mean_risk_noHF_Metabolomics
  abs_idi <- (slope_new - slope_old) * 100
  rel_idi <- abs_idi / slope_old #divided by ref model yates slope (e. g. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530349/#KWS207C26)
  return(list(abs_idi = abs_idi, rel_idi = rel_idi, slope_old = slope_old, slope_new = slope_new, mean_risk_HF = mean_risk_HF, mean_risk_noHF = mean_risk_noHF, mean_risk_HF_Metabolomics = mean_risk_HF_Metabolomics, mean_risk_noHF_Metabolomics = mean_risk_noHF_Metabolomics))
}

#define bootstrapping functions
bootstrap_Cindex <- function(data, indices, lp_column) {
  d <- data[indices, ]
  return(Cindex(d[[lp_column]], Surv(d$followup_HF, d$HF)))
}
bootstrap_Sensitivity <- function(data, indices, lp_column) {
  roc_result <- roc(data[indices, ]$HF, data[indices, ][[lp_column]][,1])
  youden <- coords(roc_result, "best", best.method="youden")
  result <- youden$sensitivity
  return(result)
}
bootstrap_Sensitivity <- function(data, indices, lp_column) {
  roc_result <- roc(data[indices, ]$HF, data[indices, ][[lp_column]][,1])
  youden <- coords(roc_result, "best", best.method="youden")
  result <- youden$sensitivity
  if(is.numeric(result) && length(result) == 1) {
    return(result)
  } else {
    return(NA)
  }
}
bootstrap_Specificity <- function(data, indices, lp_column) {
  roc_result <- roc(data[indices, ]$HF, data[indices, ][[lp_column]][,1])
  youden <- coords(roc_result, "best", best.method="youden")
  result <- youden$specificity
  return(result)
}
bootstrap_Specificity <- function(data, indices, lp_column) {
  roc_result <- roc(data[indices, ]$HF, data[indices, ][[lp_column]][,1])
  youden <- coords(roc_result, "best", best.method="youden")
  result <- youden$specificity
  if(is.numeric(result) && length(result) == 1) {
    return(result)
  } else {
    return(NA)
  }
}
bootstrap_NRI_Cases <- function(data, indices, lp_column1, lp_column2, outcome_col) {
  d <- data[indices, ]
  return(calculate_nri(d, lp_column1, lp_column2, outcome_col)$events_nri)
}
bootstrap_NRI_NonCases <- function(data, indices, lp_column1, lp_column2, outcome_col) {
  d <- data[indices, ]
  return(calculate_nri(d, lp_column1, lp_column2, outcome_col)$non_events_nri)
}
bootstrap_NRI_Overall <- function(data, indices, lp_column1, lp_column2, outcome_col) {
  d <- data[indices, ]
  return(calculate_nri(d, lp_column1, lp_column2, outcome_col)$nri)
}
bootstrap_IDI_Absolute <- function(data, indices, lp_column1, lp_column2, outcome_col) {
  d <- data[indices, ]
  return(calculate_idi(d, lp_column1, lp_column2, outcome_col)$abs_idi)
}
bootstrap_IDI_Relative <- function(data, indices, lp_column1, lp_column2, outcome_col) {
  d <- data[indices, ]
  return(calculate_idi(d, lp_column1, lp_column2, outcome_col)$rel_idi)
}
bootstrap_delta_Cindex <- function(data, indices, lp_column1, lp_column2) {
  d <- data[indices, ]
  cindex1 <- Cindex(d[[lp_column1]], Surv(d$followup_HF, d$HF))
  cindex2 <- Cindex(d[[lp_column2]], Surv(d$followup_HF, d$HF))
  return(cindex2 - cindex1)
}

#set inputs
set.seed(123)
lp_columns <- c('LP_agesex', 'LP_agesex_met', 'LP_PCP', 'LP_PCP_met')
risk_columns <- c('risk_agesex', 'risk_agesex_met', 'risk_PCP', 'risk_PCP_met')
replicates <- 1000
boot_results <- list()

#absolute metrics
for (lp_column in lp_columns) {
  boot_results[[paste0("C_index_", lp_column)]] <- boot(data = testData, statistic = function(data, indices) bootstrap_Cindex(data, indices, lp_column), R = replicates)
  boot_results[[paste0("Sensitivity_", lp_column)]] <- boot(data = testData, statistic = function(data, indices) bootstrap_Sensitivity(data, indices, lp_column), R = replicates)
  boot_results[[paste0("Specificity_", lp_column)]] <- boot(data = testData, statistic = function(data, indices) bootstrap_Specificity(data, indices, lp_column), R = replicates)
}
#relative metrics
model_pairs <- combn(lp_columns, 2, simplify = FALSE)
for (pair in model_pairs) {
  boot_results[[paste0("NRI_Cases_", pair[1], "_vs_", pair[2])]] <- boot(data = testData, statistic = function(data, indices) bootstrap_NRI_Cases(data, indices, pair[1], pair[2], "HF"), R = replicates)
  boot_results[[paste0("NRI_NonCases_", pair[1], "vs", pair[2])]] <- boot(data = testData, statistic = function(data, indices) bootstrap_NRI_NonCases(data, indices, pair[1], pair[2], "HF"), R = replicates)
  boot_results[[paste0("NRI_Overall_", pair[1], "vs", pair[2])]] <- boot(data = testData, statistic = function(data, indices) bootstrap_NRI_Overall(data, indices, pair[1], pair[2], "HF"), R = replicates)
  boot_results[[paste0("Delta_C_index_", pair[1], "_vs_", pair[2])]] <- boot(data = testData, statistic = function(data, indices) bootstrap_delta_Cindex(data, indices, pair[1], pair[2]), R = replicates)
}
model_pairs <- combn(risk_columns, 2, simplify = FALSE)
for (pair in model_pairs) {
  boot_results[[paste0("IDI_Absolute_", pair[1], "vs", pair[2])]] <- boot(data = Calibration, statistic = function(data, indices) bootstrap_IDI_Absolute(data, indices, pair[1], pair[2], "HF10y"), R = replicates)
  boot_results[[paste0("IDI_Relative_", pair[1], "vs", pair[2])]] <- boot(data = Calibration, statistic = function(data, indices) bootstrap_IDI_Relative(data, indices, pair[1], pair[2], "HF10y"), R = replicates)
}

#bootstrapping for train dataset
set.seed(123)
boot_results_train <- list()

for (lp_column in lp_columns) {
  boot_results_train[[paste0("C_index_", lp_column)]] <- boot(data = trainData, statistic = function(data, indices) bootstrap_Cindex(data, indices, lp_column), R = replicates)
  boot_results_train[[paste0("Sensitivity_", lp_column)]] <- boot(data = trainData, statistic = function(data, indices) bootstrap_Sensitivity(data, indices, lp_column), R = replicates)
  boot_results_train[[paste0("Specificity_", lp_column)]] <- boot(data = trainData, statistic = function(data, indices) bootstrap_Specificity(data, indices, lp_column), R = replicates)
}
model_pairs <- combn(lp_columns, 2, simplify = FALSE)
for (pair in model_pairs) {
  boot_results_train[[paste0("NRI_Cases_", pair[1], "_vs_", pair[2])]] <- boot(data = trainData, statistic = function(data, indices) bootstrap_NRI_Cases(data, indices, pair[1], pair[2], "HF"), R = replicates)
  boot_results_train[[paste0("NRI_NonCases_", pair[1], "vs", pair[2])]] <- boot(data = trainData, statistic = function(data, indices) bootstrap_NRI_NonCases(data, indices, pair[1], pair[2], "HF"), R = replicates)
  boot_results_train[[paste0("NRI_Overall_", pair[1], "vs", pair[2])]] <- boot(data = trainData, statistic = function(data, indices) bootstrap_NRI_Overall(data, indices, pair[1], pair[2], "HF"), R = replicates)
  boot_results_train[[paste0("Delta_C_index_", pair[1], "_vs_", pair[2])]] <- boot(data = trainData, statistic = function(data, indices) bootstrap_delta_Cindex(data, indices, pair[1], pair[2]), R = replicates)
}
model_pairs <- combn(risk_columns, 2, simplify = FALSE)
for (pair in model_pairs) {
  boot_results_train[[paste0("IDI_Absolute_", pair[1], "vs", pair[2])]] <- boot(data = Calibration_train, statistic = function(data, indices) bootstrap_IDI_Absolute(data, indices, pair[1], pair[2], "HF10y"), R = replicates)
  boot_results_train[[paste0("IDI_Relative_", pair[1], "vs", pair[2])]] <- boot(data = Calibration_train, statistic = function(data, indices) bootstrap_IDI_Relative(data, indices, pair[1], pair[2], "HF10y"), R = replicates)
}


#extract
#define function
extract_boot_results <- function(boot_obj, metric_name) {
    original_statistic <- boot_obj$t0
    ci <- boot.ci(boot_obj, type = "perc")[["percent"]][, 4:5]
    data.frame(Metric = metric_name, 
               original = round(mean(original_statistic), 3), 
               Lower_CI = round(ci[1], 3), 
               Upper_CI = round(ci[2], 3))
  }
#test
boot_summary <- data.frame()
for (metric_name in names(boot_results)) {
  boot_summary <- rbind(boot_summary, 
                        extract_boot_results(boot_results[[metric_name]], metric_name))
}
boot_summary
write.table(boot_summary, file="boot_test_170124.tsv", quote=FALSE, sep='\t')
system("dx upload 'boot_test_170124.tsv'")
#train
boot_summary_train <- data.frame()
for (metric_name in names(boot_results)) {
  boot_summary_train <- rbind(boot_summary_train, 
                              extract_boot_results(boot_results_train[[metric_name]], metric_name))
}
boot_summary_train
write.table(boot_summary_train, file="boot_train_170124.tsv", quote=FALSE, sep='\t')
system("dx upload 'boot_train_170124.tsv'")
