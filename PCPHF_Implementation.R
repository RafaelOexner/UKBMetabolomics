#Original PCP score

setwd("/Users/rafael/Library/CloudStorage/OneDrive-King'sCollegeLondon/KCL/KCL-PhD/UKBB/MetabolomicBiomarkers")

#setup
library(caret)
library(doParallel)
library(hdnom)
library(readxl)
library(dplyr)
library(survival)
library(glmnet)
library(pROC)
library(boot)

registerDoParallel(detectCores())

#fetch df
Complete <- read.delim("Final dataframes/Before_exclusion.tsv", sep = '\t')

#exclude missing values in PCP categories
PCP <- Complete[complete.cases(Complete[c(169,170,172,173,195,337,336,212,335,208,209,210)]), ]
#change colnames
annotation <- read_excel("Final dataframes/UKB metabolite field ID.xlsx")
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

#fit models on training dataset (hdnom(glmnet), 5fold CV, alpha optimisation, defaults to optimisation of deviance using partial likelihood)
#fit_agesex <- fit_enet(as.matrix(trainData[,c(169:170)]), Surv(trainData$followup_HF, trainData$HF), nfolds = 10, alphas = seq(0.1, 0.9, 0.1), rule = "lambda.1se", seed = 123, parallel = TRUE)
#fit_agesex_met <- fit_enet(as.matrix(trainData[,c(1:170)]), Surv(trainData$followup_HF, trainData$HF), nfolds = 10, alphas = seq(0.1, 0.9, 0.1), rule = "lambda.1se", seed = 123, parallel = TRUE)
#fit_PCP <- fit_enet(as.matrix(trainData[,c(169:170,172,173,195,208,209,210,212,335,336,337)]), Surv(trainData$followup_HF, trainData$HF), nfolds = 10, alphas = seq(0.1, 0.9, 0.1), rule = "lambda.1se", seed = 123, parallel = TRUE) #includes age and sex
#fit_PCP_met <- fit_enet(as.matrix(trainData[,c(1:170,172,173,195,208,209,210,212,335,336,337)]), Surv(trainData$followup_HF, trainData$HF), nfolds = 10, alphas = seq(0.1, 0.9, 0.1), rule = "lambda.1se", seed = 123, parallel = TRUE) #includes age and sex

#Use unscaled values
testData_transformed <- testData #this has scaled values inside
testData <- Complete[rownames(Complete) %in% rownames(testData_transformed), ]
colnames(testData)

#calculate risk
calculateRisk <- function() {
  sum <- calculateCoValSum()
  difference <- calculateDifference(sum)
  exp <- calculateExp(difference)
  baselineSurvival <- c(0.98752, 0.99348, 0.98295, 0.99260)
  index <- calculateIndex()
  value <- baselineSurvival[index]^exp
  risk <- 1 - value
  return(risk)
}

calculateIndex <- function() {
  i <- numeric(nrow(testData))
  i[testData$Ethnic.background...Instance.0 == "W" & testData$Sex == "M"] <- 1
  i[testData$Ethnic.background...Instance.0 == "W" & testData$Sex == "F"] <- 2
  i[testData$Ethnic.background...Instance.0 == "B" & testData$Sex == "M"] <- 3
  i[testData$Ethnic.background...Instance.0 == "B" & testData$Sex == "F"] <- 4
  return(i)
}

calculateCoValSum <- function() {
  i <- calculateIndex()
  
  lnAge <- c(41.94101, 20.54973, 2.88334, 51.75667)
  lnAgeSquared <- c(-0.88115, 0, 0, 0)
  lnTreatedSystolicBP <- c(1.03508, 12.94937, 2.31106, 28.97791)
  lnAgeXLnTreatedSystolicBP <- c(0, -2.96923, 0, -6.59777)
  lnUntreatedSystolicBP <- c(0.91252, 11.86273, 2.17229, 28.18530)
  lnAgeXlnUntreatedSystolicBP <- c(0, -2.72538, 0, -6.42425)
  currentSmoker <- c(0.73839, 11.01752, 1.65337, 0.76532)
  lnAgeXCurrentSmoker <- c(0, -2.50777, -0.24665, 0)
  lnTreatedGlucose <- c(0.90072, 1.04503, 0.64704, 0.96695)
  lnUntreatedGlucose <- c(0.77805, 0.91807, 0.57891, 0.79561)
  lnTotalCholesterol <- c(0.49323, 0, 0, 0.32646)
  lnHDLC <- c(-0.43686, -0.07455, -0.80691, 0)
  lnBMI <- c(37.21577, 1.32948, 1.16289, 21.24763)
  lnAgeXLnBMI <- c(-8.83278, 0, 0, -5.00068)
  lnQRS <- c(0.63224, 1.06089, 0.72646, 1.27475)
  
  mean(c(41.94101, 20.54973, 2.88334, 51.75667))
  mean(c(37.21577, 1.32948, 1.16289, 21.24763))
  mean(c(0.63224, 1.06089, 0.72646, 1.27475))
  
  currentSmokerVal <- ifelse(testData$Smoking.status...Instance.0 == "TRUE", 1, 0)
  
  lnAgeVal <- log(testData$Age.at.recruitment)
  lnAgeSquaredVal <- lnAgeVal^2
  lnSystolicBPVal <- log(testData$SBPmin)
  lnAgeXLnSystolicBPVal <- lnAgeVal * lnSystolicBPVal
  lnAgeXCurrentSmokerVal <- lnAgeVal * currentSmokerVal
  lnGlucoseVal <- log(testData$Glucose...Instance.0)
  lnTotalCholesterolVal <- log(testData$Cholesterol...Instance.0)
  lnHDLCVal <- log(testData$HDL.cholesterol...Instance.0)
  lnBMIVal <- log(testData$Body.mass.index..BMI....Instance.0)
  lnAgeXLnBMIVal <- lnAgeVal * lnBMIVal
  #lnQRSVal <- log(testData$QRSDuration)
  
  val3 <- ifelse(testData$BP_Meds == "1", 
                 lnTreatedSystolicBP[i] * lnSystolicBPVal, 
                 lnUntreatedSystolicBP[i] * lnSystolicBPVal)
  val4 <- ifelse(testData$BP_Meds == "1", 
                 lnAgeXLnTreatedSystolicBP[i] * lnAgeXLnSystolicBPVal, 
                 lnAgeXlnUntreatedSystolicBP[i] * lnAgeXLnSystolicBPVal)
  val7 <- ifelse(testData$DM_codes == "1", 
                 lnTreatedGlucose[i] * lnGlucoseVal, 
                 lnUntreatedGlucose[i] * lnGlucoseVal)
  
  val <- lnAge[i] * lnAgeVal
  val2 <- lnAgeSquared[i] * lnAgeSquaredVal
  val5 <- currentSmoker[i] * currentSmokerVal
  val6 <- lnAgeXCurrentSmoker[i] * lnAgeXCurrentSmokerVal
  val8 <- lnTotalCholesterol[i] * lnTotalCholesterolVal
  val9 <- lnHDLC[i] * lnHDLCVal
  val10 <- lnBMI[i] * lnBMIVal
  val11 <- lnAgeXLnBMI[i] * lnAgeXLnBMIVal
  #val12 <- lnQRS[i] * lnQRSVal

  sum <- val + val2 + val3 + val4 + val5 + val6 + val7 + val8 + val9 + val10 + val11 #+ val12
  return(sum)
}

calculateDifference <- function(sum) {
  i <- calculateIndex()
  mean <- c(171.590, 99.7321, 28.7369, 233.978)
  difference <- sum - mean[i]
  return(difference)
}

calculateExp <- function(difference) {
  return(exp(difference))
}  

testData$Risk <- calculateRisk()

#calculate C-index
C_PCPHF <- Cindex(testData$Risk, Surv(testData$followup_HF, testData$HF))
roc_PCPHF <- roc(testData$HF, testData$Risk)
youden_PCPHF <- coords(roc_PCPHF, "best", best.method="youden")

#bootstrap for CIs
bootstrap_Cindex <- function(data, indices) {
  d <- data[indices, ]
  return(Cindex(d[["Risk"]], Surv(d$followup_HF, d$HF)))
}
bootstrap_Sensitivity <- function(data, indices) {
  d <- data[indices, ]
  roc_result <- roc(d$HF, d[["Risk"]])
  youden <- coords(roc_result, "best", best.method="youden")
  return(youden$sensitivity)
}
bootstrap_Specificity <- function(data, indices) {
  d <- data[indices, ]
  roc_result <- roc(d$HF, d[["Risk"]])
  youden <- coords(roc_result, "best", best.method="youden")
  return(youden$specificity)
}

replicates = 1000
boot_results <- list()
set.seed(123)

boot_results[["boot_C"]] <- boot(data = testData, statistic = function(data, indices) bootstrap_Cindex(data, indices), R = replicates)
boot_results[["Sensitivity"]] <- boot(data = testData, statistic = function(data, indices) bootstrap_Sensitivity(data, indices), R = replicates)
boot_results[["Specificity"]] <- boot(data = testData, statistic = function(data, indices) bootstrap_Specificity(data, indices), R = replicates)

#extract
#define function
extract_boot_results <- function(boot_obj, metric_name) {
  original_statistic <- boot_obj$t0
  ci <- boot.ci(boot_obj, type = "perc")[["percent"]][, 4:5]  # Extract CI using percentile method
  data.frame(Metric = metric_name, 
             original = round(mean(original_statistic), 3), 
             Lower_CI = round(ci[1], 3), 
             Upper_CI = round(ci[2], 3))
}
boot_summary <- data.frame()
for (metric_name in names(boot_results)) {
  boot_summary <- rbind(boot_summary, 
                        extract_boot_results(boot_results[[metric_name]], metric_name))
}
boot_summary

