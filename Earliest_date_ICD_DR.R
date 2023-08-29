#Find earliest dates RAP

#packages
library(dplyr)
library(tidyr)
library(data.table)
library(splitstackshape)
library(stringr)
library(tibble)
library(lubridate)

#download & merge ICD files
  #in-patient all ICD10 codes
    #downloaded in-patient ICD codes (i.e. datafields 41270 & 41280) manually via the selecting groups of fields in the cohort browser and then exporting via table exporter
    ICD10_1 <- read.delim("/Users/rafael/Library/CloudStorage/OneDrive-King'sCollegeLondon/KCL/KCL-PhD/UKBB/MetabolicGenes/TableDownloads/HR_ICD10_1-50_participant.tsv", sep = "\t")
    ICD10_2 <- read.delim("/Users/rafael/Library/CloudStorage/OneDrive-King'sCollegeLondon/KCL/KCL-PhD/UKBB/MetabolicGenes/TableDownloads/HR_ICD10_51-100_participant.tsv", sep = "\t")
    ICD10_3 <- read.delim("/Users/rafael/Library/CloudStorage/OneDrive-King'sCollegeLondon/KCL/KCL-PhD/UKBB/MetabolicGenes/TableDownloads/HR_ICD10_101-150_participant.tsv", sep = "\t")
    ICD10_4 <- read.delim("/Users/rafael/Library/CloudStorage/OneDrive-King'sCollegeLondon/KCL/KCL-PhD/UKBB/MetabolicGenes/TableDownloads/HR_ICD10_151-200_participant.tsv", sep = "\t")
    ICD10_5 <- read.delim("/Users/rafael/Library/CloudStorage/OneDrive-King'sCollegeLondon/KCL/KCL-PhD/UKBB/MetabolicGenes/TableDownloads/HR_ICD10_201-242_participant.tsv", sep = "\t")
    #merge
    ICD10 <- merge(ICD10_1, ICD10_2, by = "eid")
    ICD10 <- merge(ICD10, ICD10_3, by = "eid")
    ICD10 <- merge(ICD10, ICD10_4, by = "eid")
    ICD10 <- merge(ICD10, ICD10_5, by = "eid")
    #reformat
    colnames(ICD10)
    dim(ICD10)
    ICD10t <- as.data.table(ICD10)
    dim(ICD10t)
    ICD10t <- cSplit(ICD10t, "X41270.0.0", "|")
    dim(ICD10t)
    colnames(ICD10t)
    head(ICD10t[1:5,c(1:5,245:247)])
    #save reformatted table
    write.table(ICD10t, file='ICD10t.tsv', quote=FALSE, sep='\t')
  #in patient all ICD9 codes  
    #downloaded in-patient ICD codes (i.e. datafields 41270 & 41280) manually via the selecting groups of fields in the cohort browser and then exporting via table exporter
    ICD9 <- read.delim("/Users/rafael/Library/CloudStorage/OneDrive-King'sCollegeLondon/KCL/KCL-PhD/UKBB/MetabolicGenes/TableDownloads/HR_ICD9_0-42_participant.tsv", sep = "\t")
    #reformat
    colnames(ICD9)
    dim(ICD9)
    ICD9t <- as.data.table(ICD9)
    dim(ICD9t)
    ICD9t <- cSplit(ICD9t, "X41271.0.0", "|")
    dim(ICD9t)
    colnames(ICD9t)
    #save reformatted table
    write.table(ICD9t, file='ICD9t.tsv', quote=FALSE, sep='\t')
  #merge to one big file
    #split ICD9t
    colnames(ICD9t)
    dim(ICD9t)
    ICD9t_diagnoses <- ICD9t[,c(1,49:95)]
    ICD9t_dates <- ICD9t[,c(1:48)]
    dim(ICD9t_dates)
    dim(ICD9t_diagnoses)
    #split ICD10t
    colnames(ICD10t)
    dim(ICD10t)
    ICD10t_diagnoses <- ICD10t[,c(1,245:487)]
    ICD10t_dates <- ICD10t[,c(1:244)]
    dim(ICD10t_dates)
    dim(ICD10t_diagnoses)
    #merge both
    ICD <- merge(ICD9t_dates, ICD10t_dates , by = "eid")
    ICD <- merge(ICD, ICD9t_diagnoses , by = "eid")
    ICD <- merge(ICD, ICD10t_diagnoses , by = "eid")
    dim(ICD)
  #save reformatted table
  write.table(ICD, file='ICD.tsv', quote=FALSE, sep='\t')

#process ICD file for earliest dates
    system ("dx download '/ICD.tsv'")
    ICD <- read.delim("ICD.tsv", sep = "\t")

#codes
    #HF
    ICD10_HF <- c("I110", "I130", "I132", "I255", "I420", "I425", "I428", "I429", "I500", "I501", "I509")
    ICD9_HF <- c("4254", "4280", "4281", "4289")
    ICD_HF <- c(ICD10_HF,ICD9_HF)
    #CAD
    ICD10_CAD <- c("I21", "I210", "I211", "I212", "I213", "I214", "I219",
                   "I22", "I220", "I221", "I228", "I229",
                   "I23", "I231", "I232", "I233", "I236", "I238",
                   "I241", "I252")
    ICD_CAD <- c(ICD10_CAD)
    #Diabetes
    ICD10_DM <- c("E10", "E100", "E101", "E102", "E103", "E104", "E105", "E106", "E107", "E108", "E109",
                  "E11", "E110", "E111", "E112", "E113", "E114", "E115", "E116", "E117", "E118", "E119",
                  "E13", "E130", "E131", "E132", "E133", "E136", "E138", "E139",
                  "E14", "E140", "E141", "E142", "E143", "E144", "E145", "E146", "E147", "E148", "E149")
    ICD_DM <- c(ICD10_DM)
    #Hypertension
    ICD10_HT <- c("I10",
                  "I11", "I110", "I119",
                  "I12", "I120", "I129",
                  "I13", "I130", "I131", "I132", "I139",
                  "I15", "I150", "I151", "I152", "I159")
    ICD_HT <- c(ICD10_HT)

#earliest_dates_HF
    earliest_dates_HF = c()
    for (x in 1:nrow(ICD)) { 
      dates_of_interest_HF = c()
      for (i in 292:581) {
        if (ICD[[i]][x] %in% ICD_HF) {  
          dates_of_interest_HF = c(dates_of_interest_HF, ICD[[i - 290]][x]) 
        }
      }
      if (length(dates_of_interest_HF) > 0) {
        min_date_HF = min(dates_of_interest_HF)
        earliest_dates_HF = c(earliest_dates_HF, min_date_HF)
      } else {
        earliest_dates_HF = c(earliest_dates_HF, NA)
      }
    }

#earliest_dates_CAD
    earliest_dates_CAD = c()
    for (x in 1:nrow(ICD)) { 
      dates_of_interest_CAD = c()
      for (i in 292:581) {
        if (ICD[[i]][x] %in% ICD_CAD) { 
          dates_of_interest_CAD = c(dates_of_interest_CAD, ICD[[i - 290]][x])  
        }
      }
      if (length(dates_of_interest_CAD) > 0) {
        min_date_CAD = min(dates_of_interest_CAD)
        earliest_dates_CAD = c(earliest_dates_CAD, min_date_CAD)
      } else {
        earliest_dates_CAD = c(earliest_dates_CAD, NA)
      }
    }

#earliest_dates_DM
    earliest_dates_DM = c()
    for (x in 1:nrow(ICD)) { 
      dates_of_interest_DM = c()
      for (i in 292:581) {
        if (ICD[[i]][x] %in% ICD_DM) {  
          dates_of_interest_DM = c(dates_of_interest_DM, ICD[[i - 290]][x]) 
        }
      }
      if (length(dates_of_interest_DM) > 0) {
        min_date_DM = min(dates_of_interest_DM)
        earliest_dates_DM = c(earliest_dates_DM, min_date_DM)
      } else {
        earliest_dates_DM = c(earliest_dates_DM, NA)
      }
    }

#earliest_dates_HT
    earliest_dates_HT = c()
    for (x in 1:nrow(ICD)) { 
      dates_of_interest_HT = c()
      for (i in 292:581) {
        if (ICD[[i]][x] %in% ICD_HT) {  
          dates_of_interest_HT = c(dates_of_interest_HT, ICD[[i - 290]][x])  
        }
      }
      if (length(dates_of_interest_HT) > 0) {
        min_date_HT = min(dates_of_interest_HT)
        earliest_dates_HT = c(earliest_dates_HT, min_date_HT)
      } else {
        earliest_dates_HT = c(earliest_dates_HT, NA)
      }
    }

# Add earliest_dates vector as a new column to ICD
    ICD$earliest_date_HF = earliest_dates_HF
    ICD$earliest_date_CAD = earliest_dates_CAD
    ICD$earliest_date_DM = earliest_dates_DM
    ICD$earliest_date_HT = earliest_dates_HT
    
#export
    write.table(ICD, file='ICD_earliest_date_sum.tsv', quote=FALSE, sep='\t')
    system ("dx upload 'ICD_earliest_date_sum.tsv'")
    
    

#download DR files
  #Death records
    #non-cancer ilness codes
    DR <- read.delim("/Users/rafael/Library/CloudStorage/OneDrive-King'sCollegeLondon/KCL/KCL-PhD/UKBB/MetabolicGenes/TableDownloads/DR_Complete_participant.tsv", sep = "\t")
    DR <- as.data.table(DR)
    colnames(DR)
    head(DR)
    #save reformatted table
    write.table(DR, file='DR.tsv', quote=FALSE, sep='\t')
    
    #DR
    DR <- read.delim("DR.tsv", sep = "\t")
    colnames(DR)
    #40000 is date of death, 40001 primary cause, 40002 secondary cause
    head(DR)
    DR[!is.na(DR$X40000.0.0)&!is.na(DR$X40000.1.0),]
    DR[ trimws(DR$X40000.1.0) != "", ]
    
  #HF
    #Death from HF
    DR$Death_HF <- "FALSE"
    for (i in 4:33) {
      if (any(DR[,i] %in% ICD_HF)) {
        DR$Death_HF[DR[,i] %in% ICD_HF] <- "TRUE"
      }
    }        
    
    table(DR$Death_HF)
    DR[DR$Death_HF == TRUE,]
    
  #CAD
    #Death from CAD
    DR$Death_CAD <- "FALSE"
    for (i in 4:33) {
      if (any(DR[,i] %in% ICD10_CAD)) {
        DR$Death_CAD[DR[,i] %in% ICD10_CAD] <- "TRUE"
      }
    }        
    
    table(DR$Death_CAD)
    DR[DR$Death_CAD == TRUE,]
    
  #DM
    #Death from DM
    DR$Death_DM <- "FALSE"
    for (i in 4:33) {
      if (any(DR[,i] %in% ICD10_DM)) {
        DR$Death_DM[DR[,i] %in% ICD10_DM] <- "TRUE"
      }
    }        
    
    table(DR$Death_DM)
    DR[DR$Death_DM == TRUE,]
    
  #HT
    #Death from HT
    DR$Death_HT <- "FALSE"
    for (i in 4:33) {
      if (any(DR[,i] %in% ICD10_HT)) {
        DR$Death_HT[DR[,i] %in% ICD10_HT] <- "TRUE"
      }
    }        
    
    table(DR$Death_HT)
    DR[DR$Death_HT == TRUE,]

    write_tsv(DR, "Death_record_status.tsv")
    