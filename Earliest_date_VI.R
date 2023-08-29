setwd('/Users/Intercalation')

#Find earliest dates RAP

#packages
library(dplyr)
library(tidyr)
library(data.table)
library(splitstackshape)
library(stringr)
library(tibble)
library(lubridate)

#Verbal interviews
    #non-cancer ilness codes
    NCIC <- read.delim("/Users/rafael/Library/CloudStorage/OneDrive-King'sCollegeLondon/KCL/KCL-PhD/UKBB/MetabolicGenes/TableDownloads/NCIC_raw_participant.tsv", sep = "\t")
    NCIC <- as.data.table(NCIC)
    colnames(NCIC)
    NCIC <- cSplit(NCIC, "X20002.0.0", "|")
    NCIC <- cSplit(NCIC, "X20002.1.0", "|")
    NCIC <- cSplit(NCIC, "X20002.2.0", "|")
    NCIC <- cSplit(NCIC, "X20002.3.0", "|")
    #operation codes
    OPC <- read.delim("/Users/rafael/Library/CloudStorage/OneDrive-King'sCollegeLondon/KCL/KCL-PhD/UKBB/MetabolicGenes/TableDownloads/Operation_code_date_v1_participant.tsv", sep = "\t")
    colnames(OPC)
    OPC <- OPC %>% rename(eid = Participant.ID)
    OPC <- OPC[,c(1:5)]
    OPC <- as.data.table(OPC)
    OPC <- cSplit(OPC, "Operation.code...Instance.0", "|")
    OPC <- cSplit(OPC, "Operation.code...Instance.1", "|")
    OPC <- cSplit(OPC, "Operation.code...Instance.2", "|")
    OPC <- cSplit(OPC, "Operation.code...Instance.3", "|")
    #merge
    VI <- merge(NCIC, OPC , by = "eid")
    #save reformatted table
    write.table(VI, file='VI.tsv', quote=FALSE, sep='\t')

#Verbal interviews
    #HF
    NCIC_HF <- c("1076", "1079")
    #CAD
    NCIC_CAD <- c("1075")
    OPC_CAD <- c("1070, 1095")
    #DM
    NCIC_DM <- c("1220", "1222", "1223")
    #HT 
    NCIC_HT <- c("1065", "1072")
    
#VI
VI <- read.delim("VI.tsv", sep = "\t")
head(VI)
colnames(VI) #20002 NCIC
  
  #HF
    #NCIC at baseline
    VI$NCIC_HF_0 <- "FALSE"
    for (i in 2:30) {
      if (any(VI[,i] %in% NCIC_HF)) {
        VI$NCIC_HF_0[VI[,i] %in% NCIC_HF] <- "TRUE"
      }
    }
    #NCIC during follow-up
    VI$NCIC_HF_1 <- "FALSE"
    for (i in 31:99) {
      if (any(VI[,i] %in% NCIC_HF)) {
        VI$NCIC_HF_1[VI[,i] %in% NCIC_HF] <- "TRUE"
      }
    }
    
    head(VI[VI$NCIC_HF_0==TRUE,],10)
    table(VI$NCIC_HF_1)

  #CAD
    #NCIC at baseline
    VI$NCIC_CAD_0 <- "FALSE"
    for (i in 2:30) {
      if (any(VI[,i] %in% NCIC_CAD)) {
        VI$NCIC_CAD_0[VI[,i] %in% NCIC_CAD] <- "TRUE"
      }
    }
    #NCIC during follow-up
    VI$NCIC_CAD_1 <- "FALSE"
    for (i in 31:99) {
      if (any(VI[,i] %in% NCIC_CAD)) {
        VI$NCIC_CAD_1[VI[,i] %in% NCIC_CAD] <- "TRUE"
      }
    }
    #OP Codes at baseline
    VI$OPC_CAD_0 <- "FALSE"
    for (i in 100:117) {
      if (any(VI[,i] %in% OPC_CAD)) {
        VI$OPC_CAD_0[VI[,i] %in% OPC_CAD] <- "TRUE"
      }
    }
    #OP Codes during follow-up
    VI$OPC_CAD_1 <- "FALSE"
    for (i in 118:160) {
      if (any(VI[,i] %in% OPC_CAD)) {
        VI$OPC_CAD_1[VI[,i] %in% OPC_CAD] <- "TRUE"
      }
    }
    
    head(VI[VI$NCIC_CAD_0==TRUE,],10)
    table(VI$NCIC_CAD_1)
    
  #DM
    #NCIC at baseline
    VI$NCIC_DM_0 <- "FALSE"
    for (i in 2:30) {
      if (any(VI[,i] %in% NCIC_DM)) {
        VI$NCIC_DM_0[VI[,i] %in% NCIC_DM] <- "TRUE"
      }
    }
    #NCIC during follow-up
    VI$NCIC_DM_1 <- "FALSE"
    for (i in 31:99) {
      if (any(VI[,i] %in% NCIC_DM)) {
        VI$NCIC_DM_1[VI[,i] %in% NCIC_DM] <- "TRUE"
      }
    }
    
    head(VI[VI$NCIC_DM_0==TRUE,],10)
    table(VI$NCIC_DM_1)
    
  #HT
    #NCIC at baseline
    VI$NCIC_HT_0 <- "FALSE"
    for (i in 2:30) {
      if (any(VI[,i] %in% NCIC_HT)) {
        VI$NCIC_HT_0[VI[,i] %in% NCIC_HT] <- "TRUE"
      }
    }
    #NCIC during follow-up
    VI$NCIC_HT_1 <- "FALSE"
    for (i in 31:99) {
      if (any(VI[,i] %in% NCIC_HT)) {
        VI$NCIC_HT_1[VI[,i] %in% NCIC_HT] <- "TRUE"
      }
    }
    
    head(VI[VI$NCIC_HT_0==TRUE,],10)
    table(VI$NCIC_HT_1)
  
#save file
colnames(VI)  
write_tsv(VI, "Verbal_interview_disease_status.tsv")