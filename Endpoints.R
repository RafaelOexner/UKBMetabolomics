#Define Endpoints UKB

setwd('/Users/Intercalation')

library(readr)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(magrittr)
library(survival)
library(survminer)
library(ROCR)
library(stats)
library(tidyverse)
library(stringr)
library(cmprsk)
library(ggsurvfit)
library(tidycmprsk)
library(survC1)
library(dcurves)
library(nricens)
library(Hmisc)
library(survIDINRI)
library(tcltk)
library(PredictABEL)
library(rmda)
library(splitstackshape)

              
#Time to event
  #Read in files
    #DR
    Death_Record <- read.delim("Death_record_status.tsv", sep = "\t")
    Death_Record <- Death_Record[,c(1,34:46)]
    
    #VR
    Verbal <- read.delim("Verbal_interview_disease_status.tsv", sep = '\t')
    Verbal <- Verbal[,c(1,161:182)]
    
    #ICD
    ICD <- read.delim("ICD_earliest_date_sum.tsv", sep = "\t")
    ICD <- ICD[,c(1,582:593)]
    
    #TTE
    Time_to_event <- read.delim("Time_to_event_participant.tsv", sep = "\t")
    Time_to_event <- rename(Time_to_event, 'Date_of_attending_assessment_centre' = 'X53.0.0', 
                            'Date_of_death_0'  = 'X40000.0.0', 'Date_of_death_1'  = 'X40000.1.0', 
                            'Date_lost_to_follow_up' = 'X191.0.0')
    
    #OPCS4
    OPCS <- read.delim("OPCS4_disease_status.tsv", sep = "\t")
    OPCS <- rename(OPCS, 'eid' = 'Participant.ID', 'OPCS_date_CAD' = 'earliest_date_CAD')

    #merge
    TF2 <- merge(Death_Record, Verbal, by="eid")
    TF2 <- merge(TF2, ICD, by="eid")
    TF2 <- merge(TF2, Time_to_event, by ="eid")
    TF2 <- merge(TF2, OPCS, by ="eid")
    colnames(TF2)
    

#reformat to date
cols_to_change <- c("earliest_date_HF", "earliest_date_CAD", "earliest_date_DM", "earliest_date_HT",
                    "Date_of_attending_assessment_centre", "Date_of_death_0", "Date_of_death_1", 
                    "Date_lost_to_follow_up", "OPCS_date_CAD")
TF2[cols_to_change] <- lapply(TF2[cols_to_change], as.Date)

#calculate TTE
    # Find earliest date among 'earliest_date', 'Date_of_death', 'Date_lost_to_follow_up' and "2021-9-30"
    TF2$sum_date_HF <- pmin(TF2$Date_of_death_0, TF2$Date_of_death_1, 
                            TF2$earliest_date_HF, TF2$Date_lost_to_follow_up, 
                            as.Date("2021-9-30"), na.rm = TRUE)
    TF2$sum_date_CAD <- pmin(TF2$Date_of_death_0, TF2$Date_of_death_1, 
                             TF2$earliest_date_CAD, TF2$OPCS_date_CAD, TF2$Date_lost_to_follow_up, 
                             as.Date("2021-9-30"), na.rm = TRUE)

    # Calculate time to event in years
    TF2$followup_HF <- as.numeric(difftime(TF2$sum_date_HF, TF2$Date_of_attending_assessment_centre, 
                                           units = "days")) / 365.25
    TF2$followup_CAD <- as.numeric(difftime(TF2$sum_date_CAD, TF2$Date_of_attending_assessment_centre, 
                                            units = "days")) / 365.25

#Yes/no/yes at bl
TF2$HF <- ifelse(!is.na(TF2$earliest_date_HF) | TF2$Death_HF == TRUE | 
                   TF2$NCIC_HF_0 == TRUE, TRUE, FALSE)
TF2$CAD <- ifelse(!is.na(TF2$earliest_date_CAD) | TF2$Death_CAD == TRUE | 
                    !is.na(TF2$OPCS_date_CAD) | TF2$NCIC_CAD_0 == TRUE | 
                    TF2$OPC_CAD_0 == TRUE, TRUE, FALSE)

TF2$HF_at_base <- ifelse(TF2$sum_date_HF <= TF2$Date_of_attending_assessment_centre | 
                           TF2$followup_HF <= 0 | TF2$NCIC_HF_0 == TRUE, TRUE, FALSE)
TF2$CAD_at_base <- ifelse(TF2$sum_date_CAD <= TF2$Date_of_attending_assessment_centre | 
                            TF2$followup_CAD <= 0 | TF2$NCIC_CAD_0 == TRUE | 
                            TF2$OPC_CAD_0 == TRUE, TRUE, FALSE)


write_tsv(TF2, "Time_to_follow_up.tsv")
                
                