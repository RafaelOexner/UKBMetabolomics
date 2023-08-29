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

#codes
    #CAD
    OPCS4_CAD <- c("K401", "K402", "K403", "K404", 
                   "K411","K412","K413","K414", 
                   "K451","K452","K453","K454","K455", 
                   "K491","K492", "K498","K499", 
                   "K502", 
                   "K751","K752","K753","K754", "K758","K759")
  
#process OPCS4 file for earliest dates
OPCS4 <- read.delim("Operation_code_date_v1_participant.tsv", sep = "\t")
OPCS4t <- as.data.table(OPCS4)
OPCS4t <- cSplit(OPCS4t, "Operative.procedures...OPCS4", "|")
write_tsv(OPCS4t, 'OPCS4_split.tsv')

#earliest_dates_CAD
    earliest_dates_CAD = c()
    for (x in 1:nrow(OPCS4t)) { 
      dates_of_interest_CAD = c()
      for (i in 147:270) {
        if (OPCS4t[[i]][x] %in% OPCS4_CAD) {  # check if current disease code is in list of codes of interest
          dates_of_interest_CAD = c(dates_of_interest_CAD, OPCS4t[[i - 140]][x])  # add corresponding date to dates_of_interest vector
        }
      }
      if (length(dates_of_interest_CAD) > 0) {
        min_date_CAD = min(dates_of_interest_CAD)
        earliest_dates_CAD = c(earliest_dates_CAD, min_date_CAD)
      } else {
        earliest_dates_CAD = c(earliest_dates_CAD, NA)
      }
    }     

# Add earliest_dates vector as a new column to OPCS4t
OPCS4t$earliest_date_CAD = earliest_dates_CAD

#export
write.table(OPCS4t, file='OPCS4_earliest_date_sum.tsv', quote=FALSE, sep='\t')
system ("dx upload 'OPCS4_earliest_date_sum.tsv'")



#Fetch age and sex, time of enrollment and time of loss of follow-up
    #check whether individuals have disease
        #OPCS4 Codes
        Time_to_event <- read.delim('Time_to_event_participant.tsv', sep = "\t")
        Time_to_event <- arrange(Time_to_event, eid)
        
        OPCS_date <- read.delim("OPCS4_earliest_date_sum.tsv", sep = "\t")
        colnames(Time_to_event)
        OPCStte <- OPCS_date[,c("Participant.ID", "earliest_date_CAD")]
        Time_to_event <- rename(Time_to_event, 'Date_of_attending_assessment_centre' = '53-0.0', 'Date_of_death_0'  = '40000-0.0', 'Date_of_death_1'  = '40000-1.0', 'Date_lost_to_follow_up' = '191-0.0')

    # Convert date columns to Date data type
    Time_to_event$Date_of_attending_assessment_centre <- 
      as.Date(Time_to_event$Date_of_attending_assessment_centre)
    Time_to_event$Date_of_death_0 <- as.Date(Time_to_event$Date_of_death_0)
    Time_to_event$Date_of_death_1 <- as.Date(Time_to_event$Date_of_death_1)
    Time_to_event$Date_lost_to_follow_up <- as.Date(Time_to_event$Date_lost_to_follow_up)
    Time_to_event$earliest_date_CAD <- as.Date(OPCStte$earliest_date_CAD)
    
    # Find earliest date among 'earliest_date', 'Date_of_death', 'Date_lost_to_follow_up' and '2021-9-30'
    OPCStte$earliest_date_CAD <- pmin(Time_to_event$earliest_date_CAD, 
                                      Time_to_event$Date_of_death_0, Time_to_event$Date_of_death_1, 
                                      Time_to_event$Date_lost_to_follow_up, as.Date("2021-9-30"), 
                                      na.rm = TRUE)

    # Calculate time to event in years
    OPCStte$time_to_event_CAD <- as.numeric(difftime(OPCStte$earliest_date_CAD, 
                                                     Time_to_event$Date_of_attending_assessment_centre, 
                                                     units = "days")) / 365.25
    
    # Exclude NAs
    TTE <- TTE[complete.cases(TTE$Risk), ]
    OPCStte$CAD <- ifelse(OPCStte$earliest_date_CAD < Time_to_event$Date_of_attending_assessment_centre | 
                          OPCStte$time_to_event_CAD < 0, "Yes at baseline", 
                          ifelse(OPCStte$earliest_date_CAD != "2021-9-30", "Yes", "No"))
    table(OPCStte$CAD)

write_tsv(OPCStte, "OPCS4_disease_status.tsv")
