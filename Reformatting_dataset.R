#NMR Metabolomics Regression

setwd("/Users/rafael/Library/CloudStorage/OneDrive-King'sCollegeLondon/KCL/KCL-PhD/UKBB/MetabolomicBiomarkers")
library(dplyr)
library(broom)
library(ggplot2)
library(gridExtra)
library(splitstackshape)

#import NMR
NMR1 <- read.delim("TableDownloads/All_metabolomics_1_v1_participant.tsv", sep = "\t")
NMR2 <- read.delim("TableDownloads/All_metabolomics_2_v1_participant.tsv", sep = "\t")
NMR <- merge(NMR1, NMR2, by="eid")

#reformat & remove non NMR
rownames(NMR) <- NMR$eid
NMR <- NMR[,-1]
NMR <- NMR[complete.cases(NMR),] #remove everything with missing values
NMR <- NMR[apply(NMR, 1, function(x) all(x > 0)),] #remove 0/negative values
write.table(NMR, file='NMR.tsv', quote=FALSE, sep='\t')

#subset for only original NMRs
NMR <- read.delim("NMR.tsv", sep = '\t')
NMR_anno <- read.delim("UKBID 168 original.tsv", sep = "\t")
original_NMR <- NMR_anno$UKB.ID
NMR <- NMR[,original_NMR]
colnames(NMR)
#output
write.table(NMR, file='NMR_168.tsv', quote=FALSE, sep='\t')



#import stratification
Strata <- read.delim("TableDownloads/Stratification.tsv", sep = "\t")
Strata <- rename(Strata, "eid" = "Participant.ID")

#import SBP
SBP <- read.delim("TableDownloads/SBP_v1_participant.tsv", sep = "\t")
SBP <- rename(SBP,
              "SBP_auto_0" = "X4080.0.0",
              "SBP_auto_1" = "X4080.0.1",
              "SBP_man_0" = "X93.0.0",
              "SBP_man_1" = "X93.0.1")
Strata <- merge(Strata, SBP, by="eid")

#reformatting
Strata$Sex <- ifelse(Strata$Sex == 1, "M", "F")
Strata$UK.Biobank.assessment.centre...Instance.0 <- as.character(Strata$UK.Biobank.assessment.centre...Instance.0)
Strata$Ethnic.background...Instance.0 <- ifelse(Strata$Ethnic.background...Instance.0 %in% 
                                                  c("1","1001","1002","1003"),"W", Strata$Ethnic.background...Instance.0)
Strata$Ethnic.background...Instance.0 <- ifelse(Strata$Ethnic.background...Instance.0 %in% 
                                                  c("4","4001","4002","4003"),"B", Strata$Ethnic.background...Instance.0)
Strata$Ethnic.background...Instance.0 <- ifelse(Strata$Ethnic.background...Instance.0 %in% 
                                                  c("2","2001","2002","2003","2004",
                                                    "3","3001","3002","3003","3004",
                                                    "5","6"),"Other", Strata$Ethnic.background...Instance.0)
Strata$Ethnic.background...Instance.0 <- ifelse(Strata$Ethnic.background...Instance.0 %in% 
                                                  c("-1","-3"),NA, Strata$Ethnic.background...Instance.0)
Strata$Smoking.status...Instance.0 <- ifelse(Strata$Smoking.status...Instance.0 %in% 
                                               c("2"), "TRUE", Strata$Smoking.status...Instance.0)
Strata$Smoking.status...Instance.0 <- ifelse(Strata$Smoking.status...Instance.0 %in% 
                                               c("0","1"), "FALSE", Strata$Smoking.status...Instance.0)
Strata$Smoking.status...Instance.0 <- ifelse(Strata$Smoking.status...Instance.0 %in% 
                                               c("-3"), NA, Strata$Smoking.status...Instance.0)
Strata$Alcohol.drinker.status...Instance.0 <- ifelse(Strata$Alcohol.drinker.status...Instance.0 %in% 
                                               c("2"), "TRUE", Strata$Alcohol.drinker.status...Instance.0)
Strata$Alcohol.drinker.status...Instance.0 <- ifelse(Strata$Alcohol.drinker.status...Instance.0 %in% 
                                               c("0","1"), "FALSE", Strata$Alcohol.drinker.status...Instance.0)
Strata$Alcohol.drinker.status...Instance.0 <- ifelse(Strata$Alcohol.drinker.status...Instance.0 %in% 
                                               c("-3"), NA, Strata$Alcohol.drinker.status...Instance.0)
Strata$WHR <- Strata$Waist.circumference...Instance.0 / Strata$Hip.circumference...Instance.0
Strata$Spectrometer...Instance.0 <- as.character(Strata$Spectrometer...Instance.0)

#summarize family history
Strata <- as.data.frame(Strata)
Strata <- cSplit(Strata, c( "Illnesses.of.father...Instance.0", 
                            "Illnesses.of.mother...Instance.0", 
                            "Illnesses.of.siblings...Instance.0" ), "|")
Strata <- as.data.frame (Strata)
Strata$Familyhistory_heartdisease <- FALSE
Strata$Familyhistory_heartdisease[apply(Strata[, 69:104] == 1, 1, any)] <- TRUE
Strata$Familyhistory_heartdisease[apply(is.na(Strata[69]), 1, any)] <- NA
Strata$Familyhistory_heartdisease[apply(Strata[, 69:104] == c(-13,-23), 1, any)] <- NA
Strata$Familyhistory_DM <- FALSE
Strata$Familyhistory_DM[apply(Strata[, 69:104] == 9, 1, any)] <- TRUE
Strata$Familyhistory_DM[apply(Strata[, 69:104] == c(-13,-23), 1, any)] <- NA
Strata$Familyhistory_DM[apply(is.na(Strata[69]), 1, any)] <- NA

#look for medication
Strata <- as.data.frame(Strata)
Strata <- cSplit(Strata, c("Medication.for.cholesterol..blood.pressure.or.diabetes...Instance.0", 
                           "Medication.for.cholesterol..blood.pressure..diabetes..or.take.exogenous.hormones...Instance.0"), "|")
Strata <- as.data.frame(Strata)
Strata <- cSplit(Strata, "Treatment.medication.code...Instance.0" , "|")
Strata <- as.data.frame(Strata)
#touchscreen
colnames(Strata)
Strata$LLMed_TS <- apply(Strata[,104:111], 1, function(x) any(x == 1))
Strata$BPMed_TS <- apply(Strata[,104:111], 1, function(x) any(x == 2))
Strata$DM_TS <- apply(Strata[,104:111], 1, function(x) any(x == 3))
#interview
LLM <- c(1141146234,1140888594,1140888648,1141192410,1140861958,1140865576,1140865576,
         1140861936,1140861944,1140862028,1141175908,1140851880,1140851882,1140861324,
         1140861868,1140861892,1141162544,1141192736,1141192740,1141157416,1140861924,
         1141157260,1140861926,1140861928,1140861922,1140861942,1140861946,1140861954,
         1140862026,1141168568,1141171548,1141201306,1140888590,1140861848,1140861856,
         1141157262,1140861858,1140926582,1140861866,1141188546,1140861876,1140861878,
         1140861884,1141181868,1141172214,1141182910,1140865752,1141157494,1141145830)
LLM <- as.character(LLM)
APsy <- c(1140867420,1140867432,1140867444,1140927956,1140927970,1140928916,1141152848,1141152860,1141153490,1141167976,1141177762,1141195974,1141202024)
APsy <- as.character(APsy)
BPMed <- c(1140860332,1140860334,1140860336,1140860338,1140860340,1140860342,1140860348,1140860352,1140860356,1140860358,1140860362,1140860380,1140860390,1140860394,1140860396,1140860398,1140860402,1140860410,1140860418,1140860422,1140860426,1140860434,1140860478,1140860492,1140860498,1140860520,1140860532,1140860552,1140860558,1140860562,1140860564,1140860580,1140860628,1140860632,1140860638,1140860654,1140860658,1140860706,1140860714,1140860728,1140860736,1140860738,1140860758,1140860764,1140860776,1140860784,1140860790,1140860828,1140860830,1140860834,1140860836,1140860838,1140860846,1140860848,1140860862,1140860878,1140860882,1140860912,1140860918,1140860938,1140860942,1140860952,1140860972,1140860976,1140860982,1140860988,1140860994,1140861008,1140861010,1140861016,1140861022,1140861024,1140861068,1140861070,1140861088,1140861090,1140861106,1140861120,1140861128,1140861130,1140861136,1140861138,1140861190,1140861194,1140861202,1140861266,1140861268,1140861326,1140861384,1140864950,1140864952,1140866072,1140866084,1140866086,1140866090,1140866092,1140866094,1140866104,1140866108,1140866110,1140866116,1140866122,1140866136,1140866138,1140866140,1140866144,1140866146,1140866162,1140866164,1140866168,1140866182,1140866192,1140866202,1140866206,1140866210,1140866212,1140866220,1140866230,1140866232,1140866236,1140866244,1140866248,1140866282,1140866306,1140866308,1140866312,1140866318,1140866330,1140866332,1140866334,1140866340,1140866352,1140866360,1140866388,1140866390,1140866396,1140866400,1140866406,1140866408,1140866410,1140866412,1140866416,1140866422,1140866426,1140866438,1140866440,1140866442,1140866448,1140866450,1140866460,1140866466,1140866484,1140866554,1140866692,1140866704,1140866712,1140866724,1140866756,1140866758,1140866764,1140866766,1140866778,1140866798,1140866800,1140866802,1140866804,1140875808,1140879762,1140879778,1140879782,1140879786,1140879794,1140879806,1140879810,1140879818,1140879822,1140879824,1140879834,1140879842,1140879854,1140879866,1140888510,1140888556,1140888560,1140888578,1140888582,1140888586,1140888760,1140888762,1140909368,1140911698,1140916356,1140923572,1140923712,1140923718,1140926778,1140926780,1141145668,1141151016,1141151018,1141151382,1141152600,1141153026,1141153032,1141153328,1141156754,1141156808,1141157252,1141157254,1141164148,1141164154,1141164276,1141165476,1141166006,1141167822,1141167832,1141171152,1141172682,1141172686,1141172698,1141173888,1141180592,1141187790,1141190160,1141192064,1141193282,1141193346,1141194804,1141194808,1141194810,1141201038,1141201040,1140860382,1140860386,1140860404,1140860406,1140860454,1140860470,1140860534,1140860544,1140860590,1140860610,1140860690,1140860696,1140860750,1140860752,1140860802,1140860806,1140860840,1140860842,1140860892,1140860904,1140860954,1140860966,1140861000,1140861002,1140861034,1140861046,1140861110,1140861114,1140861166,1140861176,1140861276,1140861282,1140866074,1140866078,1140866096,1140866102,1140866128,1140866132,1140866156,1140866158,1140866194,1140866200,1140866222,1140866226,1140866262,1140866280,1140866324,1140866328,1140866354,1140866356,1140866402,1140866404,1140866418,1140866420,1140866444,1140866446,1140866506,1140866546,1140866726,1140866738,1140866782,1140866784,1140879758,1140879760,1140879798,1140879802,1140879826,1140879830,1140888512,1140888552,1140888646,1140888686,1140916362,1140917428,1141145658,1141145660,1141152998,1141153006,1141156836,1141156846,1141164280,1141165470,1141171336,1141171344,1141180598,1141187788,1141194794,1141194800)
BPMed <- as.character(BPMed)
DMMed <- c(1140857494, 1140857496, 1140857500, 1140857502, 1140857506, 1140857584, 1140857586, 1140857590, 1140868902, 1140868908, 1140874646, 1140874650, 1140874652, 1140874658, 1140874660, 1140874664, 1140874666, 1140874674, 1140874678, 1140874680, 1140874686, 1140874690, 1140874706, 1140874712, 1140874716, 1140874718, 1140874724, 1140874726, 1140874728, 1140874732, 1140874736, 1140874740, 1140874744, 1140874746, 1140882964, 1140883066, 1140884600, 1140910564, 1140910566, 1140910818, 1140921964, 1141152590, 1141153254, 1141153262, 1141156984, 1141157284, 1141168660, 1141168668, 1141169504, 1141171508, 1141171646, 1141171652, 1141173786, 1141173882, 1141177600, 1141177606, 1141189090, 1141189094)
DMMed <- as.character(DMMed)
colnames(Strata)
Strata$LLM_codes <- FALSE
Strata$LLM_codes <- apply(Strata[,112:159], 1, function(x) any(x %in% LLM))
Strata$APsy_codes <- FALSE
Strata$APsy_codes <- apply(Strata[,112:159], 1, function(x) any(x %in% APsy))
Strata$DM_codes <- FALSE
Strata$DM_codes <- apply(Strata[,112:159], 1, function(x) any(x %in% DMMed))
Strata$BP_codes <- FALSE
Strata$BP_codes <- apply(Strata[,112:159], 1, function(x) any(x %in% BPMed))
#individuals with TRUE for either touchscreen or interview
Strata$LLM_Meds <- ifelse(Strata$LLMed_TS == TRUE | Strata$LLM_codes == TRUE, TRUE, FALSE)
Strata$LLM_Meds[is.na(Strata$LLM_Meds)] <- FALSE
Strata$APsy_Meds <- ifelse(Strata$APsy_codes == TRUE, TRUE, FALSE)
Strata$DM_Meds <- ifelse(Strata$DM_TS == TRUE | Strata$DM_codes == TRUE, TRUE, FALSE)
Strata$DM_Meds[is.na(Strata$DM_Meds)] <- FALSE
Strata$BP_Meds <- ifelse(Strata$BPMed_TS == TRUE | Strata$BP_codes == TRUE, TRUE, FALSE)
Strata$BP_Meds[is.na(Strata$BP_Meds)] <- FALSE

#reformat SBP to include minimum SBP
Strata$SBPmin <- pmin(Strata$SBP_auto_0,
                      Strata$SBP_auto_1,
                      Strata$SBP_man_0,
                      Strata$SBP_man_1,
                      na.rm=TRUE) 

#eduction
Strata <- cSplit(Strata, "Qualifications...Instance.0" , "|")
Strata$college <- FALSE
colnames(Strata)
Strata$college[apply(Strata[, 171:176] == 1, 1, any)] <- TRUE
Strata$college[apply(is.na(Strata[, 171]), 1, any)] <- NA
table(Strata$college)

#write table
Strata <- as.data.frame(Strata)
rownames(Strata) <- Strata$eid
head(Strata)
Strata <- Strata [,-1]
write.table(Strata, file='Stratification.tsv', quote=FALSE, sep='\t')


#import Endpoints
TTE <- read.delim("TableDownloads/Time_to_follow_up.tsv", sep = "\t")
row.names(TTE) <- TTE$eid
head(TTE)
TTE <- TTE[,-1]
TTE <- as.data.frame(TTE)
#add some summary diagnoses
TTE$DM <- ifelse(TTE$Death_DM==TRUE | TTE$NCIC_DM_0==TRUE |TTE$NCIC_DM_1 == TRUE | !is.na(TTE$earliest_date_DM), TRUE, FALSE)
TTE$HTN <- ifelse(TTE$Death_HT==TRUE | TTE$NCIC_HT_0 ==TRUE | TTE$NCIC_HT_1 ==TRUE | !is.na(TTE$earliest_date_HT), TRUE, FALSE)
#DM at BL
TTE$sum_date_DM <- pmin(as.Date(TTE$Date_of_death_0), as.Date(TTE$Date_of_death_1), as.Date(TTE$earliest_date_DM), as.Date(TTE$Date_lost_to_follow_up), as.Date("2021-9-30"), na.rm = TRUE)
TTE$followup_DM <- as.numeric(difftime(TTE$sum_date_DM, as.Date(TTE$Date_of_attending_assessment_centre), units = "days")) /365.25
TTE$DM_at_base <- ifelse(TTE$sum_date_DM <= TTE$Date_of_attending_assessment_centre |
                           TTE$followup_DM <= 0 |
                           TTE$NCIC_DM_0 == TRUE, TRUE, FALSE)

#merge NMR, stratification and endpoints
head(NMR)
head(Strata)
data <- merge(NMR, Strata, by = "row.names")
rownames(data) <- data$Row.names
data <- data[,-1]
data <- merge(data, TTE, by = "row.names")
rownames(data) <- data$Row.names
data <- data[,-1]
colnames(data)
head(data)
#write 
write.table(data, file='Before_exclusion.tsv', quote=FALSE, sep='\t')
