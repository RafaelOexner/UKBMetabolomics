#setup
setwd("/Users/rafael/Library/CloudStorage/OneDrive-King'sCollegeLondon/KCL/KCL-PhD/UKBB/MetabolomicBiomarkers")
library(dplyr)
library(readxl)
library(survival)
library(caret)
library(doParallel)
library(hdnom)
library(Hmisc)
library(pROC)
library(ggplot2)
library(gridExtra)
library(survminer)
library(corrplot)
library(RColorBrewer)
library(WGCNA)
library(igraph)
library(tidyverse)
library(dcurves)
library(glmnet)
library(cowplot)
library(boot)

registerDoParallel(detectCores())

#palettes
palette_models <- c("#9bd7d5", "#00999a", "#ffc79a", "#ff7400")
palette_risk_5 <- c("#d1defa", "#759bf0", "#1858e7", "#0f358a", "#05122e")
palette_risk_10 <- c("#e8eefd", "#bacdf8", "#8cacf3", "#5e8aee", "#2f69e9", "#164fd0", "#113ea1", "#0c2c73", "#071a45", "#020917")

#fetch data
Complete <- read.delim("Final dataframes/Before_exclusion.tsv", sep = '\t')
table(Complete$Sex=="F")
table(Complete$Sex=="M")

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
PCP <- PCP[(PCP$Ethnicity=="B" | PCP$Ethnicity=="W"),] 
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

fit_agesex$lambda
fit_agesex_met$lambda
fit_PCP$lambda
fit_PCP_met$lambda

fit_agesex$alpha
fit_agesex_met$alpha
fit_PCP$alpha
fit_PCP_met$alpha

#table S1 ENET coefs
coef_agesex <- as.data.frame(as.matrix(coef(fit_agesex$model)))
coef_agesex_met <- as.data.frame(as.matrix(coef(fit_agesex_met$model)))
coef_PCP <- as.data.frame(as.matrix(coef(fit_PCP$model)))
coef_PCP_met <- as.data.frame(as.matrix(coef(fit_PCP_met$model)))

write.table(coef_agesex, file="/Users/rafael/Library/CloudStorage/OneDrive-King'sCollegeLondon/KCL/KCL-PhD/Writing/UKB Serum Metabolomics to predict incident HF/Supplementary tables/coef_agesex.tsv", quote=FALSE, sep='\t')
write.table(coef_agesex_met, file="/Users/rafael/Library/CloudStorage/OneDrive-King'sCollegeLondon/KCL/KCL-PhD/Writing/UKB Serum Metabolomics to predict incident HF/Supplementary tables/coef_agesex_met.tsv", quote=FALSE, sep='\t')
write.table(coef_PCP, file="/Users/rafael/Library/CloudStorage/OneDrive-King'sCollegeLondon/KCL/KCL-PhD/Writing/UKB Serum Metabolomics to predict incident HF/Supplementary tables/coef_PCP.tsv", quote=FALSE, sep='\t')
write.table(coef_PCP_met, file="/Users/rafael/Library/CloudStorage/OneDrive-King'sCollegeLondon/KCL/KCL-PhD/Writing/UKB Serum Metabolomics to predict incident HF/Supplementary tables/coef_PCP_met.tsv", quote=FALSE, sep='\t')


#figure 1: plot individual metabolite HRs with age+sex or PCP-HF adjusted COX-PH models
#Individual metabolite COX-R
results_agesex <- data.frame(variable = character(),
                      hr = numeric(),
                      ci_low = numeric(),
                      ci_high = numeric(),
                      p_val = numeric(),
                      stringsAsFactors = FALSE)
results_PCP <- data.frame(variable = character(),
                          hr = numeric(),
                          ci_low = numeric(),
                          ci_high = numeric(),
                          p_val = numeric(),
                          stringsAsFactors = FALSE)
for (i in 1:168) {
  model <- coxph(Surv(followup_HF, HF) ~ PCP[,i] 
                 + Age + Sex,
                 data = PCP)
  hr <- exp(coef(model)[1])
  ci <- exp(confint(model)[1,])
  p_val <- summary(model)$coefficients[1, "Pr(>|z|)"]
  results_agesex <- results_agesex %>%
    add_row(variable = (colnames(PCP)[i]), hr = hr, ci_low = ci[1], ci_high = ci[2], p_val = p_val)
}
for (i in 1:168) {
  model <- coxph(Surv(followup_HF, HF) ~ PCP[,i] 
                 + Age + Sex 
                 + Ethnicity + `Smoking status` + BMI + Cholesterol + HDL + LDL + `Glucose (non-NMR)` + `DM Medication` + `HT Medication` + `SBP`,
                 data = PCP)
  hr <- exp(coef(model)[1])
  ci <- exp(confint(model)[1,])
  p_val <- summary(model)$coefficients[1, "Pr(>|z|)"]
  results_PCP <- results_PCP %>%
    add_row(variable = (colnames(PCP)[i]), hr = hr, ci_low = ci[1], ci_high = ci[2], p_val = p_val)
}

#add group label
results_agesex <- merge(results_agesex, annotation, by.x = "variable", by.y = "title", all.x = TRUE)[,c(1:5,8:9)]
results_PCP <- merge(results_PCP, annotation, by.x = "variable", by.y = "title", all.x = TRUE)[,c(1:5,8:9)]

#add adjusted pvals
results_agesex$p_adj <- p.adjust(results_agesex$p_val, "BH")
results_PCP$p_adj <- p.adjust(results_PCP$p_val, "BH")

#number of features
table(results_agesex$p_adj<0.01)
table(results_PCP$p_adj<0.01)
table(results_agesex$p_adj<0.01 & results_PCP$p_adj<0.01)

#save
write.table(results_agesex, file="/Users/rafael/Library/CloudStorage/OneDrive-King'sCollegeLondon/KCL/KCL-PhD/Writing/UKB Serum Metabolomics to predict incident HF/Supplementary tables/COXresults_agesex.tsv", quote=FALSE, sep='\t')
write.table(results_PCP, file="/Users/rafael/Library/CloudStorage/OneDrive-King'sCollegeLondon/KCL/KCL-PhD/Writing/UKB Serum Metabolomics to predict incident HF/Supplementary tables/COXresults_PCP.tsv", quote=FALSE, sep='\t')


#forrest plots
#agesex
top20_agesex <- head((results_agesex %>% arrange(p_adj)), 20)
corresponding_PCP <- results_PCP %>% filter(variable %in% top20_agesex$variable)
corresponding_PCP <- corresponding_PCP %>%
  arrange(match(variable, top20_agesex$variable))
top20_agesex$variable <- factor(top20_agesex$variable, levels = top20_agesex$variable)
corresponding_PCP$variable <- factor(corresponding_PCP$variable, levels = corresponding_PCP$variable) 
forrest_agesex <- ggplot() +
  geom_pointrange(data = corresponding_PCP, aes(x = fct_rev(variable), y = hr, ymin = ci_low, ymax = ci_high), color = "#ff7400", position = position_nudge(x = -0.1)) +
  geom_pointrange(data = top20_agesex, aes(x = variable, y = hr, ymin = ci_low, ymax = ci_high), color = "#00999a", position = position_nudge(x = 0.1)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  coord_flip() +
  xlab("") +
  ylab("Hazard Ratio") +
  ggtitle("Top 20 Metabolite Features\nAge + Sex adjusted") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0.6, 1.5)) +
  guides(color = FALSE) 
#PCP
top20_PCP <- head((results_PCP %>% arrange(p_adj)), 20)
corresponding_agesex <- results_agesex %>% filter(variable %in% top20_PCP$variable)
corresponding_agesex <- corresponding_agesex %>%
  arrange(match(variable, top20_PCP$variable))
top20_PCP$variable <- factor(top20_PCP$variable, levels = top20_PCP$variable)
corresponding_agesex$variable <- factor(corresponding_agesex$variable, levels = corresponding_agesex$variable)
forrest_PCP <- ggplot() +
  geom_pointrange(data = corresponding_agesex, aes(x = fct_rev(variable), y = hr, ymin = ci_low, ymax = ci_high), color = "#00999a", position = position_nudge(x = -0.1)) +
  geom_pointrange(data = top20_PCP, aes(x = variable, y = hr, ymin = ci_low, ymax = ci_high), color = "#ff7400", position = position_nudge(x = 0.1)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  coord_flip() +
  xlab("") +
  ylab("Hazard Ratio") +
  ggtitle("Top 20 Metabolite Features\nPCP-HF adjusted") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0.6, 1.5))
#plot
grid.arrange(forrest_agesex, forrest_PCP, nrow = 1)
#legends
legend_forrest <- data.frame(variable = c("corresponding_agesex", "top20_PCP"),
                        label = c("Adjusted for Age + Sex", "Adjusted for PCP-HF"),
                        color = c("#00999a", "#ff7400"))
legend_plot <- ggplot(legend_df, aes(x = factor(1), y = factor(1), color = variable)) +
  geom_point(size = 5) +  # Adjust the size of the dots here
  guides(color = guide_legend(title = "Legend")) +
  scale_color_manual(values = legend_df$color, labels = legend_df$label) +
  theme_void() +
  theme(legend.position = "bottom", legend.justification = "center")
print(legend_plot)

top20_PCP

#circular plots
#agesex
circ_agesex <- results_agesex
circ_agesex$Class <- ""
circ_agesex$Class[circ_agesex$Group == "Amino acids"] <- "Amino Acids"
circ_agesex$Class[circ_agesex$Group == "Ketone bodies" | circ_agesex$Group == "Fluid balance" | circ_agesex$Group == "Inflammation"] <- "Ketone bodies, Fluid balance & Inflammation"
circ_agesex$Class[circ_agesex$Group == "Apolipoproteins" | circ_agesex$Group == "Other lipids" | circ_agesex$Group == "Lipoprotein particle sizes"] <- "Apo-LP, LP sizes and other lipids"
circ_agesex$Class[circ_agesex$Group == "Glycolysis related metabolites"] <- "Glycolysis"
circ_agesex$Class[circ_agesex$Group == "Fatty acids"] <- "Fatty Acids"
circ_agesex$Class[circ_agesex$Group == "Lipoprotein subclasses" & grepl("Concentration", circ_agesex$variable)] <- "Lipoprotein particles"
circ_agesex$Class[circ_agesex$Group == "Lipoprotein particle concentrations"] <- "Lipoprotein particles"
circ_agesex$Class[circ_agesex$Group == "Lipoprotein subclasses" & grepl("Cholesterol", circ_agesex$variable)] <- "Cholesterol"
circ_agesex$Class[circ_agesex$Group == "Cholesterol"] <- "Cholesterol"
circ_agesex$Class[circ_agesex$Group == "Lipoprotein subclasses" & grepl("Cholesteryl Esters", circ_agesex$variable)] <- "Cholesteryl Esters"
circ_agesex$Class[circ_agesex$Group == "Cholesteryl esters"] <- "Cholesteryl Esters"
circ_agesex$Class[circ_agesex$Group == "Lipoprotein subclasses" & grepl("Free Cholesterol", circ_agesex$variable)] <- "Free Cholesterol"
circ_agesex$Class[circ_agesex$Group == "Free cholesterol"] <- "Free Cholesterol"
circ_agesex$Class[circ_agesex$Group == "Lipoprotein subclasses" & grepl("Phospholipids", circ_agesex$variable)] <- "Phospholipids"
circ_agesex$Class[circ_agesex$Group == "Phospholipids"] <- "Phospholipids"
circ_agesex$Class[circ_agesex$Group == "Lipoprotein subclasses" & grepl("Total Lipids", circ_agesex$variable)] <- "Total Lipids"
circ_agesex$Class[circ_agesex$Group == "Total lipids"] <- "Total Lipids"
circ_agesex$Class[circ_agesex$Group == "Lipoprotein subclasses" & grepl("Triglycerides", circ_agesex$variable)] <- "Triglycerides"
circ_agesex$Class[circ_agesex$Group == "Triglycerides"] <- "Triglycerides"
circ_agesex$Name[circ_agesex$Class == "Amino Acids"] <- circ_agesex$variable[circ_agesex$Class == "Amino Acids"]
circ_agesex$Name[circ_agesex$Name == "Total Concentration of Branched-Chain Amino Acids (Leucine + Isoleucine + Valine)"] <- "Total BCAAs"
circ_agesex$Name[circ_agesex$Class == "Apo-LP, LP sizes and other lipids"] <- c("Apo A1", "Apo B", "HDL Size", "LDL Size", "VLDL size", "Phosphatidylcholines", "Phosphoglycerides", "Sphingomyelins", "Total Cholines")
circ_agesex$Name[circ_agesex$Class == "Cholesterol"] <- c("XL-VLDL", "IDL", "L-HDL", "L-LDL", "L-VLDL", "M-HDL", "M-LDL", "M-VLDL", "S-HDL", "S-LDL", "S-VLDL", "VL-HDL", "VL-VLDL", "VS-VLDL", 
                                                    "Clinical LDL", "HDL", "LDL", "Remnant", "Total", "Total-HDLC", "VLDL")
circ_agesex$Name[circ_agesex$Class == "Cholesteryl Esters"] <- c("XL-VLDL", "HDL", "IDL", "L-HDL", "L-LDL", "L-VLDL", "LDL", "M-HDL", "M-LDL", "M-VLDL", "S-HDL", "S-LDL", "S-VLDL", "VL-HDL", "VL-VLDL", "VS-VLDL", "VLDL", "Total")
circ_agesex$Name[circ_agesex$Class == "Fatty Acids"] <- c("Unsaturation", "DHA", "LA", "MUFA", "n3-FA", "n6-FA", "PUFA", "SFA", "Total")
circ_agesex$Name[circ_agesex$Class == "Free Cholesterol"] <- c("XL-VLDL", "HDL", "IDL", "L-HDL", "L-LDL", "L-VLDL", "LDL", "M-HDL", "M-LDL", "M-VLDL", "S-HDL", "S-LDL", "S-VLDL", "VL-HDL", "VL-VLDL", "VS-VLDL", "VLDL", "Total")
circ_agesex$Name[circ_agesex$Class == "Glycolysis"] <-c("Citrate", "Glucose", "Lactate", "Pyruvate")
circ_agesex$Name[circ_agesex$Class == "Ketone bodies, Fluid balance & Inflammation"] <- c("3-Hydroxybutyrate", "Acetate", "Acetoacetate", "Acetone", "Albumin", "Creatinine", "Glycoprotein Acetyls")
circ_agesex$Name[circ_agesex$Class == "Lipoprotein particles"] <- c("XL-VLDL", "HDL", "IDL", "L-HDL", "L-LDL", "L-VLDL", "LDL", "M-HDL", "M-LDL", "M-VLDL", "S-HDL", "S-LDL", "S-VLDL", "VL-HDL", "VL-VLDL", "VS-VLDL", "VLDL", "Total")
circ_agesex$Name[circ_agesex$Class == "Phospholipids"] <- c("XL-VLDL", "HDL", "IDL", "L-HDL", "L-LDL", "L-VLDL", "LDL", "M-HDL", "M-LDL", "M-VLDL", "S-HDL", "S-LDL", "S-VLDL", "VL-HDL", "VL-VLDL", "VS-VLDL", "VLDL", "Total")
circ_agesex$Name[circ_agesex$Class == "Total Lipids"] <- c("XL-VLDL", "HDL", "IDL", "L-HDL", "L-LDL", "L-VLDL", "LDL", "Total", "M-HDL", "M-LDL", "M-VLDL", "S-HDL", "S-LDL", "S-VLDL", "VL-HDL", "VL-VLDL", "VS-VLDL", "VLDL")
circ_agesex$Name[circ_agesex$Class == "Phospholipids"] <- c("XL-VLDL", "HDL", "IDL", "L-HDL", "L-LDL", "L-VLDL", "LDL", "M-HDL", "M-LDL", "M-VLDL", "S-HDL", "S-LDL", "S-VLDL", "VL-HDL", "VL-VLDL", "VS-VLDL", "VLDL", "Total")
circ_agesex$Name[circ_agesex$Class == "Triglycerides"] <- c("Total", "XL-VLDL", "HDL", "IDL", "L-HDL", "L-LDL", "L-VLDL", "LDL", "M-HDL", "M-LDL", "M-VLDL", "S-HDL", "S-LDL", "S-VLDL", "VL-HDL", "VL-VLDL", "VS-VLDL", "VLDL")
circ_agesex$Class <- as.factor(circ_agesex$Class)
#add empty bars at the end of each group
to_add <- data.frame(matrix(NA, 2*nlevels(circ_agesex$Class), ncol(circ_agesex)))   # NA values to add at end of each groups
colnames(to_add) <- colnames(circ_agesex)   # to_add gets same colnames as results
to_add$Class <- rep(levels(circ_agesex$Class), each=2)    # empty bar added between groups
circ_agesex <- rbind(circ_agesex, to_add)   # merge
circ_agesex <- circ_agesex %>% arrange(Class)   # arrange again by group
#add some empty bars at the very end
to_add <- data.frame(matrix(NA, 10, ncol(circ_agesex)))
colnames(to_add) <- colnames(circ_agesex)
to_add$Class <- "Ketone bodies, Fluid balance & Inflammation"
circ_agesex <- rbind(circ_agesex, to_add)
#include line breaks for long class labels
circ_agesex$Class <- gsub("Apo-LP, LP sizes and other lipids", "Apo-LP, LP sizes\n& other lipids", circ_agesex$Class)
circ_agesex$Class <- gsub("Ketone bodies, Fluid balance & Inflammation", "Ketone Bodies, Fluid balance\n& Inflammation", circ_agesex$Class)
#custom reorder
custom_order <- c("Glycolysis",
                  "Amino Acids",
                  "Apo-LP, LP sizes\n& other lipids",
                  "Fatty Acids",
                  "Lipoprotein particles",
                  "Cholesterol", 
                  "Cholesteryl Esters",
                  "Free Cholesterol",
                  "Phospholipids",
                  "Total Lipids",
                  "Triglycerides",
                  "Ketone Bodies, Fluid balance\n& Inflammation")
circ_agesex <- circ_agesex[order(match(circ_agesex$Class, custom_order)), ]
custom_order_rows <- as.character(c(98,99,100,101,102,103,
  1:12,
  13:23,
  75,67,74,70,73,71,72,68,69,76:77,
  130, 114,126,116,120,123, 115, 119,117,121,124, 129,113,127,118,122,125,128, 131:132,
  42,43,41, 39,35,26,29,32, 25, 40,38,27,30,33, 44,24,36,28,31,34,37, 45:46,
  64, 48,60,50,54,57, 49, 53,51,55,58, 63,47,61,52,56,59,62, 65:66,
  95, 79,91,81,85,88, 80, 84,82,86,89, 94,78,92,83,87,90,93, 96:97,
  150, 134,146,136,140,143, 135, 139,137,141,144, 149,133,147,138,142,145,148, 151:152,
  160, 154,167,156,161,164, 155, 159,157,162,165, 170,153,168,158,163,166,169, 171:172,
  173, 175,187,177,181,184, 176, 180,178,182,185, 190,174,188,179,183,186,189, 191:192,
  104:112,
  193:202))
circ_agesex <- circ_agesex[order(match(rownames(circ_agesex), custom_order_rows)), ]
#assign ID
circ_agesex$id <- seq(1, nrow(circ_agesex))
#get name any y position of each label
label_data <- as.data.frame(circ_agesex)    # copy results to label
label_data <- label_data %>% mutate(id = row_number())    # id as row number
number_of_bar <- nrow(label_data)   # count number of bars to label
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # calculate angle of each row so it can be plotted on a circle (360 degrees)
#subtract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)    # horizontal location of the label, relative to the bar
label_data$angle <- ifelse(angle < -90, angle+180, angle)   # flip angle of label by 180 to make them readable
label_data$color <- ifelse(label_data$p_adj < 0.01, "black", "darkgrey")   # label will show significance, if not significant, biomarker will appear grey
base_data <- circ_agesex %>% 
  group_by(Class) %>% 
  dplyr::summarize(start=min(id) - 0.5, end=max(id) - 2 + 0.5) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))   # group biomarkers according to their categories, excluding empty bars. 
base_data$title[8] <- base_data$title[8]-5 #subtract half of the additional emtpy rows because its the mean
title_angle <- base_data$title*360/number_of_bar    # calculate the angle of label for each category so it changes angle as it goes around a circle
base_data$title_angle <- ifelse(90<title_angle & title_angle<270, title_angle - 180, title_angle)   #flip angle by 180 degrees if title angle is between 90 and 270 degrees
#remove final 10 rows
base_data$end[8] <- base_data$end[8]-10
#plot
circ_p_agesex <- ggplot(circ_agesex, aes(x=as.factor(id), y=hr-1, fill=hr-1)) +
  geom_bar(aes(x=as.factor(id), y=hr-1, fill=hr-1), stat="identity", alpha=0.8, na.rm = FALSE) +    # alpha changes transparency
  geom_segment(data=base_data, aes(x = start, y = 0.2, xend = end, yend = 0.2),   
               colour = "darkgrey", alpha=0.8, size=0.15, inherit.aes = FALSE) +
  geom_segment(data=base_data, aes(x = start, y = 0.1, xend = end, yend = 0.1),   
               colour = "darkgrey", alpha=0.8, size=0.15, inherit.aes = FALSE) +
  geom_segment(data=base_data, aes(x = start, y = 0, xend = end, yend = 0),
               colour = "black", alpha=0.8, size=0.3, inherit.aes = FALSE)  +
  geom_segment(data=base_data, aes(x = start, y = -0.1, xend = end, yend = -0.1),   
               colour = "darkgrey", alpha=0.8, size=0.15, inherit.aes = FALSE) +
  geom_segment(data=base_data, aes(x = start, y = -0.2, xend = end, yend = -0.2),   
               colour = "darkgrey", alpha=0.8, size=0.15, inherit.aes = FALSE) +
  annotate("text", x = rep(max(circ_agesex$id),5), y = c(-0.2, -0.1, 0, 0.1, 0.2), 
           label = c("0.8", "0.9", "1", "1.1", "1.2") , 
           color="darkgrey", size=1.5 , angle=0, fontface="bold", hjust=1) +
  scale_y_continuous(labels = function(x) x + 1, limits = c(-0.5,1)) +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0, limits=c(-0.35, 0.35)) +  
  geom_errorbar(aes(ymin = ci_low-1, ymax = ci_high-1), 
                width = 0, color = "black", alpha=0.3, na.rm = TRUE, show.legend = FALSE, size=0.25) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=0.4, label=Name, hjust=hjust),
            color=label_data$color, alpha=1, size=1.5,
            angle= label_data$angle, inherit.aes = FALSE) +
  geom_segment(data=base_data, aes(x = start, y = 0.75, xend = end, yend = 0.75),    # add grouping bracket
               colour = "black", alpha=0.8, size=0.6, inherit.aes = FALSE)  +
  geom_text(data=base_data, aes(x = title, y = 0.85, label=Class, hjust=0.5),         # Add category labels, to the centre (0.5) of each bracket
            colour = "black", alpha=0.8, size=1.8,         
            angle = -base_data$title_angle,              # Category title is angled, wraps around the graph
            fontface="bold", inherit.aes = FALSE)

#PCP
circ_PCP <- results_PCP
circ_PCP$Class <- ""
circ_PCP$Class[circ_PCP$Group == "Amino acids"] <- "Amino Acids"
circ_PCP$Class[circ_PCP$Group == "Ketone bodies" | circ_PCP$Group == "Fluid balance" | circ_PCP$Group == "Inflammation"] <- "Ketone bodies, Fluid balance & Inflammation"
circ_PCP$Class[circ_PCP$Group == "Apolipoproteins" | circ_PCP$Group == "Other lipids" | circ_PCP$Group == "Lipoprotein particle sizes"] <- "Apo-LP, LP sizes and other lipids"
circ_PCP$Class[circ_PCP$Group == "Glycolysis related metabolites"] <- "Glycolysis"
circ_PCP$Class[circ_PCP$Group == "Fatty acids"] <- "Fatty Acids"
circ_PCP$Class[circ_PCP$Group == "Lipoprotein subclasses" & grepl("Concentration", circ_PCP$variable)] <- "Lipoprotein particles"
circ_PCP$Class[circ_PCP$Group == "Lipoprotein particle concentrations"] <- "Lipoprotein particles"
circ_PCP$Class[circ_PCP$Group == "Lipoprotein subclasses" & grepl("Cholesterol", circ_PCP$variable)] <- "Cholesterol"
circ_PCP$Class[circ_PCP$Group == "Cholesterol"] <- "Cholesterol"
circ_PCP$Class[circ_PCP$Group == "Lipoprotein subclasses" & grepl("Cholesteryl Esters", circ_PCP$variable)] <- "Cholesteryl Esters"
circ_PCP$Class[circ_PCP$Group == "Cholesteryl esters"] <- "Cholesteryl Esters"
circ_PCP$Class[circ_PCP$Group == "Lipoprotein subclasses" & grepl("Free Cholesterol", circ_PCP$variable)] <- "Free Cholesterol"
circ_PCP$Class[circ_PCP$Group == "Free cholesterol"] <- "Free Cholesterol"
circ_PCP$Class[circ_PCP$Group == "Lipoprotein subclasses" & grepl("Phospholipids", circ_PCP$variable)] <- "Phospholipids"
circ_PCP$Class[circ_PCP$Group == "Phospholipids"] <- "Phospholipids"
circ_PCP$Class[circ_PCP$Group == "Lipoprotein subclasses" & grepl("Total Lipids", circ_PCP$variable)] <- "Total Lipids"
circ_PCP$Class[circ_PCP$Group == "Total lipids"] <- "Total Lipids"
circ_PCP$Class[circ_PCP$Group == "Lipoprotein subclasses" & grepl("Triglycerides", circ_PCP$variable)] <- "Triglycerides"
circ_PCP$Class[circ_PCP$Group == "Triglycerides"] <- "Triglycerides"
circ_PCP$Name[circ_PCP$Class == "Amino Acids"] <- circ_PCP$variable[circ_PCP$Class == "Amino Acids"]
circ_PCP$Name[circ_PCP$Name == "Total Concentration of Branched-Chain Amino Acids (Leucine + Isoleucine + Valine)"] <- "Total BCAAs"
circ_PCP$Name[circ_PCP$Class == "Apo-LP, LP sizes and other lipids"] <- c("Apo A1", "Apo B", "HDL Size", "LDL Size", "VLDL size", "Phosphatidylcholines", "Phosphoglycerides", "Sphingomyelins", "Total Cholines")
circ_PCP$Name[circ_PCP$Class == "Cholesterol"] <- c("XL-VLDL", "IDL", "L-HDL", "L-LDL", "L-VLDL", "M-HDL", "M-LDL", "M-VLDL", "S-HDL", "S-LDL", "S-VLDL", "VL-HDL", "VL-VLDL", "VS-VLDL", 
                                                    "Clinical LDL", "HDL", "LDL", "Remnant", "Total", "Total-HDLC", "VLDL")
circ_PCP$Name[circ_PCP$Class == "Cholesteryl Esters"] <- c("XL-VLDL", "HDL", "IDL", "L-HDL", "L-LDL", "L-VLDL", "LDL", "M-HDL", "M-LDL", "M-VLDL", "S-HDL", "S-LDL", "S-VLDL", "VL-HDL", "VL-VLDL", "VS-VLDL", "VLDL", "Total")
circ_PCP$Name[circ_PCP$Class == "Fatty Acids"] <- c("Unsaturation", "DHA", "LA", "MUFA", "n3-FA", "n6-FA", "PUFA", "SFA", "Total")
circ_PCP$Name[circ_PCP$Class == "Free Cholesterol"] <- c("XL-VLDL", "HDL", "IDL", "L-HDL", "L-LDL", "L-VLDL", "LDL", "M-HDL", "M-LDL", "M-VLDL", "S-HDL", "S-LDL", "S-VLDL", "VL-HDL", "VL-VLDL", "VS-VLDL", "VLDL", "Total")
circ_PCP$Name[circ_PCP$Class == "Glycolysis"] <-c("Citrate", "Glucose", "Lactate", "Pyruvate")
circ_PCP$Name[circ_PCP$Class == "Ketone bodies, Fluid balance & Inflammation"] <- c("3-Hydroxybutyrate", "Acetate", "Acetoacetate", "Acetone", "Albumin", "Creatinine", "Glycoprotein Acetyls")
circ_PCP$Name[circ_PCP$Class == "Lipoprotein particles"] <- c("XL-VLDL", "HDL", "IDL", "L-HDL", "L-LDL", "L-VLDL", "LDL", "M-HDL", "M-LDL", "M-VLDL", "S-HDL", "S-LDL", "S-VLDL", "VL-HDL", "VL-VLDL", "VS-VLDL", "VLDL", "Total")
circ_PCP$Name[circ_PCP$Class == "Phospholipids"] <- c("XL-VLDL", "HDL", "IDL", "L-HDL", "L-LDL", "L-VLDL", "LDL", "M-HDL", "M-LDL", "M-VLDL", "S-HDL", "S-LDL", "S-VLDL", "VL-HDL", "VL-VLDL", "VS-VLDL", "VLDL", "Total")
circ_PCP$Name[circ_PCP$Class == "Total Lipids"] <- c("XL-VLDL", "HDL", "IDL", "L-HDL", "L-LDL", "L-VLDL", "LDL", "Total", "M-HDL", "M-LDL", "M-VLDL", "S-HDL", "S-LDL", "S-VLDL", "VL-HDL", "VL-VLDL", "VS-VLDL", "VLDL")
circ_PCP$Name[circ_PCP$Class == "Phospholipids"] <- c("XL-VLDL", "HDL", "IDL", "L-HDL", "L-LDL", "L-VLDL", "LDL", "M-HDL", "M-LDL", "M-VLDL", "S-HDL", "S-LDL", "S-VLDL", "VL-HDL", "VL-VLDL", "VS-VLDL", "VLDL", "Total")
circ_PCP$Name[circ_PCP$Class == "Triglycerides"] <- c("Total", "XL-VLDL", "HDL", "IDL", "L-HDL", "L-LDL", "L-VLDL", "LDL", "M-HDL", "M-LDL", "M-VLDL", "S-HDL", "S-LDL", "S-VLDL", "VL-HDL", "VL-VLDL", "VS-VLDL", "VLDL")
circ_PCP$Class <- as.factor(circ_PCP$Class)
#add empty bars at the end of each group
to_add <- data.frame(matrix(NA, 2*nlevels(circ_PCP$Class), ncol(circ_PCP)))   # NA values to add at end of each groups
colnames(to_add) <- colnames(circ_PCP)   # to_add gets same colnames as results
to_add$Class <- rep(levels(circ_PCP$Class), each=2)    # empty bar added between groups
circ_PCP <- rbind(circ_PCP, to_add)   # merge
circ_PCP <- circ_PCP %>% arrange(Class)   # arrange again by group
#add some empty bars at the very end
to_add <- data.frame(matrix(NA, 10, ncol(circ_PCP)))
colnames(to_add) <- colnames(circ_PCP)
to_add$Class <- "Ketone bodies, Fluid balance & Inflammation"
circ_PCP <- rbind(circ_PCP, to_add)
#include line breaks for long class labels
circ_PCP$Class <- gsub("Apo-LP, LP sizes and other lipids", "Apo-LP, LP sizes\n& other lipids", circ_PCP$Class)
circ_PCP$Class <- gsub("Ketone bodies, Fluid balance & Inflammation", "Ketone Bodies, Fluid balance\n& Inflammation", circ_PCP$Class)
#custom reorder
custom_order <- c("Glycolysis",
                  "Amino Acids",
                  "Apo-LP, LP sizes\n& other lipids",
                  "Fatty Acids",
                  "Lipoprotein particles",
                  "Cholesterol", 
                  "Cholesteryl Esters",
                  "Free Cholesterol",
                  "Phospholipids",
                  "Total Lipids",
                  "Triglycerides",
                  "Ketone Bodies, Fluid balance\n& Inflammation")
circ_PCP <- circ_PCP[order(match(circ_PCP$Class, custom_order)), ]
custom_order_rows <- as.character(c(98,99,100,101,102,103,
                                    1:12,
                                    13:23,
                                    75,67,74,70,73,71,72,68,69,76:77,
                                    130, 114,126,116,120,123, 115, 119,117,121,124, 129,113,127,118,122,125,128, 131:132,
                                    42,43,41, 39,35,26,29,32, 25, 40,38,27,30,33, 44,24,36,28,31,34,37, 45:46,
                                    64, 48,60,50,54,57, 49, 53,51,55,58, 63,47,61,52,56,59,62, 65:66,
                                    95, 79,91,81,85,88, 80, 84,82,86,89, 94,78,92,83,87,90,93, 96:97,
                                    150, 134,146,136,140,143, 135, 139,137,141,144, 149,133,147,138,142,145,148, 151:152,
                                    160, 154,167,156,161,164, 155, 159,157,162,165, 170,153,168,158,163,166,169, 171:172,
                                    173, 175,187,177,181,184, 176, 180,178,182,185, 190,174,188,179,183,186,189, 191:192,
                                    104:112,
                                    193:202))
circ_PCP <- circ_PCP[order(match(rownames(circ_PCP), custom_order_rows)), ]
#assign ID
circ_PCP$id <- seq(1, nrow(circ_PCP))
#get name any y position of each label
label_data <- as.data.frame(circ_PCP)    # copy results to label
label_data <- label_data %>% mutate(id = row_number())    # id as row number
number_of_bar <- nrow(label_data)   # count number of bars to label
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # calculate angle of each row so it can be plotted on a circle (360 degrees)
#subtract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)    # horizontal location of the label, relative to the bar
label_data$angle <- ifelse(angle < -90, angle+180, angle)   # flip angle of label by 180 to make them readable
label_data$color <- ifelse(label_data$p_adj < 0.01, "black", "darkgrey")   # label will show significance, if not significant, biomarker will appear grey
base_data <- circ_PCP %>% 
  group_by(Class) %>% 
  dplyr::summarize(start=min(id) - 0.5, end=max(id) - 2 + 0.5) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))   # group biomarkers according to their categories, excluding empty bars. 
base_data$title[8] <- base_data$title[8]-5 #subtract half of the additional emtpy rows because its the mean
title_angle <- base_data$title*360/number_of_bar    # calculate the angle of label for each category so it changes angle as it goes around a circle
base_data$title_angle <- ifelse(90<title_angle & title_angle<270, title_angle - 180, title_angle)   #flip angle by 180 degrees if title angle is between 90 and 270 degrees
#remove final 10 rows
base_data$end[8] <- base_data$end[8]-10
circ_p_PCP <- ggplot(circ_PCP, aes(x=as.factor(id), y=hr-1, fill=hr-1)) +
  geom_bar(aes(x=as.factor(id), y=hr-1, fill=hr-1), stat="identity", alpha=0.8, na.rm = FALSE) +    # alpha changes transparency
  geom_segment(data=base_data, aes(x = start, y = 0.2, xend = end, yend = 0.2),   
               colour = "darkgrey", alpha=0.8, size=0.15, inherit.aes = FALSE) +
  geom_segment(data=base_data, aes(x = start, y = 0.1, xend = end, yend = 0.1),   
               colour = "darkgrey", alpha=0.8, size=0.15, inherit.aes = FALSE) +
  geom_segment(data=base_data, aes(x = start, y = 0, xend = end, yend = 0),
               colour = "black", alpha=0.8, size=0.3, inherit.aes = FALSE)  +
  geom_segment(data=base_data, aes(x = start, y = -0.1, xend = end, yend = -0.1),   
               colour = "darkgrey", alpha=0.8, size=0.15, inherit.aes = FALSE) +
  geom_segment(data=base_data, aes(x = start, y = -0.2, xend = end, yend = -0.2),   
               colour = "darkgrey", alpha=0.8, size=0.15, inherit.aes = FALSE) +
  annotate("text", x = rep(max(circ_PCP$id),5), y = c(-0.2, -0.1, 0, 0.1, 0.2), 
           label = c("0.8", "0.9", "1", "1.1", "1.2") , 
           color="darkgrey", size=1.5 , angle=0, fontface="bold", hjust=1) +
  scale_y_continuous(labels = function(x) x + 1, limits = c(-0.5,1)) +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0, limits=c(-0.35, 0.35)) +
  geom_errorbar(aes(ymin = ci_low-1, ymax = ci_high-1), 
                width = 0, color = "black", alpha=0.3, na.rm = TRUE, show.legend = FALSE, size=0.25) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=0.4, label=Name, hjust=hjust),
            color=label_data$color, alpha=1, size=1.5, 
            angle= label_data$angle, inherit.aes = FALSE) +
  geom_segment(data=base_data, aes(x = start, y = 0.75, xend = end, yend = 0.75),    # add grouping bracket
               colour = "black", alpha=0.8, size=0.6, inherit.aes = FALSE)  +
  geom_text(data=base_data, aes(x = title, y = 0.85, label=Class, hjust=0.5),         # Add category labels, to the centre (0.5) of each bracket
            colour = "black", alpha=0.8, size=1.8,         
            angle = -base_data$title_angle,              # Category title is angled, wraps around the graph
            fontface="bold", inherit.aes = FALSE)
#print plots
circ_p_agesex
circ_p_PCP

grid.arrange(circ_p_agesex, circ_p_PCP, nrow = 1)


#table 1: Participant characteristics of training and test datasets and total
#generate subsetted df
BL_characteristics <- Complete[,c("Age.at.recruitment","Sex","Ethnic.background...Instance.0", "college", 
                                  "Smoking.status...Instance.0", "Alcohol.drinker.status...Instance.0", "Summed.minutes.activity...Instance.0",
                                  "Cholesterol...Instance.0", "HDL.cholesterol...Instance.0" , "LDL.direct...Instance.0", "Triglycerides...Instance.0", "Glucose...Instance.0" , "Glycated.haemoglobin..HbA1c....Instance.0", "Creatinine...Instance.0", "C.reactive.protein...Instance.0",
                                  "SBPmin", "Body.mass.index..BMI....Instance.0", "WHR",
                                  "DM_Meds", "BP_Meds", "Familyhistory_heartdisease",
                                  "HF", "CAD", "HTN")]

BL_characteristics_PCP <- BL_characteristics[rownames(PCP),]
BL_characteristics_Train <- BL_characteristics[rownames(trainData),]
BL_characteristics_Test <- BL_characteristics[rownames(testData),]
BL_characteristics_HF <- BL_characteristics_PCP[BL_characteristics_PCP$HF==TRUE,]
BL_characteristics_noHF <- BL_characteristics_PCP[BL_characteristics_PCP$HF==FALSE,]

continuous_vars <- c("Age.at.recruitment", "Summed.minutes.activity...Instance.0", "Cholesterol...Instance.0",
                     "HDL.cholesterol...Instance.0", "LDL.direct...Instance.0", "Triglycerides...Instance.0",
                     "Glucose...Instance.0", "Glycated.haemoglobin..HbA1c....Instance.0", "Creatinine...Instance.0",
                     "C.reactive.protein...Instance.0", "SBPmin", "Body.mass.index..BMI....Instance.0", "WHR")
categorical_vars <- c("Sex", "Ethnic.background...Instance.0", "college", "Smoking.status...Instance.0",
                      "Alcohol.drinker.status...Instance.0", "DM_Meds", "BP_Meds", "Familyhistory_heartdisease",
                      "HF", "CAD", "HTN")

# Create an empty dataframe
result_df <- data.frame(Group = character(), Variable = character(),
                        BL_characteristics_PCP = numeric(),
                        BL_characteristics_Train = numeric(),
                        BL_characteristics_Test = numeric(),
                        BL_characteristics_HF = numeric(),
                        BL_characteristics_noHF = numeric(),
                        stringsAsFactors = FALSE)

# Continuous variables
result_df <- data.frame() 

for (var in continuous_vars) {
  median_val_PCP <- median(BL_characteristics_PCP[[var]], na.rm = TRUE)
  median_val_Train <- median(BL_characteristics_Train[[var]], na.rm = TRUE)
  median_val_Test <- median(BL_characteristics_Test[[var]], na.rm = TRUE)
  median_val_HF <- median(BL_characteristics_HF[[var]], na.rm = TRUE)
  median_val_noHF <- median(BL_characteristics_noHF[[var]], na.rm = TRUE)
  
  lower_q_PCP <- quantile(BL_characteristics_PCP[[var]], 0.25, na.rm = TRUE)
  upper_q_PCP <- quantile(BL_characteristics_PCP[[var]], 0.75, na.rm = TRUE)
  lower_q_Train <- quantile(BL_characteristics_Train[[var]], 0.25, na.rm = TRUE)
  upper_q_Train <- quantile(BL_characteristics_Train[[var]], 0.75, na.rm = TRUE)
  lower_q_Test <- quantile(BL_characteristics_Test[[var]], 0.25, na.rm = TRUE)
  upper_q_Test <- quantile(BL_characteristics_Test[[var]], 0.75, na.rm = TRUE)
  lower_q_HF <- quantile(BL_characteristics_HF[[var]], 0.25, na.rm = TRUE)
  upper_q_HF <- quantile(BL_characteristics_HF[[var]], 0.75, na.rm = TRUE)
  lower_q_noHF <- quantile(BL_characteristics_noHF[[var]], 0.25, na.rm = TRUE)
  upper_q_noHF <- quantile(BL_characteristics_noHF[[var]], 0.75, na.rm = TRUE)
  
  result_df <- rbind(result_df, data.frame(Group = "Continuous", Variable = var,
                                           BL_characteristics_PCP = paste(median_val_PCP, "(", lower_q_PCP, "-", upper_q_PCP, ")", sep = ""),
                                           BL_characteristics_Train = paste(median_val_Train, "(", lower_q_Train, "-", upper_q_Train, ")", sep = ""),
                                           BL_characteristics_Test = paste(median_val_Test, "(", lower_q_Test, "-", upper_q_Test, ")", sep = ""),
                                           BL_characteristics_HF = paste(median_val_HF, "(", lower_q_HF, "-", upper_q_HF, ")", sep = ""),
                                           BL_characteristics_noHF = paste(median_val_noHF, "(", lower_q_noHF, "-", upper_q_noHF, ")", sep = ""),
                                           stringsAsFactors = FALSE))
}

# Categorical variables
for (var in categorical_vars) {
  if (var == "Sex") {
    for (value in c("F", "M")) {
      count_PCP <- sum(BL_characteristics_PCP[[var]] == value, na.rm = TRUE)
      count_Train <- sum(BL_characteristics_Train[[var]] == value, na.rm = TRUE)
      count_Test <- sum(BL_characteristics_Test[[var]] == value, na.rm = TRUE)
      count_HF <- sum(BL_characteristics_HF[[var]] == value, na.rm = TRUE)
      count_noHF <- sum(BL_characteristics_noHF[[var]] == value, na.rm = TRUE)
      percentage_PCP <- count_PCP / length(BL_characteristics_PCP[[var]]) * 100
      percentage_Train <- count_Train / length(BL_characteristics_Train[[var]]) * 100
      percentage_Test <- count_Test / length(BL_characteristics_Test[[var]]) * 100
      percentage_HF <- count_HF / length(BL_characteristics_HF[[var]]) * 100
      percentage_noHF <- count_noHF / length(BL_characteristics_noHF[[var]]) * 100
      result_df <- rbind(result_df, data.frame(Group = "Categorical", Variable = paste(var, value, sep = "_"),
                                               BL_characteristics_PCP = paste(count_PCP, "(", percentage_PCP, "%)", sep = ""),
                                               BL_characteristics_Train = paste(count_Train, "(", percentage_Train, "%)", sep = ""),
                                               BL_characteristics_Test = paste(count_Test, "(", percentage_Test, "%)", sep = ""),
                                               BL_characteristics_HF = paste(count_HF, "(", percentage_HF, "%)", sep = ""),
                                               BL_characteristics_noHF = paste(count_noHF, "(", percentage_noHF, "%)", sep = ""),
                                               stringsAsFactors = FALSE))
    }
  } else if (var == "Ethnic.background...Instance.0") {
    for (value in c("B", "W")) {
      count_PCP <- sum(BL_characteristics_PCP[[var]] == value, na.rm = TRUE)
      count_Train <- sum(BL_characteristics_Train[[var]] == value, na.rm = TRUE)
      count_Test <- sum(BL_characteristics_Test[[var]] == value, na.rm = TRUE)
      count_HF <- sum(BL_characteristics_HF[[var]] == value, na.rm = TRUE)
      count_noHF <- sum(BL_characteristics_noHF[[var]] == value, na.rm = TRUE)
      percentage_PCP <- count_PCP / length(BL_characteristics_PCP[[var]]) * 100
      percentage_Train <- count_Train / length(BL_characteristics_Train[[var]]) * 100
      percentage_Test <- count_Test / length(BL_characteristics_Test[[var]]) * 100
      percentage_HF <- count_HF / length(BL_characteristics_HF[[var]]) * 100
      percentage_noHF <- count_noHF / length(BL_characteristics_noHF[[var]]) * 100
      result_df <- rbind(result_df, data.frame(Group = "Categorical", Variable = paste(var, value, sep = "_"),
                                               BL_characteristics_PCP = paste(count_PCP, "(", percentage_PCP, "%)", sep = ""),
                                               BL_characteristics_Train = paste(count_Train, "(", percentage_Train, "%)", sep = ""),
                                               BL_characteristics_Test = paste(count_Test, "(", percentage_Test, "%)", sep = ""),
                                               BL_characteristics_HF = paste(count_HF, "(", percentage_HF, "%)", sep = ""),
                                               BL_characteristics_noHF = paste(count_noHF, "(", percentage_noHF, "%)", sep = ""),
                                               stringsAsFactors = FALSE))
    }
  } else {
    count_PCP <- sum(BL_characteristics_PCP[[var]] == TRUE, na.rm = TRUE)
    count_Train <- sum(BL_characteristics_Train[[var]] == TRUE, na.rm = TRUE)
    count_Test <- sum(BL_characteristics_Test[[var]] == TRUE, na.rm = TRUE)
    count_HF <- sum(BL_characteristics_HF[[var]] == TRUE, na.rm = TRUE)
    count_noHF <- sum(BL_characteristics_noHF[[var]] == TRUE, na.rm = TRUE)
    percentage_PCP <- count_PCP / length(BL_characteristics_PCP[[var]]) * 100
    percentage_Train <- count_Train / length(BL_characteristics_Train[[var]]) * 100
    percentage_Test <- count_Test / length(BL_characteristics_Test[[var]]) * 100
    percentage_HF <- count_HF / length(BL_characteristics_HF[[var]]) * 100
    percentage_noHF <- count_noHF / length(BL_characteristics_noHF[[var]]) * 100
    result_df <- rbind(result_df, data.frame(Group = "Categorical", Variable = var,
                                             BL_characteristics_PCP = paste(count_PCP, "(", percentage_PCP, "%)", sep = ""),
                                             BL_characteristics_Train = paste(count_Train, "(", percentage_Train, "%)", sep = ""),
                                             BL_characteristics_Test = paste(count_Test, "(", percentage_Test, "%)", sep = ""),
                                             BL_characteristics_HF = paste(count_HF, "(", percentage_HF, "%)", sep = ""),
                                             BL_characteristics_noHF = paste(count_noHF, "(", percentage_noHF, "%)", sep = ""),
                                             stringsAsFactors = FALSE))
  }
}

result_df
write.table(result_df, file="/Users/rafael/Library/CloudStorage/OneDrive-King'sCollegeLondon/KCL/KCL-PhD/Writing/UKB Serum Metabolomics to predict incident HF/Supplementary tables/BLcharacteristics.tsv", quote=FALSE, sep='\t')

#stats
pvalue_df <- data.frame(Variable = character(),
                        HF_vs_noHF_pvalue = numeric(),
                        Test_vs_Train_pvalue = numeric(),
                        stringsAsFactors = FALSE)

#continuous variables - Mann-Whitney U
for (var in continuous_vars) {
  wilcox_test_HF_vs_noHF <- wilcox.test(BL_characteristics_HF[[var]], BL_characteristics_noHF[[var]])
  pvalue_HF_vs_noHF <- wilcox_test_HF_vs_noHF$p.value
  
  wilcox_test_Test_vs_Train <- wilcox.test(BL_characteristics_Test[[var]], BL_characteristics_Train[[var]])
  pvalue_Test_vs_Train <- wilcox_test_Test_vs_Train$p.value
  
  pvalue_df <- rbind(pvalue_df, data.frame(Variable = var,
                                           HF_vs_noHF_pvalue = pvalue_HF_vs_noHF,
                                           Test_vs_Train_pvalue = pvalue_Test_vs_Train,
                                           stringsAsFactors = FALSE))
}

#categorical variables - X2
for (var in categorical_vars) {
  df1 <- data.frame("1" = BL_characteristics_HF[[var]])
  df2 <- data.frame("2" = BL_characteristics_noHF[[var]])
  
  m <- cbind(table(df1), table(df2))
  
  chisq_test_HF_vs_noHF <- chisq.test(m)
  pvalue_HF_vs_noHF <- chisq_test_HF_vs_noHF$p.value
  
  df1 <- data.frame("1" = BL_characteristics_train[[var]])
  df2 <- data.frame("2" = BL_characteristics_Train[[var]])
  
  m <- cbind(table(df1), table(df2))
  
  chisq_test_Test_vs_Train <- chisq.test(m)
  pvalue_Test_vs_Train <- chisq_test_Test_vs_Train$p.value
  
  pvalue_df <- rbind(pvalue_df, data.frame(Variable = var,
                                           HF_vs_noHF_pvalue = pvalue_HF_vs_noHF[1],
                                           Test_vs_Train_pvalue = pvalue_Test_vs_Train[1],
                                           stringsAsFactors = FALSE))
}

#adjustment for multiple testing
pvalue_df$HF_vs_noHF_pvalue_adj <- p.adjust(pvalue_df$HF_vs_noHF_pvalue, method = "BH")
pvalue_df$Test_vs_Train_pvalue_adj <- p.adjust(pvalue_df$Test_vs_Train_pvalue, method = "BH")
pvalue_df
write.table(pvalue_df, file="/Users/rafael/Library/CloudStorage/OneDrive-King'sCollegeLondon/KCL/KCL-PhD/Writing/UKB Serum Metabolomics to predict incident HF/Supplementary tables/Pvalues_BLcharacteristics.tsv", quote=FALSE, sep='\t')



#model predictions
testData$LP_agesex <- predict(fit_agesex$model, newx = as.matrix(testData[,c(169:170)]), type = "link")
testData$LP_agesex_met <- predict(fit_agesex_met$model, newx = as.matrix(testData[,c(1:170)]), type = "link")
testData$LP_PCP <- predict(fit_PCP$model, newx = as.matrix(testData[,c(169:170,172,173,195,208,209,210,212,335,336,337)]), type = "link")
testData$LP_PCP_met <- predict(fit_PCP_met$model, newx = as.matrix(testData[,c(1:170,172,173,195,208,209,210,212,335,336,337)]), type = "link")

trainData$LP_agesex <- predict(fit_agesex$model, newx = as.matrix(trainData[,c(169:170)]), type = "link")
trainData$LP_agesex_met <- predict(fit_agesex_met$model, newx = as.matrix(trainData[,c(1:170)]), type = "link")
trainData$LP_PCP <- predict(fit_PCP$model, newx = as.matrix(trainData[,c(169:170,172,173,195,208,209,210,212,335,336,337)]), type = "link")
trainData$LP_PCP_met <- predict(fit_PCP_met$model, newx = as.matrix(trainData[,c(1:170,172,173,195,208,209,210,212,335,336,337)]), type = "link")

#figure 2: ROC curve
roc_agesex <- roc(testData$HF, testData$LP_agesex[,1])
roc_agesex_df <- data.frame(Sensitivity = roc_agesex$sensitivities, Specificity = 1 - roc_agesex$specificities)
roc_agesex_met <- roc(testData$HF, testData$LP_agesex_met[,1])
roc_agesex_met_df <- data.frame(Sensitivity = roc_agesex_met$sensitivities, Specificity = 1 - roc_agesex_met$specificities)
roc_PCP <- roc(testData$HF, testData$LP_PCP[,1])
roc_PCP_df <- data.frame(Sensitivity = roc_PCP$sensitivities, Specificity = 1 - roc_PCP$specificities)
roc_PCP_met <- roc(testData$HF, testData$LP_PCP_met[,1])
roc_PCP_met_df <- data.frame(Sensitivity = roc_PCP_met$sensitivities, Specificity = 1 - roc_PCP_met$specificities)

roc_agesex_train <- roc(trainData$HF, trainData$LP_agesex[,1])
roc_agesex_df_train <- data.frame(Sensitivity = roc_agesex_train$sensitivities, Specificity = 1 - roc_agesex_train$specificities)
roc_agesex_met_train <- roc(trainData$HF, trainData$LP_agesex_met[,1])
roc_agesex_met_df_train <- data.frame(Sensitivity = roc_agesex_met_train$sensitivities, Specificity = 1 - roc_agesex_met_train$specificities)
roc_PCP_train <- roc(trainData$HF, trainData$LP_PCP[,1])
roc_PCP_df_train <- data.frame(Sensitivity = roc_PCP_train$sensitivities, Specificity = 1 - roc_PCP_train$specificities)
roc_PCP_met_train <- roc(trainData$HF, trainData$LP_PCP_met[,1])
roc_PCP_met_df_train <- data.frame(Sensitivity = roc_PCP_met_train$sensitivities, Specificity = 1 - roc_PCP_met_train$specificities)

smooth_agesex <- smooth.spline(1 - roc_agesex$specificities, roc_agesex$sensitivities, spar = 0.5)
smooth_agesex_met <- smooth.spline(1 - roc_agesex_met$specificities, roc_agesex_met$sensitivities, spar = 0.5)
smooth_PCP <- smooth.spline(1 - roc_PCP$specificities, roc_PCP$sensitivities, spar = 0.5)
smooth_PCP_met <- smooth.spline(1 - roc_PCP_met$specificities, roc_PCP_met$sensitivities, spar = 0.5)
smooth_agesex_df <- data.frame(Specificity = smooth_agesex$x, Sensitivity = smooth_agesex$y)
smooth_agesex_met_df <- data.frame(Specificity = smooth_agesex_met$x, Sensitivity = smooth_agesex_met$y)
smooth_PCP_df <- data.frame(Specificity = smooth_PCP$x, Sensitivity = smooth_PCP$y)
smooth_PCP_met_df <- data.frame(Specificity = smooth_PCP_met$x, Sensitivity = smooth_PCP_met$y)

plot_ROC <- ggplot() +
  geom_line(data = roc_agesex_df, aes(x = Specificity, y = Sensitivity), color = "#9bd7d5", size = 0.5) +
  geom_line(data = roc_agesex_met_df, aes(x = Specificity, y = Sensitivity), color = "#00999a", size = 0.5) +
  geom_line(data = roc_PCP_df, aes(x = Specificity, y = Sensitivity), color = "#ffc79a", size = 0.5) +
  geom_line(data = roc_PCP_met_df, aes(x = Specificity, y = Sensitivity), color = "#ff7400", size = 0.5) +
  geom_line(data = smooth_agesex_df, aes(x = Specificity, y = Sensitivity), color = "#9bd7d5", size = 0.8) +
  geom_line(data = smooth_agesex_met_df, aes(x = Specificity, y = Sensitivity), color = "#00999a", size = 0.8) +
  geom_line(data = smooth_PCP_df, aes(x = Specificity, y = Sensitivity), color = "#ffc79a", size = 0.8) +
  geom_line(data = smooth_PCP_met_df, aes(x = Specificity, y = Sensitivity), color = "#ff7400", size = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkgrey", size = 0.8) +
  labs(x = "1 - Specificity", y = "Sensitivity", title = "Receiver Operating Characteristic Curve") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5))
plot_ROC

#DCA plot
df_dca <- Calibration
df_dca$risk_agesex <- 1-df_dca$risk_agesex
df_dca$risk_agesex_met <- 1-df_dca$risk_agesex_met
df_dca$risk_PCP <- 1-df_dca$risk_PCP
df_dca$risk_PCP_met <- 1-df_dca$risk_PCP_met

dca <- dca(Surv(followup_HF, HF) ~ risk_agesex + risk_agesex_met + risk_PCP + risk_PCP_met, 
    data = df_dca,
    time = 10,
    thresholds = 1:50 / 1000)
#dca package implementation doesn't work in those cases where TP is below min risk, where it sets NB to 0 whilst it should be equal to the corresponding "Treat All" value (all cases are TP), manually fixed here)
dca$dca$net_benefit[dca$dca$variable=="risk_agesex" & (dca$dca$threshold < min(df_dca$risk_agesex))] <- dca$dca$net_benefit[dca$dca$variable=="all" & (dca$dca$threshold < min(df_dca$risk_agesex))]
dca$dca$net_benefit[dca$dca$variable=="risk_agesex_met" & (dca$dca$threshold < min(df_dca$risk_agesex_met))] <- dca$dca$net_benefit[dca$dca$variable=="all" & (dca$dca$threshold < min(df_dca$risk_agesex_met))]
dca$dca$net_benefit[dca$dca$variable=="risk_PCP" & (dca$dca$threshold < min(df_dca$risk_PCP))] <- dca$dca$net_benefit[dca$dca$variable=="all" & (dca$dca$threshold < min(df_dca$risk_PCP))]
dca$dca$net_benefit[dca$dca$variable=="risk_PCP_met" & (dca$dca$threshold < min(df_dca$risk_PCP_met))] <- dca$dca$net_benefit[dca$dca$variable=="all" & (dca$dca$threshold < min(df_dca$risk_PCP_met))]

plot_dca <- plot(dca, smooth=TRUE, lwd = 0.1) +
  scale_color_manual(values = c("darkred", "darkgrey", "#9bd7d5", "#00999a", "#ffc79a", "#ff7400")) +
  geom_line (data = dca$dca, aes(x = threshold, y = net_benefit), lwd = 0.5) +
  labs(title = "Decicion Curve Analysis Curve") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5), legend.position = "none")

#plot
grid.arrange(plot_ROC, plot_dca, nrow = 1)
#legends
legend_dca <- data.frame(variable = c("A", "B", "C", "D", "E", "F"),
                             label = c("Treat All", "Treat None","Age & Sex", "Age & Sex + Metabolomics", "PCP-HF", "PCP-HF + Metabolomics"),
                             color = c("darkred", "darkgrey", "#9bd7d5", "#00999a", "#ffc79a", "#ff7400"))
legend_dca_plot <- ggplot(legend_dca, aes(x = factor(1), y = factor(1), color = variable)) +
  geom_line(aes(group = variable), size = 2) +  # Use geom_line() instead of geom_point() for lines
  guides(color = guide_legend(title = "Legend")) +
  scale_color_manual(values = legend_dca$color, labels = legend_dca$label) +
  theme_void() +
  theme(legend.position = c(0.5,0.5),
        legend.justification = "center", 
        legend.box.just = "center")
print(legend_dca_plot)

legend_forrest <- data.frame(variable = c("corresponding_agesex", "top20_PCP"),
                             label = c("Adjusted for Age + Sex", "Adjusted for PCP-HF"),
                             color = c("#00999a", "#ff7400"))
legend_plot <- ggplot(legend_df, aes(x = factor(1), y = factor(1), color = variable)) +
  geom_point(size = 5) +  # Adjust the size of the dots here
  guides(color = guide_legend(title = "Legend")) +
  scale_color_manual(values = legend_df$color, labels = legend_df$label) +
  theme_void() +
  theme(legend.position = "bottom", legend.justification = "center")
print(legend_plot)



#table 2: Discriminate performance
#C
model_C_agesex <- coxph(Surv(testData$followup_HF, testData$HF) ~ testData$LP_PCP_met)
summary(model_C_agesex) #alternatively getting concordance/SE

C_agesex <- Cindex(testData$LP_agesex, Surv(testData$followup_HF, testData$HF))
C_agesex_met <- Cindex(testData$LP_agesex_met, Surv(testData$followup_HF, testData$HF))
C_PCP <- Cindex(testData$LP_PCP, Surv(testData$followup_HF, testData$HF))
C_PCP_met <- Cindex(testData$LP_PCP_met, Surv(testData$followup_HF, testData$HF))
C_agesex
C_agesex_met
C_PCP
C_PCP_met

#deltaC
C_agesex_met-C_agesex
C_PCP_met-C_PCP
C_PCP-C_agesex_met

#training split
C_agesex_train <- Cindex(trainData$LP_agesex, Surv(trainData$followup_HF, trainData$HF))
C_agesex_met_train <- Cindex(trainData$LP_agesex_met, Surv(trainData$followup_HF, trainData$HF))
C_PCP_train <- Cindex(trainData$LP_PCP, Surv(trainData$followup_HF, trainData$HF))
C_PCP_met_train <- Cindex(trainData$LP_PCP_met, Surv(trainData$followup_HF, trainData$HF))
C_agesex_train
C_agesex_met_train
C_PCP_train
C_PCP_met_train

#bootstrap C
bootstrap_Cindex <- function(data, indices, lp_column) {
  d <- data[indices,]
  return(Cindex(d[[lp_column]], Surv(d$followup_HF, d$HF)))
}

set.seed(123)
boot_results_agesex <- boot(data = testData, statistic = function(data, indices) bootstrap_Cindex(data, indices, 'LP_agesex'), R = 1000)
boot_results_agesex
boot_ci_agesex <- boot.ci(boot_results_agesex, type = "perc") 
boot_ci_agesex

#Sens and Spec at Jouden
youden_agesex <- coords(roc_agesex, "best", best.method="youden")
youden_agesex_met <- coords(roc_agesex_met, "best", best.method="youden")
youden_PCP <- coords(roc_PCP, "best", best.method="youden")
youden_PCP_met <- coords(roc_PCP_met, "best", best.method="youden")
youden_agesex
youden_agesex_met
youden_PCP
youden_PCP_met

youden_agesex_train <- coords(roc_agesex_train, "best", best.method="youden")
youden_agesex_met_train <- coords(roc_agesex_met_train, "best", best.method="youden")
youden_PCP_train <- coords(roc_PCP_train, "best", best.method="youden")
youden_PCP_met_train <- coords(roc_PCP_met_train, "best", best.method="youden")
youden_agesex_train
youden_agesex_met_train
youden_PCP_train
youden_PCP_met_train

#NRI
#NRI = P(up∣event)−P(down∣event) + P(down∣nonevent)−P(up∣nonevent) 
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3341973/#FD1
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

NRI_agesex <- calculate_nri(testData, "LP_agesex", "LP_agesex_met", "HF")
NRI_nomet <- calculate_nri(testData, "LP_agesex", "LP_PCP", "HF")
NRI_agesexvspcpmet <- calculate_nri(testData, "LP_agesex", "LP_PCP_met", "HF")
NRI_agesexmetvspcp <- calculate_nri(testData, "LP_agesex_met", "LP_PCP", "HF")
NRI_met <- calculate_nri(testData, "LP_agesex_met", "LP_PCP_met", "HF")
NRI_PCP <- calculate_nri(testData, "LP_PCP", "LP_PCP_met", "HF")
NRI_agesex
NRI_nomet
NRI_agesexvspcpmet
NRI_agesexmetvspcp
NRI_met
NRI_PCP

NRI_agesex_train <- calculate_nri(trainData, "LP_agesex", "LP_agesex_met", "HF")
NRI_nomet_train <- calculate_nri(trainData, "LP_agesex", "LP_PCP", "HF")
NRI_agesexvspcpmet_train <- calculate_nri(trainData, "LP_agesex", "LP_PCP_met", "HF")
NRI_agesexmetvspcp_train <- calculate_nri(trainData, "LP_agesex_met", "LP_PCP", "HF")
NRI_met_train <- calculate_nri(trainData, "LP_agesex_met", "LP_PCP_met", "HF")
NRI_PCP_train <- calculate_nri(trainData, "LP_PCP", "LP_PCP_met", "HF")
NRI_agesex_train
NRI_nomet_train
NRI_agesexvspcpmet_train
NRI_agesexmetvspcp_train
NRI_met_train
NRI_PCP_train

#IDI
#IDI = (ave pcases − ave pcontrols)new model − (ave pcases − ave pcontrols)old model
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3341978/
#rel IDI divides by Yates slope in the reference model
#calculated on the calibration dataset with exclusion of pre-10y non-HF censoring and predicted risk at 10y

#generate calibration datasets
#calibration
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

#calculate IDI
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

IDI_agesex <- calculate_idi(Calibration, "risk_agesex", "risk_agesex_met", "HF10y")
IDI_nomet <- calculate_idi(Calibration, "risk_agesex", "risk_PCP", "HF10y")
IDI_agesexvspcpmet <- calculate_idi(Calibration, "risk_agesex", "risk_PCP_met", "HF10y")
IDI_agesexmetvspcp <- calculate_idi(Calibration, "risk_agesex_met", "risk_PCP", "HF10y")
IDI_met <- calculate_idi(Calibration, "risk_agesex_met", "risk_PCP_met", "HF10y")
IDI_PCP <- calculate_idi(Calibration, "risk_PCP", "risk_PCP_met", "HF10y")
IDI_agesex
IDI_nomet
IDI_agesexvspcpmet
IDI_agesexmetvspcp
IDI_met
IDI_PCP

IDI_agesex_train <- calculate_idi(Calibration_train, "risk_agesex", "risk_agesex_met", "HF10y")
IDI_nomet_train <- calculate_idi(Calibration_train, "risk_agesex", "risk_PCP", "HF10y")
IDI_agesexvspcpmet_train <- calculate_idi(Calibration_train, "risk_agesex", "risk_PCP_met", "HF10y")
IDI_agesexmetvspcp_train <- calculate_idi(Calibration_train, "risk_agesex_met", "risk_PCP", "HF10y")
IDI_met_train <- calculate_idi(Calibration_train, "risk_agesex_met", "risk_PCP_met", "HF10y")
IDI_PCP_train <- calculate_idi(Calibration_train, "risk_PCP", "risk_PCP_met", "HF10y")
IDI_agesex_train
IDI_nomet_train
IDI_agesexvspcpmet_train
IDI_agesexmetvspcp_train
IDI_met_train
IDI_PCP_train

#figure 3: calibration
Calibration$riskdeciles_agesex <- ntile(Calibration$risk_agesex, 5)
Calibration$riskdeciles_agesex_met <- ntile(Calibration$risk_agesex_met, 5)
Calibration$riskdeciles_PCP <- ntile(Calibration$risk_PCP, 5)
Calibration$riskdeciles_PCP_met <- ntile(Calibration$risk_PCP_met, 5)

meanobserved_agesex <- aggregate(HF10y ~ riskdeciles_agesex, data = Calibration, FUN = mean)
meanobserved_agesex_met <- aggregate(HF10y ~ riskdeciles_agesex_met, data = Calibration, FUN = mean)
meanobserved_PCP <- aggregate(HF10y ~ riskdeciles_PCP, data = Calibration, FUN = mean)
meanobserved_PCP_met <- aggregate(HF10y ~ riskdeciles_PCP_met, data = Calibration, FUN = mean)

meanpredicted_agesex <- aggregate(risk_agesex ~ riskdeciles_agesex, data = Calibration, FUN = mean)
meanpredicted_agesex_met <- aggregate(risk_agesex_met ~ riskdeciles_agesex_met, data = Calibration, FUN = mean)
meanpredicted_PCP <- aggregate(risk_PCP ~ riskdeciles_PCP, data = Calibration, FUN = mean)
meanpredicted_PCP_met <- aggregate(risk_PCP_met ~ riskdeciles_PCP_met, data = Calibration, FUN = mean)

plot_data_agesex <- merge(meanobserved_agesex, meanpredicted_agesex, by = "riskdeciles_agesex")
plot_data_agesex$'10'<- (1 - plot_data_agesex$`10`) * 100
plot_data_agesex_met <- merge(meanobserved_agesex_met, meanpredicted_agesex_met, by = "riskdeciles_agesex_met")
plot_data_agesex_met$'10'<- (1 - plot_data_agesex_met$`10`) * 100
plot_data_PCP <- merge(meanobserved_PCP, meanpredicted_PCP, by = "riskdeciles_PCP")
plot_data_PCP$'10'<- (1 - plot_data_PCP$`10`) * 100
plot_data_PCP_met <- merge(meanobserved_PCP_met, meanpredicted_PCP_met, by = "riskdeciles_PCP_met")
plot_data_PCP_met$'10'<- (1 - plot_data_PCP_met$`10`) * 100

cal_agesex <- ggplot(plot_data_agesex, aes(x = `10`, y = HF10y * 100, color = factor(riskdeciles_agesex))) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkgrey") +
  stat_smooth(method = "lm", fullrange = TRUE, color = "darkred", linewidth = 0.5, alpha = 0.075) +
  geom_point() +
  labs(x = "Predicted Risk (%)",
       y = "Observed Risk (%)") +
  scale_color_manual(values = rev(palette_risk_5)) +
  expand_limits(x = c(0, 5), y = c(0, 5)) +
  coord_cartesian(xlim = c(0,5), ylim = c(0,5)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5)) +
  guides(color = "none")
cal_agesex_met <- ggplot(plot_data_agesex_met, aes(x = `10`, y = HF10y * 100, color = factor(riskdeciles_agesex_met))) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
  stat_smooth(method = "lm", fullrange = TRUE, color = "darkred", linewidth = 0.5, alpha = 0.075) +
  geom_point() +
  labs(x = "Predicted Risk (%)",
       y = "Observed Risk (%)") +
  scale_color_manual(values = rev(palette_risk_5)) +
  expand_limits(x = c(0, 5), y = c(0, 5)) +
  coord_cartesian(xlim = c(0,5), ylim = c(0,5)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5)) +
  guides(color = "none")
cal_PCP <- ggplot(plot_data_PCP, aes(x = `10`, y = HF10y * 100, color = factor(riskdeciles_PCP))) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
  stat_smooth(method = "lm", fullrange = TRUE, color = "darkred", linewidth = 0.5, alpha = 0.075) +
  geom_point() +
  labs(x = "Predicted Risk (%)",
       y = "Observed Risk (%)") +
  scale_color_manual(values = rev(palette_risk_5)) +
  expand_limits(x = c(0, 5), y = c(0, 5)) +
  coord_cartesian(xlim = c(0,5), ylim = c(0,5)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5)) +
  guides(color = "none")
cal_PCP_met <- ggplot(plot_data_PCP_met, aes(x = `10`, y = HF10y * 100, color = factor(riskdeciles_PCP_met))) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
  stat_smooth(method = "lm", fullrange = TRUE, color = "darkred", linewidth = 0.5, alpha = 0.075) +
  geom_point() +
  labs(x = "Predicted Risk (%)",
       y = "Observed Risk (%)") +
  scale_color_manual(values = rev(palette_risk_5)) +
  expand_limits(x = c(0, 5), y = c(0, 5)) +
  coord_cartesian(xlim = c(0,5), ylim = c(0,5)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5)) +
  guides(color = "none")


legend_Calibration <- data.frame(variable = c("1", "2", "3", "4", "5"),
                             label = c("Risk quintile 1", "Risk quintile 2", "Risk quintile 3", "Risk quintile 4", "Risk quintile 5"),
                             color = c("#d1defa", "#759bf0", "#1858e7", "#0f358a", "#05122e"))
legend_Calibration_plot <- ggplot(legend_Calibration, aes(x = factor(1), y = factor(1), color = variable)) +
  geom_point(size = 5) +  # Adjust the size of the dots here
  guides(color = guide_legend(title = "Legend")) +
  scale_color_manual(values = legend_Calibration$color, labels = legend_Calibration$label) +
  theme_void() +
  theme(legend.position = "bottom", legend.justification = "center")
print(legend_Calibration_plot)

#figure 4: Incidence and survival curves
#incidence
testData$riskdeciles_LP_agesex <- ntile(testData$LP_agesex, 5)
testData$riskdeciles_LP_agesex_met <- ntile(testData$LP_agesex_met, 5)
testData$riskdeciles_LP_PCP <- ntile(testData$LP_PCP, 5)
testData$riskdeciles_LP_PCP_met <- ntile(testData$LP_PCP_met, 5)

meanobserved_LP_agesex <- aggregate(HF ~ riskdeciles_LP_agesex, data = testData, FUN = mean) * 100
meanobserved_LP_agesex_met <- aggregate(HF ~ riskdeciles_LP_agesex_met, data = testData, FUN = mean) * 100
meanobserved_LP_PCP <- aggregate(HF ~ riskdeciles_LP_PCP, data = testData, FUN = mean) * 100
meanobserved_LP_PCP_met <- aggregate(HF ~ riskdeciles_LP_PCP_met, data = testData, FUN = mean) * 100

meanobserved_LP_agesex
meanobserved_LP_agesex_met
meanobserved_LP_PCP
meanobserved_LP_PCP_met

incidence_agesex <- ggplot(meanobserved_LP_agesex, aes(x = factor(riskdeciles_LP_agesex), y = HF, color = factor(riskdeciles_LP_agesex))) +
  geom_point() +
  labs(x = "Risk quintile",
       y = "Cumulative Incidence (%)", 
       title = "Age + Sex") +
  scale_color_manual(values = palette_risk_5) +
  scale_x_discrete(labels = c("1", "2", "3", "4", "5")) +
  ylim(0, 10) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5)) +
  guides(color = "none")
incidence_agesex_met <- ggplot(meanobserved_LP_agesex_met, aes(x = factor(riskdeciles_LP_agesex_met), y = HF, color = factor(riskdeciles_LP_agesex_met))) +
  geom_point() +
  labs(x = "Risk quintile",
       y = "Cumulative Incidence (%)", 
       title = "Age + Sex + Metabolomics") +
  scale_color_manual(values = palette_risk_5) +
  scale_x_discrete(labels = c("1", "2", "3", "4", "5")) +
  ylim(0, 10) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5)) +
  guides(color = "none")
incidence_PCP <- ggplot(meanobserved_LP_PCP, aes(x = factor(riskdeciles_LP_PCP), y = HF, color = factor(riskdeciles_LP_PCP))) +
  geom_point() +
  labs(x = "Risk quintile",
       y = "Cumulative Incidence (%)", 
       title = "PCP-HF") +
  scale_color_manual(values = palette_risk_5) +
  scale_x_discrete(labels = c("1", "2", "3", "4", "5")) +
  ylim(0, 10) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5)) +
  guides(color = "none")
incidence_PCP_met <- ggplot(meanobserved_LP_PCP_met, aes(x = factor(riskdeciles_LP_PCP_met), y = HF, color = factor(riskdeciles_LP_PCP_met))) +
  geom_point() +
  labs(x = "Risk quintile",
       y = "Cumulative Incidence (%)", 
       title = "PCP-HF + Metabolomics") +
  scale_color_manual(values = palette_risk_5) +
  scale_x_discrete(labels = c("1", "2", "3", "4", "5")) +
  ylim(0, 10) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5)) +
  guides(color = "none")

#survival curves
survfit_agesex <- survfit(Surv(testData$followup_HF, testData$HF) ~ riskdeciles_LP_agesex, 
                    data = testData) 
survfit_agesex_met <- survfit(Surv(testData$followup_HF, testData$HF) ~ riskdeciles_LP_agesex_met, 
                          data = testData)
survfit_PCP <- survfit(Surv(testData$followup_HF, testData$HF) ~ riskdeciles_LP_PCP, 
                          data = testData)
survfit_PCP_met <- survfit(Surv(testData$followup_HF, testData$HF) ~ riskdeciles_LP_PCP_met, 
                          data = testData)

survplot_agesex <- ggsurvplot(survfit_agesex, ylim = c(0.9,1), censor.shape="", palette = palette_risk_5, legend = "none", ylab = "HF-free survival", xlab = "Time (years)")
survplot_agesex <- survplot_agesex$plot + theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
                                                axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))
survplot_agesex_met <- ggsurvplot(survfit_agesex_met, ylim = c(0.9,1), censor.shape="", palette = palette_risk_5, legend = "none", ylab = "HF-free survival", xlab = "Time (years)")
survplot_agesex_met <- survplot_agesex_met$plot + theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
                                                        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))
survplot_PCP <- ggsurvplot(survfit_PCP, ylim = c(0.9,1), censor.shape="", palette = palette_risk_5, legend = "none", ylab = "HF-free survival", xlab = "Time (years)")
survplot_PCP <- survplot_PCP$plot + theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
                                          axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))
survplot_PCP_met <- ggsurvplot(survfit_PCP_met, ylim = c(0.9,1), censor.shape="", palette = palette_risk_5, legend = "none", ylab = "HF-free survival", xlab = "Time (years)")
survplot_PCP_met <- survplot_PCP_met$plot + theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
                                                  axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))

grid.arrange(incidence_agesex,incidence_agesex_met,incidence_PCP,incidence_PCP_met, 
             cal_agesex, cal_agesex_met, cal_PCP, cal_PCP_met,
             survplot_agesex,survplot_agesex_met,survplot_PCP,survplot_PCP_met,
             nrow=3)


#figure 5:EN validation
coefs_agesex <- as.data.frame(as.matrix(coef(fit_agesex$model)))
coefs_agesex_met <- as.data.frame(as.matrix(coef(fit_agesex_met$model)))
coefs_PCP <- as.data.frame(as.matrix(coef(fit_PCP$model)))
coefs_PCP_met <- as.data.frame(as.matrix(coef(fit_PCP_met$model)))

coefs_agesex <- coefs_agesex[coefs_agesex$s0 != 0, , drop=FALSE]
coefs_agesex_met <- coefs_agesex_met[coefs_agesex_met$s0 != 0, , drop=FALSE]
coefs_PCP <- coefs_PCP[coefs_PCP$s0 != 0, , drop=FALSE]
coefs_PCP_met <- coefs_PCP_met[coefs_PCP_met$s0 != 0, , drop=FALSE]

met_cols <- rownames(coefs_PCP_met)[1:16]
clinical_cols <- rownames(coefs_PCP_met)[17:25]
corr_matrix <- cor(testData[, c(met_cols, clinical_cols)], method = "spearman")[met_cols, clinical_cols]
corrplot_PCP <- corrplot(corr_matrix, tl.cex = 0.8, tl.col = "black",  col.lim = c(-0.75, 0.75), is.corr=FALSE, addgrid.col = NA, cl.cex = 0.5)
#insert to have legend at the bottom , then stictch together manually type = "lower", 

#Network plot
Metabolites <- PCP[,c(1:168)]
head(Metabolites)

#soft treshold
powers = c(c(1:10), seq(from = 12, to=50, by=2))
allowWGCNAThreads() 
sft = pickSoftThreshold(Metabolites, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
softPower = 30
#adjacency matrix
adjacency = adjacency(Metabolites, power = softPower)
#create TOM object
TOM = TOMsimilarity(adjacency,
                    TOMType="unsigned")
#reintroduce labels
rownames(TOM) <- rownames(adjacency)
colnames(TOM) <- colnames(adjacency)
TOM[1:5,1:5]

#adjust which connections to keep
adj <- TOM
adj[adj > 0.2] = 1
adj[adj != 1] = 0
network <- graph.adjacency(adj)
self_edges <- which(is.loop(network))
network <- delete.edges(network, self_edges)

#color coding for metabolite groups and metabolites used in the final model
annotation2 <- circ_agesex[!is.na(circ_agesex$variable),]
head(annotation2)

annotation2$Class <- gsub("Lipoprotein particles", "Lipoprotein subclasses", annotation2$Class, fixed = TRUE)
annotation2$Class <- gsub("Cholesteryl Esters", "Lipoprotein subclasses", annotation2$Class, fixed = TRUE)
annotation2$Class <- gsub("Free Cholesterol", "Lipoprotein subclasses", annotation2$Class, fixed = TRUE)
annotation2$Class <- gsub("Cholesterol", "Lipoprotein subclasses", annotation2$Class, fixed = TRUE)
annotation2$Class <- gsub("Phospholipids", "Lipoprotein subclasses", annotation2$Class, fixed = TRUE)
annotation2$Class <- gsub("Total Lipids", "Lipoprotein subclasses", annotation2$Class, fixed = TRUE)
annotation2$Class <- gsub("Triglycerides", "Lipoprotein subclasses", annotation2$Class, fixed = TRUE)
annotation2$Class <- gsub("Apo-LP, LP sizes\n& other lipids", "Apo-LP, LP sizes & other lipids", annotation2$Class, fixed = TRUE)
annotation2$Class <- gsub("Ketone Bodies, Fluid balance\n& Inflammation", "Ketone Bodies, Fluid balance & Inflammation", annotation2$Class, fixed = TRUE)

groups <- unique(annotation2$Class)

color_palette <- c("#ffd502","#ff0100","#03ce03","#00999a", "blue", "#7309ae")
lightness_factor <- 0.15
node_colors <- rep("gray", vcount(network))  
highlighted_nodes <- met_cols

for (i in 1:length(groups)) {
  group <- groups[i]
  group_indices <- annotation2$Class == group
  metabolites <- annotation2$variable[group_indices]
  matching_indices <- match(metabolites, V(network)$name)
  color <- col2rgb(color_palette[i])
  if (length(highlighted_nodes) > 0 && any(metabolites %in% highlighted_nodes)) {
    highlighted_indices <- matching_indices[metabolites %in% highlighted_nodes]
    node_colors[highlighted_indices] <- rgb(color[1], color[2], color[3], maxColorValue = 255, names = "rgb")
    nonhighlighted_indices <- setdiff(matching_indices, highlighted_indices)
    if (length(nonhighlighted_indices) > 0) {
      lighter_color <- rgb(
        (1 - lightness_factor) * 255 + lightness_factor * color[1],
        (1 - lightness_factor) * 255 + lightness_factor * color[2],
        (1 - lightness_factor) * 255 + lightness_factor * color[3],
        maxColorValue = 255,
        names = "rgb"
      )
      node_colors[nonhighlighted_indices] <- lighter_color
    }
  } else {
    lighter_color <- rgb(
      (1 - lightness_factor) * 255 + lightness_factor * color[1],
      (1 - lightness_factor) * 255 + lightness_factor * color[2],
      (1 - lightness_factor) * 255 + lightness_factor * color[3],
      maxColorValue = 255,
      names = "rgb"
    )
    node_colors[matching_indices] <- lighter_color
  }
}
V(network)$color <- node_colors

#Add labels specifically for metabolites in met_cols
metabolite_indices <- match(met_cols, V(network)$name)
V(network)$label <- NA
V(network)$label[metabolite_indices] <- met_cols
#shrink nodes with many connections
deg <- degree(network, mode="all")

V(network)$label

#plot the graph
plot(network, edge.arrow.size = 0, 
     vertex.size = 1/(0.1 +(deg*0.005)),
     vertex.label.color = "black", 
     vertex.label.cex = 0.5,
     vertex.label.dist = 1,
     vertex.border.color = V(network)$color)
dev.off()

#plot the legend
legend_title <- "Metabolite Category"
legend_labels <- groups
legend_colors <- color_palette
plot(0, 0, xlim = c(-5, 5), ylim = c(-5, 5), type = "n", xlab = "", ylab = "", axes = FALSE)
legend_width <- 0.15
legend_height <- 0.25
legend("center", title = legend_title, legend = legend_labels, fill = legend_colors,
       inset = c(0.5, 0.5, 0.5 - legend_width, 0.5 - legend_height))
dev.off()


#CV performance within derivation cohort
#partition derivation split
partition_indices <- list()
for (i in 1:10) {
  partition_indices[[i]] <- createDataPartition(trainData$HF, p = 0.9, list = FALSE, times = 1)
}
#CV loop
absolute_results_df <- data.frame(Model = character(),
                                  Sensitivity = numeric(),
                                  Specificity = numeric(),
                                  C_Index = numeric(),
                                  stringsAsFactors = FALSE)
relative_results_df <- data.frame(ModelComparison = character(),
                                  Absolute_IDI = numeric(),
                                  Relative_IDI = numeric(),
                                  Cases_NRI = numeric(),
                                  NonCases_NRI = numeric(),
                                  Overall_NRI = numeric(),
                                  stringsAsFactors = FALSE)
for (i in 1:10) {
  fold_indices <- partition_indices[[i]]
  split <- trainData[fold_indices,]
  #LPs
  split$LP_agesex <- predict(fit_agesex$model, newx = as.matrix(split[,c(169:170)]), type = "link")
  split$LP_agesex_met <- predict(fit_agesex_met$model, newx = as.matrix(split[,c(1:170)]), type = "link")
  split$LP_PCP <- predict(fit_PCP$model, newx = as.matrix(split[,c(169:170,172,173,195,208,209,210,212,335,336,337)]), type = "link")
  split$LP_PCP_met <- predict(fit_PCP_met$model, newx = as.matrix(split[,c(1:170,172,173,195,208,209,210,212,335,336,337)]), type = "link")
  
  #C
  C_agesex <- Cindex(split$LP_agesex, Surv(split$followup_HF, split$HF))
  C_agesex_met <- Cindex(split$LP_agesex_met, Surv(split$followup_HF, split$HF))
  C_PCP <- Cindex(split$LP_PCP, Surv(split$followup_HF, split$HF))
  C_PCP_met <- Cindex(split$LP_PCP_met, Surv(split$followup_HF, split$HF))
  
  #Sens and Spec
  roc_agesex_train <- roc(split$HF, split$LP_agesex[,1])
  roc_agesex_met_train <- roc(split$HF, split$LP_agesex_met[,1])
  roc_PCP_train <- roc(split$HF, split$LP_PCP[,1])
  roc_PCP_met_train <- roc(split$HF, split$LP_PCP_met[,1])
  youden_agesex_train <- coords(roc_agesex_train, "best", best.method="youden")
  youden_agesex_met_train <- coords(roc_agesex_met_train, "best", best.method="youden")
  youden_PCP_train <- coords(roc_PCP_train, "best", best.method="youden")
  youden_PCP_met_train <- coords(roc_PCP_met_train, "best", best.method="youden")
  
  #NRI
  NRI_agesex <- calculate_nri(split, "LP_agesex", "LP_agesex_met", "HF")
  NRI_nomet <- calculate_nri(split, "LP_agesex", "LP_PCP", "HF")
  NRI_agesexvspcpmet <- calculate_nri(split, "LP_agesex", "LP_PCP_met", "HF")
  NRI_agesexmetvspcp <- calculate_nri(split, "LP_agesex_met", "LP_PCP", "HF")
  NRI_met <- calculate_nri(split, "LP_agesex_met", "LP_PCP_met", "HF")
  NRI_PCP <- calculate_nri(split, "LP_PCP", "LP_PCP_met", "HF")
  
  #idi
  Calibration_cv <- split[!(split$followup_HF < 10 & split$HF == FALSE),] 
  Calibration_cv$HF10y <- ifelse(Calibration_cv$followup_HF <= 10 & Calibration_cv$HF == TRUE, TRUE, FALSE)
  Calibration_cv$risk_agesex <- predict(fit_agesex, as.matrix(split[,c(169:170)]), Surv(split$followup_HF, split$HF), newx = as.matrix(Calibration_cv[,c(169:170)]), pred.at = 10)
  Calibration_cv$risk_agesex_met <- predict(fit_agesex_met, as.matrix(split[,c(1:170)]), Surv(split$followup_HF, split$HF), newx = as.matrix(Calibration_cv[,c(1:170)]), pred.at = 10)
  Calibration_cv$risk_PCP <- predict(fit_PCP, as.matrix(split[,c(169:170,172,173,195,208,209,210,212,335,336,337)]), Surv(split$followup_HF, split$HF), newx = as.matrix(Calibration_cv[,c(169:170,172,173,195,208,209,210,212,335,336,337)]), pred.at = 10)
  Calibration_cv$risk_PCP_met <- predict(fit_PCP_met, as.matrix(split[,c(1:170,172,173,195,208,209,210,212,335,336,337)]), Surv(split$followup_HF, split$HF), newx = as.matrix(Calibration_cv[,c(1:170,172,173,195,208,209,210,212,335,336,337)]), pred.at = 10)
  IDI_agesex <- calculate_idi(Calibration_cv, "risk_agesex", "risk_agesex_met", "HF10y")
  IDI_nomet <- calculate_idi(Calibration_cv, "risk_agesex", "risk_PCP", "HF10y")
  IDI_agesexvspcpmet <- calculate_idi(Calibration_cv, "risk_agesex", "risk_PCP_met", "HF10y")
  IDI_agesexmetvspcp <- calculate_idi(Calibration_cv, "risk_agesex_met", "risk_PCP", "HF10y")
  IDI_met <- calculate_idi(Calibration_cv, "risk_agesex_met", "risk_PCP_met", "HF10y")
  IDI_PCP <- calculate_idi(Calibration_cv, "risk_PCP", "risk_PCP_met", "HF10y")
  
  #append results
  #rows for abs
  abs_row <- data.frame(Model = c("AgeSex", "AgeSex_met", "PCP", "PCP_met"),
                        Sensitivity = c(youden_agesex_train$sensitivity, youden_agesex_met_train$sensitivity,
                                        youden_PCP_train$sensitivity, youden_PCP_met_train$sensitivity),
                        Specificity = c(youden_agesex_train$specificity, youden_agesex_met_train$specificity,
                                        youden_PCP_train$specificity, youden_PCP_met_train$specificity),
                        C_Index = c(C_agesex, C_agesex_met, C_PCP, C_PCP_met))
        
  #create rows for rel
  rel_row <- data.frame(ModelComparison = c("AgeSex vs AgeSex_met", "AgeSex vs PCP", "AgeSex vs PCP_met",
                                            "AgeSex_met vs PCP", "AgeSex_met vs PCP_met", "PCP vs PCP_met"),
                        Absolute_IDI = c(IDI_agesex$abs_idi, IDI_nomet$abs_idi, IDI_agesexvspcpmet$abs_idi, 
                                         IDI_agesexmetvspcp$abs_idi, IDI_met$abs_idi, IDI_PCP$abs_idi),
                        Relative_IDI = c(IDI_agesex$rel_idi, IDI_nomet$rel_idi, IDI_agesexvspcpmet$rel_idi, 
                                         IDI_agesexmetvspcp$rel_idi, IDI_met$rel_idi, IDI_PCP$rel_idi),
                        Cases_NRI = c(NRI_agesex$events_nri, NRI_nomet$events_nri, NRI_agesexvspcpmet$events_nri, 
                                      NRI_agesexmetvspcp$events_nri, NRI_met$events_nri, NRI_PCP$events_nri),
                        NonCases_NRI = c(NRI_agesex$non_events_nri, NRI_nomet$non_events_nri, NRI_agesexvspcpmet$non_events_nri, 
                                         NRI_agesexmetvspcp$non_events_nri, NRI_met$non_events_nri, NRI_PCP$non_events_nri),
                        Overall_NRI = c(NRI_agesex$nri, NRI_nomet$nri, NRI_agesexvspcpmet$nri, 
                                        NRI_agesexmetvspcp$nri, NRI_met$nri, NRI_PCP$nri))
  
  #append
  absolute_results_df <- rbind(absolute_results_df, abs_row)
  relative_results_df <- rbind(relative_results_df, rel_row)
}

#CI
absolute_results_df
#absolute
condensed_absolute_results <- absolute_results_df %>%
  dplyr::group_by(Model) %>%
  dplyr::summarize(
    Condensed_Sensitivity = paste(
      formatC(round(mean(Sensitivity), 3), format = "f", digits = 3),
      "[",
      formatC(round(mean(Sensitivity) - 1.96 * sd(Sensitivity) / sqrt(n()), 3), format = "f", digits = 3),
      "-",
      formatC(round(mean(Sensitivity) + 1.96 * sd(Sensitivity) / sqrt(n()), 3), format = "f", digits = 3),
      "]",
      sep = ""),
    Condensed_Specificity = paste(
      formatC(round(mean(Specificity), 3), format = "f", digits = 3),
      "[",
      formatC(round(mean(Specificity) - 1.96 * sd(Specificity) / sqrt(n()), 3), format = "f", digits = 3),
      "-",
      formatC(round(mean(Specificity) + 1.96 * sd(Specificity) / sqrt(n()), 3), format = "f", digits = 3),
      "]",
      sep = ""),
    Condensed_C_Index = paste(
      formatC(round(mean(C_Index), 3), format = "f", digits = 3),
      "[",
      formatC(round(mean(C_Index) - 1.96 * sd(C_Index) / sqrt(n()), 3), format = "f", digits = 3),
      "-",
      formatC(round(mean(C_Index) + 1.96 * sd(C_Index) / sqrt(n()), 3), format = "f", digits = 3),
      "]",
      sep = "")
  )

condensed_absolute_results

#relative
condensed_relative_results <- relative_results_df %>%
  dplyr::group_by(ModelComparison) %>%
  dplyr::summarize(
    Condensed_Absolute_IDI = paste(
      formatC(round(mean(Absolute_IDI), 3), format = "f", digits = 3),
      "[",
      formatC(round(mean(Absolute_IDI) - 1.96 * sd(Absolute_IDI) / sqrt(n()), 3), format = "f", digits = 3),
      "-",
      formatC(round(mean(Absolute_IDI) + 1.96 * sd(Absolute_IDI) / sqrt(n()), 3), format = "f", digits = 3),
      "]",
      sep = ""),
    Condensed_Relative_IDI = paste(
      formatC(round(mean(Relative_IDI), 3), format = "f", digits = 3),
      "[",
      formatC(round(mean(Relative_IDI) - 1.96 * sd(Relative_IDI) / sqrt(n()), 3), format = "f", digits = 3),
      "-",
      formatC(round(mean(Relative_IDI) + 1.96 * sd(Relative_IDI) / sqrt(n()), 3), format = "f", digits = 3),
      "]",
      sep = ""),
    Condensed_Cases_NRI = paste(
      formatC(round(mean(Cases_NRI), 3), format = "f", digits = 3),
      "[",
      formatC(round(mean(Cases_NRI) - 1.96 * sd(Cases_NRI) / sqrt(n()), 3), format = "f", digits = 3),
      "-",
      formatC(round(mean(Cases_NRI) + 1.96 * sd(Cases_NRI) / sqrt(n()), 3), format = "f", digits = 3),
      "]",
      sep = ""),
    Condensed_NonCases_NRI = paste(
      formatC(round(mean(NonCases_NRI), 3), format = "f", digits = 3),
      "[",
      formatC(round(mean(NonCases_NRI) - 1.96 * sd(NonCases_NRI) / sqrt(n()), 3), format = "f", digits = 3),
      "-",
      formatC(round(mean(NonCases_NRI) + 1.96 * sd(NonCases_NRI) / sqrt(n()), 3), format = "f", digits = 3),
      "]",
      sep = ""),
    Condensed_Overall_NRI = paste(
      formatC(round(mean(Overall_NRI), 3), format = "f", digits = 3),
      "[",
      formatC(round(mean(Overall_NRI) - 1.96 * sd(Overall_NRI) / sqrt(n()), 3), format = "f", digits = 3),
      "-",
      formatC(round(mean(Overall_NRI) + 1.96 * sd(Overall_NRI) / sqrt(n()), 3), format = "f", digits = 3),
      "]",
      sep = "")
  )
condensed_relative_results

#stats
#age sex vs agesexmet
age_sex <- subset(absolute_results_df, Model %in% c("AgeSex"))
age_sex_met <- subset(absolute_results_df, Model %in% c("AgeSex_met"))
metrics_to_compare <- c("Sensitivity", "Specificity", "C_Index")
t_test_results <- lapply(metrics_to_compare, function(metric) {
  t_test_result <- t.test(age_sex[[metric]], age_sex_met[[metric]], paired = FALSE)
  return(data.frame(Metric = metric,
                    Mean_AgeSex = mean(age_sex[[metric]]),
                    Mean_AgeSexMet = mean(age_sex_met[[metric]]),
                    p_value = t_test_result$p.value))
})
t_test_results_df <- do.call(rbind, t_test_results)
t_test_results_df

#pcp vs pcpmet
pcp_m <- subset(absolute_results_df, Model %in% c("PCP"))
pcp_m_met <- subset(absolute_results_df, Model %in% c("PCP_met"))
metrics_to_compare <- c("Sensitivity", "Specificity", "C_Index")
t_test_results <- lapply(metrics_to_compare, function(metric) {
  t_test_result <- t.test(pcp_m[[metric]], pcp_m_met[[metric]], paired = FALSE)
  return(data.frame(Metric = metric,
                    Mean_PCP = mean(pcp_m[[metric]]),
                    Mean_PCPMet = mean(pcp_m_met[[metric]]),
                    p_value = t_test_result$p.value))
})
t_test_results_df <- do.call(rbind, t_test_results)
t_test_results_df

#agesexmet vs pcp
metrics_to_compare <- c("Sensitivity", "Specificity", "C_Index")
t_test_results <- lapply(metrics_to_compare, function(metric) {
  t_test_result <- t.test(age_sex_met[[metric]], pcp_m[[metric]], paired = FALSE)
  return(data.frame(Metric = metric,
                    Mean_PCP = mean(age_sex_met[[metric]]),
                    Mean_PCPMet = mean(pcp_m[[metric]]),
                    p_value = t_test_result$p.value))
})
t_test_results_df <- do.call(rbind, t_test_results)
t_test_results_df



age_sex_age_sex_met <- subset(relative_results_df, ModelComparison %in% c("AgeSex vs AgeSex_met"))
age_sex_pcp <- subset(relative_results_df, ModelComparison %in% c("AgeSex vs PCP"))
metrics_to_compare <- c("Relative_IDI", "Absolute_IDI", "Cases_NRI", "NonCases_NRI", "Overall_NRI")
t_test_results <- lapply(metrics_to_compare, function(metric) {
  t_test_result <- t.test(age_sex_age_sex_met[[metric]], age_sex_pcp[[metric]], paired = TRUE)
  return(data.frame(Metric = metric,
                    Mean_AgeSex_AgeSex_met = mean(age_sex_age_sex_met[[metric]]),
                    Mean_AgeSex_PCP = mean(age_sex_pcp[[metric]]),
                    p_value = t_test_result$p.value))
})
t_test_results_df <- do.call(rbind, t_test_results)
t_test_results_df



















































