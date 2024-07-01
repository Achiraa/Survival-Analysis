# Survival Analysis
## This repository contains R scripts and data for performing survival analysis using the survival, ggplot2, ggpubr, and survminer packages. The analysis includes Kaplan-Meier survival curves for different cancer subtypes, gender, and race.

## Prerequisites
Ensure you have the following R packages installed:
install.packages(c("survival", "ggplot2", "ggpubr", "survminer"))

## Setup
Set the working directory to the location of your data files:
setwd("C:/Users/achir/Desktop/TCGA/Survival Analysis")

## Data
Survival Data: Survivalfile.csv
## Steps
## 1. Load Libraries

library(survival)
library(ggplot2)
library(ggpubr)
library(survminer)
## 2. Load and Prepare Data

sfile = read.csv("Survivalfile.csv")
dim(sfile)
head(sfile)
## 3. Convert Days to Years and Create Vital Status Column

sfile["OS"] = sfile["days_to_death"] * 0.00273973
sfile$VS = ifelse(sfile$vital_status == "dead", 1, 0)
dim(sfile)
sfile = sfile[1:1044, 1:10]
dim(sfile)
write.csv(sfile, file = "survivaladd.csv")
head(sfile)
View(sfile)
## 4. Set Row Names

sr = sfile[,-1]
rownames(sr) = make.names(sfile[,1], unique = TRUE)
View(sr)
head(sr)
sum(is.na(sr['OS']))
## 5. Kaplan-Meier Curve by Cancer Subtypes

data_s = survfit(Surv(as.numeric(OS), as.numeric(VS)) ~ Cancer.Subtypes, data = sfile)
data_s

ggsurvplot(data_s, data = sfile, surv.median.line = 'hv', pval = TRUE, conf.int = FALSE,
           risk.table = TRUE, tables.height = 0.2, tables.theme = theme_cleantable(), break.time.by = 2,
           legend.labs = c("HER2+", "Luminal A", "Luminal B", "TNBC"), legend.title = "Cancer Subtypes",
           title = "Kaplan-Meier Curve for BIC Survival", font.main = c(16, "bold"), ggtheme = theme_bw())
## 6. Kaplan-Meier Curve by Gender

data_s1 = survfit(Surv(as.numeric(OS), as.numeric(VS)) ~ gender, data = sfile)
data_s1

ggsurvplot(data_s1, data = sfile, surv.median.line = 'hv', pval = TRUE, conf.int = FALSE,
           risk.table = TRUE, tables.height = 0.2, tables.theme = theme_cleantable(), break.time.by = 2,
           legend.labs = c("Female", "Male"), legend.title = "Gender",
           title = "Kaplan-Meier Curve for BIC Survival", font.main = c(16, "bold"), ggtheme = theme_bw())
           
## 7. Kaplan-Meier Curve by Race

data_s2 = survfit(Surv(as.numeric(OS), as.numeric(VS)) ~ race, data = sfile)
data_s2

ggsurvplot(data_s2, data = sfile, surv.median.line = 'hv', pval = TRUE, conf.int = FALSE,
           risk.table = TRUE, tables.height = 0.2, tables.theme = theme_cleantable(), break.time.by = 2,
           legend.labs = c("Asian", "African American", "White"), legend.title = "Race",
           title = "Kaplan-Meier Curve for BIC Survival", font.main = c(16, "bold"), ggtheme = theme_bw())
           
## Output
Survival Data with Additional Columns: survivaladd.csv
Kaplan-Meier Curves:
a.By Cancer Subtypes
b.By Gender
c.By Race

## This analysis performs survival analysis on cancer patients using the Kaplan-Meier method. It includes survival curves based on different cancer subtypes, gender, and race. The data is read from a CSV file, processed to convert days to years, and Kaplan-Meier plots are generated for visualization.
