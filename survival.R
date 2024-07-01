#survival analysis
library(survival)
library(ggplot2)
library(ggpubr)
library(survminer)

#loading csv file
setwd("C:/Users/achir/Desktop/TCGA/Survival Analysis")
sfile=read.csv("Survivalfile.csv")
dim(sfile)
head(sfile)

#conversion to year
sfile["OS"]= sfile["days_to_death"]*0.00273973
sfile$VS=ifelse(sfile$vital_status=="dead",1,0)
dim(sfile)
sfile=sfile[1:1044,1:10]
dim(sfile)
sfile1=write.csv(sfile, file = "survivaladd.csv")
head(sfile)
View(sfile)

#rownames
sr=sfile[,-1]
rownames(sr)= make.names(sfile[,1], unique=TRUE) 
View(sr)
head(sr)
sum(is.na(sr['OS']))

data_s=survfit(Surv(as.numeric(OS),as.numeric(VS))~Cancer.Subtypes,data=sfile)
data_s

ggsurvplot(data_s,data=sfile,surv.median.line = 'hv',pval=T,conf.int = F,
           risk.table = T,tables.height = 0.2, tables.theme = theme_cleantable(),break.time.by = 2,
           legend.labs=c("HER2+", "Luminal A", "Luminal B", "TNBC"), legend.title="Cancer Subtypes",
           title="Kaplan-Meier Curve for BIC Survival",font.main = c(16, "bold"),ggtheme = theme_bw())

data_s1=survfit(Surv(as.numeric(OS),as.numeric(VS))~gender,data=sfile)
data_s1

ggsurvplot(data_s1,data=sfile,surv.median.line = 'hv',pval=T,conf.int = F,
           risk.table = T,tables.height = 0.2, tables.theme = theme_cleantable(),break.time.by = 2,
           legend.labs=c("Female","Male"), legend.title="Gender",
           title="Kaplan-Meier Curve for BIC Survival",font.main = c(16, "bold"),ggtheme = theme_bw())

data_s2=survfit(Surv(as.numeric(OS),as.numeric(VS))~race,data=sfile)
data_s2

ggsurvplot(data_s2,data=sfile,surv.median.line = 'hv',pval=T,conf.int = F,
           risk.table = T,tables.height = 0.2, tables.theme = theme_cleantable(),break.time.by = 2,
           legend.labs=c("Asian","African American","White"), legend.title="Race",
           title="Kaplan-Meier Curve for BIC Survival",font.main = c(16, "bold"),ggtheme = theme_bw())





















