library(haven)
library(dplyr)
library(tidyr)
library(ggplot2)
library(geepack)
library(lme4)
library(car)
library(lmtest)

alz <- read_sas("C:/Users/win11/Documents/Guillem/Erasmus/Assignatures/Longitudinal Data Analysis/Project-2-LDA-Personal/alzheimer25.sas7bdat")

head(alz)
summary(alz)

#Factor Variables
alz$trial <- as.factor(alz$trial)
alz$sex <- as.factor(alz$sex)
alz$edu <- as.factor(alz$edu)
alz$job <- as.factor(alz$job)
alz$wzc <- as.factor(alz$wzc)
alz$adl <- as.factor(alz$adl)
alz$adl_num <- as.numeric(alz$adl)
alz$n_obs_data <- rowSums(!is.na(alz[, c(18:24)]))


alz_long$CDRSBbin <- ifelse(alz_long$cdrsb > 10, 1, 0)

#Bins
alz$cdrsb_bin0 <- ifelse(alz$cdrsb0 > 10, 1, 0)
alz$cdrsb_bin1 <- ifelse(alz$cdrsb1 > 10, 1, 0)
alz$cdrsb_bin2 <- ifelse(alz$cdrsb2 > 10, 1, 0)
alz$cdrsb_bin3 <- ifelse(alz$cdrsb3 > 10, 1, 0)
alz$cdrsb_bin4 <- ifelse(alz$cdrsb4 > 10, 1, 0)
alz$cdrsb_bin5 <- ifelse(alz$cdrsb5 > 10, 1, 0)
alz$cdrsb_bin6 <- ifelse(alz$cdrsb6 > 10, 1, 0)

#Baseline
alz$ab_base <- alz$abpet0
alz$tau_base <- alz$taupet0
alz$cdrsb_base <- alz$cdrsb0
alz$bprs_base <- alz$bprs0

summary(alz)

