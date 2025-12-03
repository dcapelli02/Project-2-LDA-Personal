library(haven)
library(sas7bdat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(GGally)
library(corrplot)
alz <- read_sas("/Users/vitto/Documents/LONGITUDINAL DATA ANALYSIS/alzheimer25.sas7bdat")

alz$cdrsb_bin0 <- alz$cdrsb0
alz$cdrsb_bin0[alz$cdrsb0 > 10] <- 1
alz$cdrsb_bin0[alz$cdrsb0 <= 10] <- 0

alz$cdrsb_bin1 <- alz$cdrsb1
alz$cdrsb_bin1[alz$cdrsb1 > 10] <- 1
alz$cdrsb_bin1[alz$cdrsb1 <= 10] <- 0

alz$cdrsb_bin2 <- alz$cdrsb2
alz$cdrsb_bin2[alz$cdrsb2 > 10] <- 1
alz$cdrsb_bin2[alz$cdrsb2 <= 10] <- 0

alz$cdrsb_bin3 <- alz$cdrsb3
alz$cdrsb_bin3[alz$cdrsb3 > 10] <- 1
alz$cdrsb_bin3[alz$cdrsb3 <= 10] <- 0

alz$cdrsb_bin4 <- alz$cdrsb4
alz$cdrsb_bin4[alz$cdrsb4 > 10] <- 1
alz$cdrsb_bin4[alz$cdrsb4 <= 10] <- 0

alz$cdrsb_bin5 <- alz$cdrsb5
alz$cdrsb_bin5[alz$cdrsb5 > 10] <- 1
alz$cdrsb_bin5[alz$cdrsb5 <= 10] <- 0

alz$cdrsb_bin6 <- alz$cdrsb6
alz$cdrsb_bin6[alz$cdrsb6 > 10] <- 1
alz$cdrsb_bin6[alz$cdrsb6 <= 10] <- 0


## Create baseline values
alz$ab_base <- alz$abpet0
alz$tau_base <- alz$taupet0
alz$cdrsb_base <- alz$cdrsb0
alz$bprs_base <- alz$bprs0

summary(alz)

## Create longitudinal dataset
alz_df <- data.frame(alz)

alz_long <- alz_df %>%
  pivot_longer(
    
    cols = matches("^(bprs|cdrsb_bin|abpet|taupet)\\d+$"),
    
    
    names_to = c(".value", "year"),
    
    names_pattern = "(bprs|cdrsb_bin|abpet|taupet)(\\d+)"
  ) %>%
  mutate(
    year = as.numeric(year),                          
    sample = factor(rep(1:nrow(alz_df), each = 7)) # ID per ogni paziente
  )

alz_long$cdrsb <- alz_long$cdrsb_bin
alz_long$cdrsb_bin <- as.factor(alz_long$cdrsb_bin)

# -------------------------------------------------------------
#INFORMATIVE DROPOUT (?)
# -------------------------------------------------------------

#I find the latest reading for each patient

ultima_rilevazione <- alz_long %>%
  group_by(patid) %>%
  summarise(last_time = max(year))

#  indidcator 1 = there is untile the last one observation
dati <- alz_long %>%
  left_join(ultima_rilevazione, by = "patid") %>%
  mutate(last_visit = ifelse(year == last_time, 1, 0))

ultimo_anno <- max(alz_long$year, na.rm = TRUE)

cdrsb_baseline <- alz_long %>%
  group_by(patid) %>%
  summarise(
    # first value of BPRS
    CDRSB_start = cdrsb_bin[!is.na(cdrsb_bin) & year == min(year[!is.na(cdrsb_bin)])][1],
    # latest reading
    last_year = max(year[!is.na(cdrsb_bin)], na.rm = TRUE)
  ) %>%
  mutate(
    completed = ifelse(last_year == ultimo_anno, 1, 0)
  )

tabella <- table(cdrsb_baseline$CDRSB_start, cdrsb_baseline$completed)
rownames(tabella) <- c("CDRSB ≤ 10", "CDRSB > 10")
colnames(tabella) <- c("Abandoned", "Completed")

tabella
# Test chi-quadro
test_chi <- chisq.test(tabella)
test_chi

ggplot(cdrsb_baseline, aes(x = factor(completed, labels = c("Abandoned", "Completed")), 
                           fill = factor(CDRSB_start, labels = c("CDRSB ≤ 10", "CDRSB > 10")))) +
  geom_bar(position = "fill") +
  labs(x = "Study Completion", y = "Proportion", 
       fill = "Baseline CDRSB",
       title = "Proportion of CDRSB Categories by Completion Status") +
  theme_minimal() +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e"))



#conclusion:
#there is an informative drop out, who start the study with high levels of bprs drop out earlier than those who start with low levels of bprs



