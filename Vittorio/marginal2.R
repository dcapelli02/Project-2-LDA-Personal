library(haven)
library(sas7bdat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(GGally)
library(corrplot)
library(MuMIn)

alz <- read_sas("/Users/vitto/Documents/LONGITUDINAL DATA ANALYSIS/alzheimer25.sas7bdat")

alz$trial <- as.factor(alz$trial)
alz$sex <- as.factor(alz$sex)
alz$edu <- as.factor(alz$edu)
alz$job <- as.factor(alz$job)
alz$wzc <- as.factor(alz$wzc)

alz$ab_base <- alz$abpet0
alz$tau_base <- alz$taupet0
alz$cdrsb_base <- alz$cdrsb0

## Bins

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
alz_long$year_seq <- ave(alz_long$year, alz_long$sample,
                         FUN = function(x) as.integer(factor(x)))




#exploratory

df_plot <- alz_long %>%
  filter(!is.na(year), !is.na(cdrsb_bin)) %>%   
  mutate(
    severity = ifelse(cdrsb_bin == 1, "Severe", "Not severe"),
    
    # converti year in stringa coerente
    time = case_when(
      year == 0 ~ "Baseline",
      year == 1 ~ "year1",
      year == 2 ~ "year2",
      year == 3 ~ "year3",
      year == 4 ~ "year4",
      year == 5 ~ "year5",
      year == 6 ~ "year6",
      TRUE ~ as.character(year)
    ),
    
    # fattore con ordine desiderato
    time = factor(time, levels = c("Baseline","year1","year2","year3",
                                   "year4","year5","year6"))
  ) %>%
  group_by(time, severity) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(time) %>%
  mutate(
    total_nonmiss = sum(n),
    percent = n / total_nonmiss * 100
  ) %>%
  ungroup()

# ---- Grafico -----------------------------------------------------------------

ggplot(df_plot, aes(x = severity, y = n, fill = severity)) +
  geom_bar(stat = "identity", width = 0.75, color = "black") +
  geom_text(aes(label = sprintf("%.1f%%", percent)),
            vjust = -0.5, size = 4) +
  facet_wrap(~ time, ncol = 4) +
  scale_fill_grey(start = 0.2, end = 0.6) +
  labs(x = "", y = "# observations") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )
alz_long %>%
  group_by(year) %>%
  summarize(
    prop = mean(cdrsb, na.rm = TRUE),
    n = n(),
    n_nonmiss = sum(!is.na(cdrsb))
  )

#intresting!!


# LOOK FOR THE BEST MODEL USING QIC (Guillem procedure)

vars_to_keep <- c("cdrsb", "sex", "age", "bmi", 
                  "ab_base", "tau_base", "adl", 
                  "year", "year_seq", "sample")

alz_subset <- alz_long[, vars_to_keep]

alz_clean <- na.omit(alz_subset) #idk if it is reasonable to delate NA

# With Unstructured
full_model_un <- geeglm(
  cdrsb ~ (sex + age + bmi + ab_base + tau_base + adl) * year,
  family = binomial("logit"), 
  data = alz_subset, 
  id = sample, 
  #waves = year_seq, i dont think it is necessary
  corstr = "unstructured",
  #na.action = "na.fail"#without
)
# summary(full_model_un)
# full_model_un2 <- geeglm(
#   cdrsb ~ (sex + age + bmi + ab_base + tau_base + adl) * year,
#   family = binomial("logit"), 
#   data = alz_clean, 
#   id = sample, 
#   #waves = year_seq, i dont think it is necessary
#   corstr = "unstructured",
#   na.action = "na.fail"#without
# )
# summary(full_model_un2)

#it is the same, so we can use alz_clean <- na.omit(alz_subset) and na.action = "na.fail"
# or directly na.action = na.omit or defoult

full_model_un1 <- geeglm(
  cdrsb ~ (sex + age + bmi + ab_base + tau_base + adl) * year,
  family = binomial("logit"), 
  data = alz_clean, 
  id = sample, 
  waves = year_seq, #i dont think it is necessary
  corstr = "unstructured",
  na.action = "na.fail"
)

# Try all combinations that have sex, age, bmi and year
results_un <- dredge(
  full_model_un, 
  rank = "QIC", 
  fixed = c("sex", "age", "bmi", "year") 
)

results_un1 <- dredge(
  full_model_un1, 
  rank = "QIC", 
  fixed = c("sex", "age", "bmi", "year") 
)
# Compare using QIC
best_model_un <- get.models(results_un, subset = 1)[[1]]
summary(best_model_un)
best_model_un1 <- get.models(results_un1, subset = 1)[[1]]
summary(best_model_un)

#is the same with or without waves

#not bad as procedure, good! 




