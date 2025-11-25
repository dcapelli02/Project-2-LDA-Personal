#================================================================#

## LIBRARIES ##

library(haven)
library(dplyr)
library(tidyr)
library(ggplot2)
library(geepack)
library(lme4)
library(car)
library(lmtest)
library(MuMIn)
library(multgee)
library(ORTH.Ord)

#================================================================#

## DATASETS ##

alz <- read_sas("C:/Users/win11/Documents/Guillem/Erasmus/Assignatures/Longitudinal Data Analysis/Project-2-LDA-Personal/alzheimer25.sas7bdat")

#Factor Variables

alz$trial <- as.factor(alz$trial)
alz$sex <- as.factor(alz$sex)
alz$edu <- as.factor(alz$edu)
alz$job <- as.factor(alz$job)
alz$wzc <- as.factor(alz$wzc)
alz$adl <- as.factor(alz$adl)
alz$adl_num <- as.numeric(alz$adl)
alz$n_obs_data <- rowSums(!is.na(alz[, c(18:24)]))

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
    sample = factor(rep(1:nrow(alz_df), each = 7))
  )

alz_long$cdrsb <- alz_long$cdrsb_bin
alz_long$cdrsb_bin <- as.factor(alz_long$cdrsb_bin)

## Discretize variables

## Maybe any 5 years?
alz_long$age_disc <- (alz_long$age %/% 5) * 5
alz_long$age_disc <- as.factor(alz_long$age_disc)

## bmi any 4
alz_long$bmi_disc <- (alz_long$bmi %/% 4) * 4
alz_long$bmi_disc <- as.factor(alz_long$bmi_disc)

## inkomen any 500
alz_long$inkomen_disc <- (alz_long$inkomen %/% 500) * 500
alz_long$inkomen_disc <- as.factor(alz_long$inkomen_disc)

## adl any 5
alz_long$adl_disc <- (as.numeric(alz_long$adl) %/% 5) * 5
alz_long$adl_disc <- as.factor(alz_long$adl_disc)

## cdrsb any 5
#alz_long$cdrsb_disc <- (alz_long$cdrsb %/% 5) * 5
#alz_long$cdrsb_disc <- as.factor(alz_long$cdrsb_disc)

## abpet any 0.2
alz_long$abpet_disc <- (alz_long$abpet %/% 0.2) * 0.2
alz_long$abpet_disc <- as.factor(alz_long$abpet_disc)

## taupet any 0.2
alz_long$taupet_disc <- (alz_long$taupet %/% 0.2) * 0.2
alz_long$taupet_disc <- as.factor(alz_long$taupet_disc)

## adl any 5
alz$adl_disc <- (as.numeric(alz$adl) %/% 5) * 5
alz$adl_disc <- as.factor(alz$adl_disc)

## year discrete
alz_long$year_seq <- ave(alz_long$year, alz_long$sample,
                         FUN = function(x) as.integer(factor(x)))

summary(alz_long)

#================================================================#

## EXPLORATORY ANALYSIS ##

# General Mean Behaviour
ggplot(alz_long, aes(x = year, y = cdrsb)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
               alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior with variance")

ggplot(alz_long, aes(x = year, y = cdrsb)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  #stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
  #             alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior")

cat("Since this curve is not a straight line, we should  treat time as a factor.")
cat("We are allowing the model to estimate a different prob for each year.")

## SEX
ggplot(alz_long, aes(x = year, y = cdrsb, group = sex, color = sex)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  #stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
  #             alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior wrt age_disc without variance")

## AGE
ggplot(alz_long, aes(x = year, y = cdrsb, group = age_disc, color = age_disc)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  #stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
  #             alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior wrt age_disc without variance")

## BMI
ggplot(alz_long, aes(x = year, y = cdrsb, group = bmi_disc, color = bmi_disc)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  #stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
  #             alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior wrt age_disc without variance")

## ADL
ggplot(alz_long, aes(x = year, y = cdrsb, group = adl_disc, color = adl_disc)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  #stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
  #             alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior wrt age_disc without variance")

## AB
ggplot(alz_long, aes(x = year, y = cdrsb, group = abpet_disc, color = abpet_disc)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  #stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
  #             alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior wrt age_disc without variance")

## TAU
ggplot(alz_long, aes(x = year, y = cdrsb, group = taupet_disc, color = taupet_disc)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  #stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
  #             alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior wrt age_disc without variance")

cat("We see that the mean structure follows a logistic function.")

## CORRELATION
cor_matrix_cdrsb <- cor(alz[, c(11:17)], use = "pairwise.complete.obs")
round(cor_matrix_cdrsb, 2)

heatmap(cor_matrix_cdrsb, Rowv = NA, Colv = NA)  

cat("We see that the first 3 measurements correlated with each other but not with the other last 4.")
cat("And the same happens with the last 4 measurement that are more correlated to each other than to the first.")
cat("Seeing the mean structure behaviour in the first years is lower, and in the last years it is higher.")

## SPAGHETTI PLOT
## Random sample

set.seed(1)
casual1 <- sample(1:length(alz$patid), 50)

## Now we can start by looking at random values for the mean and see
## if we can work on the mean and so on

alz_rist1 <- alz[casual1, ]
alz_rist1_df <- data.frame(alz_rist1)

alz_long_rist <- alz_rist1_df %>%
  pivot_longer(
    
    cols = matches("^(bprs|cdrsb_bin|abpet|taupet)\\d+$"),
    
    
    names_to = c(".value", "year"),
    
    names_pattern = "(bprs|cdrsb_bin|abpet|taupet)(\\d+)"
  ) %>%
  mutate(
    year = as.numeric(year),                          
    sample = factor(rep(1:nrow(alz_rist1_df), each = 7)) # ID per ogni paziente
  )

alz_long_rist$cdrsb <- alz_long_rist$cdrsb_bin
alz_long_rist$cdrsb_bin <- as.factor(alz_long_rist$cdrsb_bin)

ggplot(alz_long_rist, aes(x = year, y = cdrsb, group = patid, 
                          color = sample, show.legend = FALSE)) + 
  geom_line(alpha = 1, show.legend = FALSE, size = 0.7) +
  theme_bw() +
  labs(title = "Time Evolution of CDRSB")

ggplot(alz_long_rist, aes(x = year, y = taupet, group = patid, 
                          color = sample, show.legend = FALSE)) + 
  geom_line(alpha = 1, show.legend = FALSE) +
  theme_bw() +
  labs(title = "TAUPET against year of measuring")

cat("I have an observation, look that all patients dropout after their taupet increases just a little bit.")
cat("Maybe we could consider something about it.")

ggplot(alz_long_rist, aes(x = year, y = abpet, group = patid, 
                          color = sample, show.legend = FALSE)) + 
  geom_line(alpha = 1, show.legend = FALSE) +
  theme_bw() +
  labs(title = "ABPET against year of measuring")

#================================================================#

## MARGINAL MODEL ##

# SUMMARY

# We are modeling the probability of cdrsb being 1. p = P(Y = 1).
# Since p c [0,1]. We need a transformation (link function) to stretch the range to cover (-inf,+inf).
# Because then we can fit a linear equation to it.

# logit(p) = ln( p/(1 - p) )



# MARGINAL MODEL:

# Y_ij  ~  Be(mu_ij)
# logit(mu_ij)  =  beta_0  +  ...

# logit( P(Y_ij = 1) )  =  Intercept  +  Main Effects  +  Time Effects  +  Interaction Effects

# logit( P(Y_ij = 1) ) = ( beta_0 ) + ( beta_S * Sex_i + beta_A * Age_i + ... ) + ( beta_1 * Year_ij ) + 
#                         + ( beta_SY * (Sex_i * Year_ij ) + beta_AY * (Age_i * Year_ij) + ... )

# This model has no random effects b_i, these will appear in the Conditional Model GLMM:

# CONDITIONAL MODEL

# (Y_ij | bi )  ~  Be(mu_ij)
# logit(mu_ij)  =  beta_0  +  ...

# logit( P(Y_ij = 1 | b_i) )  =  Intercept  +  Random Intercept  +  Main Effects  +  Time Effects  +  Interaction Effects

# logit( P(Y_ij = 1 | b_i) ) = ( beta_0 ) + (b_i) + ( beta_S * Sex_i + beta_A * Age_i + ... ) + ( beta_1 * Year_ij ) + 
#                         + ( beta_SY * (Sex_i * Year_ij ) + beta_AY * (Age_i * Year_ij) + ... )



# family = binomial gives us mean-variance relationship: Var[y] = mu * (1 - mu):
# Expect the data to be most variable when the probability is 0.5
# And least variable when the probability is near 0 or 1

model_gee_un <- geeglm(
  cdrsb ~ (sex + age + bmi + ab_base + tau_base + adl_num) * year,
  family = binomial(link = "logit"),             # Logistic link for binary data
  data = alz_long,
  id = sample,                            # Clustering by patient
  waves = year_seq,
  corstr = "unstructured"                 # Most flexible correlation structure
)

summary(model_gee_un)


model_gee_ex <- geeglm(
  cdrsb ~ (sex + age + bmi + ab_base + tau_base + adl_num) * year,
  family = binomial(link = "logit"),             # Logistic link for binary data
  data = alz_long,
  id = sample,                            # Clustering by patient
  waves = year_seq,
  corstr = "exchangeable"                 # Most flexible correlation structure
)

summary(model_gee_ex)


model_gee_ar1 <- geeglm(
  cdrsb ~ (sex + age + bmi + ab_base + tau_base + adl_num) * year,
  family = binomial(link = "logit"),             # Logistic link for binary data
  data = alz_long,
  id = sample,                            # Clustering by patient
  waves = year_seq,
  corstr = "ar1"                 # Most flexible correlation structure
)

summary(model_gee_ar1)


model_gee_ind <- geeglm(
  cdrsb ~ (sex + age + bmi + ab_base + tau_base + adl_num) * year,
  family = binomial(link = "logit"),             # Logistic link for binary data
  data = alz_long,
  id = sample,                            # Clustering by patient
  waves = year_seq,
  corstr = "independence"                 # Most flexible correlation structure
)

summary(model_gee_ind)


QIC(model_gee_un, model_gee_ex, model_gee_ar1, model_gee_ind)

cat("These models seem unuseful: they don't give much information.")








# LOOK FOR THE BEST MODEL USING QIC

vars_to_keep <- c("cdrsb", "sex", "age", "bmi", 
                  "ab_base", "tau_base", "adl_num", 
                  "year", "sample", "year_seq")

alz_subset <- alz_long[, vars_to_keep]

alz_clean <- na.omit(alz_subset)





# With AR(1)
full_model_ar1 <- geeglm(
  cdrsb ~ (sex + age + bmi + ab_base + tau_base + adl_num) * year,
  family = binomial("logit"), 
  data = alz_clean, 
  id = sample, 
  waves = year_seq, 
  corstr = "ar1",
  na.action = "na.fail"
)

# Try all combinations that have sex, age, bmi and year
results_ar1 <- dredge(
  full_model_ar1, 
  rank = "QIC", 
  fixed = c("sex", "age", "bmi", "year") 
)

# Compare using QIC
best_model_ar1 <- get.models(results_ar1, subset = 1)[[1]]
summary(best_model_ar1)

# Models with less QIC
head(results_ar1)





# With Independence
full_model_ind <- geeglm(
  cdrsb ~ (sex + age + bmi + ab_base + tau_base + adl_num) * year,
  family = binomial("logit"), 
  data = alz_clean, 
  id = sample, 
  waves = year_seq, 
  corstr = "independence",
  na.action = "na.fail"
)

# Try all combinations that have sex, age, bmi and year
results_ind <- dredge(
  full_model_ind, 
  rank = "QIC", 
  fixed = c("sex", "age", "bmi", "year") 
)

# Compare using QIC
best_model_ind <- get.models(results_ind, subset = 1)[[1]]
summary(best_model_ind)

# Models with less QIC
head(results_ind)





# With Exchangeable
full_model_ex <- geeglm(
  cdrsb ~ (sex + age + bmi + ab_base + tau_base + adl_num) * year,
  family = binomial("logit"), 
  data = alz_clean, 
  id = sample, 
  waves = year_seq, 
  corstr = "exchangeable",
  na.action = "na.fail"
)

# Try all combinations that have sex, age, bmi and year
results_ex <- dredge(
  full_model_ex, 
  rank = "QIC", 
  fixed = c("sex", "age", "bmi", "year") 
)

# Compare using QIC
best_model_ex <- get.models(results_ex, subset = 1)[[1]]
summary(best_model_ex)

# Models with less QIC
head(results_ex)





# With Unstructured
full_model_un <- geeglm(
  cdrsb ~ (sex + age + bmi + ab_base + tau_base + adl_num) * year,
  family = binomial("logit"), 
  data = alz_clean, 
  id = sample, 
  waves = year_seq, 
  corstr = "unstructured",
  na.action = "na.fail"
)

# Try all combinations that have sex, age, bmi and year
results_un <- dredge(
  full_model_un, 
  rank = "QIC", 
  fixed = c("sex", "age", "bmi", "year") 
)

# Compare using QIC
best_model_un <- get.models(results_un, subset = 1)[[1]]
summary(best_model_un)

# Models with less QIC
head(results_un)

head(results_un)
summary(best_model_un)
cat("Best Unstrucured: age, bmi, sex, year, age:year, bmi:year and sex:year. QIC = 7802")

head(results_ex)
summary(best_model_ex)
cat("Best Exchangeable: age, bmi, sex, tau_base, year, age:year and bmi:year. QIC = 7808")

head(results_ar1)
summary(best_model_ar1)
cat("Best AR(1): age, bmi, sex, tau_base, year, age:year and bmi:year. QIC = 7807")

head(results_ind)
summary(best_model_ind)
sucat("Best Independence: age, bmi, sex, tau_base, year, age:year and bmi:year. QIC = 7808")




# LOOK FOR THE BEST MODEL MANUALLY


cat("Tried eliminating the variable with the greatest p-value but got a model with high p-values from adl, ab_base and tau_base.")
cat("This means they are not rellevant because age, sex and bmi are 'stealing' explanatory power?")



# ALTERNATING LOGISTIC REGRESSION (ALR)

# I tried using the libraries ORTH.Ord and multgee but didn't work for now.
