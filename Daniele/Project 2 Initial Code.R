## Libraries
library(haven)
library(sas7bdat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(GGally)
library(corrplot)
library(nlme)
library(lmtest)
library(DescTools)
library(Matrix)
library(MASS)
library(metafor)
library(geepack)
library(car)
library(lme4)

## Import and fix the data
alz <- read_sas("C:/Users/Daniele/Desktop/2025 - 26 Primo Semestre/Longitudinal Data Analysis/Project 1 Alzheimer LDA/alzheimer25.sas7bdat")

head(alz)
summary(alz)

alz$trial <- as.factor(alz$trial)
alz$sex <- as.factor(alz$sex)
alz$edu <- as.factor(alz$edu)
alz$job <- as.factor(alz$job)
alz$wzc <- as.factor(alz$wzc)
alz$adl <- as.factor(alz$adl)
alz$adl_num <- as.numeric(alz$adl)
alz$n_obs_data <- rowSums(!is.na(alz[, c(18:24)]))

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

#####################################################################
#####################################################################

#### INITIAL EXPLORATORY DATA ANALYSIS ####

ggplot(alz_long, aes(x = year, y = cdrsb)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
               alpha = 0.2, fill = "lightblue") +
  theme_minimal() +
  labs(title = "General mean behavior")

## Sembra una funzione sigmoide --> la logistic potrebbe fittarla bene

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

### QUA SEMBRA ESSERCI DIFFERENZA TRA I VARI GRUPPI PER AB E TAU


### CORRELATION
cor_matrix_cdrsb <- cor(alz[, c(11:17)], use = "pairwise.complete.obs")
round(cor_matrix_cdrsb, 2)

heatmap(cor_matrix_cdrsb, Rowv = NA, Colv = NA)  ## --> come interpretare??


### SPAGHETTI PLOT
## For the spaghetti plot we need a random sample

set.seed(1234)
casual1 <- sample(1:length(alz$patid), 20)

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


#### 1. MARGINAL MODEL ####

### Bisogna scegliere una link function penso


### PRIMA PROVA MOLTO CASUALE
mod_prova_un <- geeglm(
  cdrsb ~ (sex + age + bmi + ab_base + tau_base + adl_num) * year,
  family = binomial(link = "logit"),
  data = alz_long,
  id = sample,
  waves = year_seq,
  corstr = "unstructured"
)
summary(mod_prova_un)  ## --> sembra tutto inutile...

mod_prova_ex <- geeglm(
  cdrsb ~ (sex + age + bmi + ab_base + tau_base + adl_num) * year,
  family = binomial(link = "logit"),
  data = alz_long,
  id = sample,
  waves = year_seq,
  corstr = "exchangeable"
)
summary(mod_prova_ex)  ## --> sembra tutto inutile...

mod_prova_ar1 <- geeglm(
  cdrsb ~ (sex + age + bmi + ab_base + tau_base + adl_num) * year,
  family = binomial(link = "logit"),
  data = alz_long,
  id = sample,
  waves = year_seq,
  corstr = "ar1"
)
summary(mod_prova_ar1)  ## --> sembra tutto inutile...

mod_prova_ind <- geeglm(
  cdrsb ~ (sex + age + bmi + ab_base + tau_base + adl_num) * year,
  family = binomial(link = "logit"),
  data = alz_long,
  id = sample,
  waves = year_seq,
  corstr = "independence"
)
summary(mod_prova_ind)  ## --> sembra tutto inutile...

## Sembra che bmi sia utile sia come baseline che nel tempo, e sex nella baseline
## Però solo prima impressione

## Penso si debba fare Wald test per vedere di ridurre il modello
## (ha senso?)


#### 1b. ALTERNATING LOGISTIC REGRESSION ####

## Cerca di capire che pacchetto usare per farla
## (forse va citato nelle references se usiamo R)


#### 2. MIXED EFFECTS MODEL ####

## I random effects vanno solo nella baseline? Anche temporal slope?
## Per ora farei solo random intercept

## Veramente prima prova a caso  --> controlla se ci sono funzioni migliori

mod_mixed <- glmer(
  cdrsb ~ (sex + age + bmi + ab_base + tau_base + adl_num) * year +
    (1 | sample),
  data = alz_long,
  family = binomial(link = "logit")
)

## Cosa fare poi? Confrontare i modelli? Confrontare la previsione?

## --> studio empirical bayes e shrinkage