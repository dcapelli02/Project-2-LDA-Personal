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
library(ORTH.Ord)

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
print(cor_matrix_cdrsb)


### SPAGHETTI PLOT
## For the spaghetti plot we need a random sample

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
  geom_line(alpha = 1, show.legend = FALSE, size = 0.5) +
  theme_bw() +
  labs(title = "Time Evolution of CDRSB")



#### 1. MARGINAL MODEL ####

### GEE RESULT

seq <- seq(0, 6, 0.01)
sex_prova <- 0
bmi_prova <- mean(stdize(alz$bmi))

coeff <- c(-0.88, 0.11, -0.018, 0.062, 0.243, -0.04)   ## data from Guillem code

const <- exp(coeff[1] + coeff[4]*(1-sex_prova) + coeff[2]*bmi_prova + coeff[3]*mean(stdize(alz$age)) +
                 seq * (coeff[5] + coeff[6]*bmi_prova))
theta <- const / (1 + const)

df_mean <- data.frame(seq = seq, theta = theta)

# Sovrapposizione con il dataset reale
ggplot() +
  stat_summary(data = alz_long, aes(x = year, y = cdrsb),
               fun = mean, geom = "line", size = 1.2, color = "black") +
  geom_line(data = df_mean, aes(x = seq, y = theta),
            size = 1.3, color = "red") +
  labs(
    title = "Dati reali + profilo medio della simulazione",
    x = "Anno",
    y = "CDRSB / Expr"
  ) +
  theme_minimal()





#### 2. RANDOM EFFECTS ####


### LOGIT RESULTS

seq <- seq(0, 6, 0.01)
sex_prova <- 0
bmi_prova <- mean(stdize(alz$bmi))

coeff <- c(-0.9851, -0.1973, 0.1212, 0.26, 0.0452, -0.04849)

calc_expr <- function(offset) {
  const <- exp(coeff[1] + coeff[2]*(1-sex_prova) + coeff[3]*bmi_prova + offset +
                 seq * (coeff[4] + coeff[5]*(1-sex_prova) + coeff[6]*bmi_prova))
  const / (1 + const)
}

n_offsets <- 500
offset_values <- rnorm(n_offsets, mean = 0, sd = 0.14)

df <- lapply(offset_values, function(off) {
  tibble(
    seq = seq,
    expr = calc_expr(off),
    offset = off
  )
}) %>% bind_rows()

# Mean profile
df_mean <- df %>%
  group_by(seq) %>%
  summarise(expr_mean = mean(expr), .groups = "drop")

# Plot con molte curve + profilo medio
ggplot() +
  geom_line(data = df, aes(x = seq, y = expr, group = offset),
            alpha = 0.1, color = "steelblue") +  # molte curve: trasparenza
  geom_line(data = df_mean, aes(x = seq, y = expr_mean),
            size = 1.3, color = "red") +         # profilo medio in rosso
  labs(
    title = "Curve logistiche con molti offset e profilo medio",
    x = "Sequenza",
    y = "Expr"
  ) +
  theme_minimal()

# Sovrapposizione con il dataset reale
ggplot() +
  stat_summary(data = alz_long, aes(x = year, y = cdrsb),
               fun = mean, geom = "line", size = 1.2, color = "black") +
  geom_line(data = df, aes(x = seq, y = expr, group = offset),
            alpha = 0.1, color = "steelblue") +  # molte curve: trasparenza
  geom_line(data = df_mean, aes(x = seq, y = expr_mean),
            size = 1.3, color = "red") +
  labs(
    title = "Dati reali + profilo medio della simulazione",
    x = "Anno",
    y = "CDRSB / Expr"
  ) +
  theme_minimal()



### PROBIT


seq <- seq(0, 6, 0.01)
sex_prova <- 0
bmi_prova <- mean(stdize(alz$bmi))

coeff <- c(-0.61, -0.11, 0.07, 0.16, 0.026, -0.029, 0.09)

calc_expr_probit <- function(offset) {
  eta <- coeff[1] + coeff[2]*(1-sex_prova) + coeff[3]*bmi_prova + offset +
    seq * (coeff[4] + coeff[5]*(1-sex_prova) + coeff[6]*bmi_prova)
  
  pnorm(eta)   # <-- probit invece della logit
}

# Numero più grande di offset
n_offsets <- 500
offset_values <- rnorm(n_offsets, mean = 0, sd = 0.14)

# Dataframe con tutte le curve
df <- lapply(offset_values, function(off) {
  tibble(
    seq = seq,
    expr = calc_expr_probit(off),
    offset = off
  )
}) %>% bind_rows()

# Calcolo del profilo medio nel tempo
df_mean <- df %>%
  group_by(seq) %>%
  summarise(expr_mean = mean(expr), .groups = "drop")

# Plot con molte curve + profilo medio
ggplot() +
  geom_line(data = df, aes(x = seq, y = expr, group = offset),
            alpha = 0.1, color = "steelblue") +  # molte curve: trasparenza
  geom_line(data = df_mean, aes(x = seq, y = expr_mean),
            size = 1.3, color = "red") +         # profilo medio in rosso
  labs(
    title = "Curve logistiche con molti offset e profilo medio",
    x = "Sequenza",
    y = "Expr"
  ) +
  theme_minimal()

# Sovrapposizione con il dataset reale
ggplot() +
  stat_summary(data = alz_long, aes(x = year, y = cdrsb),
               fun = mean, geom = "line", size = 1.2, color = "black") +
  geom_line(data = df, aes(x = seq, y = expr, group = offset),
            alpha = 0.1, color = "steelblue") +  # molte curve: trasparenza
  geom_line(data = df_mean, aes(x = seq, y = expr_mean),
            size = 1.3, color = "red") +
  labs(
    title = "Dati reali + profilo medio della simulazione",
    x = "Anno",
    y = "CDRSB / Expr"
  ) +
  theme_minimal()
