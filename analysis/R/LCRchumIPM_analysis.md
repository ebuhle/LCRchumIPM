---
title: 'Fitting Integrated Population Models to Lower Columbia River Chum Salmon Monitoring Data'
author: "Eric Buhle, Kale Bentley, Thomas Buehrens, Mark Scheuerell, Todd Hillson, and ..."
date: "June 25, 2020"
output:
  html_document:
    keep_md: true
    df_print: paged
    fig_caption: true
    toc: true
    toc_float: true
  word_document:
    toc: true
  pdf_document:
    toc: true
---





# Overview

Background on IPMs, outline of **salmonIPM**...

# Setup and data

Load the packages we'll need...


```r
options(device = ifelse(.Platform$OS.type == "windows", "windows", "quartz"))
options(mc.cores = parallel::detectCores(logical = FALSE) - 1)

library(salmonIPM)
library(rstan)
library(shinystan)
library(matrixStats)
library(tibble)
library(dplyr)
library(reshape2)
library(yarrr)
library(corrplot)
library(magicaxis)
library(zoo)
library(here)

if(file.exists(here("analysis","results","LCRchumIPM.RData")))
  load(here("analysis","results","LCRchumIPM.RData"))
```

Read in and manipulate the data...


```r
# Mapping of location to population
location_pop <- read.csv(here("data","Location.Reach_Population.csv"), 
                         header = TRUE, stringsAsFactors = TRUE) %>% 
  rename(strata = Strata, location = Location.Reach, pop1 = Population1, pop2 = Population2)

# Mapping of disposition to hatchery vs. wild (i.e., broodstock vs. natural spawner)
disposition_HW <- read.csv(here("data","Disposition_HW.csv"), 
                           header = TRUE, stringsAsFactors = TRUE) %>% 
  rename(disposition = Disposition) %>% arrange(HW)

# Start dates of hatcheries associated with populations
hatcheries <- read.csv(here("data","Hatchery_Programs.csv"), header = TRUE, stringsAsFactors = TRUE)

# Spawner abundance data
# Assumptions:
# (1) NAs in hatchery dispositions (incl. Duncan Channel) are really zeros
# (2) NAs in Duncan Creek from 2004-present are really zeros
# (3) All other NAs are real missing observations
spawner_data <- read.csv(here("data","Data_ChumSpawnerAbundance_2019-12-12.csv"), 
                         header = TRUE, stringsAsFactors = TRUE) %>% 
  rename(year = Return.Yr., strata = Strata, location = Location.Reach, 
         disposition = Disposition, method = Method, S_obs = Abund.Mean, SD = Abund.SD) %>% 
  mutate(pop = location_pop$pop2[match(location, location_pop$location)],
         disposition_HW = disposition_HW$HW[match(disposition, disposition_HW$disposition)],
         S_obs = replace(S_obs, is.na(S_obs) & disposition_HW == "H", 0),
         S_obs = replace(S_obs, is.na(S_obs) & pop == "Duncan_Creek" & year >= 2004, 0)) %>% 
  select(year:location, pop, disposition, disposition_HW, method:SD) %>% 
  arrange(strata, location, year)

names_S_obs <- disposition_HW$disposition
names_B_take_obs <- disposition_HW$disposition[disposition_HW$HW == "H"]

spawner_data_agg <- aggregate(S_obs ~ year + strata + pop + disposition,
                              data = spawner_data, FUN = sum, na.action = na.pass) %>%
  dcast(year + strata + pop ~ disposition, value.var = "S_obs", 
        fun.aggregate = identity, fill = 0) %>% 
  add_column(S_obs = rowSums(select(., all_of(names_S_obs))),
             B_take_obs = rowSums(select(., all_of(names_B_take_obs)))) %>% 
  select(-all_of(names_S_obs)) %>% arrange(strata, pop, year)

# Spawner age-, sex-, and origin-frequency (aka BioData)
bio_data <- read.csv(here("data","Data_ChumSpawnerBioData_2019-12-12.csv"), 
                     header = TRUE, stringsAsFactors = TRUE) %>% 
  rename(year = Return.Yr., strata = Strata, location = Location.Reach, 
         disposition = Disposition, origin = Origin, sex = Sex, age = Age, count = Count) %>% 
  mutate(pop = location_pop$pop2[match(location, location_pop$location)],
         origin_HW = ifelse(origin == "Natural_spawner", "W", "H"),
         count = ifelse(is.na(count), 0, count)) %>% 
  select(year:location, pop, disposition, origin, origin_HW, sex:count) %>%
  arrange(strata, location, year, origin, age, sex)

# age of wild spawners only
bio_data_age <- bio_data %>% filter(origin_HW == "W") %>%  
  dcast(year + strata + pop ~ age, value.var = "count", fun.aggregate = sum)

bio_data_origin <- bio_data %>% 
  dcast(year + strata + pop ~ origin_HW, value.var = "count", fun.aggregate = sum)

# Juvenile abundance data
# Assumptions:
# (1) Smolts from Duncan Channel represent all naturally produced offspring of spawners
#     in Duncan Creek (hence Duncan_Channel -> Duncan_Creek in location_pop)
# (2) Duncan_North + Duncan_South = Duncan_Channel, so the former two are redundant 
#     (not really an assumption, although the equality isn't perfect in all years)
juv_data <- read.csv(here("data", "Data_ChumJuvenileAbundance_2020-06-09.csv"), 
                     header = TRUE, stringsAsFactors = TRUE) %>% 
  rename(brood_year = Brood.Year, year = Outmigration.Year, strata = Strata, 
         location = Location.Reach, origin = Origin, trap_type = TrapType, 
         analysis = Analysis, partial_spawners = Partial.Spawners, raw_catch = RawCatch,
         M_obs = Abund_Mean, SD = Abund_SD, L95 = Abund_L95, U95 = Abund_U95, CV = Abund_CV,
         comments = Comments) %>% 
  mutate(pop = location_pop$pop2[match(location, location_pop$location)]) %>% 
  select(brood_year:location, pop, origin:comments) %>% arrange(strata, location, year)

# drop hatchery or redundant pops and cases with leading or trailing NAs in M_obs
head_noNA <- function(x) { cumsum(!is.na(x)) > 0 }
juv_data_incl <- juv_data %>% filter(pop %in% spawner_data$pop) %>% 
  mutate(location = factor(location), pop = factor(pop, levels = levels(spawner_data$pop))) %>% 
  group_by(pop) %>% filter(head_noNA(M_obs) & rev(head_noNA(rev(M_obs)))) %>% as.data.frame()

# Fish data formatted for salmonIPM
# Drop age-2 and age-6 samples (each is < 0.1% of aged spawners)
# Use A = 1 for now (so Rmax in units of spawners)
fish_data <- full_join(spawner_data_agg, bio_data_age, by = c("year","strata","pop")) %>% 
  full_join(bio_data_origin, by = c("year","strata","pop")) %>% 
  full_join(juv_data_incl, by = c("year","strata","pop")) %>%
  mutate(B_take_obs = replace(B_take_obs, is.na(B_take_obs), 0)) %>% 
  rename_at(vars(contains("Age-")), list(~ paste0(sub("Age-","n_age",.), "_obs"))) %>% 
  select(-c(n_age2_obs, n_age6_obs)) %>% 
  rename(n_H_obs = H, n_W_obs = W) %>% mutate(A = 1, fit_p_HOS = NA, F_rate = 0) %>% 
  mutate_at(vars(contains("n_")), list(~ replace(., is.na(.), 0))) %>% 
  select(strata, pop, year, A, S_obs, M_obs, n_age3_obs:n_W_obs, 
         fit_p_HOS, B_take_obs, F_rate) %>% arrange(strata, pop, year) 

# fill in fit_p_HOS
for(i in 1:nrow(fish_data)) {
  pop_i <- as.character(fish_data$pop[i])
  start_year <- ifelse(pop_i %in% hatcheries$pop,
                       min(hatcheries$start_brood_year[hatcheries$pop == pop_i]) + 1,
                       NA)
  fish_data$fit_p_HOS[i] <- ifelse((!is.na(start_year) & fish_data$year[i] >= start_year) |
                                     fish_data$n_H_obs[i] > 0, 1, 0)
}

# # drop cases with initial NAs in S_obs unless bio data is present
# fish_data <- fish_data %>% mutate(n_age = rowSums(select(., n_age2_obs:n_age6_obs))) %>% 
#   group_by(pop) %>%  filter(head_noNA(S_obs) | cumsum(n_age) > 0) %>% 
#   select(-n_age) %>% as.data.frame()

# subsets for models with specific stage structure
# spawner-spawner: drop cases with initial NAs in S_obs, even if bio data is present
fish_data_SS <- fish_data %>% group_by(pop) %>% filter(head_noNA(S_obs)) %>% as.data.frame()
# spawner-spawner: drop cases with initial NAs in M_obs, even if bio data is present
fish_data_SMS <- fish_data %>% group_by(pop) %>% 
  filter(head_noNA(S_obs) | head_noNA(M_obs)) %>% as.data.frame()
```


# Retrospective models

Fit two-stage spawner-smolt-spawner models and explore output...

Density-independent


```r
SMS_exp <- salmonIPM(fish_data = fish_data_SMS, ages = list(M = 1),
                     stan_model = "IPM_SMS_pp", SR_fun = "exp",
                     pars = c("mu_alpha","sigma_alpha","alpha",
                              "mu_Rmax","sigma_Rmax","Rmax","rho_alphaRmax",
                              "beta_phi_M","rho_phi_M","sigma_phi_M","phi_M","sigma_M",
                              "mu_MS","beta_phi_MS","rho_phi_MS",
                              "sigma_phi_MS","phi_MS","sigma_MS","s_MS",
                              "mu_p","sigma_gamma","R_gamma","sigma_p","R_p","p",
                              "tau_M","tau_S","p_HOS","M","S","q","LL"),
                     chains = 3, iter = 1500, warmup = 500,
                     control = list(adapt_delta = 0.99, max_treedepth = 13))
```

```r
print(SMS_exp, prob = c(0.025,0.5,0.975),
      pars = c("alpha","Rmax","phi_M","phi_MS","p","p_HOS","S","M","s_MS","q","LL"), 
      include = FALSE, use_cache = FALSE)
```

```
Inference for Stan model: IPM_SMS_pp.
3 chains, each with iter=1500; warmup=500; thin=1; 
post-warmup draws per chain=1000, total post-warmup draws=3000.

                    mean se_mean    sd      2.5%       50%     97.5% n_eff Rhat
mu_alpha            6.68    0.00  0.17      6.35      6.68      7.01  1170 1.00
sigma_alpha         0.26    0.00  0.10      0.10      0.24      0.50   902 1.00
mu_Rmax             0.01    0.15 10.17    -20.21     -0.08     19.65  4648 1.00
sigma_Rmax          1.49    0.03  0.86      0.07      1.48      2.99   984 1.00
rho_alphaRmax       0.00    0.01  0.58     -0.95      0.00      0.96  4713 1.00
rho_phi_M           0.10    0.01  0.40     -0.68      0.13      0.77  1082 1.00
sigma_phi_M         0.43    0.01  0.19      0.07      0.42      0.85   561 1.01
sigma_M             0.11    0.00  0.07      0.01      0.10      0.27   331 1.00
mu_MS               0.00    0.00  0.00      0.00      0.00      0.00  1470 1.00
rho_phi_MS          0.42    0.01  0.29     -0.27      0.47      0.82   918 1.00
sigma_phi_MS        1.02    0.01  0.25      0.64      0.99      1.61  1153 1.00
sigma_MS            0.19    0.01  0.09      0.01      0.19      0.37   243 1.01
mu_p[1]             0.19    0.00  0.01      0.16      0.18      0.22  1162 1.00
mu_p[2]             0.75    0.00  0.01      0.72      0.75      0.78  1186 1.00
mu_p[3]             0.06    0.00  0.01      0.05      0.06      0.07  1175 1.00
sigma_gamma[1]      0.23    0.01  0.15      0.02      0.21      0.57   658 1.00
sigma_gamma[2]      0.14    0.00  0.10      0.01      0.12      0.39   664 1.00
R_gamma[1,1]        1.00     NaN  0.00      1.00      1.00      1.00   NaN  NaN
R_gamma[1,2]        0.29    0.02  0.54     -0.90      0.41      0.98  1013 1.00
R_gamma[2,1]        0.29    0.02  0.54     -0.90      0.41      0.98  1013 1.00
R_gamma[2,2]        1.00    0.00  0.00      1.00      1.00      1.00   438 1.00
sigma_p[1]          0.86    0.00  0.11      0.65      0.85      1.08   701 1.01
sigma_p[2]          0.43    0.01  0.10      0.23      0.43      0.63   286 1.01
R_p[1,1]            1.00     NaN  0.00      1.00      1.00      1.00   NaN  NaN
R_p[1,2]            0.58    0.01  0.19      0.16      0.61      0.84   542 1.01
R_p[2,1]            0.58    0.01  0.19      0.16      0.61      0.84   542 1.01
R_p[2,2]            1.00    0.00  0.00      1.00      1.00      1.00   197 1.00
tau_M               0.68    0.01  0.13      0.46      0.67      0.95   483 1.01
tau_S               0.90    0.00  0.07      0.77      0.90      1.05   936 1.00
lp__           -22491.99    1.26 30.20 -22552.91 -22491.06 -22434.93   575 1.01

Samples were drawn using NUTS(diag_e) at Fri Jun 26 14:47:32 2020.
For each parameter, n_eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor on split chains (at 
convergence, Rhat=1).
```


Beverton-Holt


```r
SMS_BH <- salmonIPM(fish_data = fish_data_SMS, ages = list(M = 1),
                    stan_model = "IPM_SMS_pp", SR_fun = "BH",
                    pars = c("mu_alpha","sigma_alpha","alpha",
                             "mu_Rmax","sigma_Rmax","Rmax","rho_alphaRmax",
                             "beta_phi_M","rho_phi_M","sigma_phi_M","phi_M","sigma_M",
                             "mu_MS","beta_phi_MS","rho_phi_MS",
                             "sigma_phi_MS","phi_MS","sigma_MS","s_MS",
                             "mu_p","sigma_gamma","R_gamma","sigma_p","R_p","p",
                             "tau_M","tau_S","p_HOS","M","S","q","LL"),
                    chains = 3, iter = 1500, warmup = 500,
                    control = list(adapt_delta = 0.99, max_treedepth = 13))
```

```r
print(SMS_BH, prob = c(0.025,0.5,0.975),
      pars = c("alpha","Rmax","phi_M","phi_MS","p","p_HOS","S","M","s_MS","q","LL"), 
      include = FALSE, use_cache = FALSE)
```

```
Inference for Stan model: IPM_SMS_pp.
3 chains, each with iter=1500; warmup=500; thin=1; 
post-warmup draws per chain=1000, total post-warmup draws=3000.

                    mean se_mean    sd      2.5%       50%     97.5% n_eff Rhat
mu_alpha           10.20    0.18  1.96      8.07      9.59     15.86   121 1.05
sigma_alpha         1.31    0.03  0.73      0.18      1.17      2.84   821 1.00
mu_Rmax            12.48    0.01  0.47     11.53     12.49     13.34  1195 1.01
sigma_Rmax          1.47    0.01  0.37      0.91      1.42      2.36  1342 1.00
rho_alphaRmax       0.40    0.04  0.46     -0.67      0.50      0.98   148 1.03
rho_phi_M           0.08    0.01  0.46     -0.77      0.11      0.80  2392 1.00
sigma_phi_M         0.13    0.00  0.10      0.01      0.10      0.37   832 1.01
sigma_M             0.10    0.00  0.06      0.00      0.10      0.23   235 1.01
mu_MS               0.00    0.00  0.00      0.00      0.00      0.00  1596 1.00
rho_phi_MS          0.33    0.01  0.31     -0.34      0.36      0.81   878 1.00
sigma_phi_MS        1.14    0.01  0.24      0.77      1.11      1.70  1219 1.00
sigma_MS            0.15    0.01  0.07      0.02      0.16      0.28   180 1.01
mu_p[1]             0.18    0.00  0.02      0.15      0.18      0.21  1662 1.00
mu_p[2]             0.76    0.00  0.02      0.72      0.76      0.79  1736 1.00
mu_p[3]             0.07    0.00  0.00      0.06      0.07      0.07  1032 1.01
sigma_gamma[1]      0.28    0.00  0.14      0.04      0.27      0.61   883 1.00
sigma_gamma[2]      0.12    0.00  0.09      0.00      0.10      0.33  1012 1.00
R_gamma[1,1]        1.00     NaN  0.00      1.00      1.00      1.00   NaN  NaN
R_gamma[1,2]        0.35    0.01  0.51     -0.84      0.49      0.98  1848 1.00
R_gamma[2,1]        0.35    0.01  0.51     -0.84      0.49      0.98  1848 1.00
R_gamma[2,2]        1.00    0.00  0.00      1.00      1.00      1.00   328 1.00
sigma_p[1]          0.63    0.00  0.12      0.41      0.62      0.86   575 1.01
sigma_p[2]          0.37    0.01  0.10      0.18      0.37      0.58   193 1.05
R_p[1,1]            1.00     NaN  0.00      1.00      1.00      1.00   NaN  NaN
R_p[1,2]           -0.07    0.03  0.41     -0.90     -0.04      0.61   162 1.04
R_p[2,1]           -0.07    0.03  0.41     -0.90     -0.04      0.61   162 1.04
R_p[2,2]            1.00    0.00  0.00      1.00      1.00      1.00  2182 1.00
tau_M               0.67    0.00  0.10      0.50      0.66      0.88  1357 1.00
tau_S               0.96    0.00  0.08      0.81      0.96      1.15   736 1.01
lp__           -22509.65    1.25 30.40 -22570.93 -22509.37 -22450.78   590 1.00

Samples were drawn using NUTS(diag_e) at Fri Jun 26 17:20:34 2020.
For each parameter, n_eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor on split chains (at 
convergence, Rhat=1).
```

Ricker


```r
SMS_Ricker <- salmonIPM(fish_data = fish_data_SMS, ages = list(M = 1),
                        stan_model = "IPM_SMS_pp", SR_fun = "Ricker",
                        pars = c("mu_alpha","sigma_alpha","alpha",
                                 "mu_Rmax","sigma_Rmax","Rmax","rho_alphaRmax",
                                 "beta_phi_M","rho_phi_M","sigma_phi_M","phi_M","sigma_M",
                                 "mu_MS","beta_phi_MS","rho_phi_MS",
                                 "sigma_phi_MS","phi_MS","sigma_MS","s_MS",
                                 "mu_p","sigma_gamma","R_gamma","sigma_p","R_p","p",
                                 "tau_M","tau_S","p_HOS","M","S","q","LL"),
                        chains = 3, iter = 1500, warmup = 500,
                        control = list(adapt_delta = 0.99, max_treedepth = 13))
```

```r
print(SMS_Ricker, prob = c(0.025,0.5,0.975),
      pars = c("alpha","Rmax","phi_M","phi_MS","p","p_HOS","S","M","s_MS","q","LL"), 
      include = FALSE, use_cache = FALSE)
```

```
Inference for Stan model: IPM_SMS_pp.
3 chains, each with iter=1500; warmup=500; thin=1; 
post-warmup draws per chain=1000, total post-warmup draws=3000.

                    mean se_mean    sd      2.5%       50%     97.5% n_eff Rhat
mu_alpha            7.31    0.01  0.20      6.92      7.31      7.70   527 1.02
sigma_alpha         0.23    0.00  0.10      0.07      0.22      0.46   987 1.00
mu_Rmax            13.04    0.02  0.45     12.13     13.05     13.87   576 1.01
sigma_Rmax          1.25    0.01  0.35      0.71      1.19      2.06   859 1.00
rho_alphaRmax       0.68    0.01  0.28     -0.06      0.76      0.99   615 1.00
rho_phi_M           0.03    0.01  0.44     -0.77      0.05      0.78  1999 1.00
sigma_phi_M         0.16    0.00  0.13      0.01      0.14      0.47   775 1.01
sigma_M             0.08    0.00  0.05      0.00      0.07      0.19   377 1.01
mu_MS               0.00    0.00  0.00      0.00      0.00      0.00  1025 1.00
rho_phi_MS          0.33    0.01  0.30     -0.33      0.34      0.80   599 1.00
sigma_phi_MS        1.08    0.01  0.25      0.69      1.04      1.68   967 1.00
sigma_MS            0.11    0.00  0.06      0.01      0.11      0.24   256 1.02
mu_p[1]             0.18    0.00  0.01      0.15      0.18      0.21  1625 1.00
mu_p[2]             0.75    0.00  0.01      0.72      0.75      0.78  1651 1.00
mu_p[3]             0.06    0.00  0.01      0.05      0.06      0.07  1481 1.00
sigma_gamma[1]      0.23    0.01  0.13      0.02      0.22      0.53   681 1.00
sigma_gamma[2]      0.12    0.00  0.09      0.01      0.11      0.35   744 1.00
R_gamma[1,1]        1.00     NaN  0.00      1.00      1.00      1.00   NaN  NaN
R_gamma[1,2]        0.26    0.02  0.54     -0.88      0.36      0.97  1071 1.00
R_gamma[2,1]        0.26    0.02  0.54     -0.88      0.36      0.97  1071 1.00
R_gamma[2,2]        1.00    0.00  0.00      1.00      1.00      1.00  1549 1.00
sigma_p[1]          0.77    0.01  0.11      0.56      0.77      1.00   511 1.01
sigma_p[2]          0.43    0.01  0.10      0.23      0.43      0.61   298 1.01
R_p[1,1]            1.00     NaN  0.00      1.00      1.00      1.00   NaN  NaN
R_p[1,2]            0.37    0.01  0.24     -0.17      0.41      0.73   361 1.01
R_p[2,1]            0.37    0.01  0.24     -0.17      0.41      0.73   361 1.01
R_p[2,2]            1.00    0.00  0.00      1.00      1.00      1.00   422 1.00
tau_M               0.70    0.00  0.10      0.52      0.69      0.93  1056 1.00
tau_S               0.91    0.00  0.07      0.79      0.90      1.05   755 1.00
lp__           -22496.44    1.09 29.78 -22556.05 -22496.22 -22438.51   751 1.00

Samples were drawn using NUTS(diag_e) at Fri Jun 26 15:15:42 2020.
For each parameter, n_eff is a crude measure of effective sample size,
and Rhat is the potential scale reduction factor on split chains (at 
convergence, Rhat=1).
```

Model comparison based on LOO (unhelpful because Pareto ks are too high)


```r
LL_SMS <- lapply(list(exp = SMS_exp, BH = SMS_BH, Ricker = SMS_Ricker),
             loo::extract_log_lik, parameter_name = "LL", merge_chains = FALSE)

# Relative ESS of posterior draws of observationwise likelihood 
r_eff_SMS <- lapply(LL_SMS, function(x) relative_eff(exp(x)))

# PSIS-LOO
LOO_SMS <- lapply(1:length(LL_SMS), function(i) loo(LL_SMS[[i]], r_eff = r_eff_SMS[[i]]))
```

```
Warning: Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.

Warning: Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.

Warning: Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.
```

```r
names(LOO_SMS) <- names(LL_SMS)

## Compare all three models
loo_compare(LOO_SMS)
```

```
       elpd_diff se_diff
Ricker   0.0       0.0  
exp     -7.5       7.0  
BH     -20.1       8.8  
```

```r
## Exponential vs. Ricker
loo_compare(LOO_SMS[c("exp","Ricker")])
```

```
       elpd_diff se_diff
Ricker  0.0       0.0   
exp    -7.5       7.0   
```

```r
## Exponential vs. Beverton-Holt
loo_compare(LOO_SMS[c("exp","BH")])
```

```
    elpd_diff se_diff
exp   0.0       0.0  
BH  -12.6      11.1  
```

```r
## Beverton-Holt vs. Ricker
loo_compare(LOO_SMS[c("BH","Ricker")])
```

```
       elpd_diff se_diff
Ricker   0.0       0.0  
BH     -20.1       8.8  
```

Plot estimated spawner-smolt production curves and parameters for the Beverton-Holt model.


<img src="LCRchumIPM_analysis_files/figure-html/plot_SR_params-1.png" width="50%" style="display: block; margin: auto;" />

Now do the same for the Ricker model.


<img src="LCRchumIPM_analysis_files/figure-html/plot_SR_params-1.png" width="50%" style="display: block; margin: auto;" />

The Ricker model is more biologically plausible, so let's proceed with that model for now. Here are the fits to the spawner data:


<img src="LCRchumIPM_analysis_files/figure-html/plot_spawner_smolt_ts-1.png" width="90%" style="display: block; margin: auto;" />

And here are the fits to the much sparser smolt data:


<img src="LCRchumIPM_analysis_files/figure-html/plot_spawner_smolt_ts-1.png" width="90%" style="display: block; margin: auto;" />

We can examine how the model partitions shared interannual fluctuations between the two life-stage transitions...


<img src="LCRchumIPM_analysis_files/figure-html/plot_phi-1.png" width="50%" style="display: block; margin: auto;" />


