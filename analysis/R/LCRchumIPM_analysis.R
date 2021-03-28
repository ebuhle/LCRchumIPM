## @knitr getting_started
options(device = ifelse(.Platform$OS.type == "windows", "windows", "quartz"))
options(mc.cores = parallel::detectCores(logical = FALSE) - 1)

library(salmonIPM)
library(rstan)
library(shinystan)
library(matrixStats)
library(Hmisc)
library(tibble)
library(dplyr)
library(tidyr)
library(reshape2)
library(yarrr)
library(magicaxis)
library(viridis)
library(zoo)
library(here)

if(file.exists(here("analysis","results","LCRchumIPM.RData")))
  load(here("analysis","results","LCRchumIPM.RData"))
## @knitr


#===========================================================================
# DATA
#===========================================================================

## @knitr data

# Start dates of hatcheries associated with populations
hatcheries <- read.csv(here("data","Hatchery_Programs.csv"), header = TRUE, stringsAsFactors = TRUE)

# Spawner abundance data
# Assumptions:
# (0) Fix coding error in data that assigns some Duncan Creek rows to Cascade stratum
# (1) NAs in hatchery dispositions are really zeros
# (2) NAs in Duncan Creek and Duncan Channel are really zeros
# (3) All other NAs are real missing observations
# (4) When calculating the observation error of log(S_obs), tau_S_obs, assume
#     Abund.Mean and Abund.SD are the mean and SD of a lognormal posterior distribution
#     of spawner abundance based on the sample
spawner_data <- read.csv(here("data","Data_ChumSpawnerAbundance_2019-12-12.csv"), 
                         header = TRUE, stringsAsFactors = FALSE) %>% 
  rename(year = Return.Yr., strata = Strata, location = Location.Reach, 
         disposition = Disposition, method = Method, S_obs = Abund.Mean, SD = Abund.SD) %>% 
  mutate(strata = replace(strata, disposition == "Duncan_Channel" & strata != "Gorge", "Gorge"),
         S_obs = replace(S_obs, is.na(S_obs) & grepl("Hatchery|Duncan", disposition), 0),
         tau_S_obs = sqrt(log((SD/S_obs)^2 + 1))) %>% 
  select(year:location, disposition, method, S_obs, SD, tau_S_obs) %>% 
  arrange(strata, location, year)

# broodstock take: 
# all spawners taken from a given location to a different disposition
# (summarized by location)
broodstock_data <- spawner_data %>% group_by(strata, location, year) %>% 
  summarize(B_take_obs = sum(S_obs[disposition != location])) %>% 
  rename(pop = location) %>% as.data.frame()

# total spawners:
# all spawners with a given disposition, regardless of original return location
# (summarized by disposition)
spawner_data_agg <- spawner_data %>% group_by(strata, disposition, year) %>% 
  summarize(S_obs = sum(S_obs), tau_S_obs = unique(tau_S_obs)) %>% 
  rename(pop = disposition) %>% left_join(broodstock_data, by = c("strata","pop","year")) %>% 
  mutate(B_take_obs = replace(B_take_obs, is.na(B_take_obs), 0)) %>% as.data.frame()

# Spawner age-, sex-, and origin-frequency (aka BioData)
# Assumptions:
# (0) Fix coding error in data that assigns some Duncan Creek rows to Cascade stratum
# (1) The current generation of salmonIPM models must treat any known nonlocal-origin
#     spawner as "hatchery-origin" to avoid counting as local recruitment, so "H"
#     includes true hatchery fish based on origin *plus* any whose known origin 
#     or return location do not match their disposition *unless* they are NOR or
#     Duncan Channel fish returning to Duncan Creek and disposed to Duncan Channel,
#     which are in fact local recruitment.
bio_data <- read.csv(here("data","Data_ChumSpawnerBioData_2019-12-12.csv"), 
                     header = TRUE, stringsAsFactors = FALSE) %>% 
  rename(year = Return.Yr., strata = Strata, location = Location.Reach, 
         disposition = Disposition, origin = Origin, sex = Sex, age = Age, count = Count) %>% 
  mutate(strata = replace(strata, disposition == "Duncan_Channel" & strata != "Gorge", "Gorge"),
         count = replace(count, is.na(count), 0),
         HW = ifelse((grepl("Hatchery", origin) | location != disposition | 
                        (origin != disposition & origin != "Natural_spawner")) &
                       !(origin %in% c("Duncan_Channel","Natural_spawner") & 
                           location == "Duncan_Creek" & disposition == "Duncan_Channel"), 
                     "H", "W")) %>% 
  select(year:location, disposition, origin, HW, sex:count) %>%
  arrange(strata, location, year, origin, age, sex)

# age of wild (i.e., potentially local) spawners only
bio_data_age <- bio_data %>% filter(HW == "W") %>%
  dcast(year + strata + disposition ~ age, value.var = "count", fun.aggregate = sum) %>% 
  rename(pop = disposition)

# H/W (nonlocal/potentially local)  
bio_data_HW <- bio_data %>%
  dcast(year + strata + disposition ~ HW, value.var = "count", fun.aggregate = sum) %>% 
  rename(pop = disposition)

# Juvenile abundance data
# Assumptions:
# (1) Duncan_North + Duncan_South = Duncan_Channel, so the former two are redundant 
#     (not really an assumption, although the equality isn't perfect in all years)
# (2) When calculating the observation error of log(M_obs), tau_M_obs, assume
#     Abund_Median and Abund_SD are the median and SD of a lognormal posterior 
#     distribution of smolt abundance based on the sample
# (3) If Abund_SD == 0 (when Analysis=="Census": some years in Duncan_Channel and 
#     Hamilton_Channel) treat as NA
juv_data <- read.csv(here("data", "Data_ChumJuvenileAbundance_2020-06-09.csv"), 
                     header = TRUE, stringsAsFactors = FALSE) %>% 
  rename(brood_year = Brood.Year, year = Outmigration.Year, strata = Strata, 
         location = Location.Reach, origin = Origin, trap_type = TrapType, 
         analysis = Analysis, partial_spawners = Partial.Spawners, raw_catch = RawCatch,
         M_obs = Abund_Median, mean = Abund_Mean, SD = Abund_SD, 
         L95 = Abund_L95, U95 = Abund_U95, CV = Abund_CV, comments = Comments) %>% 
  mutate(tau_M_obs = replace(sqrt(log((SD/mean)^2 + 1)), SD==0, NA)) %>% 
  select(strata, location, year, brood_year, origin:CV, tau_M_obs, comments) %>% 
  arrange(strata, location, year)

# drop hatchery or redundant pops and cases with leading or trailing NAs in M_obs
# use only pooled Duncan Channel data for now, not North / South
head_noNA <- function(x) { cumsum(!is.na(x)) > 0 }
juv_data_incl <- juv_data %>% group_by(location) %>% 
  filter(head_noNA(M_obs) & rev(head_noNA(rev(M_obs)))) %>%
  filter(!(location %in% c("Duncan_North","Duncan_South"))) %>% 
  rename(pop = location) %>% as.data.frame()

# Fish data formatted for salmonIPM
# Drop age-2 and age-6 samples (each is < 0.1% of aged spawners)
# Use A = 1 for now (so Rmax in units of spawners)
# Drop Duncan Creek and "populations" that are actually hatcheries
# Change S_obs and tau_S_obs to NA in Hamilton Channel 2011-2012 based on
# https://github.com/mdscheuerell/chumIPM/issues/6#issuecomment-807885445
# Pad data as necessary so Grays_MS, Grays_WF, and Grays_CJ have the same set of years
# (since their estimated smolts will be summed)  
fish_data <- full_join(spawner_data_agg, bio_data_age, by = c("strata","pop","year")) %>% 
  full_join(bio_data_HW, by = c("strata","pop","year")) %>% 
  full_join(juv_data_incl, by = c("strata","pop","year")) %>%
  full_join({ complete(select(filter(., strata == "Coastal"), c(strata,pop,year)), 
                       pop, year, fill = list(strata = "Coastal")) },
            by = c("strata","pop","year")) %>%
  rename_at(vars(contains("Age-")), list(~ paste0(sub("Age-","n_age",.), "_obs"))) %>% 
  select(-c(n_age2_obs, n_age6_obs)) %>% 
  rename(n_H_obs = H, n_W_obs = W) %>% mutate(A = 1, fit_p_HOS = NA, F_rate = 0) %>% 
  mutate_at(vars(contains("n_")), ~ replace(., is.na(.), 0)) %>% 
  filter(!grepl("Hatchery|Duncan_Creek", pop)) %>% 
  mutate(S_obs = replace(S_obs, pop == "Hamilton_Channel" & year %in% 2011:2012, NA),
         tau_S_obs = replace(tau_S_obs, pop == "Hamilton_Channel" & year %in% 2011:2012, NA),
         strata = factor(strata), pop = factor(pop), 
         B_take_obs = replace(B_take_obs, is.na(B_take_obs), 0)) %>%
  select(strata, pop, year, A, S_obs, tau_S_obs, M_obs, tau_M_obs, n_age3_obs:n_W_obs, 
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

# subsets for models with specific stage structure
# spawner-spawner: drop cases with initial NAs in S_obs, even if bio data is present
fish_data_SS <- fish_data %>% group_by(pop) %>% filter(head_noNA(S_obs)) %>% as.data.frame()
# spawner-smolt-spawner: drop cases with initial NAs in both S_obs and M_obs
fish_data_SMS <- fish_data %>% group_by(pop) %>% filter(head_noNA(S_obs) | head_noNA(M_obs)) %>% 
  add_column(downstream_trap = NA, .after = "tau_M_obs")  %>% as.data.frame()

# fish_data_SMS: 
# assign Grays_WF and Grays_CJ smolts to the downstream trap in Grays_MS
# where they will be counted (or double-counted, in the case of Grays_CJ),
# assuming no mortality between the upstream tributary and the downstream trap
fish_data_SMS <- fish_data_SMS %>%
  mutate(downstream_trap = replace(downstream_trap, pop %in% c("Grays_WF","Grays_CJ"),
                                   which(pop == "Grays_MS")))

# pad data with future years to generate forecasts
# use 5-year (1-generation) time horizon
fish_data_SMS_fore <- fish_data_SMS %>% group_by(pop) %>% 
  slice(rep(n(), max(fish_data_SMS$year) + 5 - max(year))) %>% 
  mutate(year = (unique(year) + 1):(max(fish_data_SMS$year) + 5),
         S_obs = NA, tau_S_obs = NA, M_obs = NA, tau_M_obs = NA,
         fit_p_HOS = 0, B_take_obs = 0, F_rate = 0) %>% 
  mutate_at(vars(starts_with("n_")), ~ 0) %>% 
  full_join(fish_data_SMS) %>% arrange(pop, year) %>% 
  mutate(forecast = year > max(fish_data_SMS$year), downstream_trap = NA) %>% 
  select(strata:year, forecast, A:F_rate) %>% as.data.frame()

# fish_data_SMS_fore: 
# assign Grays_WF and Grays_CJ smolts to the downstream trap in Grays_MS
# (need to do this again b/c row indices have changed)
fish_data_SMS_fore <- fish_data_SMS_fore %>%
  mutate(downstream_trap = replace(downstream_trap, pop %in% c("Grays_WF","Grays_CJ"),
                                   which(pop == "Grays_MS")))

# Fecundity data
# Note that L95% and U95% are reversed
fecundity <- read.csv(here("data","Data_ChumFecundity_fromHatcheryPrograms_2017-01-25.csv"),
                      header = TRUE, stringsAsFactors = FALSE) %>% 
  rename(stock = Stock, year = BY, ID = Female.., age_E = Age, L95 = U95., U95 = L95.,
         reproductive_effort = Reproductive.Effort, E_obs = Estimated.Fecundity,
         mean_mass = Green.egg.avg.weight, comments = Comments) %>% 
  mutate(stock = factor(stock))

# drop cases with age not in c(3,4,5), with estimated fecundity missing, 
# or with reproductive effort <= 16%
# add strata based on stock: Grays -> Coastal, I-205 -> Cascade, Lower Gorge -> Gorge
fecundity_data <- fecundity %>% 
  filter(age_E %in% 3:5 & !is.na(E_obs) & !is.na(reproductive_effort) & reproductive_effort > 16) %>% 
  mutate(strata = recode(stock, Grays = "Coastal", `I-205` = "Cascade", `Lower Gorge` = "Gorge")) %>% 
  select(strata, year, ID, age_E, E_obs) %>% arrange(strata, year, age_E) 
## @knitr

#--------------------------------------------------------------
# Data exploration
#--------------------------------------------------------------

# Time series of sex ratio by population
windows()
bio_data %>% filter(HW=="W" & !grepl("Hatchery|Duncan_Creek", disposition)) %>% 
  group_by(disposition, year, sex) %>% 
  summarize(n = sum(count)) %>% 
  dcast(disposition + year ~ sex, value.var = "n", fun.aggregate = sum) %>% 
  mutate(total = Female + Male) %>% 
  data.frame(., with(., binconf(x = Female, n = total))) %>%
  mutate(prop_female = PointEst) %>% 
  ggplot(aes(x = year, y = prop_female, ymin = Lower, ymax = Upper)) + 
  geom_abline(intercept = 0.5, slope = 0, color = "gray") + 
  geom_point(size = 2) + geom_line() + geom_errorbar(width = 0) +
  facet_wrap(vars(disposition), nrow = 3, ncol = 4) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
# Histogram of spawner age by sex and H/W origin
windows()
bio_data %>% mutate(age = substring(age,5,5)) %>% group_by(HW, sex, age) %>% 
  summarize(n = sum(count)) %>% mutate(prop = n/sum(n)) %>% 
  ggplot(aes(x = age, y = prop)) + geom_bar(stat = "identity") + 
  facet_wrap(vars(HW, sex), nrow = 2, ncol = 2) + theme_bw()

# Boxplots of fecundity by age
windows()
fecundity_data %>% mutate(age_E = factor(age_E)) %>% 
  ggplot(aes(x = age_E, y = E_obs)) + geom_boxplot() + theme_bw()

# Boxplots of fecundity by age, grouped by strata
windows()
fecundity_data %>% mutate(age_E = factor(age_E)) %>% 
  ggplot(aes(x = age_E, y = E_obs)) + geom_boxplot() + 
  facet_wrap(vars(strata), nrow = 1, ncol = 3) + theme_bw()

# Histograms of fecundity, grouped by age and strata
windows()
fecundity_data %>% mutate(age_E = factor(age_E)) %>%  ggplot(aes(x = E_obs)) +
  geom_histogram(aes(y = stat(density)), bins = 15, color = "white", fill = "darkgray") + 
  facet_grid(rows = vars(strata), cols = vars(age_E)) + theme_bw()

# Normal QQ plots of fecundity, grouped by age and strata
windows()
fecundity_data %>% mutate(age_E = factor(age_E)) %>%  ggplot(aes(sample = E_obs)) +
  geom_qq(distribution = qnorm) + geom_qq_line(distribution = qnorm) +
  facet_grid(rows = vars(strata), cols = vars(age_E)) + theme_bw()

# Smolts/spawner vs. spawners, grouped by population
windows()
fish_data_SMS %>% group_by(pop) %>% mutate(M0_obs = lead(M_obs), MperS = M0_obs/S_obs) %>% 
  ggplot(aes(x = S_obs, y = MperS)) + scale_x_log10() + scale_y_log10() +
  geom_point(size = 2) + labs(x = "spawners", y = "smolts / spawner") +
  facet_wrap(vars(pop), nrow = 3, ncol = 4) + theme_bw() +
  theme(panel.grid = element_blank())


#===========================================================================
# FIT RETROSPECTIVE MODELS
#===========================================================================

#--------------------------------------------------------------
# Spawner-to-spawner IPM
#--------------------------------------------------------------

# Density-independent
## @knitr fit_SS_exp
SS_exp <- salmonIPM(fish_data = fish_data_SS, stan_model = "IPM_SS_pp", SR_fun = "exp",
                    pars = c("B_rate_all","mu_Rmax","sigma_Rmax","Rmax"), include = FALSE, 
                    log_lik = TRUE, chains = 3, iter = 1500, warmup = 500,
                    control = list(adapt_delta = 0.99, max_treedepth = 13))

## @knitr print_SS_exp
print(SS_exp, prob = c(0.025,0.5,0.975),
      pars = c("alpha","phi","p_HOS","q","gamma","p","S","R","LL"), 
      include = FALSE, use_cache = FALSE)
## @knitr

launch_shinystan(SS_exp)

# Beverton-Holt
## @knitr fit_SS_BH
SS_BH <- salmonIPM(fish_data = fish_data_SS, stan_model = "IPM_SS_pp", SR_fun = "BH",
                   pars = "B_rate_all", include = FALSE, log_lik = TRUE, 
                   chains = 3, iter = 1500, warmup = 500,
                   control = list(adapt_delta = 0.99, max_treedepth = 13))

## @knitr print_SS_BH
print(SS_BH, prob = c(0.025,0.5,0.975),
      pars = c("alpha","Rmax","phi","p_HOS","q","gamma","p","S","R","LL"), 
      include = FALSE, use_cache = FALSE)
## @knitr

launch_shinystan(SS_BH)

# Ricker
## @knitr fit_SS_Ricker
SS_Ricker <- salmonIPM(fish_data = fish_data_SS, stan_model = "IPM_SS_pp", SR_fun = "Ricker",
                       pars = "B_rate_all", include = FALSE, log_lik = TRUE, 
                       chains = 3, iter = 1500, warmup = 500,
                       control = list(adapt_delta = 0.99, max_treedepth = 13))

## @knitr print_SS_Ricker
print(SS_Ricker, prob = c(0.025,0.5,0.975),
      pars = c("alpha","Rmax","phi","p_HOS","q","gamma","p","S","R","LL"), 
      include = FALSE, use_cache = FALSE)
## @knitr

launch_shinystan(SS_Ricker)


#--------------------------------------------------------------
# Spawner-smolt-spawner IPM
#
# NOTE: Deprecated because certain features of the data
# are incompatible with the model (e.g., pooling of smolts
# from Grays_MS and Grays_WF and double-counting of Grays_CJ)
#--------------------------------------------------------------

# # Density-independent
# ## @knitr fit_SMS_exp
# SMS_exp <- salmonIPM(fish_data = fish_data_SMS, ages = list(M = 1),
#                      stan_model = "IPM_SMS_pp", SR_fun = "exp",
#                      pars = c("B_rate_all","mu_Rmax","sigma_Rmax","Rmax"), include = FALSE, 
#                      log_lik = TRUE, chains = 3, iter = 1500, warmup = 500,
#                      control = list(adapt_delta = 0.99, max_treedepth = 13))
# 
# ## @knitr print_SMS_exp
# print(SMS_exp, prob = c(0.025,0.5,0.975),
#       pars = c("alpha","phi_M","phi_MS","gamma","p","p_HOS","S","M","s_MS","q","LL"), 
#       include = FALSE, use_cache = FALSE)
# ## @knitr
# 
# launch_shinystan(SMS_exp)
# 
# # Beverton-Holt
# ## @knitr fit_SMS_BH
# SMS_BH <- salmonIPM(fish_data = fish_data_SMS, ages = list(M = 1),
#                     stan_model = "IPM_SMS_pp", SR_fun = "BH",
#                     pars = "B_rate_all", include = FALSE, log_lik = TRUE, 
#                     chains = 3, iter = 1500, warmup = 500,
#                     control = list(adapt_delta = 0.99, max_treedepth = 13))
# 
# ## @knitr print_SMS_BH
# print(SMS_BH, prob = c(0.025,0.5,0.975),
#       pars = c("alpha","Rmax","phi_M","phi_MS","gamma","p","p_HOS","S","M","s_MS","q","LL"), 
#       include = FALSE, use_cache = FALSE)
# ## @knitr
# 
# launch_shinystan(SMS_BH)
# 
# # Ricker
# ## @knitr fit_SMS_Ricker
# SMS_Ricker <- salmonIPM(fish_data = fish_data_SMS, ages = list(M = 1),
#                         stan_model = "IPM_SMS_pp", SR_fun = "Ricker",
#                         pars = "B_rate_all", include = FALSE, log_lik = TRUE, 
#                         chains = 3, iter = 1500, warmup = 500,
#                         control = list(adapt_delta = 0.99, max_treedepth = 13))
# 
# ## @knitr print_SMS_Ricker
# print(SMS_Ricker, prob = c(0.025,0.5,0.975),
#       pars = c("alpha","Rmax","phi_M","phi_MS","gamma","p","p_HOS","S","M","s_MS","q","LL"), 
#       include = FALSE, use_cache = FALSE)
# ## @knitr
# 
# launch_shinystan(SMS_Ricker)


#--------------------------------------------------------------
# Lower Columbia chum spawner-smolt-spawner IPM
#--------------------------------------------------------------

# Density-independent
## @knitr fit_LCRchum_exp
LCRchum_exp <- salmonIPM(fish_data = fish_data_SMS,  fecundity_data = fecundity_data,
                         ages = list(M = 1), stan_model = "IPM_LCRchum_pp", SR_fun = "exp",
                         pars = c("B_rate_all","mu_Emax","sigma_Emax","Emax"), 
                         include = FALSE, log_lik = TRUE, 
                         chains = 3, iter = 1500, warmup = 500,
                         control = list(adapt_delta = 0.99, max_treedepth = 15))

## @knitr print_LCRchum_exp
print(LCRchum_exp, prob = c(0.025,0.5,0.975),
      pars = c("eta_pop_EM","eta_year_EM","eta_year_MS","eta_pop_p","p",
               "tau_M","tau_S","p_HOS","E","S","M","s_EM","s_MS","q","LL"), 
      include = FALSE, use_cache = FALSE)
## @knitr

launch_shinystan(LCRchum_exp)

# Beverton-Holt
## @knitr fit_LCRchum_BH
LCRchum_BH <- salmonIPM(fish_data = fish_data_SMS, fecundity_data = fecundity_data,
                        ages = list(M = 1), stan_model = "IPM_LCRchum_pp", SR_fun = "BH",
                        pars = "B_rate_all", include = FALSE, log_lik = TRUE, 
                        chains = 3, iter = 1500, warmup = 500,
                        control = list(adapt_delta = 0.95, max_treedepth = 13))

## @knitr print_LCRchum_BH
print(LCRchum_BH, prob = c(0.025,0.5,0.975),
      pars = c("psi","Mmax","eta_year_M","eta_year_MS","eta_pop_p",
               "p","tau_M","tau_S","p_HOS","E_hat","M","S","s_MS","q","LL"), 
      include = FALSE, use_cache = FALSE)
## @knitr

launch_shinystan(LCRchum_BH)

# Ricker
## @knitr fit_LCRchum_Ricker
LCRchum_Ricker <- salmonIPM(fish_data = fish_data_SMS, fecundity_data = fecundity_data,
                            ages = list(M = 1), stan_model = "IPM_LCRchum_pp", SR_fun = "Ricker",
                            pars = "B_rate_all", include = FALSE, log_lik = TRUE, 
                            chains = 3, iter = 1500, warmup = 500,
                            control = list(adapt_delta = 0.99, max_treedepth = 14))

## @knitr print_LCRchum_Ricker
print(LCRchum_Ricker, prob = c(0.025,0.5,0.975),
      pars = c("psi","Mmax","eta_year_M","eta_year_MS","eta_pop_p",
               "p","tau_M","tau_S","p_HOS","E_hat","M","S","s_MS","q","LL"), 
      include = FALSE, use_cache = FALSE)
## @knitr

launch_shinystan(LCRchum_Ricker)


#--------------------------------------------------------------
# Model selection using LOO
#--------------------------------------------------------------

# Observationwise log-likelihood of each fitted model
# Here an observation is a row of fish_data, and the total likelihood includes 
# components for smolt abundance, spawner abundance, age-frequency, H/W-frequency,
# and smolt and spawner observation error SDs
## @knitr loo_LCRchum
LL_LCRchum <- lapply(list(exp = LCRchum_exp, BH = LCRchum_BH, Ricker = LCRchum_Ricker),
                 loo::extract_log_lik, parameter_name = "LL", merge_chains = FALSE)

# Relative ESS of posterior draws of observationwise likelihood 
r_eff_LCRchum <- lapply(LL_LCRchum, function(x) relative_eff(exp(x)))

# PSIS-LOO
LOO_LCRchum <- lapply(1:length(LL_LCRchum), 
                      function(i) loo(LL_LCRchum[[i]], r_eff = r_eff_LCRchum[[i]]))
names(LOO_LCRchum) <- names(LL_LCRchum)

## Compare all three models
loo_compare(LOO_LCRchum)

## Exponential vs. Ricker
loo_compare(LOO_LCRchum[c("exp","Ricker")])

## Exponential vs. Beverton-Holt
loo_compare(LOO_LCRchum[c("exp","BH")])

## Beverton-Holt vs. Ricker
loo_compare(LOO_LCRchum[c("BH","Ricker")])
## @knitr


#===========================================================================
# FIT PROSPECTIVE FORECASTING MODELS
#===========================================================================

#--------------------------------------------------------------
# Lower Columbia chum spawner-smolt-spawner IPM
#--------------------------------------------------------------

# Ricker
## @knitr fit_LCRchum_BH_fore
LCRchum_BH_fore <- salmonIPM(fish_data = fish_data_SMS_fore, fecundity_data = fecundity_data,
                             ages = list(M = 1), stan_model = "IPM_LCRchum_pp", SR_fun = "BH",
                             pars = c("Emax","p","p_HOS","B_rate_all"), 
                             include = FALSE, chains = 3, iter = 1500, warmup = 500,
                             control = list(adapt_delta = 0.99, max_treedepth = 15))

## @knitr
print(LCRchum_BH_fore, prob = c(0.025,0.5,0.975),
      pars = c("eta_pop_EM","eta_year_EM","eta_year_MS","eta_pop_p",
               "E","S","M","s_EM","s_MS","q","tau_M","tau_S"), 
      include = FALSE, use_cache = FALSE)


#--------------------------------------------------------------
# Save stanfit objects
#--------------------------------------------------------------

save(list = ls()[sapply(ls(), function(x) do.call(class, list(as.name(x)))) == "stanfit"], 
     file = here("analysis","results","LCRchumIPM.RData"))



#===========================================================================
# FIGURES
#===========================================================================

#--------------------------------------------------------------------------#
############# FIGURES THAT WORK WITH LOWER COLUMBIA CHUM IPM ###############
#--------------------------------------------------------------------------#

#--------------------------------------------------------------------
# Lower Columbia chum spawner-egg-smolt-spawner #
# S-R curves (spawners to egg deposition)
# Posterior distributions of fecundity and Emax parameters
# Time series of egg-smolt survival and SAR
#--------------------------------------------------------------------

mod_name <- "LCRchum_BH"

dev.new(width = 7, height = 8)
# png(filename=here("analysis","results",paste0("SR_",mod_name,".png")),
#     width=7, height=8, units="in", res=200, type="cairo-png")

## @knitr plot_LCM_params
life_cycle <- "LCRchum"
dat <- fish_data_SMS
SR_fun <- unlist(strsplit(mod_name, "_"))[2]

SR_eval <- function(alpha, Rmax = NULL, S, SR_fun) 
{
  switch(SR_fun,
         exp = alpha*S,
         BH = alpha*S/(1 + alpha*S/Rmax),
         Ricker = alpha*S*exp(-alpha*S/(exp(1)*Rmax)))
}

# fecundity
mu_E <- do.call(extract1, list(as.name(mod_name), "mu_E"))
ages <- substring(names(select(dat, starts_with("n_age"))), 6, 6)
# egg deposition
q <- do.call(extract1, list(as.name(mod_name), "q"))
alpha <- apply(sweep(q, c(1,3), mu_E, "*"), 1:2, sum)*0.5  # age-weighted fecundity / 2
mu_pop_alpha <- t(as.matrix(aggregate(t(log(alpha)), list(pop = dat$pop), mean)[,-1]))
mu_alpha <- rowMeans2(log(alpha))
mu_Emax <- as.vector(do.call(extract1, list(as.name(mod_name), "mu_Emax")))
Emax <- do.call(extract1, list(as.name(mod_name), "Emax"))
S <- colMedians(do.call(extract1, list(as.name(mod_name), "S")))
S_grid <- matrix(seq(0, quantile(S/dat$A, 0.9, na.rm = TRUE), length = 100),
                 nrow = length(mu_alpha), ncol = 100, byrow = TRUE)
E_ESU <- SR_eval(alpha = exp(mu_alpha), Rmax = 1e6*exp(mu_Emax), S = S_grid, SR_fun = SR_fun)
E_pop <- sapply(1:ncol(Emax), function(i) {
  colMedians(SR_eval(alpha = alpha[,i], Rmax = 1e6*Emax[,i], S = S_grid, SR_fun = SR_fun))
})
# egg-smolt survival
y <- sort(unique(dat$year))
eta_year_EM <- do.call(extract1, list(as.name(mod_name), "eta_year_EM"))
mu_EM <- do.call(extract1, list(as.name(mod_name), "mu_EM"))
s_hat_EM <- plogis(sweep(eta_year_EM, 1, qlogis(mu_EM), "+"))
s_EM <- do.call(extract1, list(as.name(mod_name), "s_EM"))
# SAR
eta_year_MS <- do.call(extract1, list(as.name(mod_name), "eta_year_MS"))
mu_MS <- do.call(extract1, list(as.name(mod_name), "mu_MS"))
s_hat_MS <- plogis(sweep(eta_year_MS, 1, qlogis(mu_MS), "+"))
s_MS <- do.call(extract1, list(as.name(mod_name), "s_MS"))
# colors
c1 <- "slategray4"
c1t <- transparent(c1, trans.val = 0.5)
c1tt <- transparent(c1, trans.val = 0.7)
ac <- viridis(length(ages), end = 0.8, alpha = 0.7) 

par(mfrow = c(3,2), mar = c(5.1,5.1,1,1))

# Egg deposition vs. spawners
plot(S_grid[1,], colMedians(E_ESU)*1e-6, type = "l", lwd=3, col = c1, las = 1,
     cex.axis = 1.2, cex.lab = 1.5, xaxs = "i", yaxs = "i", #yaxt = "n", 
     ylim = range(E_pop*1e-6), xlab = "Spawners", ylab = "")
for(i in 1:ncol(E_pop)) 
  lines(S_grid[1,], E_pop[,i]*1e-6, col = c1t)
polygon(c(S_grid[1,], rev(S_grid[1,])), 
        c(colQuantiles(E_ESU, probs = 0.05)*1e-6, 
          rev(colQuantiles(E_ESU, probs = 0.95)*1e-6)), 
        col = c1tt, border = NA)
rug(S, col = c1)
title(ylab = "Egg deposition (millions)", line = 3.5, cex.lab = 1.5)
text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "A", cex = 1.5)

# log(recruits/spawner) vs. spawners
plot(S_grid[1,], colMedians(log(E_ESU/S_grid)), type = "l", lwd=3, col = c1, las = 1,
     cex.axis = 1.2, cex.lab = 1.5, xaxs = "i", yaxs = "i", 
     xlab="Spawners", ylab = "log(eggs/spawner)",
     ylim = range(colQuantiles(log(E_ESU/S_grid), probs = c(0.05,0.95)), na.rm = TRUE))
for(i in 1:ncol(E_pop))
  lines(S_grid[1,], log(E_pop[,i]/S_grid[1,]), col = c1t)
polygon(c(S_grid[1,], rev(S_grid[1,])),
        c(colQuantiles(log(E_ESU/S_grid), probs = 0.05),
          rev(colQuantiles(log(E_ESU/S_grid), probs = 0.95))),
        col = c1tt, border = NA)
rug(S, col = c1)
text(par("usr")[2], par("usr")[4], adj = c(1.5,1.5), "B", cex = 1.5)

# Posterior distributions of fecundity by age
dd_age <- vector("list", length(ages))
for(a in 1:length(dd_age))
  dd_age[[a]] <- density(mu_E[,a])

plot(dd_age[[1]]$x, dd_age[[1]]$y, pch = "", las = 1, cex.axis = 1.2, cex.lab = 1.5,
     xlab = "Mean fecundity", ylab = "", xaxs = "i", yaxt = "n",
     xlim = range(sapply(dd_age, function(m) m$x)),
     ylim = range(sapply(dd_age, function(m) m$y)))
for(a in 1:length(ages))
  lines(dd_age[[a]]$x, dd_age[[a]]$y, col = ac[a], lwd = 2)
title(ylab = "Probability density", line = 1, cex.lab = 1.5)
text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "C", cex = 1.5)
legend("topright", paste("age", ages, "  "), x.intersp = 0.5, col = ac, lwd = 2, 
       xjust = 0.5, bty = "n")

# Posterior densities of log(Emax)
dd_ESU <- density(mu_Emax)
dd_pop <- vector("list", length(levels(dat$pop)))
for(i in 1:length(dd_pop))
  dd_pop[[i]] <- density(log(Emax[,i]))

plot(dd_ESU$x, dd_ESU$y, type = "l", lwd = 3, col = c1, las = 1, xaxs = "i", yaxt = "n",
     cex.axis = 1.2, cex.lab = 1.5, xlab = bquote(log(italic(E)[max])), ylab = "",
     xlim = range(c(dd_ESU$x, sapply(dd_pop, function(m) m$x))),
     ylim = range(c(dd_ESU$y, sapply(dd_pop, function(m) m$y))))
for(i in 1:length(dd_pop))
  lines(dd_pop[[i]]$x, dd_pop[[i]]$y, col = c1t)
title(ylab = "Probability density", line = 1, cex.lab = 1.5)
text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "D", cex = 1.5)

# Egg-smolt survival
plot(y, colMedians(s_hat_EM), type = "n", las = 1, cex.axis = 1.2, cex.lab = 1.5,
     ylim = range(0, colQuantiles(s_hat_EM, probs = 0.95), 
                  colQuantiles(s_EM, probs = 0.95)), 
     xaxs = "i", xaxt = "n", xlab = "Outmigration year", ylab = "")
mtext("Egg-to-smolt survival", side = 2, line = 3.7, cex = par("cex")*1.5)
polygon(c(y, rev(y)), 
        c(colQuantiles(s_hat_EM, probs = 0.05), rev(colQuantiles(s_hat_EM, probs = 0.95))),
        col = c1tt, border = NA)
lines(y, colMedians(s_hat_EM), col = c1t, lwd = 3)
for(j in levels(dat$pop))
  lines(dat$year[dat$pop == j], colMedians(s_EM[,dat$pop == j]), col = c1t)
axis(side = 1, at = y[y %% 5 == 0], cex.axis = 1.2)
rug(y[y %% 5 != 0], ticksize = -0.02)
text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "E", cex = 1.5)

# SAR
plot(y, colMedians(s_hat_MS), type = "n", las = 1, cex.axis = 1.2, cex.lab = 1.5,
     ylim = range(0, colQuantiles(s_hat_MS, probs = 0.95), 
                  colQuantiles(s_MS, probs = 0.95)), 
     xaxs = "i", xaxt = "n", xlab = "Outmigration year", ylab = "")
mtext("SAR", side = 2, line = 3.7, cex = par("cex")*1.5)
polygon(c(y, rev(y)), 
        c(colQuantiles(s_hat_MS, probs = 0.05), rev(colQuantiles(s_hat_MS, probs = 0.95))),
        col = c1tt, border = NA)
lines(y, colMedians(s_hat_MS), col = c1t, lwd = 3)
for(j in levels(dat$pop))
  lines(dat$year[dat$pop == j], colMedians(s_MS[,dat$pop == j]), col = c1t)
axis(side = 1, at = y[y %% 5 == 0], cex.axis = 1.2)
rug(y[y %% 5 != 0], ticksize = -0.02)
text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "F", cex = 1.5)

rm(list=c("mod_name","life_cycle","SR_fun","mu_alpha","mu_Emax","S","S_grid","E_ESU",
          "c1","c1t","c1tt","dd_ESU","dd_pop","dd_age","alpha","Emax","E_pop","ac","ages",
          "y","eta_year_EM","mu_EM","s_hat_EM","eta_year_MS","mu_MS","s_hat_MS","dat"))
## @knitr
# dev.off()


#--------------------------------------------------------------------------------
# Lower Columbia chum spawner-egg-smolt-spawner # 
# Estimated S-R curves for each pop, overlaid with states and observed data
#--------------------------------------------------------------------------------

mod_name <- "LCRchum_Ricker"
life_stage <- "M"   # "S" = spawners, "M" = smolts

dev.new(width=14.5,height=8)
# png(filename=here("analysis", "results", paste0("SR_fit_", mod_name, ".png")),
#     width=14.5*0.9, height=8*0.9, units="in", res=200, type="cairo-png")

## @knitr plot_SR_fits
SR_eval <- function(alpha, Rmax = NULL, S, SR_fun) 
{
  switch(SR_fun,
         exp = alpha*S,
         BH = alpha*S/(1 + alpha*S/Rmax),
         Ricker = alpha*S*exp(-alpha*S/(exp(1)*Rmax)))
}

scl <- switch(life_stage, M = 1e-6, S = 1)
life_cycle <- unlist(strsplit(mod_name, "_"))[1]
dat <- switch(life_cycle, SS = fish_data, SMS = fish_data_SMS, LCRchum = fish_data_SMS)
S_obs <- dat$S_obs
M_obs <- dat$M_obs*scl
S_IPM <- do.call(extract1, list(as.name(mod_name), "S"))
M_IPM <- do.call(extract1, list(as.name(mod_name), "M"))*scl
SR_fun <- unlist(strsplit(mod_name, "_"))[2]

# S-R params
if(life_cycle == "SMS") {
  alpha <- do.call(extract1, list(as.name(mod_name), "alpha"))
  Rmax <- do.call(extract1, list(as.name(mod_name), "Rmax"))
} else if(life_cycle == "LCRchum") {
  # fecundity
  mu_E <- do.call(extract1, list(as.name(mod_name), "mu_E"))
  ages <- substring(names(select(dat, starts_with("n_age"))), 6, 6)
  # calculate alpha and egg deposition
  q <- do.call(extract1, list(as.name(mod_name), "q"))
  alpha <- apply(sweep(q, c(1,3), mu_E, "*"), 1:2, sum)*0.5  # age-weighted fecundity / 2
  alpha <- exp(t(as.matrix(aggregate(t(log(alpha)), list(pop = dat$pop), mean)[,-1])))
  Rmax <- do.call(extract1, list(as.name(mod_name), "Emax"))*1e6
  # egg-smolt survival
  eta_pop_EM <- do.call(extract1, list(as.name(mod_name), "eta_pop_EM"))
  mu_EM <- do.call(extract1, list(as.name(mod_name), "mu_EM"))
  s_pop_EM <- plogis(sweep(eta_pop_EM, 1, qlogis(mu_EM), "+"))
}
#colors
c1 <- "slategray4"
c1t <- transparent(c1, trans.val = 0.7)
y <- sort(unique(dat$year))
cy <- plasma(length(y), end = 0.9, direction = -1, alpha = 0.8)

op <- par(mfrow=c(3,4), mar=c(1,3,4.1,1), oma=c(4.1,3.1,0,13))

for(i in levels(dat$pop))
{
  indx <- dat$pop == i
  ii <- levels(dat$pop) == i
  yi <- dat$year[indx]
  cyi <- cy[match(yi, y)]
  S_upper <- max(S_obs[indx], colQuantiles(S_IPM[,indx], probs = 0.9), na.rm = TRUE)
  S <- matrix(seq(0, S_upper, length = 1000), nrow = nrow(alpha), ncol = 1000, byrow = TRUE)
  M_fit_IPM <- SR_eval(alpha = alpha[,ii], Rmax = Rmax[,ii], S = S, SR_fun = SR_fun)*scl
  if (life_cycle == "LCRchum") M_fit_IPM <- M_fit_IPM * s_pop_EM[,ii]
  
  plot(colMedians(S_IPM[,indx]), colMedians(M_IPM[,indx]), pch = "",
       xlim = range(0, S_obs[indx], colQuantiles(S_IPM[,indx], probs = 0.6), na.rm = TRUE),
       ylim = range(0, M_obs[indx], colQuantiles(M_IPM[,indx], probs = 0.6), na.rm = TRUE), 
       las = 1, xlab = "", ylab = "", cex.axis = 1.2)
  mtext(i, side = 3, line = 0.5, cex = par("cex")*1.5)
  if(par("mfg")[2] == 1) 
    mtext("Smolts (millions)", side = 2, line = 3.5, cex = par("cex")*1.5)
  if(par("mfg")[1] == par("mfg")[3]) 
    mtext("Spawners", side = 1, line = 3, cex = par("cex")*1.5)
  polygon(c(S[1,], rev(S[1,])), 
          c(colQuantiles(M_fit_IPM, probs = 0.05), 
            rev(colQuantiles(M_fit_IPM, probs = 0.95))), 
          col = c1t, border = NA)
  lines(S[1,], colMedians(M_fit_IPM), lwd = 3, col = c1)
  points(colMedians(S_IPM[,indx]), colMedians(M_IPM[,indx]), pch = 3, cex = 1.8, col = cyi)
  segments(x0 = colMedians(S_IPM[,indx]), 
           y0 = colQuantiles(M_IPM[,indx], probs = 0.05), 
           y1 = colQuantiles(M_IPM[,indx], probs = 0.95),
           col = cyi)
  segments(x0 = colQuantiles(S_IPM[,indx], probs = 0.05), 
           x1 = colQuantiles(S_IPM[,indx], probs = 0.95),
           y0 = colMedians(M_IPM[,indx]), 
           col = cyi)
  points(S_obs[indx], M_obs[indx], pch = 16, col = cyi, cex = 1.8)
}
par(op)
par(usr = c(0,1,0,1))
legend(0.9, 0.5, c("observation", "state", "", y), fill = c(rep(NA,3), cy), border = NA,
       col = c(rep("black",3), cy), pch = c(16, 3, rep(NA, length(y) + 1)), pt.cex = 1.2, 
       yjust = 0.5, xpd = NA, box.lwd = 0.5)

if(life_cycle == "LRCchum") rm(list = c("mu_E","agess","q","eta_pop_M","mu_EM","s_pop_EM"))
rm(list = c("mod_name","SR_fun","S_IPM","M_IPM","S_obs","M_obs","alpha","Rmax",
            "M_fit_IPM","c1","c1t","yi","y","cy","cyi","indx","life_stage","scl",
            "ii","S_upper","op","dat"))
## @knitr
# dev.off()


#--------------------------------------------------------------------------------
# Time series of observed and fitted total spawners or smolts for each pop
#--------------------------------------------------------------------------------

mod_name <- "LCRchum_BH"
life_stage <- "M"   # "S" = spawners, "M" = smolts

dev.new(width=13,height=8)
# png(filename=here("analysis", "results", paste0(life_stage, "_fit_", mod_name, ".png")),
#     width=13*0.9, height=8*0.9, units="in", res=200, type="cairo-png")

## @knitr plot_spawner_smolt_ts
life_cycle <- unlist(strsplit(mod_name, "_"))[1]
forecasting <- ifelse(identical(unlist(strsplit(mod_name, "_"))[3], "fore"), "yes", "no")
dat <- switch(life_cycle, SS = fish_data_SS, 
              SMS = switch(forecasting, no = fish_data_SMS, yes = fish_data_SMS_fore),
              LCRchum = switch(forecasting, no = fish_data_SMS, yes = fish_data_SMS_fore))

N_obs <- dat[,paste0(life_stage, "_obs")]
N_IPM <- do.call(extract1, list(as.name(mod_name), life_stage))
tau <- do.call(extract1, list(as.name(mod_name),
                              switch(life_cycle, SS = "tau",
                                     SMS = switch(life_stage, M = "tau_M", S = "tau_S"),
                                     LCRchum = switch(life_stage, M = "tau_M", S = "tau_S"))))
if(life_cycle == "LCRchum" & life_stage == "M")
  N_IPM[,na.omit(dat$downstream_trap)] <- N_IPM[,na.omit(dat$downstream_trap)] + 
  N_IPM[,which(!is.na(dat$downstream_trap))]
N_obs_IPM <- N_IPM * rlnorm(length(N_IPM), 0, tau)

c1 <- "slategray4"
c1t <- transparent(c1, trans.val = 0.5)
c1tt <- transparent(c1, trans.val = 0.7)
pc <- ifelse(is.na(dat[,paste0("tau_", life_stage, "_obs")]), 1, 16)

par(mfrow=c(3,4), mar=c(1,3,4.1,1), oma=c(4.1,3.1,0,0))

for(i in levels(dat$pop))
{
  yi <- dat$year[dat$pop==i]
  plot(yi, N_obs[dat$pop==i], pch = "",
       xlim = range(dat$year),
       ylim = range(pmax(N_obs[dat$pop==i], 1),
                    colQuantiles(N_obs_IPM[,dat$pop==i], probs = c(0.05,0.95)), na.rm = T), 
       las = 1, xaxt = "n", yaxt = "n", xlab = "", ylab = "", log = "y")
  axis(side = 1, at = yi[yi %% 5 == 0], cex.axis = 1.2)
  rug(yi[yi %% 5 != 0], ticksize = -0.02)
  # tck <- maglab(10^par("usr")[3:4], log = TRUE)
  # axis(2, at = tck$labat, labels = tck$exp, cex.axis = 1.2, las = 1)
  magaxis(2, minorn = 0, crunch = FALSE, tcl = -0.4, las = 1, cex.axis = 1.2)
  mtext(i, side = 3, line = 0.5, cex = par("cex")*1.5)
  if(par("mfg")[2] == 1) 
    mtext(switch(life_stage, M = "Smolts", S = "Spawners"), 
          side = 2, line = 3.5, cex = par("cex")*1.5)
  if(par("mfg")[1] == par("mfg")[3]) 
    mtext("Year", side = 1, line = 3, cex = par("cex")*1.5)
  lines(yi, colMedians(N_IPM[,dat$pop==i]), col = c1, lwd = 3)
  polygon(c(yi, rev(yi)), 
          c(colQuantiles(N_IPM[,dat$pop==i], probs = 0.05), 
            rev(colQuantiles(N_IPM[,dat$pop==i], probs = 0.95))),
          col = c1t, border = NA)
  polygon(c(yi, rev(yi)), 
          c(colQuantiles(N_obs_IPM[,dat$pop==i], probs = 0.05), 
            rev(colQuantiles(N_obs_IPM[,dat$pop==i], probs = 0.95))),
          col = c1tt, border = NA)
  points(yi, N_obs[dat$pop==i], pch = pc[dat$pop==i], lwd = 2, cex = 1.8)
}

rm(list = c("mod_name","forecasting","life_stage","life_cycle","dat",
            "N_IPM","N_obs_IPM","N_obs","c1","c1t","c1tt","yi","tau"))
## @knitr
# dev.off()


#--------------------------------------------------------------------------------
# Time series of observed and fitted spawner age structure for each pop
#--------------------------------------------------------------------------------

mod_name <- "LCRchum_BH"

# dev.new(width=13,height=8.5)
png(filename=here("analysis", "results", paste0("q_fit_", mod_name, ".png")),
    width=13*0.9, height=8.5*0.9, units="in", res=200, type="cairo-png")

## @knitr plot_spawner_age_ts
life_cycle <- unlist(strsplit(mod_name, "_"))[1]
dat <- switch(life_cycle, SS = fish_data_SS, SMS = fish_data_SMS, LCRchum = fish_data_SMS)

n_age_obs <- select(dat, starts_with("n_age")) 
q_IPM <- do.call(extract1, list(as.name(mod_name), "q"))

c1 <- viridis(ncol(n_age_obs), end = 0.8) 
c1t <- transparent(c1, trans.val = 0.2)
c1tt <- transparent(c1, trans.val = 0.7)

op <- par(mfrow=c(3,4), mar=c(1,3,4.1,1), oma=c(4.1,3.1,4,0))

for(i in levels(dat$pop))
{
  yi <- dat$year[dat$pop==i]
  plot(yi, rep(0.5, length(yi)), pch = "", xlim = range(dat$year), ylim = c(0,1), 
       las = 1, cex.axis = 1.2, xaxt = "n", xlab = "", ylab = "")
  axis(side = 1, at = yi[yi %% 5 == 0], cex.axis = 1.2)
  rug(yi[yi %% 5 != 0], ticksize = -0.02)
  mtext(i, side = 3, line = 0.5, cex = par("cex")*1.5)
  if(par("mfg")[2] == 1) 
    mtext("Proportion at age", side = 2, line = 3.5, cex = par("cex")*1.5)
  if(par("mfg")[1] == par("mfg")[3]) 
    mtext("Year", side = 1, line = 3, cex = par("cex")*1.5)
  
  for(a in 1:ncol(n_age_obs))
  {
    q_obs <- binconf(n_age_obs[dat$pop==i,a], rowSums(n_age_obs[dat$pop==i,]), alpha = 0.1)
    lines(yi, colMedians(q_IPM[,dat$pop==i,a]), col = c1t[a], lwd = 2)
    polygon(c(yi, rev(yi)),
            c(colQuantiles(q_IPM[,dat$pop==i,a], probs = 0.05),
              rev(colQuantiles(q_IPM[,dat$pop==i,a], probs = 0.95))),
            col = c1tt[a], border = NA)
    points(yi, q_obs[,"PointEst"], pch = 16, col = c1t[a], cex = 1.8)
    segments(x0 = yi, y0 = q_obs[,"Lower"], y1 = q_obs[,"Upper"], col = c1t[a])
  }
}
par(op)
par(usr = c(0,1,0,1))
legend(0.5, 1.1, paste("age", substring(names(n_age_obs), 6, 6), "  "), x.intersp = 0.5,
       col = c1, pch = 16, pt.cex = 1.2, lwd = 2, horiz = TRUE, xjust = 0.5, 
       xpd = NA, box.lwd = 0.5)
       
rm(list = c("mod_name","life_cycle","dat","q_IPM","n_age_obs","q_obs","op",
            "c1","c1t","c1tt","yi"))
## @knitr
dev.off()


#--------------------------------------------------------------------------------
# Time series of observed and fitted p_HOS for each pop
#--------------------------------------------------------------------------------

mod_name <- "LCRchum_BH"

# dev.new(width=13,height=8.5)
png(filename=here("analysis", "results", paste0("p_HOS_fit_", mod_name, ".png")),
    width=13*0.9, height=8.5*0.9, units="in", res=200, type="cairo-png")

## @knitr plot_p_HOS_ts
life_cycle <- unlist(strsplit(mod_name, "_"))[1]
dat <- switch(life_cycle, SS = fish_data_SS, SMS = fish_data_SMS, LCRchum = fish_data_SMS)
p_HOS_IPM <- matrix(0, nrow(do.call(as.matrix, list(as.name(mod_name)))), nrow(dat))
p_HOS_IPM[,dat$fit_p_HOS==1] <- do.call(extract1, list(as.name(mod_name), "p_HOS"))

c1 <- "slategray4"
c1t <- transparent(c1, trans.val = 0.5)

par(mfrow=c(3,4), mar=c(1,3,4.1,1), oma=c(4.1,3.1,0,0))

for(i in levels(dat$pop))
{
  yi <- dat$year[dat$pop==i]
  plot(yi, rep(0.5, length(yi)), pch = "", xlim = range(dat$year), ylim = c(0,1), 
       las = 1, cex.axis = 1.2, xaxt = "n", xlab = "", ylab = "")
  axis(side = 1, at = yi[yi %% 5 == 0], cex.axis = 1.2)
  rug(yi[yi %% 5 != 0], ticksize = -0.02)
  mtext(i, side = 3, line = 0.5, cex = par("cex")*1.5)
  if(par("mfg")[2] == 1) 
    mtext(bquote(italic(p)[HOS]), side = 2, line = 3.5, cex = par("cex")*1.5)
  if(par("mfg")[1] == par("mfg")[3]) 
    mtext("Year", side = 1, line = 3, cex = par("cex")*1.5)
  p_HOS_obs <- binconf(dat$n_H_obs[dat$pop==i], 
                       dat$n_H_obs[dat$pop==i] + dat$n_W_obs[dat$pop==i], 
                       alpha = 0.1)
  lines(yi, colMedians(p_HOS_IPM[,dat$pop==i]), col = c1, lwd = 2)
  polygon(c(yi, rev(yi)),
          c(colQuantiles(p_HOS_IPM[,dat$pop==i], probs = 0.05),
            rev(colQuantiles(p_HOS_IPM[,dat$pop==i], probs = 0.95))),
          col = c1t, border = NA)
  points(yi, p_HOS_obs[,"PointEst"], pch = ifelse(dat$fit_p_HOS[dat$pop==i], 16, 1), cex = 1.8)
  segments(x0 = yi, y0 = p_HOS_obs[,"Lower"], y1 = p_HOS_obs[,"Upper"])
}

rm(list = c("mod_name","life_cycle","dat","p_HOS_IPM","p_HOS_obs","c1","c1t","yi"))
## @knitr
dev.off()


#--------------------------------------------------------------------------------
# Lower Columbia chum spawner-egg-smolt-spawner #
# Observed and fitted distributions of fecundity by age
#--------------------------------------------------------------------------------

mod_name <- "LCRchum_BH"

# dev.new(width=7,height=7)
png(filename=here("analysis","results",paste0("fecundity_fit_", mod_name, ".png")),
    width=7, height=7, units="in", res=200, type="cairo-png")

## @knitr plot_fecundity_fit
ages <- substring(names(select(fish_data_SMS, starts_with("n_age"))), 6, 6)
E_obs <- fecundity_data$E_obs
E_seq <- seq(min(E_obs, na.rm = TRUE), max(E_obs, na.rm = TRUE), length = 500)
mu_E <- do.call(extract1, list(as.name(mod_name), "mu_E"))
sigma_E <- do.call(extract1, list(as.name(mod_name), "sigma_E"))
E_fit <- array(NA, c(nrow(mu_E), length(E_seq), ncol(mu_E)))
for(a in 1:length(ages))
  E_fit[,,a] <- sapply(E_seq, function(x) dnorm(x, mu_E[,a], sigma_E[,a]))

c1 <- viridis(length(ages), end = 0.8) 
c1t <- transparent(c1, trans.val = 0.5)
c1tt <- transparent(c1, trans.val = 0.7)

par(mfrow = c(3,1), mar = c(3,2,0,2), oma = c(2,2,0,0))

for(a in 1:length(ages))
{
  hist(E_obs[fecundity_data$age_E == ages[a]], 20, prob = TRUE, 
       col = c1tt[a], border = "white", las = 1, cex.axis = 1.2, cex.lab = 1.5,
       xlim = range(E_seq), ylim = range(0, apply(E_fit, 2:3, quantile, 0.99)),
       xlab = "", ylab = "", main = "", xaxs = "i", yaxt = "n", bty = "n")
  lines(E_seq, colMedians(E_fit[,,a]), col = c1[a], lwd = 3)
  polygon(c(E_seq, rev(E_seq)),
          c(colQuantiles(E_fit[,,a], probs = 0.05), 
            rev(colQuantiles(E_fit[,,a], probs = 0.95))),
          col = c1t[a], border = NA)
  text(par("usr")[1] + 0.8*diff(par("usr")[1:2]), par("usr")[4]*0.5, 
       labels = paste("age", ages[a]), cex = 1.5, col = c1[a], adj = 1)
}
title(xlab = "Fecundity", ylab = "Probability density", cex.lab = 1.5, line = 0, outer = TRUE)

rm(list = c("mod_name","c1","c1t","c1tt","ages","E_obs","E_seq","mu_E","sigma_E","E_fit"))
## @knitr
dev.off()


#--------------------------------------------------------------------------------
# Lower Columbia chum spawner-egg-smolt-spawner #
# Observed and fitted distributions of "known" smolt and spawner 
# observation error SDs
#--------------------------------------------------------------------------------

mod_name <- "LCRchum_BH"

# dev.new(width=6,height=8)
png(filename=here("analysis","results",paste0("tau_fit_", mod_name, ".png")),
    width=6, height=8, units="in", res=200, type="cairo-png")

## @knitr plot_obs_error_fit
tau_M_obs <- fish_data_SMS$tau_M_obs
tau_M_seq <- seq(min(tau_M_obs, na.rm = TRUE), max(tau_M_obs, na.rm = TRUE), length = 500)
mu_tau_M <- do.call(extract1, list(as.name(mod_name), "mu_tau_M"))
sigma_tau_M <- do.call(extract1, list(as.name(mod_name), "sigma_tau_M"))
tau_M_fit <- sapply(tau_M_seq, function(x) dlnorm(x, log(mu_tau_M), sigma_tau_M))

tau_S_obs <- fish_data_SMS$tau_S_obs
tau_S_seq <- seq(min(tau_S_obs, na.rm = TRUE), max(tau_S_obs, na.rm = TRUE), length = 500)
mu_tau_S <- do.call(extract1, list(as.name(mod_name), "mu_tau_S"))
sigma_tau_S <- do.call(extract1, list(as.name(mod_name), "sigma_tau_S"))
tau_S_fit <- sapply(tau_S_seq, function(x) dlnorm(x, log(mu_tau_S), sigma_tau_S))

c1 <- "slategray4"
c1t <- transparent(c1, trans.val = 0.6)

par(mfcol = c(2,1), mar = c(5.1,5.1,2,2),  oma = c(0,0.1,0,0))

# smolt observation error SD
hist(tau_M_obs, 10, prob = TRUE, las = 1, cex.axis = 1.2, cex.lab = 1.5, 
     col = "lightgray", border = "white",
     ylim = c(0, max(colQuantiles(tau_M_fit, probs = 0.95))),
     xlab = bquote(tau[italic(M)]), ylab = "Probability density", 
     main = "Smolt observation error")
polygon(c(tau_M_seq, rev(tau_M_seq)),
        c(colQuantiles(tau_M_fit, probs = 0.05), rev(colQuantiles(tau_M_fit, probs = 0.95))),
        col = c1t, border = NA)
lines(tau_M_seq, colMedians(tau_M_fit), col = c1, lwd = 3)

# spawner observation error SD
hist(tau_S_obs, 20, prob = TRUE, las = 1, cex.axis = 1.2, cex.lab = 1.5, 
     col = "lightgray", border = "white",
     ylim = c(0, max(colQuantiles(tau_S_fit, probs = 0.95))),
     xlab = bquote(tau[italic(S)]), ylab = "Probability density", main = "Spawner observation error")
polygon(c(tau_S_seq, rev(tau_S_seq)),
        c(colQuantiles(tau_S_fit, probs = 0.05), rev(colQuantiles(tau_S_fit, probs = 0.95))),
        col = c1t, border = NA)
lines(tau_S_seq, colMedians(tau_S_fit), col = c1, lwd = 3)

rm(list = c("mod_name","tau_M_obs","tau_M_seq","mu_tau_M","sigma_tau_M","tau_M_fit",
            "tau_S_obs","tau_S_seq","mu_tau_S","sigma_tau_S","tau_S_fit","c1","c1t"))
## @knitr
dev.off()


#--------------------------------------------------------------------------------
# Tabular summary of 1-year-ahead forecasts by pop
#--------------------------------------------------------------------------------

mod_name <- "LCRchum_BH_fore"
life_stage <- "S"   # "S" = spawners, "M" = smolts

## @knitr forecast_df
life_cycle <- unlist(strsplit(mod_name, "_"))[1]
dat <- switch(life_cycle, SS = fish_data_SS_fore, SMS = fish_data_SMS_fore,
              LCRchum = fish_data_SMS_fore)
y <- dat$year
y1 <- min(dat$year[dat$forecast])
pop <- c(as.character(dat$pop[dat$year==y1]), "Total")
N_IPM <- do.call(extract1, list(as.name(mod_name), life_stage))
N1_IPM <- cbind(N_IPM[,dat$year==y1], rowSums(N_IPM[,dat$year==y1]))
N1_IPM_median <- round(colMedians(N1_IPM), 0)
N1_IPM_lo <- round(colQuantiles(N1_IPM, probs = 0.05), 0)
N1_IPM_up <- round(colQuantiles(N1_IPM, probs = 0.95), 0)

forecast_df <- data.frame(Population = pop, 
                          Estimate = paste0(N1_IPM_median, " (", N1_IPM_lo, ", ",
                                            N1_IPM_up, ")"))

rm("mod_name","life_stage","life_cycle","dat","y","y1","pop","N_IPM","N1_IPM",
   "N1_IPM_median","N1_IPM_lo","N1_IPM_up")
## @knitr


#--------------------------------------------------------------------------#
################# FIGURES FOR GENERIC SS OR SMS IPMs #######################
#--------------------------------------------------------------------------#

#--------------------------------------------------------------------
# Spawner-spawner / spawner-smolt-spawner #
# S-R curves and posterior distributions of parameters
#--------------------------------------------------------------------

mod_name <- "SMS_BH"

dev.new(width = 7, height = 7)
# png(filename=here("analysis","results",paste0("SR_",mod_name,".png")),
#     width=7, height=7, units="in", res=200, type="cairo-png")

## @knitr plot_SR_params
life_cycle <- unlist(strsplit(mod_name, "_"))[1]
dat <- switch(life_cycle, SS = fish_data_SS, SMS = fish_data_SMS)
SR_fun <- unlist(strsplit(mod_name, "_"))[2]

SR_eval <- function(alpha, Rmax = NULL, S, SR_fun) 
{
  switch(SR_fun,
         exp = alpha*S,
         BH = alpha*S/(1 + alpha*S/Rmax),
         Ricker = alpha*S*exp(-alpha*S/(exp(1)*Rmax)))
}

S <- colMedians(do.call(extract1, list(as.name(mod_name), "S")))
mu_alpha <- as.vector(do.call(extract1, list(as.name(mod_name), "mu_alpha")))
mu_Rmax <- as.vector(do.call(extract1, list(as.name(mod_name), "mu_Rmax")))
# S <- matrix(seq(0, quantile(dat$S_obs/dat$A, 0.9, na.rm = TRUE), length = 100),
#             nrow = length(mu_alpha), ncol = 100, byrow = TRUE)
scl <- ifelse(life_cycle == "SS", 1, 1e-6)
Smat <- matrix(seq(0, quantile(S/dat$A, 0.9, na.rm = TRUE), length = 100),
            nrow = length(mu_alpha), ncol = 100, byrow = TRUE)
R_ESU <- SR_eval(alpha = exp(mu_alpha), Rmax = exp(mu_Rmax), S = Smat, SR_fun = SR_fun)
alpha <- do.call(extract1, list(as.name(mod_name), "alpha"))
Rmax <- do.call(extract1, list(as.name(mod_name), "Rmax"))
R_pop <- sapply(1:ncol(alpha), function(i) {
  colMedians(SR_eval(alpha = alpha[,i], Rmax = Rmax[,i], S = Smat, SR_fun = SR_fun))
})

c1 <- "slategray4"
c1t <- transparent(c1, trans.val = 0.5)
c1tt <- transparent(c1, trans.val = 0.7)

par(mfrow = c(2,2), mar = c(5.1,5.1,1,1))

# Recruits vs. spawners
plot(Smat[1,], colMedians(R_ESU*scl), type = "l", lwd=3, col = c1, las = 1,
     cex.axis = 1.2, cex.lab = 1.5, xaxs = "i", yaxs = "i", #yaxt = "n", 
     ylim = range(R_pop*scl), xlab="Spawners", ylab = "")
for(i in 1:ncol(R_pop))
  lines(Smat[1,], R_pop[,i]*scl, col = c1t)
polygon(c(Smat[1,], rev(Smat[1,])), 
        c(colQuantiles(R_ESU*scl, probs = 0.05), 
          rev(colQuantiles(R_ESU*scl, probs = 0.95))), 
        col = c1tt, border = NA)
rug(S, col = c1)
title(ylab = ifelse(life_cycle=="SS", "Recruits", "Smolts (millions)"), 
      line = 3.5, cex.lab = 1.5)
text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "A", cex = 1.5)

# log(recruits/spawner) vs. spawners
plot(Smat[1,], colMedians(log(R_ESU/Smat)), type = "l", lwd=3, col = c1, las = 1,
     cex.axis = 1.2, cex.lab = 1.5, xaxs = "i", yaxs = "i", xlab="Spawners", 
     ylab = ifelse(life_cycle=="SS", "log(recruits/spawner)", "log(smolts/spawner)"),
     ylim = range(colQuantiles(log(R_ESU/Smat), probs = c(0.05,0.95)), na.rm = TRUE))
for(i in 1:ncol(R_pop))
  lines(Smat[1,], log(R_pop[,i]/Smat[1,]), col = c1t)
polygon(c(Smat[1,], rev(Smat[1,])),
        c(colQuantiles(log(R_ESU/Smat), probs = 0.05),
          rev(colQuantiles(log(R_ESU/Smat), probs = 0.95))),
        col = c1tt, border = NA)
rug(S, col = c1)
text(par("usr")[2], par("usr")[4], adj = c(1.5,1.5), "B", cex = 1.5)

# Posterior densities of log(alpha)
dd_ESU <- density(mu_alpha)
dd_pop <- vector("list", length(levels(dat$pop)))
for(i in 1:length(dd_pop))
  dd_pop[[i]] <- density(log(alpha[,i]))

plot(dd_ESU$x, dd_ESU$y, type = "l", lwd = 3, col = c1, las = 1, 
     cex.axis = 1.2, xlab = "", ylab = "", xaxs = "i",
     xlim = range(c(dd_ESU$x, sapply(dd_pop, function(m) m$x))),
     ylim = range(c(dd_ESU$y, sapply(dd_pop, function(m) m$y))))
for(i in 1:length(dd_pop))
  lines(dd_pop[[i]]$x, dd_pop[[i]]$y, col = c1t)
title(xlab = bquote(log(alpha)), ylab = "Probability density", line = 3.5, cex.lab = 1.5)
text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "C", cex = 1.5)

# Posterior densities of log(Rmax)
dd_ESU <- density(mu_Rmax)
dd_pop <- vector("list", length(levels(dat$pop)))
for(i in 1:length(dd_pop))
  dd_pop[[i]] <- density(log(Rmax[,i]))

plot(dd_ESU$x, dd_ESU$y, type = "l", lwd = 3, col = c1, las = 1, 
     cex.axis = 1.2, xlab = "", ylab = "", xaxs = "i",
     xlim = range(c(dd_ESU$x, sapply(dd_pop, function(m) m$x))),
     ylim = range(c(dd_ESU$y, sapply(dd_pop, function(m) m$y))))
for(i in 1:length(dd_pop))
  lines(dd_pop[[i]]$x, dd_pop[[i]]$y, col = c1t)
title(xlab = bquote(log(italic(R)[max])), ylab = "Probability density", 
      line = 3.5, cex.lab = 1.5)
text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "D", cex = 1.5)

rm(list=c("mod_name","life_cycle","SR_fun","scl","mu_alpha","mu_Rmax","S","R_ESU",
          "c1","c1t","c1tt","dd_ESU","dd_pop","alpha","Rmax","R_pop","Smat"))
## @knitr
# dev.off()


#--------------------------------------------------------------------------------
# Spawner-Spawner #
# Shared recruitment process errors (brood year productivity anomalies)
#--------------------------------------------------------------------------------

mod_name <- "SS_BH"

# dev.new(width=7,height=5)
png(filename=here("analysis","results",paste0("phi_", mod_name, ".png")),
    width=7, height=5, units="in", res=200, type="cairo-png")

y <- sort(unique(fish_data_SS$year))
phi <- do.call(extract1, list(as.name(mod_name), "phi"))

c1 <- "slategray4"
c1t <- transparent(c1, trans.val = 0.5)

par(oma = c(0,0.1,0,0))

plot(y, colMedians(phi), type = "n", las = 1, cex.axis = 1.2, cex.lab = 1.5,
     ylim = range(colQuantiles(phi, probs = c(0.05,0.95))), xaxs = "i", xaxt = "n",
     xlab = "Brood year", ylab = "Productivity anomaly")
abline(h = 0, col = "gray")
polygon(c(y, rev(y)), 
        c(colQuantiles(phi, probs = 0.05), rev(colQuantiles(phi, probs = 0.95))),
        col = c1t, border = NA)
lines(y, colMedians(phi), lwd = 3)
axis(side = 1, at = y[y %% 5 == 0], cex.axis = 1.2)
rug(y[y %% 5 != 0], ticksize = -0.02)

dev.off()
rm(list = c("mod_name","y","phi","c1","c1t"))


#--------------------------------------------------------------------------------
# Spawner-Smolt-Spawner #
# Shared recruitment and SAR process errors (brood year productivity anomalies)
#--------------------------------------------------------------------------------

mod_name <- "SMS_BH"

dev.new(width=6,height=8)
# png(filename=here("analysis","results",paste0("phi_", mod_name, ".png")),
#     width=6, height=8, units="in", res=200, type="cairo-png")

## @knitr plot_phi
y <- sort(unique(fish_data_SMS$year))

phi_M <- do.call(extract1, list(as.name(mod_name), "phi_M"))
sigma_M <- do.call(extract1, list(as.name(mod_name), "sigma_M"))
zeta_M <- do.call(stan_mean, list(as.name(mod_name), "zeta_M"))
epsilon_M <- t(sapply(sigma_M, function(x) x*zeta_M))

phi_MS <- do.call(extract1, list(as.name(mod_name), "phi_MS"))
mu_MS <- do.call(extract1, list(as.name(mod_name), "mu_MS"))
s_hat_MS <- plogis(sweep(phi_MS, 1, qlogis(mu_MS), "+"))
s_MS <- do.call(extract1, list(as.name(mod_name), "s_MS"))

c1 <- "slategray4"
c1t <- transparent(c1, trans.val = 0.5)
c1tt <- transparent(c1, trans.val = 0.7)

par(mfcol = c(2,1), mar = c(5.1,5.1,2,2),  oma = c(0,0.1,0,0))

# Smolt recruitment
plot(y, colMedians(phi_M), type = "n", las = 1, cex.axis = 1.2, cex.lab = 1.5,
     ylim = range(colQuantiles(phi_M, probs = c(0.05,0.95)),
                  colMedians(phi_M[,match(fish_data_SMS$year, y)] + epsilon_M)), 
     xlab = "Brood year", ylab = "Productivity anomaly", main = "Smolt recruitment",
     xaxs = "i", xaxt = "n")
abline(h = 0, col = "darkgray")
polygon(c(y, rev(y)), 
        c(colQuantiles(phi_M, probs = 0.05), rev(colQuantiles(phi_M, probs = 0.95))),
        col = c1tt, border = NA)
lines(y, colMedians(phi_M), col = c1t, lwd = 4)
for(j in levels(fish_data_SMS$pop))
{
  indx1 <- fish_data_SMS$pop == j
  indx2 <- y %in% fish_data_SMS$year[indx1]
  lines(y[indx2], colMedians(phi_M[,indx2] + epsilon_M[,indx1]), col = c1t)
}
axis(side = 1, at = y[y %% 5 == 0], cex.axis = 1.2)
rug(y[y %% 5 != 0], ticksize = -0.02)

# SAR
plot(y, colMedians(s_hat_MS), type = "n", las = 1, cex.axis = 1.2, cex.lab = 1.5,
     ylim = range(0, colQuantiles(s_hat_MS, probs = 0.95), colQuantiles(s_MS, probs = 0.95)), 
     xaxs = "i", xaxt = "n",
     xlab = "Outmigration year", ylab = "", main = "SAR")
mtext("Survival", side = 2, line = 3.7, cex = par("cex")*1.5)
polygon(c(y, rev(y)), 
        c(colQuantiles(s_hat_MS, probs = 0.05), rev(colQuantiles(s_hat_MS, probs = 0.95))),
        col = c1tt, border = NA)
lines(y, colMedians(s_hat_MS), col = c1t, lwd = 4)
for(j in levels(fish_data_SMS$pop))
  lines(fish_data_SMS$year[fish_data_SMS$pop==j], colMedians(s_MS[,fish_data_SMS$pop==j]), 
        col = c1t)
axis(side = 1, at = y[y %% 5 == 0], cex.axis = 1.2)
rug(y[y %% 5 != 0], ticksize = -0.02)

rm(list = c("mod_name","y","phi_M","phi_MS","s_hat_MS","c1","c1t",
            "sigma_M","zeta_M","epsilon_M"))
## @knitr
# dev.off()


#--------------------------------------------------------------------------------
# Spawner-Spawner #
# Posterior distributions of shared and unique recruitment process errors and 
# spawner observation error
#--------------------------------------------------------------------------------

mod_name <- "SS_Ricker"

# variance components
sigma <- do.call(extract1, list(as.name(mod_name), "sigma"))
sigma_phi <- do.call(extract1, list(as.name(mod_name), "sigma_phi"))
rho_phi <- do.call(extract1, list(as.name(mod_name), "rho_phi"))
sd_phi <- sqrt(sigma_phi^2 / sqrt(1 - rho_phi^2))
sigma_tot <- sqrt(sigma^2 + sd_phi^2)
tau <- do.call(extract1, list(as.name(mod_name), "tau"))

# densities
dd_sigma <- density(sigma)
dd_sd_phi <- density(sd_phi)
dd_sigma_tot <- density(sigma_tot)
dd_tau <- density(tau)

cols <- c("darkblue", "darkblue", "gray")

# dev.new(width=7,height=7)
png(filename=here("analysis","results",paste0("process_obs_SDs_", mod_name, ".png")),
    width=7, height=7, units="in", res=200, type="cairo-png")

par(oma = c(0,0.1,0,0))

plot(dd_sigma$x, dd_sigma$y, type = "l", col = cols[1], 
     las = 1, cex.axis = 1.2, cex.lab = 1.5, yaxs = "i",
     xlim = range(c(dd_sigma$x, dd_sd_phi$x, dd_sigma_tot$x, dd_tau$x)),
     ylim = range(c(dd_sigma$y, dd_sd_phi$y, dd_sigma_tot$y, dd_tau$y)),
     xlab = "Error SD", ylab = "Probability density")
# lines(dd_sd_phi$x, dd_sd_phi$y, col = cols[2])
lines(dd_sigma_tot$x, dd_sigma_tot$y, lwd = 3, col = cols[2])
lines(dd_tau$x, dd_tau$y, lwd = 3, col = cols[3])
legend("topright", lwd = c(1,3,3), col = cols, 
       legend = c(expression("unique " * italic(R) * " process error (" * sigma * ")"), 
                  #expression(sigma[phi]), 
                  expression("total " * italic(R) * " process error (" * sigma[tot] * ")"), 
                  expression(italic(S) * " observation error (" * tau * ")")))

dev.off()
rm(list = c("mod_name","sigma","sigma_phi","rho_phi","sd_phi","sigma_tot","tau",
            "dd_sigma","dd_sd_phi","dd_sigma_tot","dd_tau"))


