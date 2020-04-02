options(device = windows)

library(salmonIPM)
library(rstan)
library(shinystan)
library(matrixStats)
library(dplyr)
library(reshape2)
library(yarrr)
library(corrplot)
library(magicaxis)
library(zoo)
library(here)


#===========================================================================
# DATA
#===========================================================================

# Mapping of location to population
location_pop <- read.csv(here("data","Location.Reach_Population.csv"), header = TRUE) %>% 
  rename(strata = Strata, location = Location.Reach, pop = Population)

# Mapping of disposition to hatchery vs. wild (i.e., broodstock vs. natural spawner)
disposition_HW <- read.csv(here("data","Disposition_HW.csv"), header = TRUE) %>% 
  rename(disposition = Disposition) %>% arrange(HW)

# Start dates of hatcheries associated with populations
hatchery_start_dates <- read.csv(here("data","Hatchery_Start_Dates.csv"), header = TRUE)

# Spawner abundance data
spawner_data <- read.csv(here("data","Data_ChumSpawnerAbundance_2019-12-12.csv"), header = TRUE) %>% 
  rename(species = Species, stage = LifeStage, year = Return.Yr., strata = Strata,
         location = Location.Reach, disposition = Disposition, method = Method,
         S = Abund.Mean, SD = Abund.SD) %>% 
  mutate(pop = location_pop$pop[match(location, location_pop$location)],
         disposition_HW = disposition_HW$HW[match(disposition, disposition_HW$disposition)]) %>% 
  select(species:location, pop, disposition, disposition_HW, method:SD) %>% 
  arrange(strata, location, year)

spawner_data_agg <- aggregate(S ~ species + stage + year + strata + location + pop + disposition_HW,
                              data = spawner_data, FUN = sum, na.rm = TRUE) %>% 
  dcast(species + stage + year + strata + pop ~ disposition_HW, value.var = "S",
        fun.aggregate = sum, na.rm = TRUE) %>% rename(S_obs = W, B_take_obs = H) %>% 
  arrange(strata, pop, year)

spawner_data_agg$S_obs <- spawner_data_agg$S_obs + spawner_data_agg$B_take_obs


# Spawner age-, sex-, and origin-frequency (aka BioData)
bio_data <- read.csv(here("data","Data_ChumSpawnerBioData_2019-12-12.csv"), header = TRUE) %>% 
  rename(species = Species, stage = LifeStage, year = Return.Yr., strata = Strata,
         location = Location.Reach, disposition = Disposition, origin = Origin,
         sex = Sex, age = Age, count = Count) %>% 
  mutate(pop = location_pop$pop[match(location, location_pop$location)],
         origin_HW = ifelse(origin == "Natural_spawner", "W", "H"),
         count = ifelse(is.na(count), 0, count)) %>% 
  select(species:location, pop, disposition, origin, origin_HW, sex:count) %>%
  arrange(strata, location, year, origin, age, sex)

bio_data_age <- bio_data %>% 
  dcast(species + stage + year + strata + pop ~ age, value.var = "count", 
        fun.aggregate = sum)

bio_data_origin <- bio_data %>% 
  dcast(species + stage + year + strata + pop ~ origin_HW, value.var = "count", 
        fun.aggregate = sum)

# Fish data formatted for salmonIPM
# Use A = 1 for now (so Rmax in units of spawners)
fish_data <- inner_join(spawner_data_agg, bio_data_age, 
                        by = c("species","stage","year","strata","pop")) %>% 
  inner_join(bio_data_origin, by = c("species","stage","year","strata","pop")) %>% 
  rename_at(vars(contains("Age-")), funs(paste0(sub("Age-","n_age",.), "_obs"))) %>% 
  rename(n_H_obs = H, n_W_obs = W) %>% mutate(A = 1, fit_p_HOS = NA, F_rate = 0) %>% 
  select(species, stage, strata, pop, year, A, S_obs, n_age2_obs:n_W_obs, 
         fit_p_HOS, B_take_obs, F_rate) %>% arrange(strata, pop, year) 

for(i in 1:nrow(fish_data)) {
  indx <- match(fish_data$pop[i], hatchery_start_dates$pop)
  fish_data$fit_p_HOS[i] <- ifelse((!is.na(indx) & 
                                     fish_data$year[i] >= hatchery_start_dates$start_brood_year[indx] + 2) |
                                     fish_data$n_H_obs[i] > 0, 1, 0)
}


#----------------------
# Data Exploration
#----------------------




#===========================================================================
# FIT RETROSPECTIVE MODELS
#===========================================================================

#--------------------------------------------------------------
# Spawner-to-spawner IPM
#--------------------------------------------------------------

SS_BH <- salmonIPM(fish_data = fish_data, stan_model = "IPM_SS_pp", 
                    pars = c("mu_alpha","sigma_alpha","alpha",
                             "mu_Rmax","sigma_Rmax","Rmax","rho_alphaRmax",
                             "sigma_phi","rho_phi","phi",
                             "mu_p","sigma_gamma","R_gamma","gamma",
                             "sigma_p","R_p","p",
                             "p_HOS","sigma","tau","S","R","q","LL"),
                    chains = 3, iter = 1500, warmup = 500,
                    control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))

print(SS_BH, prob = c(0.025,0.5,0.975),
      pars = c("alpha","Rmax","phi","p_HOS","q","gamma","p","S","R","LL"), 
      include = FALSE, use_cache = FALSE)

launch_shinystan(SS_BH)


