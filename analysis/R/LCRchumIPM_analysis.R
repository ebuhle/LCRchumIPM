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
  rename(n_H_obs = H, n_W_obs = W) %>% mutate(A = 1, F_rate = 0) %>% 
  select(species, stage, strata, pop, year, A, S_obs, n_age2_obs:n_W_obs, B_take_obs, F_rate)


#----------------------
# Data Exploration
#----------------------




