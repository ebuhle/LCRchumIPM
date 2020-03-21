options(device = windows)

library(salmonIPM)
library(rstan)
library(shinystan)
library(matrixStats)
library(dplyr)
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
         HW = disposition_HW$HW[match(disposition, disposition_HW$disposition)]) %>% 
  select(species:location, pop, disposition, HW, method:SD) %>% 
  arrange(strata, location, year)

# Spawner age-, sex-, and origin-frequency (aka BioData)
bio_data <- read.csv(here("data","Data_ChumSpawnerBioData_2019-12-12.csv"), header = TRUE) %>% 
  rename(species = Species, stage = LifeStage, year = Return.Yr., strata = Strata,
         location = Location.Reach, disposition = Disposition, origin = Origin,
         sex = Sex, age = Age, count = Count) %>% 
  mutate(pop = location_pop$pop[match(location, location_pop$location)],
         HW = disposition_HW$HW[match(disposition, disposition_HW$disposition)],
         count = ifelse(is.na(count), 0, count)) %>% 
  select(species:location, pop, disposition, HW, origin:count) %>% 
  arrange(strata, location, year, origin, age, sex)

bio_data_age <- bio_data %>% select(-HW, -sex) %>% 
  dcast(species + stage + year + strata + location + pop + disposition ~ age, 
        value.var = "count", fun.aggregate = sum)

bio_data_HW <- bio_data %>% select(-age, -sex) %>% 
  dcast(species + stage + year + strata + location + disposition ~ HW, value.var = "count",
        fun.aggregate = sum)


#----------------------
# Data Exploration
#----------------------



