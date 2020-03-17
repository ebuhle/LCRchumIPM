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

# Spawner abundance data
spawner_data <- read.csv(here("data","Data_ChumSpawnerAbundance_2019-12-12.csv"), header = TRUE) %>% 
  rename(species = Species, stage = LifeStage, year = Return.Yr., strata = Strata,
         location = Location.Reach, disposition = Disposition, method = Method,
         S = Abund.Mean, SD = Abund.SD) %>% arrange(strata, location, year)

# Spawner age-, sex-, and origin-frequency (aka BioData)
bio_data <- read.csv(here("data","Data_ChumSpawnerBioData_2019-12-12.csv"), header = TRUE) %>% 
  rename(species = Species, stage = LifeStage, year = Return.Yr., strata = Strata,
         location = Location.Reach, disposition = Disposition, origin = Origin,
         sex = Sex, age = Age, count = Count) %>% 
  arrange(strata, location, year, origin, age, sex)


#----------------------
# Data Exploration
#----------------------



