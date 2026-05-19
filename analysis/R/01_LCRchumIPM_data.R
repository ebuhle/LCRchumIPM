#===========================================================================
# SETUP
#===========================================================================

library(matrixStats)
library(Hmisc)
library(tibble)
library(dplyr)
library(tidyr)
library(forcats)
library(lubridate)
library(reshape2)
library(here)

#===========================================================================
# DATA PROCESSING
#===========================================================================

## @knitr data

#---------------------------------------------------------------------------
# Basic population info 
#---------------------------------------------------------------------------

# Population names
pop_names <- read.csv(here("data","pop_names.csv"), header = TRUE, stringsAsFactors = FALSE) %>% 
  arrange(east_to_west) %>% 
  mutate(recovery_pop = factor(recovery_pop, levels = unique(recovery_pop)))

# Recovery targets
recovery_targets <- read.csv(here("data","recovery_targets.csv"), header = TRUE, stringsAsFactors = FALSE)

# List of "bad" or questionable observations
bad_data <- read.csv(here("data","Data_Chum_estimates_to_censor.csv"), 
                     header = TRUE, stringsAsFactors = FALSE) %>% 
  rename(pop = Population) %>% filter(pop != "") %>%  # there are bad data in the bad data
  mutate(pop = gsub("_", " ", pop), Censored_data_file = Censored_data_file == "Y", 
         Censored_R_Script = Censored_R_Script == "Y")

# Start dates of hatcheries associated with populations
hatcheries <- read.csv(here("data","Hatchery_Programs.csv"), header = TRUE, stringsAsFactors = TRUE)

# Habitat size data (currently stream length)
# Assumptions:
# (1) Duncan Creek should really be Duncan Channel
# (2) most recent years are missing; assume same as last available
habitat_data <- read.csv(here("data","Data_Habitat_Spawning_Linear.csv"),
                         header = TRUE, stringsAsFactors = FALSE) %>% 
  rename(year = Return.Yr., strata = Strata, pop = Location.Reach, mi = Habitat_Length_miles) %>% 
  complete(nesting(strata, pop), year = min(year):year(Sys.Date())) %>% fill(mi) %>% 
  mutate(pop = gsub("Duncan Creek", "Duncan Channel", gsub("I205", "I-205", gsub("_", " ", pop))), 
         m = mi*1609) %>%  # convert to m
  select(strata, pop, year, m) %>% arrange(strata, pop, year)

# Proportion of "green" females in Duncan Channel
# Non-green (ripe or partial) females are assumed to have lower fecundity
# Proportion green females outside Duncan Channel assumed to = 1
# https://github.com/ebuhle/chumIPM/issues/5
green_female_data <- read.csv(here("data","Data_Duncan_Females_by_Condition.csv"),
                              header = TRUE, stringsAsFactors = FALSE) %>% 
  rename(year = BY, disposition = Channel_Disposition, sex = Sex, condition = Condition,
         N = Qty, comment = Comment) %>% 
  dcast(year ~ condition, value.var = "N", fun.aggregate = sum) %>% 
  mutate(pop = "Duncan Channel", p_G_obs = Green / (Green + Ripe + Partial))

#---------------------------------------------------------------------------
# Spawner abundance data
#---------------------------------------------------------------------------

# Assumptions:
# (1) NAs in hatchery dispositions are really zeros
# (2) NAs in Duncan Creek and Duncan Channel are really zeros
# (3) All other NAs are real missing observations
# (4) When calculating tau_S_obs, the observation error SD of log(S_obs), assume
#     Abund.Mean and Abund.SD are the mean and SD of a lognormal posterior distribution
#     of spawner abundance based on the sample
spawner_data <- read.csv(here("data","Data_Abundance_Spawners_Chum.csv"), 
                         header = TRUE, stringsAsFactors = FALSE) %>% 
  rename(year = Return.Yr., strata = Strata, location = Location.Reach, 
         disposition = Disposition, method = Method, 
         S_obs = Abund.Mean, SD = Abund.SD) %>% 
  mutate(disposition = gsub("I205", "I-205", gsub("_", " ", disposition)),
         location = gsub("I205", "I-205", gsub("_", " ", location)),
         S_obs = replace(S_obs, is.na(S_obs) & grepl("Hatchery|Duncan", disposition), 0),
         tau_S_obs = sqrt(log((SD/S_obs)^2 + 1))) %>% 
  select(year:location, disposition, method, S_obs, SD, tau_S_obs) %>% 
  arrange(strata, location, year)

# total broodstock translocations: 
# total spawners taken from a given return location to any other disposition
# unless location is Duncan Creek and disposition is Duncan Channel,
# which is considered natural recruitment
# (summarized by *location*)
broodstock_data <- spawner_data %>% group_by(location, year) %>% 
  mutate(location = replace(location, location == "Duncan Creek", "Duncan Channel"),
         disposition = replace(disposition, disposition == "Duncan Creek", "Duncan Channel")) %>% 
  summarize(B_take_obs = sum(S_obs[disposition != location])) %>% 
  rename(pop = location) %>% as.data.frame()

# distribution of translocated spawners:
# spawners taken from a given return location to each specific disposition
# unless location is Duncan Creek and disposition is Duncan Channel,
# or the disposition is not a population included in the model
# (summarized by *location*)
translocation_data <- spawner_data %>% 
  mutate(location = replace(location, location == "Duncan Creek", "Duncan Channel"),
         disposition = replace(disposition, disposition == "Duncan Creek", "Duncan Channel")) %>% 
  filter(disposition != location) %>% 
  group_by(location, disposition, year) %>% summarize(S_obs = sum(S_obs)) %>% 
  dcast(location + year ~ disposition, value.var = "S_obs", fun.aggregate = sum) %>% 
  rename(pop = location) %>% 
  select(pop, year, `Duncan Channel`, `Duncan Hatchery`, `Lewis Hatchery`, `Grays Hatchery`)

# total spawners:
# all spawners with a given disposition, regardless of original return location
# (summarized by *disposition*)
# drop Duncan Creek, Big Creek Hatchery, Peterson RSI, Sea Resources, Skamokawa
spawner_data_agg <- spawner_data %>% 
  rename(pop = disposition) %>% 
  filter(!pop %in% c("Duncan Creek","Big Creek Hatchery","Peterson RSI","Sea Resources","Skamokawa")) %>% 
  group_by(pop, year) %>% 
  summarize(S_obs = sum(S_obs), tau_S_obs = unique(tau_S_obs)) %>% 
  full_join(broodstock_data, by = c("pop","year")) %>%
  full_join(translocation_data, by = c("pop","year")) %>%
  mutate(B_take_obs = replace_na(B_take_obs, 0)) %>%
  mutate(across(matches("Channel|Hatchery"), ~replace_na(.x, 0))) %>%
  rename_at(vars(matches("Channel|Hatchery")), list(~paste0("n_B", .x, "_obs"))) %>%
  arrange(pop, year) %>% as.data.frame()

#---------------------------------------------------------------------------
# Spawner age, sex, and origin frequencies (BioData)
#---------------------------------------------------------------------------

# Spawners from Duncan Channel or unknown natural origin are considered "W", all others "H"
bio_data <- read.csv(here("data","Data_BioData_Spawners_Chum.csv"), 
                     header = TRUE, stringsAsFactors = FALSE) %>% 
  rename(year = Return.Yr., strata = Strata, location = Location.Reach, 
         disposition = Disposition, origin = Origin, sex = Sex, age = Age, count = Count) %>% 
  mutate(disposition = gsub("I205", "I-205", gsub("_", " ", disposition)),
         location = gsub("I205", "I-205", gsub("_", " ", location)),
         origin = gsub("_", " ", origin),
         count = replace_na(count, 0), sex = substring(sex,1,1),
         HW = ifelse((grepl("Natural spawner|Duncan Channel", origin)), "W", "H")) %>% 
  select(year:location, disposition, origin, HW, sex:count) %>%
  arrange(strata, location, year, origin, age, sex)

# spawner age structure
bio_data_age <- bio_data %>% 
  dcast(year + disposition ~ age, value.var = "count", fun.aggregate = sum) %>% 
  rename(pop = disposition)

# spawner sex composition
bio_data_sex <- bio_data %>% 
  dcast(year + disposition ~ sex, value.var = "count", fun.aggregate = sum) %>% 
  rename(pop = disposition) %>% select(year:pop, M, `F`)

# spawner origin composition
# Duncan Creek escapement is considered as pop == "Duncan Channel"
# ignore Big Creek Hatchery-origin spawners for now
bio_data_origin <- bio_data %>% 
  mutate(pop = replace(disposition, disposition == "Duncan Creek", "Duncan Channel")) %>% 
  dcast(year + pop ~ origin, value.var = "count", fun.aggregate = sum) %>% 
  select(year, pop, `Natural spawner`, `Duncan Channel`, `Duncan Hatchery`, 
         `Lewis Hatchery`, `Grays Hatchery`) %>% 
  rename_at(vars(matches("Channel|Hatchery")), list(~paste0("n_O", .x, "_obs")))

#---------------------------------------------------------------------------
# Juvenile abundance data
#---------------------------------------------------------------------------

# Assumptions:
# (1) Duncan North + Duncan South = Duncan Channel, so the former two are redundant 
#     (not really an assumption, although the equality isn't perfect in all years)
# (2) Drop Big Creek (note that this is coded "Big Creek Hatchery" in other datasets)
#     and Peterson RSI
# (3) When calculating the observation error of log(M_obs), tau_M_obs, assume
#     Abund_Median and Abund_SD are the median and SD of a lognormal posterior 
#     distribution of smolt abundance based on the sample
# (4) If Abund_SD == 0 (when Analysis=="Census": some years in Duncan_Channel and 
#     Hamilton_Channel) treat as NA
juv_data <- read.csv(here("data","Data_Abundance_Juveniles_Chum.csv"), 
                     header = TRUE, stringsAsFactors = FALSE) %>% 
  rename(brood_year = Brood.Year, year = Outmigration.Year, strata = Strata, 
         location = Location.Reach, origin = Origin, trap_type = TrapType, 
         analysis = Analysis, partial_spawners = Partial.Spawners, raw_catch = RawCatch,
         M_obs = Abund_Median, mean = Abund_Mean, SD = Abund_SD, 
         L95 = Abund_L95, U95 = Abund_U95, CV = Abund_CV, comments = Comments) %>% 
  filter(!grepl("Big_Creek|Peterson_RSI", location) & !grepl("Sea_Resources_Hatchery", origin)) %>% 
  mutate(location = gsub("I205", "I-205", gsub("_", " ", location)),
         origin = sapply(gsub("_", " ", origin), function(x)
           paste(rev(strsplit(x, " ")[[1]]), collapse = " ")), # names inconsistent w/ bio_data 
         tau_M_obs = replace(sqrt(log((SD/mean)^2 + 1)), SD == 0, NA)) %>% 
  select(strata, location, year, brood_year, origin:CV, tau_M_obs, comments) %>% 
  arrange(strata, location, year)

# drop redundant pops and cases with leading or trailing NAs in M_obs
# change NAs in hatchery pops to 0
# use only pooled Duncan Channel data for now, not North / South
# drop final year if it contains only hatchery releases
#   (kludge to avoid rows with only hatchery M_obs in fish_data)
head_noNA <- function(x) { cumsum(!is.na(x) & x > 0) > 0 }
juv_data_incl <- juv_data %>% 
  mutate(pop = ifelse(grepl("Hatchery", origin), origin, location), .after = strata,
         M_obs = replace(M_obs, is.na(M_obs) & grepl("Hatchery", pop), 0)) %>% 
  filter(!pop %in% c("Duncan North","Duncan South","Big Creek Hatchery")) %>% 
  select(pop, year, M_obs, tau_M_obs) %>% arrange(pop, year) %>% 
  group_by(pop, year) %>% summarize(M_obs = sum(M_obs), tau_M_obs = unique(tau_M_obs)) %>%
  ungroup() %>% group_by(pop) %>% filter(head_noNA(M_obs) & rev(head_noNA(rev(M_obs)))) %>%
  group_by(year) %>% mutate(all_H = all(grepl("Hatchery", pop))) %>% ungroup() %>% 
  filter(!(year == max(year) & all_H)) %>% select(-all_H) %>% as.data.frame()

#---------------------------------------------------------------------------
# Fish data formatted for salmonIPM
#---------------------------------------------------------------------------

# Drop age-2 and age-6 samples (each is < 0.1% of aged spawners)
# Drop Duncan Creek
# Change S_obs and tau_S_obs to NA in Hamilton Channel 2011-2012 based on
# https://github.com/ebuhle/chumIPM/issues/6#issuecomment-807885445
# Replace unknown tau_S_obs in hatchery pops and Duncan Channel with 0.01 
#  (nearly smallest value from natural pops)
# Replace unknown tau_M_obs in hatchery pops with 0.05 
#  (37th percentile of values from smolt traps)
# Pad data as necessary so Grays MS, Grays WF, and Grays CJ have the same set of years
#  (since their estimated smolts will be summed)  
# Censor Grays Hatchery origin frequencies pre-2007
fish_data_all <- full_join(spawner_data_agg, bio_data_age, by = c("pop","year")) %>% 
  full_join(bio_data_origin, by = c("pop","year")) %>% 
  full_join(bio_data_sex, by = c("pop","year")) %>% 
  full_join(juv_data_incl, by = c("pop","year")) %>%
  full_join({ complete(select(filter(., grepl("Grays", pop)), c(pop,year)), pop, year) },
            by = c("pop","year")) %>%
  left_join(habitat_data, by = c("pop","year")) %>% 
  left_join(green_female_data, by = c("pop","year")) %>% 
  rename(A = m, n_O0_obs = `Natural spawner`, n_M_obs = M, n_F_obs= `F`) %>% 
  rename_at(vars(contains("Age-")), list(~paste0(sub("Age-","n_age",.), "_obs"))) %>% 
  select(-c(n_age2_obs, n_age6_obs)) %>% 
  filter(!pop %in% c("Duncan Creek", "Sea Resources Hatchery")) %>% 
  mutate(pop = droplevels(factor(pop, levels = pop_names$pop)), # order E-W
         A = replace(A, grepl("Hatchery", pop), 1),
         S_obs = replace(S_obs, pop == "Hamilton Channel" & year %in% 2011:2012, NA),
         tau_S_obs = replace(tau_S_obs, pop == "Hamilton Channel" & year %in% 2011:2012, NA),
         tau_S_obs = replace(tau_S_obs, grepl("Hatchery|Duncan", pop), 0.01), # kludge
         # tau_M_obs = replace(tau_M_obs, grepl("Hatchery", pop), 0.05), # kludge
         across(starts_with("n_O"), ~ifelse(pop == "Grays Hatchery" & year < 2007, 0, .)),
         B_take_obs = replace_na(B_take_obs, 0),
         p_G_obs = replace_na(p_G_obs, 1), F_rate = 0) %>%
  mutate(n_H_obs = rowSums(across(matches("n_O.*Hatchery_obs"))),
         n_W_obs = rowSums(across(starts_with("n_O"))) - n_H_obs, 
         .before = n_O0_obs) %>% 
  mutate_at(vars(contains("n_")), ~replace_na(., 0)) %>%
  do({ 
    lev <- levels(.$pop)
    .cols <- grepl("Hatchery|Channel", names(.))
    .fn <- function(.x) {
      rgx <- "n_[O|B](.*)_obs"
      streplace <- sub(rgx, "\\1", grep(rgx, .x, value = TRUE))
      sub(streplace, match(streplace, lev), .x)
    }
    setNames(., replace(names(.), .cols, sapply(names(.)[.cols], .fn)))
  }) %>% 
  select(pop, year, A, S_obs, tau_S_obs, M_obs, tau_M_obs, n_age3_obs:n_F_obs, 
         p_G_obs, B_take_obs, starts_with("n_B"), F_rate) %>% 
  arrange(pop, year) 

# drop cases with initial NAs in S_obs and M_obs
fish_data <- fish_data_all %>% group_by(pop) %>% 
  filter(head_noNA(S_obs) | head_noNA(M_obs)) %>%
  as.data.frame()

# assign Grays WF and Grays CJ smolts to the downstream trap in Grays MS
# where they will be counted (or double-counted, in the case of Grays CJ),
# assuming no mortality between the upstream tributary and the downstream trap
fish_data <- fish_data %>%
  add_column(downstream_trap = NA, .after = "tau_M_obs") %>% 
  mutate(downstream_trap = replace(downstream_trap, pop %in% c("Grays WF","Grays CJ"),
                                   which(pop == "Grays MS" & year > 2000))) %>% 
  mutate(pop_type = factor(ifelse(grepl("Hatchery", pop), "hatchery", "natural"), 
                           levels = c("natural","hatchery")),
         .after = pop)

# pad data with future years to generate forecasts
N_year_fore <- 2
fish_data_fore <- fish_data %>% group_by(pop) %>%
  slice(rep(n(), max(fish_data$year) + N_year_fore - max(year))) %>%
  reframe(year = (unique(year) + 1):(max(fish_data$year) + N_year_fore), forecast = TRUE,
          S_obs = NA, tau_S_obs = NA, M_obs = NA, tau_M_obs = NA, downstream_trap = NA, 
          p_G_obs = 1, B_take_obs = 0, F_rate = 0) %>%
  full_join(mutate(fish_data, forecast = FALSE, downstream_trap = NA)) %>% 
  mutate_at(vars(starts_with("n_")), ~ replace_na(., 0)) %>%
  arrange(pop, year) %>% fill(A, pop_type, .direction = "down") %>%
  select(pop, pop_type, year, forecast, A, S_obs, tau_S_obs, M_obs, tau_M_obs,
         downstream_trap, n_age3_obs:n_F_obs, p_G_obs, B_take_obs, 
         starts_with("n_B"), F_rate) %>% 
  as.data.frame()

# assign Grays_WF and Grays_CJ smolts to the downstream trap in Grays_MS
# (need to do this again b/c row indices have changed)
fish_data_fore <- fish_data_fore %>%
  mutate(downstream_trap = replace(downstream_trap, pop %in% c("Grays WF","Grays CJ"),
                                   which(pop == "Grays MS" & year > 2000)))

# forecast scenario 1:
# no broodstock removals or hatchery smolt releases
fish_data_foreH0 <- fish_data_fore %>% group_by(pop) %>% 
  mutate(M_obs = replace(M_obs, pop_type == "hatchery" & forecast, 0)) %>% 
  ungroup() %>% as.data.frame()

# forecast scenario 2:
# broodstock removal rates and hatchery smolt releases at maximum observed
fish_data_foreHmax <- fish_data_fore %>% group_by(pop) %>% 
  mutate(B_rate_obs = B_take_obs / (S_obs + B_take_obs), .after = B_take_obs) %>% 
  mutate(F_rate = replace(F_rate, forecast, max(B_rate_obs, na.rm = TRUE)),
         M_obs = replace(M_obs, pop_type == "hatchery" & forecast, 
                         ifelse(all(is.na(M_obs)), NA, max(M_obs, na.rm = TRUE)))) %>% 
  ungroup() %>% select(-B_rate_obs) %>% as.data.frame()

#---------------------------------------------------------------------------
# Fecundity data
#---------------------------------------------------------------------------

# Note that L95% and U95% are reversed
fecundity <- read.csv(here("data","Data_ChumFecundity_fromHatcheryPrograms.csv"),
                      header = TRUE, stringsAsFactors = FALSE) %>% 
  rename(stock = Stock, year = BY, ID = Female.., age_E = Age, L95 = U95., U95 = L95.,
         reproductive_effort = Reproductive.Effort, E_obs = Estimated.Fecundity,
         mean_mass = Green.egg.avg.weight, comments = Comments) %>% 
  mutate(stock = factor(stock))

# drop cases with age not in c(3,4,5), with estimated fecundity missing, 
#   or with reproductive effort <= 16%
# add strata based on stock: Grays -> Coastal, I-205 -> Cascade, Lower Gorge -> Gorge
fecundity_data <- fecundity %>% filter(age_E %in% 3:5 & !is.na(E_obs) & 
           !is.na(reproductive_effort) & reproductive_effort > 16) %>% 
  mutate(strata = recode(stock, Grays = "Coastal", `I-205` = "Cascade", `Lower Gorge` = "Gorge")) %>% 
  select(strata, year, ID, age_E, E_obs) %>% arrange(strata, year, age_E) 

#---------------------------------------------------------------------------
# Pairwise distance data
#---------------------------------------------------------------------------

# convert ft to km
# add missing populations: 
#  Duncan Channel, Duncan Hatchery (same location as Duncan Creek)
#  Grays Hatchery (same location as Grays MS)
# extract lower triangle to a tidy data frame
pairwise_dist <- read.csv(here("data","Chum_Pairwise_Points_Mat_Final.csv"), 
                          header = TRUE, row.names = 1) * 0.0003048
rownames(pairwise_dist) <- pop_names$pop[match(rownames(pairwise_dist), pop_names$pairwise_dist_name)]
colnames(pairwise_dist) <- rownames(pairwise_dist)
newcols <- pairwise_dist[,c("Duncan Creek", "Duncan Creek", "Grays MS")]
colnames(newcols) <- c("Duncan Channel", "Duncan Hatchery", "Grays Hatchery")
pairwise_dist <- cbind(pairwise_dist, newcols)
newrows <- pairwise_dist[c("Duncan Creek", "Duncan Creek", "Grays MS"),]
rownames(newrows) <- c("Duncan Channel", "Duncan Hatchery", "Grays Hatchery")
pairwise_dist <- rbind(pairwise_dist, newrows)
indx <- which(lower.tri(pairwise_dist, diag = FALSE), arr.ind = TRUE)
pairwise_data <- data.frame(pop1 = droplevels(factor(rownames(pairwise_dist)[indx[,2]],
                                                     levels = pop_names$pop)),
                            pop2 = droplevels(factor(rownames(pairwise_dist)[indx[,1]],
                                                     levels = pop_names$pop)),
                            dist = pairwise_dist[indx])
dist_mouth_data <- pairwise_data %>% filter(pop1 == "Columbia Mouth") %>% select(-pop1) %>% 
  rename(pop = pop2, dist_mouth = dist) %>% mutate(dist_mouth_std = scale(dist_mouth))

# #---------------------------------------------------------------------------
# # Environmental covariates
# #---------------------------------------------------------------------------
#
# # PDO
# PDO_data <- read.csv(file = "https://www.ncdc.noaa.gov/teleconnections/pdo/data.csv",
#                      header = TRUE, skip = 1, stringsAsFactors = FALSE) %>% 
#   separate(Date, into = c("year","month"), sep = 4, convert = TRUE) %>% 
#   rename(PDO = Value) %>% group_by(year) %>% summarize(PDO = mean(PDO)) %>% 
#   as.data.frame() %>% arrange(year)
# 
# # NPGO
# NPGO_data <- read.csv(file = "http://www.o3d.org/npgo/npgo.php", comment.char = "#", 
#                       stringsAsFactors = FALSE) %>% 
#   tail(nrow(.) - 2) %>% head(nrow(.) - 3) %>% # drop html tags in first 2 & last 3 rows
#   rename(X = X.html.) %>% mutate(X = trimws(X)) %>% 
#   separate(col = X, into = c("year","month","NPGO"), sep = "  ", convert = TRUE) %>% 
#   group_by(year) %>% summarize(NPGO = mean(NPGO)) %>% as.data.frame() %>% arrange(year)
# 
# # Covariate data
# env_data <- PDO_data %>% left_join(NPGO_data, by = "year") %>% 
#   filter(year %in% fish_data$year) %>% arrange(year)
## @knitr
