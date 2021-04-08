#===========================================================================
# SETUP
#===========================================================================

options(device = ifelse(.Platform$OS.type == "windows", "windows", "quartz"))

library(matrixStats)
library(Hmisc)
library(tibble)
library(dplyr)
library(tidyr)
library(reshape2)
library(here)

#===========================================================================
# DATA PROCESSING
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
spawner_data <- read.csv(here("data","Data_Abundance_Spawners_Chum_2021-04-06.csv"), 
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
bio_data <- read.csv(here("data","Data_BioData_Spawners_Chum_2021-04-06.csv"), 
                     header = TRUE, stringsAsFactors = FALSE) %>% 
  rename(year = Return.Yr., strata = Strata, location = Location.Reach, 
         disposition = Disposition, origin = Origin, sex = Sex, age = Age, count = Count) %>% 
  mutate(strata = replace(strata, disposition == "Duncan_Channel" & strata != "Gorge", "Gorge"),
         count = replace(count, is.na(count), 0), sex = substring(sex,1,1),
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

# sex composition of spawners, regardless of origin
bio_data_sex <- bio_data %>% 
  dcast(year + strata + disposition ~ sex, value.var = "count", fun.aggregate = sum) %>% 
  rename(pop = disposition) %>% select(year:pop, M, `F`)

# Juvenile abundance data
# Assumptions:
# (1) Duncan_North + Duncan_South = Duncan_Channel, so the former two are redundant 
#     (not really an assumption, although the equality isn't perfect in all years)
# (2) When calculating the observation error of log(M_obs), tau_M_obs, assume
#     Abund_Median and Abund_SD are the median and SD of a lognormal posterior 
#     distribution of smolt abundance based on the sample
# (3) If Abund_SD == 0 (when Analysis=="Census": some years in Duncan_Channel and 
#     Hamilton_Channel) treat as NA
juv_data <- read.csv(here("data", "Data_Abundance_Juveniles_Chum_2021-04-07.csv"), 
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
  full_join(bio_data_sex, by = c("strata","pop","year")) %>% 
  full_join(juv_data_incl, by = c("strata","pop","year")) %>%
  full_join({ complete(select(filter(., strata == "Coastal"), c(strata,pop,year)), 
                       pop, year, fill = list(strata = "Coastal")) },
            by = c("strata","pop","year")) %>%
  rename_at(vars(contains("Age-")), list(~ paste0(sub("Age-","n_age",.), "_obs"))) %>% 
  select(-c(n_age2_obs, n_age6_obs)) %>% 
  rename(n_H_obs = H, n_W_obs = W, n_M_obs = M, n_F_obs= `F`) %>% 
  mutate(A = 1, fit_p_HOS = NA, F_rate = 0) %>% 
  mutate_at(vars(contains("n_")), ~ replace(., is.na(.), 0)) %>% 
  filter(!grepl("Hatchery|Duncan_Creek", pop)) %>% 
  mutate(S_obs = replace(S_obs, pop == "Hamilton_Channel" & year %in% 2011:2012, NA),
         tau_S_obs = replace(tau_S_obs, pop == "Hamilton_Channel" & year %in% 2011:2012, NA),
         strata = factor(strata), pop = factor(pop), 
         B_take_obs = replace(B_take_obs, is.na(B_take_obs), 0)) %>%
  select(strata, pop, year, A, S_obs, tau_S_obs, M_obs, tau_M_obs, n_age3_obs:n_F_obs, 
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


#===========================================================================
# DATA EXPLORATION
#===========================================================================

EDA <- FALSE

if(EDA)
{
  
  # Time series of sex ratio by population
  windows()
  bio_data %>% filter(HW=="W" & !grepl("Hatchery|Duncan_Creek", disposition)) %>% 
    group_by(disposition, year, sex) %>% 
    summarize(n = sum(count)) %>% 
    dcast(disposition + year ~ sex, value.var = "n", fun.aggregate = sum) %>% 
    mutate(total = `F` + M) %>% 
    data.frame(., with(., binconf(x = `F`, n = total))) %>%
    rename(prop_female = PointEst) %>% 
    ggplot(aes(x = year, y = prop_female, ymin = Lower, ymax = Upper)) + 
    geom_abline(intercept = 0.5, slope = 0, color = "gray") + 
    geom_point(size = 2) + geom_line() + geom_errorbar(width = 0) +
    facet_wrap(vars(disposition), ncol = 4) + theme_bw() +
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
  
}
