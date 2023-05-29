#===========================================================================
# SETUP
#===========================================================================

options(device = ifelse(.Platform$OS.type == "windows", "windows", "quartz"))

library(matrixStats)
library(Hmisc)
library(tibble)
library(dplyr)
library(tidyr)
library(forcats)
library(reshape2)
library(here)

#===========================================================================
# DATA PROCESSING
#===========================================================================

## @knitr data
# Population names
pop_names <- read.csv(here("data","pop_names.csv"), header = TRUE, stringsAsFactors = FALSE) %>% 
  arrange(east_to_west)

# List of "bad" or questionable observations
bad_data <- read.csv(here("data","Data_Chum_estimates_to_censor_2021-04-16.csv"), 
                     header = TRUE, stringsAsFactors = FALSE) %>% 
  rename(pop = Population) %>% filter(pop != "") %>%  # there are bad data in the bad data
  mutate(pop = gsub("_", " ", pop), Censored_data_file = Censored_data_file == "Y", 
         Censored_R_Script = Censored_R_Script == "Y")

# Start dates of hatcheries associated with populations
hatcheries <- read.csv(here("data","Hatchery_Programs.csv"), header = TRUE, stringsAsFactors = TRUE)

# Habitat size data (currently stream length)
# Assumptions:
# (1) Duncan Creek should really be Duncan Channel
# (2) 2020 is missing; assume same as 2019
habitat_data <- read.csv(here("data","Data_Habitat_Spawning_Linear_2021-04-15.csv"),
                         header = TRUE, stringsAsFactors = FALSE) %>% 
  rename(year = Return.Yr., strata = Strata, pop = Location.Reach, mi = Habitat_Length_miles) %>% 
  complete(nesting(strata, pop), year = min(year):2021) %>% fill(mi) %>% 
  mutate(pop = gsub("Duncan Creek", "Duncan Channel", gsub("I205", "I-205", gsub("_", " ", pop))), 
         m = mi*1609) %>%  # convert to m
  select(strata, pop, year, m) %>% arrange(strata, pop, year)

# Spawner abundance data
# Assumptions:
# (0) Fix coding error in data that assigns some Duncan Creek rows to Cascade stratum
# (1) NAs in hatchery dispositions are really zeros
# (2) NAs in Duncan Creek and Duncan Channel are really zeros
# (3) All other NAs are real missing observations
# (4) When calculating the observation error of log(S_obs), tau_S_obs, assume
#     Abund.Mean and Abund.SD are the mean and SD of a lognormal posterior distribution
#     of spawner abundance based on the sample
spawner_data <- read.csv(here("data","Data_Abundance_Spawners_Chum_2022-02-22.csv"), 
                         header = TRUE, stringsAsFactors = FALSE) %>% 
  rename(year = Return.Yr., strata = Strata, location = Location.Reach, 
         disposition = Disposition, method = Method, S_obs = Abund.Mean, SD = Abund.SD) %>% 
  mutate(disposition = gsub("I205", "I-205", gsub("_", " ", disposition)),
         location = gsub("I205", "I-205", gsub("_", " ", location)),
         strata = replace(strata, disposition == "Duncan Channel" & strata != "Gorge", "Gorge"),
         S_obs = replace(S_obs, is.na(S_obs) & grepl("Hatchery|Duncan", disposition), 0),
         tau_S_obs = sqrt(log((SD/S_obs)^2 + 1))) %>% 
  select(year:location, disposition, method, S_obs, SD, tau_S_obs) %>% 
  arrange(strata, location, year)

# broodstock take: 
# all spawners taken from a given location to a different disposition
# unless location is Duncan Creek and disposition is Duncan Channel,
# which is considered natural recruitment
# (summarized by *location*)
broodstock_data <- spawner_data %>% group_by(strata, location, year) %>% 
  summarize(B_take_obs = sum(S_obs[disposition != location & 
                                     !(location == "Duncan Creek" & disposition == "Duncan Channel")])) %>% 
  rename(pop = location) %>% as.data.frame()

# added (translocated) spawners:
# all spawners brought to a given disposition from a different location
# unless disposition is Duncan Channel and location is Duncan Creek,
# which is considered natural recruitment
# (summarized by *disposition*)
translocation_data <- spawner_data %>% group_by(strata, disposition, year) %>% 
  summarize(S_add_obs = sum(S_obs[disposition != location & 
                                    !(location == "Duncan Creek" & disposition == "Duncan Channel")])) %>% 
  rename(pop = disposition) %>% as.data.frame()

# total spawners:
# all spawners with a given disposition, regardless of original return location
# (summarized by disposition)
spawner_data_agg <- spawner_data %>% group_by(strata, disposition, year) %>% 
  summarize(S_obs = sum(S_obs), tau_S_obs = unique(tau_S_obs)) %>% 
  rename(pop = disposition) %>% 
  left_join(broodstock_data, by = c("strata","pop","year")) %>% 
  left_join(translocation_data, by = c("strata","pop","year")) %>% 
  mutate(B_take_obs = replace(B_take_obs, is.na(B_take_obs), 0),
         S_add_obs = replace(S_add_obs, is.na(S_add_obs), 0)) %>% as.data.frame()

# Spawner age-, sex-, and origin-frequency (aka BioData)
# Assumptions:
# (1) The first generation of salmonIPM models treated any known nonlocal-origin
#     spawner as "hatchery-origin" to avoid counting as local recruitment, so "H"
#     included true hatchery fish based on origin *plus* any whose known origin 
#     or return location do not match their disposition *unless* they are NOR or
#     Duncan Channel fish returning to Duncan Creek and disposed to Duncan Channel,
#     which are in fact local recruitment.
# (1b) The current version of IPM_LCRchum_pp allows S_add_obs to represent 
#     nonlocal natural spawners (of known or unknown origin) translocated from
#     return location to another disposition (i.e., Duncan Channel).
#     So "H" now includes true hatchery fish based on origin *plus* any others whose
#     known origin (i.e., Duncan Channel) does not match the disposition.
bio_data <- read.csv(here("data","Data_BioData_Spawners_Chum_2022-02-22.csv"), 
                     header = TRUE, stringsAsFactors = FALSE) %>% 
  rename(year = Return.Yr., strata = Strata, location = Location.Reach, 
         disposition = Disposition, origin = Origin, sex = Sex, age = Age, count = Count) %>% 
  mutate(disposition = gsub("I205", "I-205", gsub("_", " ", disposition)),
         location = gsub("I205", "I-205", gsub("_", " ", location)),
         origin = gsub("_", " ", origin),
         count = replace(count, is.na(count), 0), sex = substring(sex,1,1),
         HW = ifelse((grepl("Hatchery", origin) | !(origin %in% c(disposition, "Natural spawner"))), 
                     "H", "W")) %>% 
  select(year:location, disposition, origin, HW, sex:count) %>%
  arrange(strata, location, year, origin, age, sex)

# age of wild (i.e., potentially local) spawners only
bio_data_age <- bio_data %>% filter(HW == "W") %>%
  dcast(year + disposition ~ age, value.var = "count", fun.aggregate = sum) %>% 
  rename(pop = disposition)

# H/W (nonlocal/potentially local)  
bio_data_HW <- bio_data %>%
  dcast(year + disposition ~ HW, value.var = "count", fun.aggregate = sum) %>% 
  rename(pop = disposition)

# origin of spawners returning to a given location, regardless of final disposition
# recruits to Duncan Creek are considered as pop == "Duncan Channel"
# ignore Duncan Channel- and Big Creek Hatchery-origin spawners for now
bio_data_origin <- bio_data %>% 
  mutate(pop = replace(location, location == "Duncan Creek", "Duncan Channel")) %>% 
  dcast(year + pop ~ origin, value.var = "count", fun.aggregate = sum) %>% 
  select(year, pop, `Natural spawner`, `Duncan Hatchery`, `Lewis Hatchery`, `Grays Hatchery`)

# sex composition of spawners, regardless of origin
bio_data_sex <- bio_data %>% 
  dcast(year + disposition ~ sex, value.var = "count", fun.aggregate = sum) %>% 
  rename(pop = disposition) %>% select(year:pop, M, `F`)

# Proportion of "green" females in Duncan Channel
# Non-green (ripe or partial) females are assumed to have lower fecundity
# Proportion green females outside Duncan Channel assumed to = 1
# https://github.com/mdscheuerell/chumIPM/issues/5
green_female_data <- read.csv(here("data","Data_Duncan_Females_by_Condition_2022-02-22.csv"),
                              header = TRUE, stringsAsFactors = FALSE) %>% 
  rename(year = BY, disposition = Channel_Disposition, sex = Sex, condition = Condition,
         N = Qty, comment = Comment) %>% 
  dcast(year ~ condition, value.var = "N", fun.aggregate = sum) %>% 
  mutate(pop = "Duncan Channel", p_G_obs = Green / (Green + Ripe + Partial))

# Juvenile abundance data
# Assumptions:
# (1) Duncan_North + Duncan_South = Duncan_Channel, so the former two are redundant 
#     (not really an assumption, although the equality isn't perfect in all years)
# (2) When calculating the observation error of log(M_obs), tau_M_obs, assume
#     Abund_Median and Abund_SD are the median and SD of a lognormal posterior 
#     distribution of smolt abundance based on the sample
# (3) If Abund_SD == 0 (when Analysis=="Census": some years in Duncan_Channel and 
#     Hamilton_Channel) treat as NA
juv_data <- read.csv(here("data", "Data_Abundance_Juveniles_Chum_2022-02-22.csv"), 
                     header = TRUE, stringsAsFactors = FALSE) %>% 
  rename(brood_year = Brood.Year, year = Outmigration.Year, strata = Strata, 
         location = Location.Reach, origin = Origin, trap_type = TrapType, 
         analysis = Analysis, partial_spawners = Partial.Spawners, raw_catch = RawCatch,
         M_obs = Abund_Median, mean = Abund_Mean, SD = Abund_SD, 
         L95 = Abund_L95, U95 = Abund_U95, CV = Abund_CV, comments = Comments) %>% 
  mutate(location = gsub("I205", "I-205", gsub("_", " ", location)),
         origin = sapply(gsub("_", " ", origin), function(x)
           paste(rev(strsplit(x, " ")[[1]]), collapse = " ")), # names inconsistent w/ bio_data 
         tau_M_obs = replace(sqrt(log((SD/mean)^2 + 1)), SD==0, NA)) %>% 
  select(strata, location, year, brood_year, origin:CV, tau_M_obs, comments) %>% 
  arrange(strata, location, year)

# drop hatchery or redundant pops, hatchery-origin releases,
# and cases with leading or trailing NAs in M_obs
# use only pooled Duncan Channel data for now, not North / South
head_noNA <- function(x) { cumsum(!is.na(x)) > 0 }
juv_data_incl <- juv_data %>% 
  mutate(pop = ifelse(grepl("Hatchery", origin), origin, location), .after = strata) %>% 
  filter(!(pop %in% c("Duncan North","Duncan South","Big Creek Hatchery"))) %>% 
  select(pop, year, M_obs, tau_M_obs) %>% arrange(pop, year) %>% 
  group_by(pop, year) %>% summarize(M_obs = sum(M_obs), tau_M_obs = unique(tau_M_obs)) %>%
  ungroup() %>% group_by(pop) %>% filter(head_noNA(M_obs) & rev(head_noNA(rev(M_obs)))) %>%
  as.data.frame()
# juv_data_incl <- juv_data %>% group_by(location) %>%
#   filter(head_noNA(M_obs) & rev(head_noNA(rev(M_obs)))) %>%
#   filter(!(location %in% c("Duncan North","Duncan South")) & !grepl("Hatchery", origin)) %>%
#   rename(pop = location) %>% as.data.frame()

# Pairwise distance data
# convert ft to km
# add missing populations: 
#  Duncan Channel, Duncan Hatchery (same location as Duncan Creek)
#  Grays Hatchery (same location as Grays MS)
# extract lower triangle to a tidy data frame
pairwise_dist <- read.csv(here("data","Chum_Pairwise_Points_2021_Mat_Final.csv"), 
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

# Fish data formatted for salmonIPM
# Drop age-2 and age-6 samples (each is < 0.1% of aged spawners)
# Drop Duncan Creek
# Change S_obs and tau_S_obs to NA in Hamilton Channel 2011-2012 based on
# https://github.com/mdscheuerell/chumIPM/issues/6#issuecomment-807885445
# Pad data as necessary so Grays_MS, Grays_WF, and Grays_CJ have the same set of years
# (since their estimated smolts will be summed)  
# Pad data as necessary so hatchery populations are represented in all years
fish_data_all <- full_join(spawner_data_agg, bio_data_age, by = c("pop","year")) %>% 
  # full_join(bio_data_HW, by = c("pop","year")) %>% 
  full_join(bio_data_origin, by = c("pop","year")) %>% 
  full_join(bio_data_sex, by = c("pop","year")) %>% 
  full_join(juv_data_incl, by = c("pop","year")) %>%
  full_join({ complete(select(filter(., grepl("Grays", pop)), c(pop,year)), pop, year) },
            by = c("pop","year")) %>%
  full_join({ complete(select(filter(., grepl("Hatchery", pop)), c(pop,year)), pop, year) },
            by = c("pop","year")) %>% 
  left_join(habitat_data, by = c("pop","year")) %>% 
  left_join(green_female_data, by = c("pop","year")) %>% 
  rename_at(vars(contains("Age-")), list(~ paste0(sub("Age-","n_age",.), "_obs"))) %>% 
  select(-c(n_age2_obs, n_age6_obs)) %>% 
  # filter(!grepl("Hatchery|Duncan Creek", pop)) %>% 
  filter(pop != "Duncan Creek") %>% 
  mutate(pop = droplevels(factor(pop, levels = pop_names$pop)), # order E-W
         m = replace(m, grepl("Hatchery", pop), 1),
         S_obs = replace(S_obs, pop == "Hamilton Channel" & year %in% 2011:2012, NA),
         M_obs = replace(M_obs, grepl("Hatchery", pop) & is.na(M_obs), 0), 
         tau_S_obs = replace(tau_S_obs, pop == "Hamilton Channel" & year %in% 2011:2012, NA),
         B_take_obs = replace(B_take_obs, is.na(B_take_obs), 0),
         S_add_obs = replace(S_add_obs, is.na(S_add_obs), 0),
         p_G_obs = replace(p_G_obs, is.na(p_G_obs), 1), fit_p_HOS = 0, F_rate = 0) %>%
  # rename(A = m, n_H_obs = H, n_W_obs = W, n_M_obs = M, n_F_obs= `F`) %>% ## put w/ data sets
  rename(A = m, n_origin0_obs = `Natural spawner`, n_M_obs = M, n_F_obs= `F`) %>% 
  do({ 
    lev <- levels(.$pop)
    rename_with(., .cols = contains("Hatchery"), 
                .fn = ~ paste0("n_origin", match(.x, lev), "_obs"))
  }) %>% 
  mutate_at(vars(contains("n_")), ~ replace(., is.na(.), 0)) %>%
  mutate(n_W_obs = n_origin0_obs, 
         n_H_obs = rowSums(across(contains("origin"))) - n_origin0_obs,
         .before = n_origin0_obs) %>% 
  select(pop, year, A, S_obs, tau_S_obs, M_obs, tau_M_obs, n_age3_obs:n_F_obs, 
         p_G_obs, fit_p_HOS, B_take_obs, S_add_obs, F_rate) %>% 
  arrange(pop, year) 

# fill in fit_p_HOS
for(i in 1:nrow(fish_data_all)) {
  pop_i <- as.character(fish_data_all$pop[i])
  start_year <- ifelse(pop_i %in% hatcheries$pop,
                       min(hatcheries$start_brood_year[hatcheries$pop == pop_i]) + 1,
                       NA)
  fish_data_all$fit_p_HOS[i] <- ifelse((!is.na(start_year) & fish_data_all$year[i] >= start_year) |
                                         fish_data_all$n_H_obs[i] > 0, 1, 0)
}

# drop cases with initial NAs in both S_obs and M_obs
# (except hatchery populations, which must be present in every year)
fish_data <- fish_data_all %>% group_by(pop) %>% 
  filter(head_noNA(S_obs) | head_noNA(M_obs) | grepl("Hatchery", pop)) %>% 
  add_column(downstream_trap = NA, .after = "tau_M_obs")  %>% as.data.frame()

# assign Grays WF and Grays CJ smolts to the downstream trap in Grays MS
# where they will be counted (or double-counted, in the case of Grays CJ),
# assuming no mortality between the upstream tributary and the downstream trap
fish_data <- fish_data %>%
  mutate(downstream_trap = replace(downstream_trap, pop %in% c("Grays WF","Grays CJ"),
                                   which(pop == "Grays MS"))) %>% 
  mutate(pop_type = factor(ifelse(grepl("Hatchery", pop), "hatchery", "natural"), 
                           levels = c("natural","hatchery")),
         .after = pop)

# pad data with future years to generate forecasts
# use 1-year time horizon
fish_data_fore <- fish_data %>% group_by(pop) %>%
  slice(rep(n(), max(fish_data$year) + 1 - max(year))) %>%
  reframe(year = (unique(year) + 1):(max(fish_data$year) + 1),
            S_obs = NA, tau_S_obs = NA, M_obs = NA, tau_M_obs = NA, downstream_trap = NA, 
            p_G_obs = 1, S_add_obs = 0, fit_p_HOS = 0, B_take_obs = 0, F_rate = 0) %>%
  full_join(fish_data) %>% arrange(pop, year) %>%
  mutate_at(vars(starts_with("n_")), ~ replace_na(., 0)) %>%
  mutate(forecast = year > max(fish_data$year), downstream_trap = NA) %>% 
  fill(A, pop_type, .direction = "down") %>%
  select(pop, pop_type, year, forecast, A, S_obs, tau_S_obs, M_obs, tau_M_obs,
         downstream_trap, n_age3_obs:n_F_obs, p_G_obs, 
         fit_p_HOS, B_take_obs, S_add_obs, F_rate) %>% 
  as.data.frame()

# assign Grays_WF and Grays_CJ smolts to the downstream trap in Grays_MS
# (need to do this again b/c row indices have changed)
fish_data_fore <- fish_data_fore %>%
  mutate(downstream_trap = replace(downstream_trap, pop %in% c("Grays WF","Grays CJ"),
                                   which(pop == "Grays MS")))

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

# # Environmental covariates
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


#===========================================================================
# DATA EXPLORATION
#===========================================================================

EDA <- FALSE

if(EDA)
{
  
  # Time series of sex ratio by population
  windows()
  bio_data %>% filter(HW=="W" & !grepl("Hatchery|Duncan Creek", disposition)) %>% 
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
  fish_data %>% group_by(pop) %>% mutate(M0_obs = lead(M_obs), MperS = M0_obs/S_obs) %>% 
    ggplot(aes(x = S_obs, y = MperS)) + scale_x_log10() + scale_y_log10() +
    geom_point(size = 2) + labs(x = "spawners", y = "smolts / spawner") +
    facet_wrap(vars(pop), nrow = 3, ncol = 4) + theme_bw() +
    theme(panel.grid = element_blank())
  
  # Distribution of hatchery adults across recipient populations
  bdat <- bio_data %>% group_by(location, year, origin) %>% summarize(count = sum(count)) 
  sdat <- spawner_data %>% group_by(location, year) %>% summarize(S_obs = sum(S_obs, na.rm = TRUE))
  dat <- inner_join(bdat, sdat, by = c("location", "year")) %>% 
    mutate(prop = count/sum(count), S_origin = S_obs*prop) %>% 
    group_by(origin, location) %>% 
    summarize(count = sum(count), S_location = sum(S_origin, na.rm = TRUE)) %>% 
    mutate(prop = S_location/sum(S_location)) %>% 
    filter(origin != "Natural spawner") %>% ungroup() %>% 
    mutate(origin = droplevels(factor(origin, levels = c(pop_names$pop, "Big Creek Hatchery"))),
           location = droplevels(factor(location, levels = pop_names$pop))) %>% 
    complete(origin, location, fill = list(count = 0, S_location = 0, prop = 0))
  pdat <- pairwise_data %>% rename(origin = pop1, location = pop2) %>% 
    full_join(rename(., location = origin, origin = location)) %>% 
    right_join(dat, by = c("origin","location")) %>%
    filter(origin != "Big Creek Hatchery")
  
  windows(width = 5, height = 9)
  dat %>% ggplot(aes(x = location, y = prop)) + geom_col() +
    labs(x = "", y = "Proportion of adults returning to location") + 
    facet_wrap(vars(origin), ncol = 1) + theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          strip.text = element_text(margin = margin(b = 2, t = 2)),
          panel.grid = element_blank()) 
  
  windows(width = 5, height = 9)
  pdat %>% ggplot(aes(x = dist, y = prop)) + 
    geom_point(size = 2.5, pch = 16, alpha = 0.4) +
    labs(x = "Distance from origin (km)", y = "Proportion of adults returning to location") + 
    facet_wrap(vars(origin), ncol = 1) + theme_bw(base_size = 14) +
    theme(strip.text = element_text(margin = margin(b = 2, t = 2)),
          panel.grid = element_blank()) 
  
  # Fit cheesy beta regression to pooled dispersal by distance data
  library(rstanarm)
  betadat <- pdat %>% mutate(prop = pmax(prop, 0.001))
  betafit <- stan_betareg(prop ~ dist | dist, data = betadat, link = "logit", link.phi = "log")
  print(betafit, 2)
  ndat <- data.frame(dist = 0:round(max(pdat$dist)))
  epd <- posterior_epred(betafit, newdata = ndat) %>% 
    colQuantiles(probs = c(0.025, 0.5, 0.975)) %>% as.data.frame
  ppd <- posterior_predict(betafit, newdata = betadat) %>% 
    colQuantiles(probs = c(0.025, 0.975)) %>% as.data.frame %>% 
    rename(L = `2.5%`, U = `97.5%`) %>% cbind(dist = betadat$dist)
  betapred <- cbind(ndat, epd)
  
  windows()
  pdat %>% ggplot(aes(x = dist, y = prop)) + 
    geom_line(aes(y = `50%`), data = betapred, col = "darkgray", size = 2) +
    geom_ribbon(aes(y = L, ymin = L, ymax = U), data = ppd, col = "gray", alpha = 0.3) +
    geom_ribbon(aes(y = `50%`, ymin = `2.5%`, ymax = `97.5%`), data = betapred,
                col = "darkgray", alpha = 0.4) +
    geom_point(size = 2.5, pch = 16, alpha = 0.4) +
    labs(x = "Distance from origin (km)", y = "Proportion of adults returning to location") + 
    theme_bw(base_size = 16) +
    theme(strip.text = element_text(margin = margin(b = 2, t = 2)),
          panel.grid = element_blank()) 
}
