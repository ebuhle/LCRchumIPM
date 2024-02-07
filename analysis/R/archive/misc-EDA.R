#===========================================================================
# DATA EXPLORATION
#===========================================================================

options(device = windows)
library(ggplot2)
theme_set(theme_bw(base_size = 16))
library(here)
source(here("analysis","R","01_LCRchumIPM_data.R"))

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
  facet_wrap(vars(disposition), ncol = 4) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Histogram of spawner age by sex and H/W origin
windows()
bio_data %>% mutate(age = substring(age,5,5)) %>% group_by(HW, sex, age) %>% 
  summarize(n = sum(count)) %>% mutate(prop = n/sum(n)) %>% 
  ggplot(aes(x = age, y = prop)) + geom_bar(stat = "identity") + 
  facet_wrap(vars(HW, sex), nrow = 2, ncol = 2)

# Histogram of age by broodstock status for each pop
windows(width = 10, height = 7)
bio_data %>% 
  mutate(status = ifelse(location == disposition, "natural", "broodstock"), 
         age = substring(age,5,5)) %>% 
  filter(age %in% 3:5 & location != "Duncan Creek") %>% 
  group_by(location, status, age) %>% summarize(n = sum(count)) %>% 
  mutate(N = sum(n)) %>% data.frame(., with(., binconf(x = n, n = N))) %>%
  rename(proportion = PointEst) %>% 
  group_by(location) %>% mutate(any_broodstock = any(status == "broodstock")) %>% 
  filter(any_broodstock) %>% 
  ggplot(aes(x = age, y = proportion, ymin = Lower, ymax = Upper, fill = status)) +
  geom_col(position = "dodge") +
  geom_errorbar(position = position_dodge(width = 0.9), width = 0) +
  facet_wrap(vars(location), nrow = 2) + theme(legend.position = c(0.9, 0.3))

# Histogram of sex by broodstock status for each pop
windows(width = 10, height = 7)
bio_data %>% 
  mutate(status = ifelse(location == disposition, "natural", "broodstock"), 
         age = substring(age,5,5)) %>% 
  filter(age %in% 3:5 & location != "Duncan Creek") %>% 
  group_by(location, status, sex) %>% summarize(n = sum(count)) %>% 
  dcast(location + status ~ sex, value.var = "n", fun.aggregate = sum) %>% 
  mutate(total = `F` + M) %>% 
  data.frame(., with(., binconf(x = `F`, n = total))) %>%
  rename(`proportion female` = PointEst) %>% 
  group_by(location) %>% mutate(any_broodstock = any(status == "broodstock")) %>% 
  filter(any_broodstock) %>% 
  ggplot(aes(x = status, y = `proportion female`, ymin = Lower, ymax = Upper, fill = status)) +
  geom_col(position = "dodge") +
  geom_errorbar(position = position_dodge(width = 0.9), width = 0) +
  facet_wrap(vars(location), nrow = 2) + theme(legend.position = c(0.9, 0.3))

# Histogram of origin (known vs. unknown) by broodstock status for each pop
windows(width = 10, height = 7)
bio_data %>% 
  mutate(status = ifelse(location == disposition, "natural", "broodstock"), 
         origin = ifelse(origin == "Natural spawner", "unknown", "known"),
         age = substring(age,5,5)) %>% 
  filter(age %in% 3:5 & location != "Duncan Creek") %>% 
  group_by(location, status, origin) %>% summarize(n = sum(count)) %>% 
  dcast(location + status ~ origin, value.var = "n", fun.aggregate = sum) %>% 
  mutate(total = known + unknown) %>% 
  data.frame(., with(., binconf(x = known, n = total))) %>%
  rename(`proportion known origin` = PointEst) %>% 
  group_by(location) %>% mutate(any_broodstock = any(status == "broodstock")) %>% 
  filter(any_broodstock) %>% 
  ggplot(aes(x = status, y = `proportion known origin`, ymin = Lower, ymax = Upper, fill = status)) +
  geom_col(position = "dodge") +
  geom_errorbar(position = position_dodge(width = 0.9), width = 0) +
  facet_wrap(vars(location), nrow = 2) + theme(legend.position = c(0.9, 0.3))

# Boxplots of fecundity by age
windows()
fecundity_data %>% mutate(age_E = factor(age_E)) %>% 
  ggplot(aes(x = age_E, y = E_obs)) + geom_boxplot()

# Boxplots of fecundity by age, grouped by strata
windows()
fecundity_data %>% mutate(age_E = factor(age_E)) %>% 
  ggplot(aes(x = age_E, y = E_obs)) + geom_boxplot() + 
  facet_wrap(vars(strata), nrow = 1, ncol = 3)

# Histograms of fecundity, grouped by age and strata
windows()
fecundity_data %>% mutate(age_E = factor(age_E)) %>%  ggplot(aes(x = E_obs)) +
  geom_histogram(aes(y = stat(density)), bins = 15, color = "white", fill = "darkgray") + 
  facet_grid(rows = vars(strata), cols = vars(age_E))

# Normal QQ plots of fecundity, grouped by age and strata
windows()
fecundity_data %>% mutate(age_E = factor(age_E)) %>%  ggplot(aes(sample = E_obs)) +
  geom_qq(distribution = qnorm) + geom_qq_line(distribution = qnorm) +
  facet_grid(rows = vars(strata), cols = vars(age_E))

# Smolts/spawner vs. spawners, grouped by population
windows(width = 10, height = 5)
fish_data %>% group_by(pop) %>% 
  mutate(M0_obs = lead(M_obs), MperS = M0_obs/S_obs) %>% 
  filter(any(!is.na(MperS)) & is.finite(MperS)) %>% 
  ggplot(aes(x = S_obs, y = MperS))  +
  geom_point(size = 2) + labs(x = "spawners", y = "smolts / spawner") +
  scale_y_log10() + facet_wrap(vars(pop), nrow = 2, scales = "free") +
  theme(panel.grid = element_blank())

# Distribution of hatchery adults across recipient populations
bdat <- bio_data %>% group_by(location, year, origin) %>% summarize(count = sum(count)) 
sdat <- spawner_data %>% group_by(location, year) %>% summarize(S_obs = sum(S_obs, na.rm = TRUE))
dat <- inner_join(bdat, sdat, by = c("location", "year")) %>% 
  mutate(prop = count/sum(count), S_O = S_obs*prop) %>% 
  group_by(origin, location) %>% 
  summarize(count = sum(count), S_location = sum(S_O, na.rm = TRUE)) %>% 
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

