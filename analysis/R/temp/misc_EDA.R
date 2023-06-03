#===========================================================================
# DATA EXPLORATION
#===========================================================================

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

