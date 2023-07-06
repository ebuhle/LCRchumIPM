## Explore non-IID patterns in observation-level SAR residuals zeta_MS
## Is there enough among-pop variance (vs. among-year) to warrant pop-level random effects?

# zeta_MS timeseries by pop
windows(width = 10, height = 7)
fish_data %>% mutate(zeta_MS = stan_mean(LCRchum_Ricker,'zeta_MS')) %>% 
  ggplot(aes(x = year, y = zeta_MS)) + 
  geom_hline(yintercept = 0) + 
  geom_line() +
  geom_point(pch = 1, size = 2) + 
  facet_wrap(vars(pop)) +
  theme(panel.grid.minor = element_blank(), 
        strip.background = element_rect(fill = NA),
        strip.text = element_text(margin = margin(b = 3, t = 3)))

# zeta_MS boxplots by pop
windows()
fish_data %>% mutate(zeta_MS = stan_mean(LCRchum_Ricker,'zeta_MS')) %>% 
  ggplot(aes(y = pop, x = zeta_MS)) + 
  geom_vline(xintercept = 0, lty = 2) + 
  geom_boxplot(fill = NA) +
  scale_y_discrete(limits = rev) + ylab("")

# zeta_MS boxplots by H/W
windows()
fish_data %>% mutate(zeta_MS = stan_mean(LCRchum_Ricker,'zeta_MS')) %>% 
  ggplot(aes(y = pop_type, x = zeta_MS)) + 
  geom_vline(xintercept = 0, lty = 2) + 
  geom_boxplot(fill = NA) +
  scale_y_discrete(limits = rev) + ylab("")

# Hierarchical model of posterior mean obs-level residuals
# Compare pop-level and residual SD
require(rstanarm)
m1 <- stan_lmer(zeta_MS ~ (1|pop),
                data = data.frame(fish_data, zeta_MS = stan_mean(LCRchum_Ricker,'zeta_MS')))
summary(m1)


## Does distance from Columbia mouth explain pop-level variation in SAR residuals?

dist_to_mouth <- pairwise_data %>% filter(pop1 == "Columbia Mouth") %>% 
  rename(pop = pop2) %>% select(-pop1)
dat <- fish_data %>% left_join(dist_to_mouth) %>% 
  mutate(zeta_MS = stan_mean(LCRchum_Ricker,'zeta_MS')) %>% 
  group_by(pop) %>% mutate(has_M_obs = any(!is.na(M_obs)))

# Scatterplot of obs-level residuals vs. distance to mouth
windows()
dat %>%
  ggplot(aes(x = dist, y = zeta_MS, group = pop, color = pop_type, shape = has_M_obs)) +
  geom_point(size = 2.5) +
  scale_color_discrete(type = c(natural = alpha("slategray4", 0.6), 
                                hatchery = alpha("salmon", 0.6))) +
  scale_shape_manual(values = c(1,16)) +
  labs(x = "Distance from Columbia mouth (km)", color = "Origin", shape = "Smolts") +
  theme(legend.position = "top")

# Scatterplot of pop-level mean residuals vs. distance to mouth
windows()
dat %>% group_by(pop, pop_type) %>% 
  summarize(dist = mean(dist), zeta_MS = mean(zeta_MS), has_M_obs = unique(has_M_obs)) %>% 
  ggplot(aes(x = dist, y = zeta_MS, color = pop_type, shape = has_M_obs)) +
  geom_point(size = 3) +
  geom_smooth(aes(x = dist, y = zeta_MS), method = "lm", inherit.aes = FALSE) +
  scale_color_discrete(type = c(natural = alpha("slategray4", 0.8), 
                                hatchery = alpha("salmon", 0.8))) +
  scale_shape_manual(values = c(1,16)) +
  labs(x = "Distance from Columbia mouth (km)", color = "Origin", shape = "Smolts") +
  theme(legend.position = "top")

# Regression model of obs-level residuals vs. distance to mouth
require(rstanarm)
m2 <- stan_glm(zeta_MS ~ dist_std, data = mutate(dat, dist_std = scale(dist)))
summary(m2)

windows()
dat %>% mutate(resid_MS = resid(m2)) %>% 
  ggplot(aes(y = pop, x = resid_MS)) + 
  geom_vline(xintercept = 0, lty = 2) + 
  geom_boxplot(fill = NA) +
  scale_y_discrete(limits = rev) + ylab("")

# Scatterplot of pop-level mean residuals vs. pop-level psi (natural pops only)
windows()
dat %>% group_by(pop, pop_type) %>% 
  summarize(zeta_MS = mean(zeta_MS), has_M_obs = unique(has_M_obs)) %>% 
  ungroup() %>% mutate(psi = stan_mean(LCRchum_Ricker, 'psi')) %>% 
  filter(pop_type == 'natural') %>% 
  ggplot(aes(x = psi, y = zeta_MS, shape = has_M_obs)) +
  geom_point(color = 'slategray4', size = 3) +
  scale_shape_manual(values = c(1,16)) +
  labs(shape = "Smolts") +
  theme(legend.position = "top")


