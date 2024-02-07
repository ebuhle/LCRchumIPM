options(device = windows)
library(ggplot2)
theme_set(theme_bw(base_size = 16))
library(here)
source(here("analysis","R","01_LCRchumIPM_data.R"))

agedat <- bio_data %>% 
  mutate(pop = factor(replace(location, location == "Duncan Creek", "Duncan Channel"),
                      levels = levels(fish_data$pop))) %>%
  group_by(pop, year, origin, age) %>% 
  summarize(n_age_obs = sum(count)) %>% 
  mutate(q_obs = n_age_obs / sum(n_age_obs))

agedat_agg <- agedat %>% group_by(pop, year, age) %>% 
  summarize(n_age_obs = sum(n_age_obs)) %>% 
  mutate(q_obs = n_age_obs / sum(n_age_obs)) %>% 
  filter(!grepl("Hatchery", pop) & !age %in% c("Age-2","Age-6"))
  
agedat %>% filter(!grepl("Hatchery", pop) & !age %in% c("Age-2","Age-6")) %>% 
  ggplot(aes(x = year, y = q_obs, group = origin, color = origin)) +
  geom_line() +
  geom_point() +
  geom_line(aes(x = year, y = q_obs), data = agedat_agg, inherit.aes = FALSE, lty = 2) +
  geom_point(aes(x = year, y = q_obs), data = agedat_agg, inherit.aes = FALSE,
            pch = 1) +
  facet_grid(rows = vars(age), cols = vars(pop)) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), 
        strip.background = element_rect(fill = NA),
        strip.text = element_text(margin = margin(b = 3, t = 3)),
        legend.position = "top")

sexdat <- bio_data %>% 
  mutate(pop = factor(replace(location, location == "Duncan Creek", "Duncan Channel"),
                      levels = levels(fish_data$pop))) %>%
  group_by(pop, year, origin, sex) %>% 
  summarize(n_MF_obs = sum(count)) %>% 
  mutate(q_F_obs = n_MF_obs / sum(n_MF_obs)) %>% 
  filter(sex == "F")

sexdat %>% filter(!grepl("Hatchery", pop)) %>% 
  ggplot(aes(x = year, y = q_F_obs, group = origin, color = origin)) +
  geom_line() +
  geom_point() +
  facet_wrap(vars(pop)) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), 
        strip.background = element_rect(fill = NA),
        strip.text = element_text(margin = margin(b = 3, t = 3)),
        legend.position = "top")



