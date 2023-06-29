windows(width = 5, height = 7)
fish_data %>% filter(pop_type == 'hatchery') %>% 
  ggplot(aes(x = year, y = M_obs)) +
  geom_line() + geom_point() +
  facet_wrap(vars(pop), ncol = 1) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(margin = margin(b = 3, t = 3)))

windows(width = 11, height = 7)
fish_data %>% filter(pop_type == 'natural') %>% 
  ggplot(aes(x = year, y = B_take_obs)) +
  geom_line() + geom_point() +
  facet_wrap(vars(pop), ncol = 4) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(margin = margin(b = 3, t = 3)))

windows(width = 11, height = 7)
fish_data %>% filter(pop_type == 'natural') %>% 
  ggplot(aes(x = S_obs + B_take_obs, y = B_take_obs)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(size = 2.5, color = 'gray50', alpha = 0.8) +
  scale_x_log10() + scale_y_log10() +
  facet_wrap(vars(pop), ncol = 4) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(margin = margin(b = 3, t = 3)))



