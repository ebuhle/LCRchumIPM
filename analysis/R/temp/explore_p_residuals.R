# Explore non-IID patterns in observation-level age-structure residuals zeta_p
# Is there enough among-pop synchrony to warrant an ESU-level trend (vs. static mu_p)?

zeta_p <- matrix(stan_mean(LCRchum_Ricker, 'zeta_p'), ncol = 2, byrow = TRUE,
                 dimnames = list(NULL, c('zeta_p1', 'zeta_p2')))
dat <- fish_data %>% select(pop, pop_type, year) %>% cbind(zeta_p) %>% 
  complete(pop, year) %>% filter(pop_type == 'natural') %>% select(-pop_type)

dat %>% 
  pivot_longer(cols = starts_with('zeta'), names_to = 'age',
               names_pattern = 'zeta_p(.)', values_to = 'zeta_p') %>% 
  ggplot(aes(x = year, y = zeta_p, group = pop)) +
  geom_line(col = 'slategray4', alpha = 0.7) +
  facet_grid(vars(age))
  
  
  