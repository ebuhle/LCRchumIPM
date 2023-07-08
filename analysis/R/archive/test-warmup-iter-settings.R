#### TEMP: Compare wall time, n_eff, and n_eff/time for shorter warmup / iter

load(here("analysis","results","archive","test-warmup-iter-settings.RData"))
pars <- stan_pars('IPM_LCRchum_pp', 'all')
names(pars) <- gsub('\\d+', '', names(pars))

## rstan default settings currently used in @knitr fit_Ricker

# fit_Ricker1 <- salmonIPM(stan_model = "IPM_LCRchum_pp", SR_fun = "Ricker", ages = list(M = 1), 
#                         par_models = list(s_MS ~ pop_type), center = FALSE, scale = FALSE, 
#                         fish_data = fish_data, fecundity_data = fecundity_data,
#                         chains = 4, iter = 2000, warmup = 1000,
#                         control = list(max_treedepth = 15))

wall_time1 <- rowSums(get_elapsed_time(fit_Ricker1))/3600
draws1 <- as_draws_rvars(as.matrix(fit_Ricker1, stan_pars('IPM_LCRchum_pp', 'all')))
sdraws1 <- summarize_draws(draws1) %>% mutate(mod = 'fit_Ricker1', .before = variable) %>% 
  mutate(par_name = sapply(strsplit(variable, "\\["), function(x) x[1]), 
         level = names(pars)[match(par_name, pars)], .after = variable)

## short warmup, same post-warmup sample size

# fit_Ricker2 <- salmonIPM(stan_model = "IPM_LCRchum_pp", SR_fun = "Ricker", ages = list(M = 1), 
#                          par_models = list(s_MS ~ pop_type), center = FALSE, scale = FALSE, 
#                          fish_data = fish_data, fecundity_data = fecundity_data,
#                          chains = 4, iter = 1500, warmup = 500,
#                          control = list(max_treedepth = 15))

wall_time2 <- rowSums(get_elapsed_time(fit_Ricker2))/3600
draws2 <- as_draws_rvars(as.matrix(fit_Ricker2, stan_pars('IPM_LCRchum_pp', 'all')))
sdraws2 <- summarize_draws(draws2) %>% mutate(mod = 'fit_Ricker2', .before = variable) %>% 
  mutate(par_name = sapply(strsplit(variable, "\\["), function(x) x[1]), 
         level = names(pars)[match(par_name, pars)], .after = variable)

## short warmup, same iter

# fit_Ricker3 <- salmonIPM(stan_model = "IPM_LCRchum_pp", SR_fun = "Ricker", ages = list(M = 1), 
#                          par_models = list(s_MS ~ pop_type), center = FALSE, scale = FALSE, 
#                          fish_data = fish_data, fecundity_data = fecundity_data,
#                          chains = 4, iter = 2000, warmup = 500,
#                          control = list(max_treedepth = 15))

wall_time3 <- rowSums(get_elapsed_time(fit_Ricker3))/3600
draws3 <- as_draws_rvars(as.matrix(fit_Ricker3, stan_pars('IPM_LCRchum_pp', 'all')))
sdraws3 <- summarize_draws(draws3) %>% mutate(mod = 'fit_Ricker3', .before = variable) %>% 
  mutate(par_name = sapply(strsplit(variable, "\\["), function(x) x[1]), 
         level = names(pars)[match(par_name, pars)], .after = variable)

## short warmup, same post-warmup sample size, higher adapt_delta

# fit_Ricker4 <- salmonIPM(stan_model = "IPM_LCRchum_pp", SR_fun = "Ricker", ages = list(M = 1),
#                          par_models = list(s_MS ~ pop_type), center = FALSE, scale = FALSE,
#                          fish_data = fish_data, fecundity_data = fecundity_data,
#                          chains = 4, iter = 1500, warmup = 500,
#                          control = list(adapt_delta = 0.99, max_treedepth = 15))

wall_time4 <- rowSums(get_elapsed_time(fit_Ricker4))/3600
draws4 <- as_draws_rvars(as.matrix(fit_Ricker4, stan_pars('IPM_LCRchum_pp', 'all')))
sdraws4 <- summarize_draws(draws4) %>% mutate(mod = 'fit_Ricker4', .before = variable) %>% 
  mutate(par_name = sapply(strsplit(variable, "\\["), function(x) x[1]), 
         level = names(pars)[match(par_name, pars)], .after = variable)

# save stanfits
# save(list = c('fit_Ricker1','fit_Ricker2','fit_Ricker3','fit_Ricker4'), 
#      file = here('analysis','results','archive','test-warmup-iter-settings.RData'))

## Compare fits

# wall time summaries
wall_time <- data.frame(mod = c('fit_Ricker1','fit_Ricker2','fit_Ricker3','fit_Ricker4'),
                   warmup = c(1000, 500, 500, 500), iter = c(2000, 1500, 2000, 1500)) %>% 
  mutate(saved = iter - warmup, chains = 4, adapt_delta = c(rep(0.95, 3), 0.99),
         hr_total = c(sum(wall_time1), sum(wall_time2), sum(wall_time3), sum(wall_time4)), 
         hr_mean = hr_total/chains, 
         hr_sd = c(sd(wall_time1), sd(wall_time2), sd(wall_time3), sd(wall_time4)), 
         hr_cv = c(hr_sd/hr_mean),
         num_divergent = c(get_num_divergent(fit_Ricker1),
                           get_num_divergent(fit_Ricker2),
                           get_num_divergent(fit_Ricker3),
                           get_num_divergent(fit_Ricker4))) %>% 
  arrange(warmup, saved) %>% 
  mutate(settings = factor(mod, levels = mod, 
                           labels = paste0("warmup = ", warmup, 
                                           "\nsaved = ", saved,
                                           "\nadapt_delta = ", adapt_delta)),
         .after = mod)



# ess/hr
ess_hr <- rbind(sdraws1, sdraws2, sdraws3, sdraws4) %>% left_join(wall_time, by = 'mod') %>% 
  select(settings, warmup = warmup, saved = saved, hr_total, variable, level, ess_bulk, ess_tail) %>% 
  mutate(ess_hr_bulk = ess_bulk/hr_total, ess_hr_tail = ess_tail/hr_total)

ess_hr_long <- ess_hr %>% select(-c(ess_bulk, ess_tail)) %>% 
  pivot_longer(cols = starts_with('ess_hr'), names_to = "type", names_pattern = "ess_hr_(.*)",
               values_to = "ess_hr") 

# boxplots of ess / hr
# (note "rows containing non-finite values" are constants such as
#  diagonals of R_pop_p and R_p, hatchery S-R params and M states, etc.)
windows()
ess_hr_long %>% 
  ggplot(aes(x = type, y = ess_hr, color = settings)) + 
  geom_boxplot(lwd = 0.8) + 
  labs(x = NULL, y = "ess / hr", color = NULL) +
  facet_wrap(~ level, nrow = 1) +
  theme(panel.grid = element_blank())

# which pars have low ess?
low_ess <- ess_hr %>% group_by(variable) %>% 
  summarize(ess_bulk = mean(ess_bulk), ess_draw_bulk = mean(ess_bulk/(4*saved)), 
            ess_tail = mean(ess_tail), ess_draw_tail = mean(ess_tail/(4*saved))) %>% 
  arrange(ess_draw_bulk)

# histograms of ess / draw: is there a distinct "bad" group?
windows()
low_ess %>% select(variable, ess_draw_bulk, ess_draw_tail) %>% 
  pivot_longer(cols = starts_with('ess_draw'), names_to = 'type', 
               names_pattern = 'ess_draw_(.*)', values_to = 'ess_draw') %>% 
  ggplot(aes(x = ess_draw)) +
  geom_histogram(col = 'white') +
  geom_vline(xintercept = 1, col = 'blue') +
  geom_vline(xintercept = 0.1, col = 'red') +
  geom_vline(xintercept = 0.5, col = 'red') +
  xlab('ess per draw') +
  facet_wrap(~ type, ncol = 1) +
  theme(panel.grid = element_blank())

# use 10% ess per draw threshold
low_ess %>% filter(ess_draw_bulk < 0.1) %>% arrange(variable) %>% print(n = 50)
low_ess %>% filter(ess_draw_tail < 0.1) %>% arrange(variable) %>% print(n = 50)

# examples of "low" (still in the hundreds) bulk or tail ess:
#  - eta_year_MS[2:4], i.e. 1999:2002 before any (?) smolt data exist
#  - s_MS for Duncan, Lewis and Grays Hatcheries in 2000; first two have M_obs == 0
#  - s_MS for Grays Hatchery in years 2000:2003, i.e. before wild Grays time series start
#  - q_O[, c(1,4)] (and subsequently p_HOS) for four Lower Gorge pops in 2003
#    (q_O[,4] is Grays Hatchery, and these cases have only q_O[,1] natural recoveries)

# compare selected hyperparameter estimates to see if divergences are causing bias
#  iter = 2000, warmup = 500, adapt_delta = 0.95 
#  iter = 2000, warmup = 1000, adapt_delta = 0.99

adapt_delta_compare <- rbind(sdraws2, sdraws4) %>% filter(level == 'hyper') %>% 
  left_join(wall_time, by = "mod") %>% 
  select(mod, warmup, iter, adapt_delta, variable, q5, median, q95, ess_bulk, ess_tail)

adapt_delta_compare %>% mutate(adapt_delta = factor(adapt_delta)) %>% 
  filter(!grepl('mu_E|sigma_E|mu_Mmax|mu_p|R_p|P_D|mu_tau|sigma_tau', variable)) %>% 
  ggplot(aes(x = median, y = variable, group = adapt_delta, color = adapt_delta)) +
  # geom_point(size = 3) +
  # geom_segment(aes(x = q5, xend = q95, yend = variable), linewidth = 1) + 
  geom_pointrange(aes(xmin = q5, xmax = q95), size = 0.8, linewidth = 0.8, 
                  position = position_dodge(width = 1)) +
  labs(x = NULL, y = NULL)



