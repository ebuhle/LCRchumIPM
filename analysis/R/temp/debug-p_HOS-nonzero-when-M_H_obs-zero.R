## Conclusion from this deep dive:
##   The reason p_origin / p_HOS are nonzero even when all hatchery M_obs == 0
##   is that stan_data() converts 0 abundances to 1. 
##   (It does not flag them as missing, though it also changes NA to 1 and then ignores them.)
##   This produces tiny trace amounts of S_H from the M_obs == 1. Hatchery M_obs is 
##   considered known w/o error in the current version of the model, but that will
##   change once we build the broodstock-to-smolt stage (true hatchery smolts will
##   be an unknown state M), so it's not worth worrying about now.

dat <- fish_data %>% mutate(M_obs = replace(M_obs, pop_type == 'hatchery', NA))
fd <- stan_data('IPM_LCRchum_pp', ages = list(M = 1),
                fish_data = dat, fecundity_data = fecundity_data)

test <- salmonIPM(stan_model = 'IPM_LCRchum_pp', SR_fun = 'Ricker', 
                  par_models = list(s_MS ~ pop_type), 
                  center = FALSE, scale = FALSE, ages = list(M = 1), 
                  fish_data = dat, fecundity_data = fecundity_data,
                  pars = c(stan_pars('IPM_LCRchum_pp'), 'S_origin'),
                  chains = 1, iter = 10, warmup = 0)

dr <- as_draws_rvars(as.matrix(test, c('S_origin','p_HOS'))) #%>% 
dr

S_origin = mean(dr$S_origin) %>% 
  matrix(ncol = 15, byrow = TRUE) %>% as.data.frame() %>% setNames(1:ncol(.)) %>% 
  cbind(data.frame(pop = rep(levels(dat$pop), length(unique(dat$year))),
                   year = rep(sort(unique(dat$year)), each = 15)), .) %>% 
  group_by(pop) %>% mutate(pop_year = year - min(year) + 1, .after = year)

S_origin %>% select(pop, year, pop_year, as.character(fd$which_H_pop)) %>% 
  filter(year > min(year) + 4 & !grepl('Hatchery', pop)) %>% print(n = 200)

