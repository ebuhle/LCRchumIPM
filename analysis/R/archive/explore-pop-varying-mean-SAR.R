logit <- rfun(qlogis)
ilogit <- rfun(plogis)

dd <- stan_data('IPM_LCRchum_pp', ages = list(M = 1), fish_data = fish_data, 
                fecundity_data = fecundity_data)
N_pop <- max(dd$pop)
which_H_pop <- dd$which_H_pop

draws <- as.matrix(fit_Ricker, c("psi","beta_psi","eta_year_M","sigma_M","mu_MS",
                                 "beta_MS","mu_pop_MS","eta_year_MS","s_MS")) %>% 
  as_draws_rvars() %>% 
  mutate_variables(logit_psi = logit(psi),
                   Xbeta_psi = replace(rep(rvar(0), N_pop), which_H_pop, beta_psi),
                   psi = ilogit(logit_psi + Xbeta_psi),
                   zeta_M = as_rvar(stan_mean(fit_Ricker,"zeta_M")), # not monitored: use mean
                   epsilon_M = sigma_M*zeta_M,
                   error_M = eta_year_M[dd$year] + epsilon_M,
                   exp_eta_year_M = exp(eta_year_M), 
                   exp_error_M = exp(error_M),
                   Xbeta_pop_MS = replace(rep(rvar(0), N_pop), which_H_pop, beta_MS),
                   logit_pop_MS_hat = logit(mu_MS) + Xbeta_pop_MS,
                   eta_pop_MS = logit(mu_pop_MS) - logit_pop_MS_hat,
                   SAR = 100*s_MS,
                   .value = c(psi, eta_pop_MS, logit(mu_pop_MS)))

dat <- fish_data %>% group_by(pop) %>% 
  summarize(pop = unique(pop),
            has_M_obs = any(!is.na(M_obs) | !is.na(downstream_trap))) %>%
  rbind(.,.,.) %>% 
  mutate(pars = rep(c("psi", "eta_pop_MS", "logit(mu_pop_MS)"), each = N_pop),
         pars = factor(pars, levels = unique(pars)),
         .value = draws$.value) %>% 
  as.data.frame()

gg <- dat %>% 
  ggplot(aes(xdist = .value, y = pop, fill = has_M_obs)) +
  stat_eye(.width = c(0.5, 0.9), normalize = "groups", 
           color = "slategray4", slab_color = "slategray4", slab_linewidth = 0.5) +
  scale_y_discrete(limits = rev) + labs(x = NULL, y = NULL) + 
  scale_fill_manual(values = alpha("slategray4", c(0, 0.15)), guide = "none") +
  facet_wrap(~ pars, scales = "free_x", strip.position = "bottom") + 
  theme(panel.grid.minor = element_blank(), strip.background = element_blank(),
        strip.text = element_text(size = 14, margin = margin(b = 3, t = 3)),
        strip.placement = "outside")

windows(width = 10, height = 7)
gg
ggsave(filename = here("analysis","results","archive","pop-varying-SAR.png"),
       width=10, height=7, units="in", dpi=300)

######

psi <- as_draws_df(draws$psi) %>% 
  pivot_longer(col = starts_with("x"), names_to = "pop", 
               names_pattern = "x\\[(.+)\\]", values_to = "psi")
mu_pop_MS <- as_draws_df(draws$mu_pop_MS) %>% 
  pivot_longer(col = starts_with("x"), names_to = "pop", 
               names_pattern = "x\\[(.+)\\]", values_to = "mu_pop_MS")
df <- left_join(psi, mu_pop_MS, by = c(".chain",".iteration",".draw","pop")) %>% 
  mutate(pop = sort(unique(fish_data$pop))[as.numeric(pop)])

gg <- df %>% 
  ggplot(aes(x = qlogis(psi), y = qlogis(mu_pop_MS))) +
  geom_point(pch = 1, alpha = 0.5) +
  facet_wrap(~ pop, ncol = 5) +
  labs(x = "logit(psi)", y = "logit(mu_pop_MS)") +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = NA),
        strip.text = element_text(size = 14, margin = margin(b = 3, t = 3)))
  
windows(width = 12, height = 7)
gg
ggsave(filename=here("analysis","results","archive","posterior-corr-pop-psi-SAR.png"),
       width=12, height=7, units="in", dpi=300)

