# smolt obs error

m <- as.data.frame(summary(LCRchum_BH, prob = c(0.05,0.5,0.95), 
                           pars = "tau_M")$summary) %>% 
  select(-c(se_mean,sd,n_eff,Rhat)) %>% round(2) %>% 
  mutate(pop = fish_data_SMS$pop, year = fish_data_SMS$year) %>% 
  filter(pop %in% c("Duncan_Creek","Hamilton_Channel")) %>% 
  left_join(filter(juv_data, pop %in% c("Duncan_Creek","Hamilton_Channel")) %>% 
              select(pop,year,trap_type,analysis,CV) %>% mutate(CV = round(CV,2))) %>% 
  select(pop, year, trap_type, analysis, CV, everything())
  
m

windows()
m %>% ggplot(aes(x = factor(CV==0), y = `50%`)) + geom_boxplot() + theme_bw()

windows()
juv_data %>% ggplot(aes(x = trap_type, y = CV)) + geom_boxplot() + theme_bw()

windows()
juv_data %>% ggplot(aes(x = analysis, y = CV)) + geom_boxplot() + theme_bw()


# implied lognormal posteriors from smolt observation model 
# are they really lognormal??

# (1) compare implied 95% quantiles to empirically reported quantiles
jd <- juv_data %>% filter(origin=="Wild") %>% select(-comments) %>% 
  filter(tau_M_obs %in% c(range(tau_M_obs, na.rm = TRUE),
                          quantile(tau_M_obs, c(0.25, 0.5, 0.75, 0.95), 
                                   type = 1, na.rm = TRUE)))

windows()
par(mfrow=c(3,2))
for(i in 1:6) {
  curve(dlnorm(x, log(jd$M_obs[i]), jd$tau_M_obs[i]), 
        from = jd$L95[i]*0.9, to = jd$U95[i]*1.1, n = 500, 
        xlab = "M", ylab = "Probability", main = paste(jd$pop[i], jd$year[i]),
        cex.axis = 1.2, cex.lab = 1.5)
  abline(v = jd$L95[i], col = "blue", lty = 2)
  abline(v = qlnorm(0.025, log(jd$M_obs[i]), jd$tau_M_obs[i]), col = "blue")
  abline(v = jd$U95[i], col = "blue", lty = 2)
  abline(v = qlnorm(0.975, log(jd$M_obs[i]), jd$tau_M_obs[i]), col = "blue")
}

# (2) compute median based on mean and tau (tau already computed from mean and SD),
#     compare to empirically reported median
juv_data %>% filter(origin=="Wild") %>% 
  mutate(median_solve = exp(log(mean) - 0.5*tau_M_obs^2)) %>% 
  ggplot(aes(x = median_solve, y = M_obs)) + geom_point(size = 2) + 
  geom_abline(intercept = 0, slope = 1) + scale_x_log10() + scale_y_log10() +
  theme_bw()


# spawner obs error

windows()
spawner_data %>% filter(disposition_HW=="W") %>% mutate(CV = SD/S_obs) %>% 
  ggplot(aes(x = method, y = CV)) + geom_boxplot() + theme_bw()


