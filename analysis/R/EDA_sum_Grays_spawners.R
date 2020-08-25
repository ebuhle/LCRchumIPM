dat <- spawner_data_agg %>% filter(pop %in% c('Grays_MS','Grays_WF')) %>% 
  select(pop, year, S_obs, tau_S_obs) %>% arrange(pop,year) %>% na.omit() %>% 
  pivot_wider(id_cols = year, names_from = pop, values_from = c(S_obs, tau_S_obs)) %>% 
  arrange(year)
  
windows(width = 10, height = 7)
par(mfrow=c(3,5), mar = c(4,4,2,1))
for(i in 1:nrow(dat)) {
  with(slice(dat, i), {
       ms <- rlnorm(10000, log(S_obs_Grays_MS), tau_S_obs_Grays_MS)
       wf <- rlnorm(10000, log(S_obs_Grays_WF), tau_S_obs_Grays_WF) 
       mean_ms <- S_obs_Grays_MS * exp(0.5*tau_S_obs_Grays_MS^2)
       var_ms <- (exp(tau_S_obs_Grays_MS^2) - 1)*mean_ms^2
       mean_wf <- S_obs_Grays_WF * exp(0.5*tau_S_obs_Grays_WF^2)
       var_wf <- (exp(tau_S_obs_Grays_WF^2) - 1)*mean_wf^2
       mean_mswf <- mean_ms + mean_wf
       var_mswf <- var_ms + var_wf
       sigma_mswf <- sqrt(log(var_mswf/mean_mswf^2 + 1))
       mu_mswf <- log(mean_mswf) - sigma_mswf/2

       if(any(!is.na(ms + wf))) {
         hist(log(ms + wf), 20, prob = TRUE, main = year)
         curve(dnorm(x, mu_mswf, sigma_mswf), col = "blue", add = TRUE)
       }
  })
}


windows(width = 10, height = 7)
par(mfrow=c(3,5), mar = c(4,4,2,1))
for(i in 1:nrow(dat)) {
  with(slice(dat, i), {
    ms <- rlnorm(10000, log(S_obs_Grays_MS), tau_S_obs_Grays_MS)
    wf <- rlnorm(10000, log(S_obs_Grays_WF), tau_S_obs_Grays_WF) 
    
    if(any(!is.na(ms + wf))) {
      qqnorm(log(ms + wf), main = year)
      qqline(log(ms + wf))
    }
  })
}

