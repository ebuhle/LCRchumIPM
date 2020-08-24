dat <- fish_data %>% filter(pop %in% c('Grays_MSWF','Grays_CJ')) %>% 
  select(pop, year, M_obs, tau_M_obs) %>% arrange(pop,year) %>% na.omit() %>% 
  pivot_wider(id_cols = year, names_from = pop, values_from = c(M_obs, tau_M_obs)) %>% 
  arrange(year)

windows(width=10, height=7)

par(mfrow=c(3,5), mar = c(4,4,2,1))
for(i in 1:nrow(dat)) {
  with(slice(dat, i), {
    mswf <- rlnorm(10000, log(M_obs_Grays_MSWF), tau_M_obs_Grays_MSWF)
    cj <- rlnorm(10000, log(M_obs_Grays_CJ), tau_M_obs_Grays_CJ) 
    mswf_only <- pmax(mswf - cj, 0)

    if(any(!is.na(mswf_only)))
      hist(log(mswf_only), 20, prob = TRUE, main = year)
  })
}


par(mfrow=c(3,5), mar = c(4,4,2,1))
for(i in 1:nrow(dat)) {
  with(slice(dat, i), {
    mswf <- rlnorm(10000, log(M_obs_Grays_MSWF), tau_M_obs_Grays_MSWF)
    cj <- rlnorm(10000, log(M_obs_Grays_CJ), tau_M_obs_Grays_CJ) 
    mswf_only <- pmax(mswf - cj, 0)
    
    if(any(!is.na(mswf_only))) {
      qqnorm(log(mswf_only[mswf_only > 0]), main = year)
      qqline(log(mswf_only[mswf_only > 0]))
    }
  })
}

