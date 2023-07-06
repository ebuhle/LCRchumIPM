eta_year_MS = as.matrix(LCRchum_Ricker_foreHmax,'eta_year_MS')
eta_mean_MS = rowMeans(eta_year_MS[, tail(1:43,20)])
S = as.matrix(LCRchum_Ricker_foreHmax,'S')
# p_HOS = as.matrix(LCRchum_Ricker_foreHmax,'p_HOS')
# S_W = S * (1 - p_HOS)
S = S[, fish_data_fore$forecast & fish_data_fore$pop=='Grays CJ']
log_S_mean = rowMeans(log(S))
plot(eta_mean_MS, log_S_mean)


eta_year_MS = as.matrix(LCRchum_Ricker_foreHmax,'eta_year_MS')
eta_mean_MS = rowMeans(eta_year_MS[, tail(1:43,20)])
S = as.matrix(LCRchum_Ricker_foreHmax,'S')
S = S[, fish_data_fore$forecast & fish_data_fore$pop=='Grays CJ']
quasi_extinct = rowAnys(S < 50)
plot(eta_mean_MS, jitter(as.integer(quasi_extinct)))


eta_year_MS_H0 = as.matrix(LCRchum_Ricker_foreH0,'eta_year_MS')
eta_mean_MS_H0 = rowMeans(eta_year_MS_H0[, tail(1:43,20)])
eta_year_MS_Hmax = as.matrix(LCRchum_Ricker_foreHmax,'eta_year_MS')
eta_mean_MS_Hmax = rowMeans(eta_year_MS_Hmax[, tail(1:43,20)])
boxplot(eta_mean_MS_H0, eta_mean_MS_Hmax)


eta_year_MS = as_draws_rvars(as.matrix(LCRchum_Ricker,'eta_year_MS'))
eta_year_MS_foreH0 = as_draws_rvars(as.matrix(LCRchum_Ricker_foreH0,'eta_year_MS'))
eta_year_MS_foreHmax = as_draws_rvars(as.matrix(LCRchum_Ricker_foreHmax,'eta_year_MS'))
dat = data.frame(year = c(sort(unique(fish_data$year)), rep(sort(unique(fish_data_fore$year)), 2)),
                 model = rep(c('Ricker','Ricker_foreH0','Ricker_foreHmax'), times = c(23,43,43)),
                 eta_year_MS = c(eta_year_MS$eta_year_MS, 
                                 eta_year_MS_foreH0 = eta_year_MS_foreH0$eta_year_MS,
                                 eta_year_MS_foreHmax = eta_year_MS_foreHmax$eta_year_MS))

windows(width = 11, height = 7)
dat %>% 
  ggplot(aes(x = year, y = mean(eta_year_MS), color = model, fill = model)) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = as.vector(quantile(eta_year_MS, 0.05)), 
                  ymax = as.vector(quantile(eta_year_MS, 0.95))),
              alpha = 0.5) + 
  scale_color_manual(values = c('black','skyblue','salmon'), aesthetics = c('color','fill')) +
  ylab('eta_year_MS')


windows()
plot(stan_mean(LCRchum_Ricker,'B_rate_all'), 
     stan_mean(LCRchum_Ricker_foreHmax,'B_rate_all')[!fish_data_fore$forecast])
abline(0,1)


B_rate_all = as_draws_rvars(as.matrix(LCRchum_Ricker_foreHmax, 'B_rate_all'))
windows(width = 11, height = 7)
fish_data_fore %>% mutate(B_rate_all = B_rate_all$B_rate_all) %>% filter(pop_type == 'natural') %>% 
  ggplot(aes(x = year, y = median(B_rate_all))) +
  geom_line(linewidth = 1, color = 'slategray4') +
  geom_ribbon(aes(ymin = as.vector(quantile(B_rate_all, 0.05)), 
                  ymax = as.vector(quantile(B_rate_all, 0.95))),
              fill = 'slategray4', alpha = 0.5) +
  ylim(c(0,1)) + ylab('B_rate') +
  facet_wrap(vars(pop), ncol = 4) +
  theme(panel.grid.minor = element_blank(), 
        strip.background = element_rect(fill = NA),
        strip.text = element_text(margin = margin(b = 3, t = 3)))

B_draws = as_draws_rvars(as.matrix(LCRchum_Ricker_foreHmax, c('S','B_rate_all'))) %>% 
  mutate_variables(B_take = B_rate_all * S / (1 - B_rate_all))
windows(width = 11, height = 7)
fish_data_foreHmax %>% mutate(B_take = B_draws$B_take) %>% filter(pop_type == 'natural') %>% 
  ggplot(aes(x = year, y = median(B_take))) +
  geom_line(linewidth = 1, color = 'slategray4') +
  geom_ribbon(aes(ymin = as.vector(quantile(B_take, 0.05)), 
                  ymax = as.vector(quantile(B_take, 0.95))),
              fill = 'slategray4', alpha = 0.5) +
  geom_point(aes(y = B_take_obs), size = 2) + 
  scale_y_log10() +
  ylab('Broodstock') +
  facet_wrap(vars(pop), ncol = 4) +
  theme(panel.grid.minor = element_blank(), 
        strip.background = element_rect(fill = NA),
        strip.text = element_text(margin = margin(b = 3, t = 3)))

windows()
fish_data_fore %>% 
  mutate(B_rate_all = B_rate_all$B_rate_all) %>% filter(pop == 'Ives' & forecast) %>% 
  ggplot(aes(x = year, ydist = B_rate_all)) +
  stat_eye(.width = c(0.5, 0.9), normalize = "groups", 
           point_size = 2, point_color = 'slategray4', interval_color = 'slategray4',
           slab_fill = 'slategray4', slab_alpha = 0.5) +
  ylim(c(0,1)) + ylab('B_rate')


  
  


                 
                 
                 
                 