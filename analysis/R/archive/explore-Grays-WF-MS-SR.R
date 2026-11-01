library(scatterplot3d)

dr <- as.matrix(fit_Ricker, c("S", "M")) %>% as_draws_rvars()

dat <- fish_data %>% 
  mutate(S = dr$S, M = dr$M) %>% 
  filter(grepl("Grays", pop) & pop_type == "natural") %>% 
  select(pop, year, S_obs, S, M_obs, M) %>% 
  mutate(pop = substring(pop, 7, 8)) %>% 
  pivot_wider(names_from = "pop", values_from = c("S_obs","S","M_obs","M")) %>% 
  as.data.frame()

# 3D scatterplot: MS and WF spawners vs MS + WF smolts (removing CJ)
windows()
plt <- scatterplot3d(dat$S_obs_MS/1e3, dat$S_obs_WF/1e3, 
                     (dat$M_obs_MS - dat$M_obs_CJ)/1e6, 
                     type = "h", pch = 16, angle = 60, y.margin.add = 0.5,
                     xlab = bquote("Spawners MS (" * 10^3 * ")"), ylab = "",
                     zlab = bquote("Smolts MS + WF (" * 10^6 * ")"))

plt$points3d(median(dat$S_MS)/1e3, median(dat$S_WF)/1e3, 
             median(dat$M_MS + dat$M_WF)/1e6,
             type = "h", pch = 23, col = "slategray4", bg = "white")

text(9.5, 1, bquote("Spawners WF (" * 10^3 * ")"), srt = 65)

# states only: proportion smolts MS / (MS + WF) vs. total spawners MS + WF
windows()
plot(median(dat$S_MS + dat$S_WF)/1e3, median(dat$M_MS / (dat$M_MS + dat$M_WF)),
     xlab = bquote("Spawners MS + WF (" * 10^3 * ")"), ylab = "Smolts MS / (MS + WF)",
     las = 1, pch = 23, col = "slategray4", bg = "white")

# bivariate posterior of psi for Grays MS and WF
windows()
pairs(qlogis(as.matrix(fit_Ricker, c('psi[13]','psi[14]'))))



