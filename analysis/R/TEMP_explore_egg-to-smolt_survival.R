# Why does egg-to-smolt survival (psi) want to be so high?
# Upper end of hyper-mean and pop-level estimates bump up against 1,
# which I suspect is the cause of the severe numerical difficulties
# and slow sampling / large treedepth in this version of the model.
# Is the model somehow underestimating eggs (E) based on spawner
# abundance (S) and weighted mean fecundity (mu_E, q), leading to
# an overestimate of phi?

mod <- LCRchum_Ricker

# data
pop <- fish_data_SMS$pop
year <- fish_data_SMS$year
S_obs <- fish_data_SMS$S_obs
M_obs <- fish_data_SMS$M_obs
# states
S <- extract1(mod,"S")
mu_E <- extract1(mod,"mu_E")
p_F <- extract1(mod,"p_F")
q <- extract1(mod,"q")
mu_psi <- extract1(mod,"mu_psi")
psi <- extract1(mod,"psi")
eta_year_M <- extract1(mod,"eta_year_M")
zeta_M <- stan_mean(mod,"zeta_M")   # only means available
epsilon_M <- outer(extract1(mod,"sigma_M"), zeta_M)
error_M <- eta_year_M[,as.numeric(factor(year))] + epsilon_M
M <- extract1(mod,"M")

# weighted fecundity and eggs (assuming 50:50 sex ratio)
# NB: This assumes egg production is density-independent and deterministic
#     (or "potential eggs", if you like)
f <- apply(sweep(q, c(1,3), mu_E, "*"), 1:2, sum)
E_hat <- p_F*S*f # states
E1 <- sweep(f, 2, colMedians(p_F)*S_obs, "*") # use observed spawners

# back-calculate egg-to-smolt survival based on E and either M or M_obs
# note 1-yr lag between eggs and smolts
# NB: This somewhat misleadingly assumes all the process error is in survival,
#     since eggs are calculated deterministically -- it could just as well be reversed
psi1 <- array(NA, dim(E_hat))
psi2 <- array(NA, dim(E_hat))
for(i in levels(pop)) {
  psi1[,pop==i][,-sum(pop==i)] <- M[,pop==i][,-1] / E_hat[,pop==i][,-sum(pop==i)]
  psi2[,pop==i][,-sum(pop==i)] <- sweep(1/E_hat[,pop==i][,-sum(pop==i)], 2, M_obs[pop==i][-1], "*")
}

# posterior distribution of mu_psi
hist(mu_psi, 20)
hist(qlogis(mu_psi), 20)

# posterior distributions of model-fitted pop-level psi
par(mfrow = c(4,3), mar = c(4.5,3,2,1))
for(i in 1:ncol(psi))
  hist(psi[,i], 20, prob = TRUE, xlab = "psi", ylab = "", main = levels(pop)[i])

# time series of back-calculated egg-to-smolt survival 
# many / most are > 1 !?!
data.frame(pop = pop, year = year, psi1_l = colQuantiles(psi1, probs = 0.05), 
           psi1_m = colMedians(psi1), psi1_u = colQuantiles(psi1, probs = 0.95),
           psi2_l = colQuantiles(psi2, probs = 0.05), psi2_m = colMedians(psi2), 
           psi2_u = colQuantiles(psi2, probs = 0.95)) %>% 
  ggplot(aes(x = year, y = psi1_m, ymin = psi1_l, ymax = psi1_u)) +   # edit to plot psi1 or psi2
  geom_ribbon(fill = "darkblue", alpha = 0.5) +
  geom_line(lwd = 1, col = "darkblue") + 
  geom_abline(intercept = 1, slope = 0) +
  ylab("M / E_hat") + facet_wrap(vars(pop), ncol = 3) + theme_bw()
  
# So if M frequently exceeds E, the annual and/or unique log-recruitment process errors
# must be disproportionately positive, right?
# => Not really => The problem is the mean relationship (psi too close to 1)
par(mar = c(8,4,1,1))
boxplot(split(zeta_M, pop), ylab = "zeta_M", las = 2)
abline(h = 0)

data.frame(year = sort(unique(year)), l = colQuantiles(eta_year_M, probs = 0.05), 
           m = colMedians(eta_year_M), u = colQuantiles(eta_year_M, probs = 0.95)) %>% 
  ggplot(aes(x = year, y = m, ymin = l, ymax = u)) +
  geom_ribbon(fill = "darkblue", alpha = 0.5) +
  geom_line(lwd = 1, col = "darkblue") + 
  geom_abline(intercept = 0, slope = 0) + 
  ylab("eta_year_M") + theme_bw()

# Maybe the positive process errors coincide with M > E ?
# => Sort of
plot(zeta_M, colMedians(psi1))
abline(h = 1)

plot(colMedians(eta_year_M)[as.numeric(factor(year))], colMedians(psi1), xlab = "eta_year_M")
abline(h = 1)

# Maybe positive smolt recruitment process errors are associated with high F:M ratios?
# => Even with sex ratio in the model, relationship is weak
bio_data %>% rename(pop = disposition) %>% filter(!grepl("Hatchery|Duncan_Creek", pop)) %>% 
  group_by(pop, year, sex) %>% summarize(n = sum(count)) %>% 
  dcast(pop + year ~ sex, value.var = "n", fun.aggregate = sum) %>% 
  mutate(total = `F` + M, prop_female = `F`/total) %>% 
  right_join(fish_data_SMS, by = c("pop","year")) %>% 
  mutate(error_M = colMeans(error_M), M0_div_E_hat = colMedians(psi1)) %>% 
  ggplot(aes(x = prop_female, y = M0_div_E_hat)) + geom_point(size = 1) + 
  facet_wrap(vars(pop), ncol = 3) + theme_bw()

# Now that we've added sex composition to the model,
# plot the estimates and data
cbind(fish_data_SMS, colQuantiles(p_F, probs = c(0.05,0.95))) %>%
  mutate(total = n_M_obs + n_F_obs) %>% 
  cbind(., with(., binconf(x = n_F_obs, n = total))) %>%
  rename(prop_female = PointEst) %>% 
  ggplot(aes(x = year, y = prop_female, ymin = Lower, ymax = Upper)) + 
  geom_abline(intercept = 0.5, slope = 0, color = "gray") + 
  geom_ribbon(aes(ymin = `5%`, ymax = `95%`), fill = "darkblue", alpha = 0.5) +
  geom_point(size = 2) + geom_line() + geom_errorbar(width = 0) +
  facet_wrap(vars(pop), ncol = 3) + theme_bw()

# Let's look at the brood years that give rise to these estimates of M >= E_hat
dat <- fish_data_SMS %>% 
  mutate(S = colMedians(S), f = colMedians(f), E_hat = colMedians(E_hat), M0 = stan_mean(mod,"M0"), 
         M_downstream = stan_mean(mod,'M_downstream'), M0_div_E_hat = colMedians(psi1)) %>%
  group_by(pop) %>% 
  mutate(M0_obs = lead(M_obs), tau_M0_obs = lead(tau_M_obs), M0_downstream = lead(M_downstream)) %>% 
  ungroup() %>% as.data.frame() %>% 
  select(pop, year, S_obs, tau_S_obs, S, f, E_hat, M0_obs, tau_M0_obs, M0, 
         downstream_trap, M0_downstream, M0_div_E_hat)

dat %>% filter(M0_div_E_hat >= 1) %>% select(-downstream_trap) %>% format(digits=2)
# just the pops with smolt data
dat %>% select(-downstream_trap) %>% 
  filter(M0_div_E_hat >= 1 & pop %in% c("Duncan_Channel","Grays_CJ","Grays_MS",
                                        "Grays_WF","Hamilton_Channel","Hardy_Creek")) %>% 
  format(digits=2)

# write out this data frame
write.csv(dat, file = here("analysis","results","TEMP_fish_data_egg-to-smolt.csv"), row.names = FALSE)

# Spawner-recruit pairs that produced M > E_hat
# filled circles indicate pairs that have a corresponding (S_obs, M_obs) pair
dat %>% 
  ggplot(aes(x = E_hat/1e3, y = M0/E_hat, 
         pch = !is.na(S_obs) & (!is.na(M0_obs) | !is.na(M0_obs[downstream_trap])))) +
  geom_point(size = 2.5, color = "darkblue", alpha = 0.7) + scale_shape_manual(values = c(1,16)) +
  geom_hline(yintercept = 1) + labs(x = "E_hat (thousands)", y = "M / E_hat") +
  facet_wrap(vars(pop), ncol = 3, scales = "free") + theme_bw() + theme(legend.position = "none")

# clean up
rm(list = c('mod','S','mu_E','q','M','S_obs','M_obs','f','E_hat','pop','year',
            'psi','mu_psi','psi1','psi2'))


