# Why does egg-to-smolt survival (psi) want to be so high?
# Upper end of hyper-mean and pop-level estimates bump up against 1,
# which I suspect is the cause of the severe numerical difficulties
# and slow sampling / large treedepth in this version of the model.
# Is the model somehow underestimating eggs (E) based on spawner
# abundance (S) and weighted mean fecundity (mu_E, q), leading to
# an overestimate of phi?

mod <- LCRchum_BH

# states
S <- extract1(mod,"S")
mu_E <- extract1(mod,"mu_E")
q <- extract1(mod,"q")
mu_psi <- extract1(mod,"mu_psi")
psi <- extract1(mod,"psi")
eta_year_M <- extract1(mod,"eta_year_M")
zeta_M <- stan_mean(mod,"zeta_M")   # only means available
M <- extract1(mod,"M")
# data
pop <- fish_data_SMS$pop
year <- fish_data_SMS$year
S_obs <- fish_data_SMS$S_obs
M_obs <- fish_data_SMS$M_obs

# weighted fecundity and eggs (assuming 50:50 sex ratio)
# NB: This assumes egg production is density-independent and deterministic
#     (or "potential eggs", if you like)
f <- apply(sweep(q, c(1,3), mu_E, "*"), 1:2, sum)
E <- 0.5*S*f # states
E1 <- sweep(f, 2, 0.5*S_obs, "*") # use observed spawners

# back-calculate egg-to-smolt survival based on E and either M or M_obs
# note 1-yr lag between eggs and smolts
# NB: This somewhat misleadingly assumes all the process error is in survival,
#     since eggs are calculated deterministically -- it could just as well be reversed
psi1 <- array(NA, dim(E))
psi2 <- array(NA, dim(E))
for(i in levels(pop)) {
  psi1[,pop==i][,-sum(pop==i)] <- M[,pop==i][,-1] / E[,pop==i][,-sum(pop==i)]
  psi2[,pop==i][,-sum(pop==i)] <- sweep(1/E[,pop==i][,-sum(pop==i)], 2, M_obs[pop==i][-1], "*")
}

# posterior distribution of mu_psi
hist(mu_psi, 20)
hist(qlogis(mu_psi), 20)

# posterior distributions of model-fitted pop-level psi
par(mfrow = c(3,4), mar = c(4.5,3,2,1))
for(i in 1:ncol(psi))
  hist(psi[,i], 20, prob = TRUE, xlab = "psi", ylab = "", main = levels(pop)[i])

for(i in 1:ncol(psi))
  hist(qlogis(psi[,i]), 20, prob = TRUE, xlab = "logit(psi)", ylab = "", main = levels(pop)[i])

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
  facet_wrap(vars(pop), ncol = 4) + theme_bw()
  
# So if M frequently exceeds E, the annual and/or unique log-recruitment process errors
# must be disproportionately positive, right?
# => Not really
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
# => Not even a hint of a relationship
bio_data %>% rename(pop = disposition) %>% filter(HW=="W" & !grepl("Hatchery|Duncan_Creek", pop)) %>% 
  group_by(pop, year, sex) %>% summarize(n = sum(count)) %>% 
  dcast(pop + year ~ sex, value.var = "n", fun.aggregate = sum) %>% 
  mutate(total = Female + Male, prop_female = Female/total) %>% 
  right_join(fish_data_SMS, by = c("pop","year")) %>% mutate(zeta_M = zeta_M) %>% 
  ggplot(aes(x = prop_female, y = zeta_M)) + geom_point(size = 1) + 
  facet_wrap(vars(pop), ncol = 4) + theme_bw()


# clean up
rm(list = c('mod','S','mu_E','q','M','S_obs','M_obs','f','E','pop','year',
            'psi','mu_psi','psi1','psi2'))


