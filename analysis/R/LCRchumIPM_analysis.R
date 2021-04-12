#===========================================================================
# SETUP
#===========================================================================

## @knitr getting_started
options(device = ifelse(.Platform$OS.type == "windows", "windows", "quartz"))
options(mc.cores = parallel::detectCores(logical = FALSE) - 1)

library(salmonIPM)
library(rstan)
library(shinystan)
library(matrixStats)
library(Hmisc)
library(tibble)
library(dplyr)
library(tidyr)
library(reshape2)
library(yarrr)
library(magicaxis)
library(viridis)
library(zoo)
library(here)

# load data
source(here("analysis","R","LCRchumIPM_data-processing.R"))
# load saved stanfit objects
if(file.exists(here("analysis","results","LCRchumIPM.RData")))
  load(here("analysis","results","LCRchumIPM.RData"))
## @knitr


#===========================================================================
# FIT RETROSPECTIVE MODELS
#===========================================================================

#--------------------------------------------------------------
# Spawner-to-spawner IPM
#--------------------------------------------------------------

# Density-independent
## @knitr fit_SS_exp
SS_exp <- salmonIPM(fish_data = fish_data_SS, stan_model = "IPM_SS_pp", SR_fun = "exp",
                    pars = c("B_rate_all","mu_Rmax","sigma_Rmax","Rmax"), include = FALSE, 
                    log_lik = TRUE, chains = 3, iter = 1500, warmup = 500,
                    control = list(adapt_delta = 0.99, max_treedepth = 13))

## @knitr print_SS_exp
print(SS_exp, prob = c(0.025,0.5,0.975),
      pars = c("alpha","phi","p_HOS","q","gamma","p","S","R","LL"), 
      include = FALSE, use_cache = FALSE)
## @knitr

# Beverton-Holt
## @knitr fit_SS_BH
SS_BH <- salmonIPM(fish_data = fish_data_SS, stan_model = "IPM_SS_pp", SR_fun = "BH",
                   pars = "B_rate_all", include = FALSE, log_lik = TRUE, 
                   chains = 3, iter = 1500, warmup = 500,
                   control = list(adapt_delta = 0.99, max_treedepth = 13))

## @knitr print_SS_BH
print(SS_BH, prob = c(0.025,0.5,0.975),
      pars = c("alpha","Rmax","phi","p_HOS","q","gamma","p","S","R","LL"), 
      include = FALSE, use_cache = FALSE)
## @knitr

# Ricker
## @knitr fit_SS_Ricker
SS_Ricker <- salmonIPM(fish_data = fish_data_SS, stan_model = "IPM_SS_pp", SR_fun = "Ricker",
                       pars = "B_rate_all", include = FALSE, log_lik = TRUE, 
                       chains = 3, iter = 1500, warmup = 500,
                       control = list(adapt_delta = 0.99, max_treedepth = 13))

## @knitr print_SS_Ricker
print(SS_Ricker, prob = c(0.025,0.5,0.975),
      pars = c("alpha","Rmax","phi","p_HOS","q","gamma","p","S","R","LL"), 
      include = FALSE, use_cache = FALSE)
## @knitr


#--------------------------------------------------------------
# Lower Columbia chum spawner-smolt-spawner IPM
#--------------------------------------------------------------

# Density-independent
## @knitr fit_LCRchum_exp
LCRchum_exp <- salmonIPM(fish_data = fish_data_SMS,  fecundity_data = fecundity_data,
                         ages = list(M = 1), stan_model = "IPM_LCRchum_pp", SR_fun = "exp",
                         pars = c("B_rate_all","mu_Emax","sigma_Emax","Emax"), 
                         include = FALSE, log_lik = TRUE, 
                         chains = 3, iter = 1500, warmup = 500,
                         control = list(adapt_delta = 0.99, max_treedepth = 15))

## @knitr print_LCRchum_exp
print(LCRchum_exp, prob = c(0.025,0.5,0.975),
      pars = c("eta_pop_EM","eta_year_EM","eta_year_MS","eta_pop_p","p",
               "tau_M","tau_S","p_HOS","E","S","M","s_EM","s_MS","q","LL"), 
      include = FALSE, use_cache = FALSE)
## @knitr

# Beverton-Holt
## @knitr fit_LCRchum_BH
LCRchum_BH <- salmonIPM(fish_data = fish_data_SMS, fecundity_data = fecundity_data,
                        ages = list(M = 1), stan_model = "IPM_LCRchum_pp", SR_fun = "BH",
                        pars = "B_rate_all", include = FALSE, log_lik = TRUE, 
                        chains = 3, iter = 1500, warmup = 500,
                        control = list(adapt_delta = 0.95, max_treedepth = 13))

## @knitr print_LCRchum_BH
print(LCRchum_BH, prob = c(0.025,0.5,0.975),
      pars = c("psi","Mmax","eta_year_M","eta_year_MS","eta_pop_p",
               "p","tau_M","tau_S","p_HOS","E_hat","M","S","s_MS","q","LL"), 
      include = FALSE, use_cache = FALSE)
## @knitr

# Ricker
## @knitr fit_LCRchum_Ricker
LCRchum_Ricker <- salmonIPM(fish_data = fish_data_SMS, fecundity_data = fecundity_data,
                            ages = list(M = 1), stan_model = "IPM_LCRchum_pp", SR_fun = "Ricker",
                            pars = "B_rate_all", include = FALSE, log_lik = TRUE, 
                            chains = 3, iter = 1500, warmup = 500,
                            control = list(adapt_delta = 0.95, max_treedepth = 14))

## @knitr print_LCRchum_Ricker
print(LCRchum_Ricker, prob = c(0.025,0.5,0.975),
      pars = c("psi","Mmax","eta_year_M","eta_year_MS","eta_pop_p","mu_pop_alr_p",
               "p","p_F","tau_M","tau_S","p_HOS","E_hat","M","S","s_MS","q","q_F","LL"), 
      include = FALSE, use_cache = FALSE)
## @knitr


#--------------------------------------------------------------
# Model selection using LOO
#--------------------------------------------------------------

# Observationwise log-likelihood of each fitted model
# Here an observation is a row of fish_data, and the total likelihood includes 
# components for smolt abundance, spawner abundance, age-frequency, H/W-frequency,
# and smolt and spawner observation error SDs
## @knitr loo_LCRchum
LL_LCRchum <- lapply(list(exp = LCRchum_exp, BH = LCRchum_BH, Ricker = LCRchum_Ricker),
                 loo::extract_log_lik, parameter_name = "LL", merge_chains = FALSE)

# Relative ESS of posterior draws of observationwise likelihood 
r_eff_LCRchum <- lapply(LL_LCRchum, function(x) relative_eff(exp(x)))

# PSIS-LOO
LOO_LCRchum <- lapply(1:length(LL_LCRchum), 
                      function(i) loo(LL_LCRchum[[i]], r_eff = r_eff_LCRchum[[i]]))
names(LOO_LCRchum) <- names(LL_LCRchum)

## Compare all three models
loo_compare(LOO_LCRchum)

## Exponential vs. Ricker
loo_compare(LOO_LCRchum[c("exp","Ricker")])

## Exponential vs. Beverton-Holt
loo_compare(LOO_LCRchum[c("exp","BH")])

## Beverton-Holt vs. Ricker
loo_compare(LOO_LCRchum[c("BH","Ricker")])
## @knitr


#===========================================================================
# FIT PROSPECTIVE FORECASTING MODELS
#===========================================================================

#--------------------------------------------------------------
# Lower Columbia chum spawner-smolt-spawner IPM
#--------------------------------------------------------------

# Ricker
## @knitr fit_LCRchum_BH_fore
LCRchum_BH_fore <- salmonIPM(fish_data = fish_data_SMS_fore, fecundity_data = fecundity_data,
                             ages = list(M = 1), stan_model = "IPM_LCRchum_pp", SR_fun = "BH",
                             pars = c("Emax","p","p_HOS","B_rate_all"), 
                             include = FALSE, chains = 3, iter = 1500, warmup = 500,
                             control = list(adapt_delta = 0.99, max_treedepth = 15))

## @knitr
print(LCRchum_BH_fore, prob = c(0.025,0.5,0.975),
      pars = c("eta_pop_EM","eta_year_EM","eta_year_MS","eta_pop_p",
               "E","S","M","s_EM","s_MS","q","tau_M","tau_S"), 
      include = FALSE, use_cache = FALSE)


#--------------------------------------------------------------
# Save stanfit objects
#--------------------------------------------------------------

save(list = ls()[sapply(ls(), function(x) do.call(class, list(as.name(x)))) == "stanfit"], 
     file = here("analysis","results","LCRchumIPM.RData"))



#===========================================================================
# FIGURES
#===========================================================================

#--------------------------------------------------------------------------#
############# FIGURES THAT WORK WITH LOWER COLUMBIA CHUM IPM ###############
#--------------------------------------------------------------------------#

#--------------------------------------------------------------------
# Lower Columbia chum spawner-egg-smolt-spawner #
# S-R curves (spawners to smolts)
# Posterior distributions of fecundity and Mmax parameters
# Time series of smolt productivity process errors and SAR
#--------------------------------------------------------------------

mod_name <- "LCRchum_Ricker"
save_plot <- TRUE

if(save_plot) {
  png(filename=here("analysis","results",paste0("SR_",mod_name,".png")),
      width=7*0.9, height=8*0.9, units="in", res=300, type="cairo-png")
} else dev.new(width = 7, height = 8)

## @knitr plot_LCM_params
dat <- fish_data_SMS
SR_fun <- unlist(strsplit(mod_name, "_"))[2]

SR_eval <- function(alpha, Rmax = NULL, S, SR_fun) 
{
  switch(SR_fun,
         exp = alpha*S,
         BH = alpha*S/(1 + alpha*S/Rmax),
         Ricker = alpha*S*exp(-alpha*S/(exp(1)*Rmax)))
}

# fecundity
mu_E <- do.call(extract1, list(as.name(mod_name), "mu_E"))
ages <- substring(names(select(dat, starts_with("n_age"))), 6, 6)
# egg deposition
q <- do.call(extract1, list(as.name(mod_name), "q"))
q_F <- do.call(extract1, list(as.name(mod_name), "q_F"))
mu_psi <- do.call(extract1, list(as.name(mod_name), "mu_psi"))
psi <- do.call(extract1, list(as.name(mod_name), "psi"))
alpha <- apply(sweep(q, c(1,3), mu_E, "*"), 1:2, sum) * q_F * psi[,as.numeric(dat$pop)]
alpha <- t(as.matrix(aggregate(t(alpha), list(pop = dat$pop), median)[,-1]))
mu_alpha <- rowMeans(log(alpha))
mu_Mmax <- as.vector(do.call(extract1, list(as.name(mod_name), "mu_Mmax")))
Mmax <- do.call(extract1, list(as.name(mod_name), "Mmax"))
S <- colMedians(do.call(extract1, list(as.name(mod_name), "S")))
S_grid <- matrix(seq(0, quantile(S/dat$A, 0.9, na.rm = TRUE), length = 100),
                 nrow = length(mu_alpha), ncol = 100, byrow = TRUE)
M_ESU <- SR_eval(alpha = exp(mu_alpha), Rmax = exp(mu_Mmax)*1e3, S = S_grid, SR_fun = SR_fun)
M_pop <- sapply(1:ncol(Mmax), function(i) {
  colMedians(SR_eval(alpha = alpha[,i], Rmax = Mmax[,i]*1e3, S = S_grid, SR_fun = SR_fun))
})
# smolt recruitment process errors
y <- sort(unique(dat$year))
eta_year_M <- do.call(extract1, list(as.name(mod_name), "eta_year_M"))
sigma_M <- do.call(extract1, list(as.name(mod_name), "sigma_M"))
zeta_M <- do.call(stan_mean, list(as.name(mod_name), "zeta_M"))
epsilon_M <- outer(sigma_M, zeta_M, "*")
# SAR
eta_year_MS <- do.call(extract1, list(as.name(mod_name), "eta_year_MS"))
mu_MS <- do.call(extract1, list(as.name(mod_name), "mu_MS"))
s_hat_MS <- plogis(sweep(eta_year_MS, 1, qlogis(mu_MS), "+"))
s_MS <- do.call(extract1, list(as.name(mod_name), "s_MS"))
# colors
c1 <- "slategray4"
c1t <- transparent(c1, trans.val = 0.3)
c1tt <- transparent(c1, trans.val = 0.5)
ac <- viridis(length(ages), end = 0.8, direction = -1, alpha = 0.5) 

par(mfcol = c(3,2), mar = c(5.1,5.1,1,1.5))

# Posterior distributions of fecundity by age
dd_age <- lapply(as.data.frame(mu_E), density)

plot(dd_age[[1]]$x, dd_age[[1]]$y, pch = "", las = 1, cex.axis = 1.2, cex.lab = 1.5,
     xlab = bquote("Mean fecundity (" * mu[italic(E)] * ")"), ylab = "", 
     xlim = range(sapply(dd_age, function(m) m$x)),
     ylim = range(sapply(dd_age, function(m) m$y)),
     xaxs = "i", yaxt = "n")
for(a in 1:length(ages))
  polygon(dd_age[[a]], col = ac[a], border = NA)
title(ylab = "Probability density", line = 1, cex.lab = 1.5)
text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "A", cex = 1.5)
legend("topright", paste("age", ages, "  "), cex = 1.2, fill = ac, border = NA,
       bty = "n", inset = c(-0.05,0))

# Posterior densities of psi
dd_ESU <- density(mu_psi)
dd_pop <- lapply(as.data.frame(psi), density)

plot(dd_ESU$x, dd_ESU$y, pch = "", lwd = 3, col = c1, las = 1, 
     xaxs = "i", yaxt = "n", cex.axis = 1.2, cex.lab = 1.5, 
     xlab = bquote("Maximum egg-to-smolt survival (" * psi * ")"), ylab = "",
     xlim = c(0,1), ylim = range(dd_ESU$y, sapply(dd_pop, function(m) m$y)))
polygon(dd_ESU, col = c1tt, border = NA)
for(i in 1:length(dd_pop))
  lines(dd_pop[[i]]$x, dd_pop[[i]]$y, col = c1t)
title(ylab = "Probability density", line = 1, cex.lab = 1.5)
text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "B", cex = 1.5)

# Posterior densities of log(Mmax)
dd_ESU <- density(mu_Mmax * log10(exp(1)))  # convert to base 10
dd_pop <- lapply(as.data.frame(log10(Mmax)), density)

plot(dd_ESU$x, dd_ESU$y, pch = "", lwd = 3, col = c1, las = 1, 
     xaxt = "n", yaxt = "n", cex.axis = 1.2, cex.lab = 1.5, 
     xlab = bquote("Maximum smolt production (" * italic(M)[max] * ")"), ylab = "",
     xlim = range(dd_ESU$x[dd_ESU$y > 0.01], sapply(dd_pop, function(m) range(m$x[m$y > 0.01]))),
     ylim = range(dd_ESU$y, sapply(dd_pop, function(m) m$y)))
polygon(dd_ESU, col = c1tt, border = NA)
for(i in 1:length(dd_pop))
  lines(dd_pop[[i]]$x, dd_pop[[i]]$y, col = c1t)
tck <- maglab(10^par("usr")[1:2], log = TRUE)
axis(1, at = log10(tck$labat), labels = tck$labat, cex.axis = 1.2)
title(ylab = "Probability density", line = 1, cex.lab = 1.5)
text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "C", cex = 1.5)

# Smolts vs. egg deposition
plot(S_grid[1,], colMedians(M_ESU)*1e-6, type = "l", lwd=3, col = c1, las = 1,
     cex.axis = 1.2, cex.lab = 1.5, xaxs = "i", yaxs = "i",
     ylim = range(M_pop*1e-6), xlab = "Spawners", ylab = "Smolts (millions)")
for(i in 1:ncol(M_pop)) 
  lines(S_grid[1,], M_pop[,i]*1e-6, col = c1t)
polygon(c(S_grid[1,], rev(S_grid[1,])), 
        c(colQuantiles(M_ESU, probs = 0.05), rev(colQuantiles(M_ESU, probs = 0.95)))*1e-6, 
        col = c1tt, border = NA)
rug(S, col = c1)
text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "D", cex = 1.5)

# Smolt recruitment process errors
plot(y, colMedians(eta_year_M), type = "n", las = 1, cex.axis = 1.2, cex.lab = 1.5,
     ylim = range(colQuantiles(eta_year_M, probs = c(0.05, 0.95)), 
                  colQuantiles(epsilon_M, probs = c(0.05, 0.95))), 
     xaxs = "i", xaxt = "n", xlab = "Brood year", 
     ylab = "Smolt recruitment anomaly", xpd = NA)
polygon(c(y, rev(y)), 
        c(colQuantiles(eta_year_M, probs = 0.05), 
          rev(colQuantiles(eta_year_M, probs = 0.95))),
        col = c1tt, border = NA)
lines(y, colMedians(eta_year_M), col = c1t, lwd = 3)
for(j in levels(dat$pop))
  lines(dat$year[dat$pop == j], colMedians(epsilon_M[,dat$pop == j]), col = c1t)
axis(side = 1, at = y[y %% 5 == 0], cex.axis = 1.2)
rug(y[y %% 5 != 0], ticksize = -0.02)
text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "E", cex = 1.5)

# SAR
plot(y, colMedians(s_hat_MS), type = "n", las = 1, cex.axis = 1.2, cex.lab = 1.5,
     ylim = range(0, colQuantiles(s_hat_MS, probs = 0.95), colMedians(s_MS)), 
     xaxs = "i", xaxt = "n", xlab = "Outmigration year", ylab = "")
mtext("SAR", side = 2, line = 3.7, cex = par("cex")*1.5)
polygon(c(y, rev(y)), 
        c(colQuantiles(s_hat_MS, probs = 0.05), rev(colQuantiles(s_hat_MS, probs = 0.95))),
        col = c1tt, border = NA)
lines(y, colMedians(s_hat_MS), col = c1t, lwd = 3)
for(j in levels(dat$pop))
  lines(dat$year[dat$pop == j], colMedians(s_MS[,dat$pop == j]), col = c1t)
axis(side = 1, at = y[y %% 5 == 0], cex.axis = 1.2)
rug(y[y %% 5 != 0], ticksize = -0.02)
text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "F", cex = 1.5)

rm(list=c("mod_name","SR_fun","mu_alpha","mu_Mmax","S","S_grid","M_ESU","q_F",
          "c1","c1t","c1tt","dd_ESU","dd_pop","dd_age","alpha","Mmax","M_pop","ac","ages",
          "y","eta_year_M","sigma_M","zeta_M","epsilon_M","eta_year_MS","mu_MS","s_hat_MS","dat"))
## @knitr
if(save_plot) dev.off()


#--------------------------------------------------------------------------------
# Time series of observed and fitted total spawners or smolts for each pop
#--------------------------------------------------------------------------------

mod_name <- "LCRchum_Ricker"
life_stage <- "S"   # "S" = spawners, "M" = smolts
save_plot <- TRUE

dev.new(width=11,height=7)

## @knitr plot_spawner_smolt_ts
life_cycle <- unlist(strsplit(mod_name, "_"))[1]
forecasting <- ifelse(identical(unlist(strsplit(mod_name, "_"))[3], "fore"), "yes", "no")
N <- do.call(extract1, list(as.name(mod_name), life_stage))
tau <- do.call(extract1, list(as.name(mod_name),
                              switch(life_cycle, SS = "tau",
                                     LCRchum = switch(life_stage, M = "tau_M", S = "tau_S"))))
if(life_cycle == "LCRchum" & life_stage == "M")
  N[,na.omit(fish_data_SMS$downstream_trap)] <- N[,na.omit(fish_data_SMS$downstream_trap)] + 
                                                N[,which(!is.na(fish_data_SMS$downstream_trap))]
N_obs <- N * rlnorm(length(N), 0, tau)

switch(life_cycle, SS = fish_data_SS, 
       LCRchum = switch(forecasting, no = fish_data_SMS, yes = fish_data_SMS_fore)) %>% 
  mutate(N_L = colQuantiles(N, probs = 0.05),
         N_U = colQuantiles(N, probs = 0.95),
         N_obs_L = colQuantiles(N_obs, probs = 0.05),
         N_obs_U = colQuantiles(N_obs, probs = 0.95),
         pch = switch(life_stage, M = ifelse(is.na(tau_M_obs), 1, 16),
                      S = ifelse(is.na(tau_S_obs), 1, 16))) %>% 
  ggplot(aes(x = year, y = !!sym(paste0(life_stage, "_obs")))) +
  geom_ribbon(aes(ymin = N_L, ymax = N_U), fill = "slategray4", alpha = 0.5) +
  geom_ribbon(aes(ymin = N_obs_L, ymax = N_obs_U), fill = "slategray4", alpha = 0.3) +
  geom_point(aes(shape = pch), size = 2.5) + scale_shape_identity() +
  labs(x = "Year", y = switch(life_stage, M = "Smolts (thousands)", S = "Spawners")) + 
  scale_x_continuous(minor_breaks = function(v) min(v):max(v)) + 
  scale_y_log10(labels = function(y) y*switch(life_stage, M = 1e-3, S = 1)) + 
  facet_wrap(vars(pop), ncol = 4) + theme_bw(base_size = 16) +
  theme(panel.grid.minor.y = element_blank(),
        strip.background = element_rect(fill = NA))

## @knitr

if(save_plot)   
  ggsave(filename=here("analysis", "results", paste0(life_stage, "_fit_", mod_name, ".png")),
         width=11, height=7, units="in", dpi=300, type="cairo-png")
rm(list = c("mod_name","forecasting","life_stage","life_cycle","N","N_obs","tau"))


#--------------------------------------------------------------------------------
# Time series of observed and fitted spawner age structure for each pop
#--------------------------------------------------------------------------------

mod_name <- "LCRchum_Ricker"
save_plot <- FALSE

dev.new(width=12,height=7)

## @knitr plot_spawner_age_ts
life_cycle <- unlist(strsplit(mod_name, "_"))[1]
q <- do.call(extract1, list(as.name(mod_name), "q"))

switch(life_cycle, SS = fish_data_SS, LCRchum = fish_data_SMS) %>% 
  select(pop, year, starts_with("n_age")) %>% 
  mutate(total = rowSums(across(starts_with("n_age"))),
         across(starts_with("n_age"), ~ binconf(.x, total, alpha = 0.1))) %>% 
  do.call(data.frame, .) %>% # unpack cols of nested data frames
  pivot_longer(cols = -c(pop, year, total), names_to = c("age",".value"),
               names_pattern = "n_age(.)_obs.(.*)") %>% 
  cbind(array(aperm(sapply(1:3, function(k) colQuantiles(q[,,k], probs = c(0.05, 0.5, 0.95)), 
                           simplify = "array"), c(3,1,2)), dim = c(nrow(.), 3), 
              dimnames = list(NULL, paste0("q_age_", c("L","m","U"))))) %>%
  ggplot(aes(x = year, group = age, color = age, fill = age)) +
  geom_line(aes(y = q_age_m), lwd = 1, alpha = 0.8) +
  geom_ribbon(aes(ymin = q_age_L, ymax = q_age_U), color = NA, alpha = 0.3) +
  geom_point(aes(y = PointEst), pch = 16, size = 2.5, alpha = 0.8) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0, alpha = 0.8) +
  scale_color_manual(values = viridis(3, end = 0.8, direction = -1)) +
  scale_fill_manual(values = viridis(3, end = 0.8, direction = -1)) +
  scale_x_continuous(minor_breaks = function(v) min(v):max(v)) +
  labs(x = "Year", y = "Proportion at age") + 
  facet_wrap(vars(pop), ncol = 4) + theme_bw(base_size = 16) + 
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        strip.background = element_rect(fill = NA), legend.box.margin = margin(0,-10,0,-15))

## @knitr

if(save_plot)
  ggsave(filename=here("analysis", "results", paste0("q_fit_", mod_name, ".png")),
         width=12, height=7, units="in", dpi=300, type="cairo-png")
rm(list = c("mod_name","life_cycle","q"))


#--------------------------------------------------------------------------------
# Time series of observed and fitted sex ratio for each pop
#--------------------------------------------------------------------------------

mod_name <- "LCRchum_Ricker"
save_plot <- TRUE

dev.new(width=11,height=7)

## @knitr plot_sex_ratio_ts
life_cycle <- unlist(strsplit(mod_name, "_"))[1]
q_F <- do.call(extract1, list(as.name(mod_name), "q_F"))

cbind(fish_data_SMS, colQuantiles(q_F, probs = c(0.05, 0.5, 0.95))) %>%
  mutate(n_MF_obs = n_M_obs + n_F_obs) %>% 
  cbind(., with(., binconf(x = n_F_obs, n = n_MF_obs))) %>%
  ggplot(aes(x = year, y = PointEst, ymin = Lower, ymax = Upper)) + 
  geom_abline(intercept = 0.5, slope = 0, color = "gray") + 
  geom_ribbon(aes(ymin = `5%`, ymax = `95%`), fill = "slategray4", alpha = 0.5) +
  geom_line(aes(y = `50%`), col = "slategray4", lwd = 1) +
  geom_point(pch = 16, size = 2.5) + geom_errorbar(width = 0) +
  labs(x = "Year", y = "Proportion female") +
  facet_wrap(vars(pop), ncol = 4) + theme_bw(base_size = 16) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        strip.background = element_rect(fill = NA))
## @knitr

if(save_plot) 
  ggsave(filename=here("analysis", "results", paste0("q_F_fit_", mod_name, ".png")),
         width=11, height=7, units="in", dpi=300, type="cairo-png")
rm(list = c("mod_name","life_cycle","q_F"))


#--------------------------------------------------------------------------------
# Time series of observed and fitted p_HOS for each pop
#--------------------------------------------------------------------------------

mod_name <- "LCRchum_Ricker"
save_plot <- FALSE

dev.new(width=11,height=7)

## @knitr plot_p_HOS_ts
life_cycle <- unlist(strsplit(mod_name, "_"))[1]
p_HOS <- do.call(extract1, list(as.name(mod_name), "p_HOS"))

switch(life_cycle, SS = fish_data_SS, LCRchum = fish_data_SMS) %>% 
  mutate(zeros = 0, fit_p_HOS = as.logical(fit_p_HOS),
         p_HOS_obs = binconf(n_H_obs, n_H_obs + n_W_obs, alpha = 0.1)) %>% 
  do.call(data.frame, .) %>% # unpack col with nested data frame
  mutate(p_HOS_L = replace(zeros, fit_p_HOS, colQuantiles(p_HOS, probs = 0.05)),
         p_HOS_m = replace(zeros, fit_p_HOS, colMedians(p_HOS)),
         p_HOS_U = replace(zeros, fit_p_HOS, colQuantiles(p_HOS, probs = 0.95))) %>% 
  ggplot(aes(x = year)) +
  geom_ribbon(aes(ymin = p_HOS_L, ymax = p_HOS_U), fill = "slategray4", alpha = 0.5) +
  geom_line(aes(y = p_HOS_m), col = "slategray4", lwd = 1) +
  geom_point(aes(y = p_HOS_obs.PointEst), pch = 16, size = 2.5) +
  geom_errorbar(aes(ymin = p_HOS_obs.Lower, ymax = p_HOS_obs.Upper), width = 0) +
  labs(x = "Year", y = bquote(italic(p)[HOS])) +
  scale_x_continuous(minor_breaks = function(v) min(v):max(v)) +
  facet_wrap(vars(pop), ncol = 4) + theme_bw(base_size = 16) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        strip.background = element_rect(fill = NA))
## @knitr

if(save_plot) 
  ggsave(filename=here("analysis", "results", paste0("p_HOS_fit_", mod_name, ".png")),
      width=11, height=7, units="in", dpi=300, type="cairo-png")
rm(list = c("mod_name","life_cycle","p_HOS"))


#--------------------------------------------------------------------------------
# Lower Columbia chum spawner-egg-smolt-spawner #
# Observed and fitted distributions of fecundity by age
#--------------------------------------------------------------------------------

mod_name <- "LCRchum_Ricker"
save_plot <- FALSE

if(save_plot) {
png(filename=here("analysis","results",paste0("fecundity_fit_", mod_name, ".png")),
    width=7, height=7, units="in", res=200, type="cairo-png")
} else dev.new(width=7,height=7)

## @knitr plot_fecundity_fit
ages <- substring(names(select(fish_data_SMS, starts_with("n_age"))), 6, 6)
E_obs <- fecundity_data$E_obs
E_seq <- seq(min(E_obs, na.rm = TRUE), max(E_obs, na.rm = TRUE), length = 500)
mu_E <- do.call(extract1, list(as.name(mod_name), "mu_E"))
sigma_E <- do.call(extract1, list(as.name(mod_name), "sigma_E"))
E_fit <- array(NA, c(nrow(mu_E), length(E_seq), ncol(mu_E)))
for(a in 1:length(ages))
  E_fit[,,a] <- sapply(E_seq, function(x) dnorm(x, mu_E[,a], sigma_E[,a]))

c1 <- viridis(length(ages), end = 0.8, direction = -1) 
c1t <- transparent(c1, trans.val = 0.5)
c1tt <- transparent(c1, trans.val = 0.7)

par(mfrow = c(3,1), mar = c(3,2,0,2), oma = c(2,2,0,0))

for(a in 1:length(ages))
{
  hist(E_obs[fecundity_data$age_E == ages[a]], 20, prob = TRUE, 
       col = c1tt[a], border = "white", las = 1, cex.axis = 1.5, cex.lab = 1.8,
       xlim = range(E_seq), ylim = range(0, apply(E_fit, 2:3, quantile, 0.99)),
       xlab = "", ylab = "", main = "", xaxs = "i", yaxt = "n", bty = "n")
  lines(E_seq, colMedians(E_fit[,,a]), col = c1[a], lwd = 3)
  polygon(c(E_seq, rev(E_seq)),
          c(colQuantiles(E_fit[,,a], probs = 0.05), 
            rev(colQuantiles(E_fit[,,a], probs = 0.95))),
          col = c1t[a], border = NA)
  text(par("usr")[1] + 0.8*diff(par("usr")[1:2]), par("usr")[4]*0.5, 
       labels = paste("age", ages[a]), cex = 1.5, col = c1[a], adj = 1)
}
title(xlab = "Fecundity", ylab = "Probability density", cex.lab = 1.9, line = 0, outer = TRUE)

rm(list = c("mod_name","c1","c1t","c1tt","ages","E_obs","E_seq","mu_E","sigma_E","E_fit"))
## @knitr
if(save_plot) dev.off()


#--------------------------------------------------------------------------------
# Lower Columbia chum spawner-egg-smolt-spawner #
# Observed and fitted distributions of "known" smolt and spawner 
# observation error SDs
#--------------------------------------------------------------------------------

mod_name <- "LCRchum_Ricker"
save_plot <- FALSE

if(save_plot) {
png(filename=here("analysis","results",paste0("tau_fit_", mod_name, ".png")),
    width=6, height=8, units="in", res=200, type="cairo-png")
} else dev.new(width=6,height=8)

## @knitr plot_obs_error_fit
tau_M_obs <- fish_data_SMS$tau_M_obs
tau_M_seq <- seq(min(tau_M_obs, na.rm = TRUE), max(tau_M_obs, na.rm = TRUE), length = 500)
mu_tau_M <- do.call(extract1, list(as.name(mod_name), "mu_tau_M"))
sigma_tau_M <- do.call(extract1, list(as.name(mod_name), "sigma_tau_M"))
tau_M_fit <- sapply(tau_M_seq, function(x) dlnorm(x, log(mu_tau_M), sigma_tau_M))

tau_S_obs <- fish_data_SMS$tau_S_obs
tau_S_seq <- seq(min(tau_S_obs, na.rm = TRUE), max(tau_S_obs, na.rm = TRUE), length = 500)
mu_tau_S <- do.call(extract1, list(as.name(mod_name), "mu_tau_S"))
sigma_tau_S <- do.call(extract1, list(as.name(mod_name), "sigma_tau_S"))
tau_S_fit <- sapply(tau_S_seq, function(x) dlnorm(x, log(mu_tau_S), sigma_tau_S))

c1 <- "slategray4"
c1t <- transparent(c1, trans.val = 0.6)

par(mfcol = c(2,1), mar = c(5,5,0,1),  oma = c(0,0,0,0))

# smolt observation error SD
hist(tau_M_obs, 10, prob = TRUE, las = 1, cex.axis = 1.2, cex.lab = 1.5, 
     col = "lightgray", border = "white",
     ylim = c(0, max(colQuantiles(tau_M_fit, probs = 0.95))),
     xlab = bquote("Smolt observation error (" * tau[italic(M)] * ")"), 
     ylab = "Probability density", main = "")
polygon(c(tau_M_seq, rev(tau_M_seq)),
        c(colQuantiles(tau_M_fit, probs = 0.05), rev(colQuantiles(tau_M_fit, probs = 0.95))),
        col = c1t, border = NA)
lines(tau_M_seq, colMedians(tau_M_fit), col = c1, lwd = 3)

# spawner observation error SD
hist(tau_S_obs, 20, prob = TRUE, las = 1, cex.axis = 1.2, cex.lab = 1.5, 
     col = "lightgray", border = "white",
     ylim = c(0, max(colQuantiles(tau_S_fit, probs = 0.95))),
     xlab = bquote("Spawner observation error (" * tau[italic(S)] * ")"), 
     ylab = "Probability density", main = "")
polygon(c(tau_S_seq, rev(tau_S_seq)),
        c(colQuantiles(tau_S_fit, probs = 0.05), rev(colQuantiles(tau_S_fit, probs = 0.95))),
        col = c1t, border = NA)
lines(tau_S_seq, colMedians(tau_S_fit), col = c1, lwd = 3)

rm(list = c("mod_name","tau_M_obs","tau_M_seq","mu_tau_M","sigma_tau_M","tau_M_fit",
            "tau_S_obs","tau_S_seq","mu_tau_S","sigma_tau_S","tau_S_fit","c1","c1t"))
## @knitr
if(save_plot) dev.off()


#--------------------------------------------------------------------------------
# Tabular summary of 1-year-ahead forecasts by pop
#--------------------------------------------------------------------------------

mod_name <- "LCRchum_BH_fore"
life_stage <- "S"   # "S" = spawners, "M" = smolts

## @knitr forecast_df
life_cycle <- unlist(strsplit(mod_name, "_"))[1]
dat <- switch(life_cycle, SS = fish_data_SS_fore, SMS = fish_data_SMS_fore,
              LCRchum = fish_data_SMS_fore)
y <- dat$year
y1 <- min(dat$year[dat$forecast])
pop <- c(as.character(dat$pop[dat$year==y1]), "Total")
N_IPM <- do.call(extract1, list(as.name(mod_name), life_stage))
N1_IPM <- cbind(N_IPM[,dat$year==y1], rowSums(N_IPM[,dat$year==y1]))
N1_IPM_median <- round(colMedians(N1_IPM), 0)
N1_IPM_lo <- round(colQuantiles(N1_IPM, probs = 0.05), 0)
N1_IPM_up <- round(colQuantiles(N1_IPM, probs = 0.95), 0)

forecast_df <- data.frame(Population = pop, 
                          Estimate = paste0(N1_IPM_median, " (", N1_IPM_lo, ", ",
                                            N1_IPM_up, ")"))

rm("mod_name","life_stage","life_cycle","dat","y","y1","pop","N_IPM","N1_IPM",
   "N1_IPM_median","N1_IPM_lo","N1_IPM_up")
## @knitr


#--------------------------------------------------------------------------#
################# FIGURES FOR GENERIC SS OR SMS IPMs #######################
#--------------------------------------------------------------------------#

#--------------------------------------------------------------------
# Spawner-spawner / spawner-smolt-spawner #
# S-R curves and posterior distributions of parameters
#--------------------------------------------------------------------

mod_name <- "SMS_BH"

dev.new(width = 7, height = 7)
# png(filename=here("analysis","results",paste0("SR_",mod_name,".png")),
#     width=7, height=7, units="in", res=200, type="cairo-png")

## @knitr plot_SR_params
life_cycle <- unlist(strsplit(mod_name, "_"))[1]
dat <- switch(life_cycle, SS = fish_data_SS, SMS = fish_data_SMS)
SR_fun <- unlist(strsplit(mod_name, "_"))[2]

SR_eval <- function(alpha, Rmax = NULL, S, SR_fun) 
{
  switch(SR_fun,
         exp = alpha*S,
         BH = alpha*S/(1 + alpha*S/Rmax),
         Ricker = alpha*S*exp(-alpha*S/(exp(1)*Rmax)))
}

S <- colMedians(do.call(extract1, list(as.name(mod_name), "S")))
mu_alpha <- as.vector(do.call(extract1, list(as.name(mod_name), "mu_alpha")))
mu_Rmax <- as.vector(do.call(extract1, list(as.name(mod_name), "mu_Rmax")))
# S <- matrix(seq(0, quantile(dat$S_obs/dat$A, 0.9, na.rm = TRUE), length = 100),
#             nrow = length(mu_alpha), ncol = 100, byrow = TRUE)
scl <- ifelse(life_cycle == "SS", 1, 1e-6)
Smat <- matrix(seq(0, quantile(S/dat$A, 0.9, na.rm = TRUE), length = 100),
            nrow = length(mu_alpha), ncol = 100, byrow = TRUE)
R_ESU <- SR_eval(alpha = exp(mu_alpha), Rmax = exp(mu_Rmax), S = Smat, SR_fun = SR_fun)
alpha <- do.call(extract1, list(as.name(mod_name), "alpha"))
Rmax <- do.call(extract1, list(as.name(mod_name), "Rmax"))
R_pop <- sapply(1:ncol(alpha), function(i) {
  colMedians(SR_eval(alpha = alpha[,i], Rmax = Rmax[,i], S = Smat, SR_fun = SR_fun))
})

c1 <- "slategray4"
c1t <- transparent(c1, trans.val = 0.5)
c1tt <- transparent(c1, trans.val = 0.7)

par(mfrow = c(2,2), mar = c(5.1,5.1,1,1))

# Recruits vs. spawners
plot(Smat[1,], colMedians(R_ESU*scl), type = "l", lwd=3, col = c1, las = 1,
     cex.axis = 1.2, cex.lab = 1.5, xaxs = "i", yaxs = "i", #yaxt = "n", 
     ylim = range(R_pop*scl), xlab="Spawners", ylab = "")
for(i in 1:ncol(R_pop))
  lines(Smat[1,], R_pop[,i]*scl, col = c1t)
polygon(c(Smat[1,], rev(Smat[1,])), 
        c(colQuantiles(R_ESU*scl, probs = 0.05), 
          rev(colQuantiles(R_ESU*scl, probs = 0.95))), 
        col = c1tt, border = NA)
rug(S, col = c1)
title(ylab = ifelse(life_cycle=="SS", "Recruits", "Smolts (millions)"), 
      line = 3.5, cex.lab = 1.5)
text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "A", cex = 1.5)

# log(recruits/spawner) vs. spawners
plot(Smat[1,], colMedians(log(R_ESU/Smat)), type = "l", lwd=3, col = c1, las = 1,
     cex.axis = 1.2, cex.lab = 1.5, xaxs = "i", yaxs = "i", xlab="Spawners", 
     ylab = ifelse(life_cycle=="SS", "log(recruits/spawner)", "log(smolts/spawner)"),
     ylim = range(colQuantiles(log(R_ESU/Smat), probs = c(0.05,0.95)), na.rm = TRUE))
for(i in 1:ncol(R_pop))
  lines(Smat[1,], log(R_pop[,i]/Smat[1,]), col = c1t)
polygon(c(Smat[1,], rev(Smat[1,])),
        c(colQuantiles(log(R_ESU/Smat), probs = 0.05),
          rev(colQuantiles(log(R_ESU/Smat), probs = 0.95))),
        col = c1tt, border = NA)
rug(S, col = c1)
text(par("usr")[2], par("usr")[4], adj = c(1.5,1.5), "B", cex = 1.5)

# Posterior densities of log(alpha)
dd_ESU <- density(mu_alpha)
dd_pop <- vector("list", length(levels(dat$pop)))
for(i in 1:length(dd_pop))
  dd_pop[[i]] <- density(log(alpha[,i]))

plot(dd_ESU$x, dd_ESU$y, type = "l", lwd = 3, col = c1, las = 1, 
     cex.axis = 1.2, xlab = "", ylab = "", xaxs = "i",
     xlim = range(c(dd_ESU$x, sapply(dd_pop, function(m) m$x))),
     ylim = range(c(dd_ESU$y, sapply(dd_pop, function(m) m$y))))
for(i in 1:length(dd_pop))
  lines(dd_pop[[i]]$x, dd_pop[[i]]$y, col = c1t)
title(xlab = bquote(log(alpha)), ylab = "Probability density", line = 3.5, cex.lab = 1.5)
text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "C", cex = 1.5)

# Posterior densities of log(Rmax)
dd_ESU <- density(mu_Rmax)
dd_pop <- vector("list", length(levels(dat$pop)))
for(i in 1:length(dd_pop))
  dd_pop[[i]] <- density(log(Rmax[,i]))

plot(dd_ESU$x, dd_ESU$y, type = "l", lwd = 3, col = c1, las = 1, 
     cex.axis = 1.2, xlab = "", ylab = "", xaxs = "i",
     xlim = range(c(dd_ESU$x, sapply(dd_pop, function(m) m$x))),
     ylim = range(c(dd_ESU$y, sapply(dd_pop, function(m) m$y))))
for(i in 1:length(dd_pop))
  lines(dd_pop[[i]]$x, dd_pop[[i]]$y, col = c1t)
title(xlab = bquote(log(italic(R)[max])), ylab = "Probability density", 
      line = 3.5, cex.lab = 1.5)
text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "D", cex = 1.5)

rm(list=c("mod_name","life_cycle","SR_fun","scl","mu_alpha","mu_Rmax","S","R_ESU",
          "c1","c1t","c1tt","dd_ESU","dd_pop","alpha","Rmax","R_pop","Smat"))
## @knitr
# dev.off()


#--------------------------------------------------------------------------------
# Spawner-Spawner #
# Shared recruitment process errors (brood year productivity anomalies)
#--------------------------------------------------------------------------------

mod_name <- "SS_BH"

# dev.new(width=7,height=5)
png(filename=here("analysis","results",paste0("phi_", mod_name, ".png")),
    width=7, height=5, units="in", res=200, type="cairo-png")

y <- sort(unique(fish_data_SS$year))
phi <- do.call(extract1, list(as.name(mod_name), "phi"))

c1 <- "slategray4"
c1t <- transparent(c1, trans.val = 0.5)

par(oma = c(0,0.1,0,0))

plot(y, colMedians(phi), type = "n", las = 1, cex.axis = 1.2, cex.lab = 1.5,
     ylim = range(colQuantiles(phi, probs = c(0.05,0.95))), xaxs = "i", xaxt = "n",
     xlab = "Brood year", ylab = "Productivity anomaly")
abline(h = 0, col = "gray")
polygon(c(y, rev(y)), 
        c(colQuantiles(phi, probs = 0.05), rev(colQuantiles(phi, probs = 0.95))),
        col = c1t, border = NA)
lines(y, colMedians(phi), lwd = 3)
axis(side = 1, at = y[y %% 5 == 0], cex.axis = 1.2)
rug(y[y %% 5 != 0], ticksize = -0.02)

dev.off()
rm(list = c("mod_name","y","phi","c1","c1t"))


#--------------------------------------------------------------------------------
# Spawner-Smolt-Spawner #
# Shared recruitment and SAR process errors (brood year productivity anomalies)
#--------------------------------------------------------------------------------

mod_name <- "SMS_BH"

dev.new(width=6,height=8)
# png(filename=here("analysis","results",paste0("phi_", mod_name, ".png")),
#     width=6, height=8, units="in", res=200, type="cairo-png")

## @knitr plot_phi
y <- sort(unique(fish_data_SMS$year))

phi_M <- do.call(extract1, list(as.name(mod_name), "phi_M"))
sigma_M <- do.call(extract1, list(as.name(mod_name), "sigma_M"))
zeta_M <- do.call(stan_mean, list(as.name(mod_name), "zeta_M"))
epsilon_M <- t(sapply(sigma_M, function(x) x*zeta_M))

phi_MS <- do.call(extract1, list(as.name(mod_name), "phi_MS"))
mu_MS <- do.call(extract1, list(as.name(mod_name), "mu_MS"))
s_hat_MS <- plogis(sweep(phi_MS, 1, qlogis(mu_MS), "+"))
s_MS <- do.call(extract1, list(as.name(mod_name), "s_MS"))

c1 <- "slategray4"
c1t <- transparent(c1, trans.val = 0.5)
c1tt <- transparent(c1, trans.val = 0.7)

par(mfcol = c(2,1), mar = c(5.1,5.1,2,2),  oma = c(0,0.1,0,0))

# Smolt recruitment
plot(y, colMedians(phi_M), type = "n", las = 1, cex.axis = 1.2, cex.lab = 1.5,
     ylim = range(colQuantiles(phi_M, probs = c(0.05,0.95)),
                  colMedians(phi_M[,match(fish_data_SMS$year, y)] + epsilon_M)), 
     xlab = "Brood year", ylab = "Productivity anomaly", main = "Smolt recruitment",
     xaxs = "i", xaxt = "n")
abline(h = 0, col = "darkgray")
polygon(c(y, rev(y)), 
        c(colQuantiles(phi_M, probs = 0.05), rev(colQuantiles(phi_M, probs = 0.95))),
        col = c1tt, border = NA)
lines(y, colMedians(phi_M), col = c1t, lwd = 4)
for(j in levels(fish_data_SMS$pop))
{
  indx1 <- fish_data_SMS$pop == j
  indx2 <- y %in% fish_data_SMS$year[indx1]
  lines(y[indx2], colMedians(phi_M[,indx2] + epsilon_M[,indx1]), col = c1t)
}
axis(side = 1, at = y[y %% 5 == 0], cex.axis = 1.2)
rug(y[y %% 5 != 0], ticksize = -0.02)

# SAR
plot(y, colMedians(s_hat_MS), type = "n", las = 1, cex.axis = 1.2, cex.lab = 1.5,
     ylim = range(0, colQuantiles(s_hat_MS, probs = 0.95), colQuantiles(s_MS, probs = 0.95)), 
     xaxs = "i", xaxt = "n",
     xlab = "Outmigration year", ylab = "", main = "SAR")
mtext("Survival", side = 2, line = 3.7, cex = par("cex")*1.5)
polygon(c(y, rev(y)), 
        c(colQuantiles(s_hat_MS, probs = 0.05), rev(colQuantiles(s_hat_MS, probs = 0.95))),
        col = c1tt, border = NA)
lines(y, colMedians(s_hat_MS), col = c1t, lwd = 4)
for(j in levels(fish_data_SMS$pop))
  lines(fish_data_SMS$year[fish_data_SMS$pop==j], colMedians(s_MS[,fish_data_SMS$pop==j]), 
        col = c1t)
axis(side = 1, at = y[y %% 5 == 0], cex.axis = 1.2)
rug(y[y %% 5 != 0], ticksize = -0.02)

rm(list = c("mod_name","y","phi_M","phi_MS","s_hat_MS","c1","c1t",
            "sigma_M","zeta_M","epsilon_M"))
## @knitr
# dev.off()


#--------------------------------------------------------------------------------
# Spawner-Spawner #
# Posterior distributions of shared and unique recruitment process errors and 
# spawner observation error
#--------------------------------------------------------------------------------

mod_name <- "SS_Ricker"

# variance components
sigma <- do.call(extract1, list(as.name(mod_name), "sigma"))
sigma_phi <- do.call(extract1, list(as.name(mod_name), "sigma_phi"))
rho_phi <- do.call(extract1, list(as.name(mod_name), "rho_phi"))
sd_phi <- sqrt(sigma_phi^2 / sqrt(1 - rho_phi^2))
sigma_tot <- sqrt(sigma^2 + sd_phi^2)
tau <- do.call(extract1, list(as.name(mod_name), "tau"))

# densities
dd_sigma <- density(sigma)
dd_sd_phi <- density(sd_phi)
dd_sigma_tot <- density(sigma_tot)
dd_tau <- density(tau)

cols <- c("darkblue", "darkblue", "gray")

# dev.new(width=7,height=7)
png(filename=here("analysis","results",paste0("process_obs_SDs_", mod_name, ".png")),
    width=7, height=7, units="in", res=200, type="cairo-png")

par(oma = c(0,0.1,0,0))

plot(dd_sigma$x, dd_sigma$y, type = "l", col = cols[1], 
     las = 1, cex.axis = 1.2, cex.lab = 1.5, yaxs = "i",
     xlim = range(c(dd_sigma$x, dd_sd_phi$x, dd_sigma_tot$x, dd_tau$x)),
     ylim = range(c(dd_sigma$y, dd_sd_phi$y, dd_sigma_tot$y, dd_tau$y)),
     xlab = "Error SD", ylab = "Probability density")
# lines(dd_sd_phi$x, dd_sd_phi$y, col = cols[2])
lines(dd_sigma_tot$x, dd_sigma_tot$y, lwd = 3, col = cols[2])
lines(dd_tau$x, dd_tau$y, lwd = 3, col = cols[3])
legend("topright", lwd = c(1,3,3), col = cols, 
       legend = c(expression("unique " * italic(R) * " process error (" * sigma * ")"), 
                  #expression(sigma[phi]), 
                  expression("total " * italic(R) * " process error (" * sigma[tot] * ")"), 
                  expression(italic(S) * " observation error (" * tau * ")")))

dev.off()
rm(list = c("mod_name","sigma","sigma_phi","rho_phi","sd_phi","sigma_tot","tau",
            "dd_sigma","dd_sd_phi","dd_sigma_tot","dd_tau"))


