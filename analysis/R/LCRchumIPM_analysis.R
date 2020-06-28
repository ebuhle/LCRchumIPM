## @knitr getting_started
options(device = ifelse(.Platform$OS.type == "windows", "windows", "quartz"))
options(mc.cores = parallel::detectCores(logical = FALSE) - 1)

library(salmonIPM)
library(rstan)
library(shinystan)
library(matrixStats)
library(tibble)
library(dplyr)
library(reshape2)
library(yarrr)
library(corrplot)
library(magicaxis)
library(zoo)
library(here)

if(file.exists(here("analysis","results","LCRchumIPM.RData")))
  load(here("analysis","results","LCRchumIPM.RData"))
## @knitr


#===========================================================================
# DATA
#===========================================================================

## @knitr data

# Mapping of location to population
location_pop <- read.csv(here("data","Location.Reach_Population.csv"), 
                         header = TRUE, stringsAsFactors = TRUE) %>% 
  rename(strata = Strata, location = Location.Reach, pop1 = Population1, pop2 = Population2)

# Mapping of disposition to hatchery vs. wild (i.e., broodstock vs. natural spawner)
disposition_HW <- read.csv(here("data","Disposition_HW.csv"), 
                           header = TRUE, stringsAsFactors = TRUE) %>% 
  rename(disposition = Disposition) %>% arrange(HW)

# Start dates of hatcheries associated with populations
hatcheries <- read.csv(here("data","Hatchery_Programs.csv"), header = TRUE, stringsAsFactors = TRUE)

# Spawner abundance data
# Assumptions:
# (1) NAs in hatchery dispositions (incl. Duncan Channel) are really zeros
# (2) NAs in Duncan Creek from 2004-present are really zeros
# (3) All other NAs are real missing observations
spawner_data <- read.csv(here("data","Data_ChumSpawnerAbundance_2019-12-12.csv"), 
                         header = TRUE, stringsAsFactors = TRUE) %>% 
  rename(year = Return.Yr., strata = Strata, location = Location.Reach, 
         disposition = Disposition, method = Method, S_obs = Abund.Mean, SD = Abund.SD) %>% 
  mutate(pop = location_pop$pop2[match(location, location_pop$location)],
         disposition_HW = disposition_HW$HW[match(disposition, disposition_HW$disposition)],
         S_obs = replace(S_obs, is.na(S_obs) & disposition_HW == "H", 0),
         S_obs = replace(S_obs, is.na(S_obs) & pop == "Duncan_Creek" & year >= 2004, 0)) %>% 
  select(year:location, pop, disposition, disposition_HW, method:SD) %>% 
  arrange(strata, location, year)

names_S_obs <- disposition_HW$disposition
names_B_take_obs <- disposition_HW$disposition[disposition_HW$HW == "H"]

spawner_data_agg <- aggregate(S_obs ~ year + strata + pop + disposition,
                              data = spawner_data, FUN = sum, na.action = na.pass) %>%
  dcast(year + strata + pop ~ disposition, value.var = "S_obs", 
        fun.aggregate = identity, fill = 0) %>% 
  add_column(S_obs = rowSums(select(., all_of(names_S_obs))),
             B_take_obs = rowSums(select(., all_of(names_B_take_obs)))) %>% 
  select(-all_of(names_S_obs)) %>% arrange(strata, pop, year)

# Spawner age-, sex-, and origin-frequency (aka BioData)
bio_data <- read.csv(here("data","Data_ChumSpawnerBioData_2019-12-12.csv"), 
                     header = TRUE, stringsAsFactors = TRUE) %>% 
  rename(year = Return.Yr., strata = Strata, location = Location.Reach, 
         disposition = Disposition, origin = Origin, sex = Sex, age = Age, count = Count) %>% 
  mutate(pop = location_pop$pop2[match(location, location_pop$location)],
         origin_HW = ifelse(origin == "Natural_spawner", "W", "H"),
         count = ifelse(is.na(count), 0, count)) %>% 
  select(year:location, pop, disposition, origin, origin_HW, sex:count) %>%
  arrange(strata, location, year, origin, age, sex)

# age of wild spawners only
bio_data_age <- bio_data %>% filter(origin_HW == "W") %>%  
  dcast(year + strata + pop ~ age, value.var = "count", fun.aggregate = sum)

bio_data_origin <- bio_data %>% 
  dcast(year + strata + pop ~ origin_HW, value.var = "count", fun.aggregate = sum)

# Juvenile abundance data
# Assumptions:
# (1) Smolts from Duncan Channel represent all naturally produced offspring of spawners
#     in Duncan Creek (hence Duncan_Channel -> Duncan_Creek in location_pop)
# (2) Duncan_North + Duncan_South = Duncan_Channel, so the former two are redundant 
#     (not really an assumption, although the equality isn't perfect in all years)
juv_data <- read.csv(here("data", "Data_ChumJuvenileAbundance_2020-06-09.csv"), 
                     header = TRUE, stringsAsFactors = TRUE) %>% 
  rename(brood_year = Brood.Year, year = Outmigration.Year, strata = Strata, 
         location = Location.Reach, origin = Origin, trap_type = TrapType, 
         analysis = Analysis, partial_spawners = Partial.Spawners, raw_catch = RawCatch,
         M_obs = Abund_Mean, SD = Abund_SD, L95 = Abund_L95, U95 = Abund_U95, CV = Abund_CV,
         comments = Comments) %>% 
  mutate(pop = location_pop$pop2[match(location, location_pop$location)]) %>% 
  select(brood_year:location, pop, origin:comments) %>% arrange(strata, location, year)

# drop hatchery or redundant pops and cases with leading or trailing NAs in M_obs
head_noNA <- function(x) { cumsum(!is.na(x)) > 0 }
juv_data_incl <- juv_data %>% filter(pop %in% spawner_data$pop) %>% 
  mutate(location = factor(location), pop = factor(pop, levels = levels(spawner_data$pop))) %>% 
  group_by(pop) %>% filter(head_noNA(M_obs) & rev(head_noNA(rev(M_obs)))) %>% as.data.frame()

# Fish data formatted for salmonIPM
# Drop age-2 and age-6 samples (each is < 0.1% of aged spawners)
# Use A = 1 for now (so Rmax in units of spawners)
fish_data <- full_join(spawner_data_agg, bio_data_age, by = c("year","strata","pop")) %>% 
  full_join(bio_data_origin, by = c("year","strata","pop")) %>% 
  full_join(juv_data_incl, by = c("year","strata","pop")) %>%
  mutate(B_take_obs = replace(B_take_obs, is.na(B_take_obs), 0)) %>% 
  rename_at(vars(contains("Age-")), list(~ paste0(sub("Age-","n_age",.), "_obs"))) %>% 
  select(-c(n_age2_obs, n_age6_obs)) %>% 
  rename(n_H_obs = H, n_W_obs = W) %>% mutate(A = 1, fit_p_HOS = NA, F_rate = 0) %>% 
  mutate_at(vars(contains("n_")), list(~ replace(., is.na(.), 0))) %>% 
  select(strata, pop, year, A, S_obs, M_obs, n_age3_obs:n_W_obs, 
         fit_p_HOS, B_take_obs, F_rate) %>% arrange(strata, pop, year) 

# fill in fit_p_HOS
for(i in 1:nrow(fish_data)) {
  pop_i <- as.character(fish_data$pop[i])
  start_year <- ifelse(pop_i %in% hatcheries$pop,
                       min(hatcheries$start_brood_year[hatcheries$pop == pop_i]) + 1,
                       NA)
  fish_data$fit_p_HOS[i] <- ifelse((!is.na(start_year) & fish_data$year[i] >= start_year) |
                                     fish_data$n_H_obs[i] > 0, 1, 0)
}

# # drop cases with initial NAs in S_obs unless bio data is present
# fish_data <- fish_data %>% mutate(n_age = rowSums(select(., n_age2_obs:n_age6_obs))) %>% 
#   group_by(pop) %>%  filter(head_noNA(S_obs) | cumsum(n_age) > 0) %>% 
#   select(-n_age) %>% as.data.frame()

# subsets for models with specific stage structure
# spawner-spawner: drop cases with initial NAs in S_obs, even if bio data is present
fish_data_SS <- fish_data %>% group_by(pop) %>% filter(head_noNA(S_obs)) %>% as.data.frame()
# spawner-spawner: drop cases with initial NAs in M_obs, even if bio data is present
fish_data_SMS <- fish_data %>% group_by(pop) %>% 
  filter(head_noNA(S_obs) | head_noNA(M_obs)) %>% as.data.frame()

## @knitr

#--------------------------------------------------------------
# Data exploration
#--------------------------------------------------------------

# Histogram of spawner age by sex and H/W origin
windows()
bio_data %>% mutate(age = substring(age,5,5)) %>% group_by(origin_HW, sex, age) %>% 
  summarize(n = sum(count)) %>% mutate(prop = n/sum(n)) %>% 
  ggplot(aes(x = age, y = prop)) + geom_bar(stat = "identity") + 
  facet_wrap(vars(origin_HW, sex), nrow = 2, ncol = 2) + theme_bw()


#===========================================================================
# FIT RETROSPECTIVE MODELS
#===========================================================================

#--------------------------------------------------------------
# Spawner-to-spawner IPM
#--------------------------------------------------------------

# Density-independent
## @knitr fit_SS_exp
SS_exp <- salmonIPM(fish_data = fish_data_SS, stan_model = "IPM_SS_pp", SR_fun = "exp",
                   pars = c("mu_alpha","sigma_alpha","alpha",
                            "sigma_phi","rho_phi","phi",
                            "mu_p","sigma_gamma","R_gamma","gamma",
                            "sigma_p","R_p","p",
                            "p_HOS","sigma","tau","S","R","q","LL"),
                   chains = 3, iter = 1500, warmup = 500,
                   control = list(adapt_delta = 0.99, max_treedepth = 13))

## @knitr print_SS_exp
print(SS_exp, prob = c(0.025,0.5,0.975),
      pars = c("alpha","phi","p_HOS","q","gamma","p","S","R","LL"), 
      include = FALSE, use_cache = FALSE)
## @knitr

launch_shinystan(SS_exp)

# Beverton-Holt
## @knitr fit_SS_BH
SS_BH <- salmonIPM(fish_data = fish_data_SS, stan_model = "IPM_SS_pp", SR_fun = "BH",
                    pars = c("mu_alpha","sigma_alpha","alpha",
                             "mu_Rmax","sigma_Rmax","Rmax","rho_alphaRmax",
                             "sigma_phi","rho_phi","phi",
                             "mu_p","sigma_gamma","R_gamma","gamma",
                             "sigma_p","R_p","p",
                             "p_HOS","sigma","tau","S","R","q","LL"),
                    chains = 3, iter = 1500, warmup = 500,
                    control = list(adapt_delta = 0.99, max_treedepth = 13))

## @knitr print_SS_BH
print(SS_BH, prob = c(0.025,0.5,0.975),
      pars = c("alpha","Rmax","phi","p_HOS","q","gamma","p","S","R","LL"), 
      include = FALSE, use_cache = FALSE)
## @knitr

launch_shinystan(SS_BH)

# Ricker
## @knitr fit_SS_Ricker
SS_Ricker <- salmonIPM(fish_data = fish_data_SS, stan_model = "IPM_SS_pp", SR_fun = "Ricker",
                       pars = c("mu_alpha","sigma_alpha","alpha",
                                "mu_Rmax","sigma_Rmax","Rmax","rho_alphaRmax",
                                "sigma_phi","rho_phi","phi",
                                "mu_p","sigma_gamma","R_gamma","gamma",
                                "sigma_p","R_p","p",
                                "p_HOS","sigma","tau","S","R","q","LL"),
                       chains = 3, iter = 1500, warmup = 500,
                       control = list(adapt_delta = 0.99, max_treedepth = 13))

## @knitr print_SS_Ricker
print(SS_Ricker, prob = c(0.025,0.5,0.975),
      pars = c("alpha","Rmax","phi","p_HOS","q","gamma","p","S","R","LL"), 
      include = FALSE, use_cache = FALSE)
## @knitr

launch_shinystan(SS_Ricker)

#--------------------------------------------------------------
# Model selection using LOO
#--------------------------------------------------------------

# Observationwise log-likelihood of each fitted model
# Here an observation is a row of fish_data, and the total likelihood includes 
# components for spawner abundance, age-frequency, and hatchery/wild-frequency
## @knitr loo_SS
LL_SS <- lapply(list(exp = SS_exp, BH = SS_BH, Ricker = SS_Ricker),
                loo::extract_log_lik, parameter_name = "LL", merge_chains = FALSE)

# Relative ESS of posterior draws of observationwise likelihood 
r_eff_SS <- lapply(LL_SS, function(x) relative_eff(exp(x)))

# PSIS-LOO
LOO_SS <- lapply(1:length(LL_SS), function(i) loo(LL_SS[[i]], r_eff = r_eff_SS[[i]]))
names(LOO_SS) <- names(LL_SS)

## Compare all three models
loo_compare(LOO_SS)

## Exponential vs. Ricker
loo_compare(LOO_SS[c("exp","Ricker")])

## Exponential vs. Beverton-Holt
loo_compare(LOO_SS[c("exp","BH")])

## Beverton-Holt vs. Ricker
loo_compare(LOO_SS[c("BH","Ricker")])
## @knitr


#--------------------------------------------------------------
# Spawner-smolt-spawner IPM
#--------------------------------------------------------------

# Density-independent
## @knitr fit_SMS_exp
SMS_exp <- salmonIPM(fish_data = fish_data_SMS, ages = list(M = 1),
                     stan_model = "IPM_SMS_pp", SR_fun = "exp",
                     pars = c("mu_alpha","sigma_alpha","alpha",
                              "mu_Rmax","sigma_Rmax","Rmax","rho_alphaRmax",
                              "beta_phi_M","rho_phi_M","sigma_phi_M","phi_M","sigma_M",
                              "mu_MS","beta_phi_MS","rho_phi_MS",
                              "sigma_phi_MS","phi_MS","sigma_MS","s_MS",
                              "mu_p","sigma_gamma","R_gamma","sigma_p","R_p","p",
                              "tau_M","tau_S","p_HOS","M","S","q","LL"),
                     chains = 3, iter = 1500, warmup = 500,
                     control = list(adapt_delta = 0.99, max_treedepth = 13))

## @knitr print_SMS_exp
print(SMS_exp, prob = c(0.025,0.5,0.975),
      pars = c("alpha","Rmax","phi_M","phi_MS","p","p_HOS","S","M","s_MS","q","LL"), 
      include = FALSE, use_cache = FALSE)
## @knitr

launch_shinystan(SMS_exp)

# Beverton-Holt
## @knitr fit_SMS_BH
SMS_BH <- salmonIPM(fish_data = fish_data_SMS, ages = list(M = 1),
                    stan_model = "IPM_SMS_pp", SR_fun = "BH",
                    pars = c("mu_alpha","sigma_alpha","alpha",
                             "mu_Rmax","sigma_Rmax","Rmax","rho_alphaRmax",
                             "beta_phi_M","rho_phi_M","sigma_phi_M","phi_M","sigma_M",
                             "mu_MS","beta_phi_MS","rho_phi_MS",
                             "sigma_phi_MS","phi_MS","sigma_MS","s_MS",
                             "mu_p","sigma_gamma","R_gamma","sigma_p","R_p","p",
                             "tau_M","tau_S","p_HOS","M","S","q","LL"),
                    chains = 3, iter = 1500, warmup = 500,
                    control = list(adapt_delta = 0.99, max_treedepth = 13))

## @knitr print_SMS_BH
print(SMS_BH, prob = c(0.025,0.5,0.975),
      pars = c("alpha","Rmax","phi_M","phi_MS","p","p_HOS","S","M","s_MS","q","LL"), 
      include = FALSE, use_cache = FALSE)
## @knitr

launch_shinystan(SMS_BH)

# Ricker
## @knitr fit_SMS_Ricker
SMS_Ricker <- salmonIPM(fish_data = fish_data_SMS, ages = list(M = 1),
                        stan_model = "IPM_SMS_pp", SR_fun = "Ricker",
                        pars = c("mu_alpha","sigma_alpha","alpha",
                                 "mu_Rmax","sigma_Rmax","Rmax","rho_alphaRmax",
                                 "beta_phi_M","rho_phi_M","sigma_phi_M","phi_M","sigma_M",
                                 "mu_MS","beta_phi_MS","rho_phi_MS",
                                 "sigma_phi_MS","phi_MS","sigma_MS","s_MS",
                                 "mu_p","sigma_gamma","R_gamma","sigma_p","R_p","p",
                                 "tau_M","tau_S","p_HOS","M","S","q","LL"),
                        chains = 3, iter = 1500, warmup = 500,
                        control = list(adapt_delta = 0.99, max_treedepth = 13))

## @knitr print_SMS_Ricker
print(SMS_Ricker, prob = c(0.025,0.5,0.975),
      pars = c("alpha","Rmax","phi_M","phi_MS","p","p_HOS","S","M","s_MS","q","LL"), 
      include = FALSE, use_cache = FALSE)
## @knitr

launch_shinystan(SMS_Ricker)


#--------------------------------------------------------------
# Model selection using LOO
#--------------------------------------------------------------

# Observationwise log-likelihood of each fitted model
# Here an observation is a row of fish_data, and the total likelihood includes 
# components for smolt abundance, spawner abundance, age-frequency, and H/W-frequency
## @knitr loo_SMS
LL_SMS <- lapply(list(exp = SMS_exp, BH = SMS_BH, Ricker = SMS_Ricker),
             loo::extract_log_lik, parameter_name = "LL", merge_chains = FALSE)

# Relative ESS of posterior draws of observationwise likelihood 
r_eff_SMS <- lapply(LL_SMS, function(x) relative_eff(exp(x)))

# PSIS-LOO
LOO_SMS <- lapply(1:length(LL_SMS), function(i) loo(LL_SMS[[i]], r_eff = r_eff_SMS[[i]]))
names(LOO_SMS) <- names(LL_SMS)

## Compare all three models
loo_compare(LOO_SMS)

## Exponential vs. Ricker
loo_compare(LOO_SMS[c("exp","Ricker")])

## Exponential vs. Beverton-Holt
loo_compare(LOO_SMS[c("exp","BH")])

## Beverton-Holt vs. Ricker
loo_compare(LOO_SMS[c("BH","Ricker")])
## @knitr


#--------------------------------------------------------------
# Save stanfit objects
#--------------------------------------------------------------

save(list = ls()[sapply(ls(), function(x) do.call(class, list(as.name(x)))) == "stanfit"], 
     file = here("analysis","results","LCRchumIPM.RData"))



#===========================================================================
# FIGURES
#===========================================================================

#--------------------------------------------------------------------
# S-R curves and posterior distributions of parameters
#--------------------------------------------------------------------

mod_name <- "SMS_Ricker"
# dev.new(width = 7, height = 7)
png(filename=here("analysis","results",paste0("SR_",mod_name,".png")),
    width=7, height=7, units="in", res=200, type="cairo-png")

## @knitr plot_SR_params
life_cycle <- unlist(strsplit(mod_name, "_"))[1]
dat <- switch(life_cycle, SS = fish_data_SS, SMS = fish_data_SMS)
SR_fun <- unlist(strsplit(mod_name, "_"))[2]

SR_eval <- function(alpha, Rmax = NULL, S, SR_fun = "BH") 
{
  switch(SR_fun,
         exp = alpha*S,
         BH = alpha*S/(1 + alpha*S/Rmax),
         Ricker = alpha*S*exp(-alpha*S/(exp(1)*Rmax)))
}

S_IPM <- colMedians(do.call(extract1, list(as.name(mod_name), "S")))
mu_alpha <- as.vector(do.call(extract1, list(as.name(mod_name), "mu_alpha")))
mu_Rmax <- as.vector(do.call(extract1, list(as.name(mod_name), "mu_Rmax")))
# S <- matrix(seq(0, quantile(dat$S_obs/dat$A, 0.9, na.rm = TRUE), length = 100),
#             nrow = length(mu_alpha), ncol = 100, byrow = TRUE)
S <- matrix(seq(0, quantile(S_IPM/dat$A, 0.9, na.rm = TRUE), length = 100),
            nrow = length(mu_alpha), ncol = 100, byrow = TRUE)
R_ESU_IPM <- SR_eval(alpha = exp(mu_alpha), Rmax = exp(mu_Rmax), S = S, SR_fun = SR_fun)
alpha <- do.call(extract1, list(as.name(mod_name), "alpha"))
Rmax <- do.call(extract1, list(as.name(mod_name), "Rmax"))
R_pop_IPM <- sapply(1:ncol(alpha), function(i) {
  colMedians(SR_eval(alpha = alpha[,i], Rmax = Rmax[,i], S = S, SR_fun = SR_fun))
  })

c1 <- "slategray4"
c1t <- transparent(c1, trans.val = 0.5)
c1tt <- transparent(c1, trans.val = 0.7)

par(mfrow = c(2,2), mar = c(5.1,5.1,1,1))

# Recruits vs. spawners
plot(S[1,], colMedians(R_ESU_IPM), type = "l", lwd=3, col = c1, las = 1,
     cex.axis = 1.2, cex.lab = 1.8, xaxs = "i", yaxs = "i", yaxt = "n", 
     ylim = range(R_pop_IPM), xlab="Spawners", ylab = "")
tck <- axTicks(side = 2)
axis(side = 2, at = tck, labels = tck * switch(life_cycle, SS = 1, SMS = 1/1000), 
     las = 1, cex.axis = 1.2)
for(i in 1:ncol(R_pop_IPM))
  lines(S[1,], R_pop_IPM[,i], col = c1t)
polygon(c(S[1,], rev(S[1,])), 
        c(colQuantiles(R_ESU_IPM, probs = 0.025), rev(colQuantiles(R_ESU_IPM, probs = 0.975))), 
        col = c1tt, border = NA)
rug(S_IPM, col = c1)
title(ylab = switch(life_cycle, SS = "Recruits", SMS = "Smolts (thousands)"), 
      line = 3.5, cex.lab = 1.8)
# text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "A", cex = 2)

# log(recruits/spawner) vs. spawners
plot(S[1,], colMedians(log(R_ESU_IPM/S)), type = "l", lwd=3, col = c1, las = 1,
     cex.axis = 1.2, cex.lab = 1.8, xaxs = "i", yaxs = "i", xlab="Spawners", 
     ylab = switch(life_cycle, SS = "log(recruits/spawner)", SMS = "log(smolts/spawner)"),
     ylim = range(colQuantiles(log(R_ESU_IPM/S), probs = c(0.025,0.975)), na.rm = TRUE))
for(i in 1:ncol(R_pop_IPM))
  lines(S[1,], log(R_pop_IPM[,i]/S[1,]), col = c1t)
polygon(c(S[1,], rev(S[1,])),
        c(colQuantiles(log(R_ESU_IPM/S), probs = 0.025),
          rev(colQuantiles(log(R_ESU_IPM/S), probs = 0.975))),
        col = c1tt, border = NA)
rug(S_IPM, col = c1)
# text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "B", cex = 2)

# Posterior densities of log(alpha)
dd_IPM_ESU <- density(mu_alpha)
dd_IPM_pop <- vector("list", length(levels(dat$pop)))
for(i in 1:length(dd_IPM_pop))
  dd_IPM_pop[[i]] <- density(log(alpha[,i]))

plot(dd_IPM_ESU$x, dd_IPM_ESU$y, type = "l", lwd = 3, col = c1, las = 1, 
     cex.axis = 1.2, xlab = "", ylab = "", xaxs = "i",
     xlim = range(c(dd_IPM_ESU$x, sapply(dd_IPM_pop, function(m) m$x))),
     ylim = range(c(dd_IPM_ESU$y, sapply(dd_IPM_pop, function(m) m$y))))
for(i in 1:length(dd_IPM_pop))
  lines(dd_IPM_pop[[i]]$x, dd_IPM_pop[[i]]$y, col = c1t)
title(xlab = bquote(log(alpha)), ylab = "Probability density", line = 3.5, cex.lab = 1.8)
# text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "C", cex = 2)

# Posterior densities of log(Rmax)
dd_IPM_ESU <- density(mu_Rmax)
dd_IPM_pop <- vector("list", length(levels(dat$pop)))
for(i in 1:length(dd_IPM_pop))
  dd_IPM_pop[[i]] <- density(log(Rmax[,i]))

plot(dd_IPM_ESU$x, dd_IPM_ESU$y, type = "l", lwd = 3, col = c1, las = 1, 
     cex.axis = 1.2, xlab = "", ylab = "", xaxs = "i",
     xlim = range(c(dd_IPM_ESU$x, sapply(dd_IPM_pop, function(m) m$x))),
     ylim = range(c(dd_IPM_ESU$y, sapply(dd_IPM_pop, function(m) m$y))))
for(i in 1:length(dd_IPM_pop))
  lines(dd_IPM_pop[[i]]$x, dd_IPM_pop[[i]]$y, col = c1t)
title(xlab = bquote(log(italic(R)[max])), ylab = "Probability density", 
      line = 3.5, cex.lab = 1.8)
# text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "D", cex = 2)

rm(list=c("mod_name","life_cycle","SR_fun","mu_alpha","mu_Rmax","S","R_ESU_IPM","tck",
          "c1","c1t","c1tt","dd_IPM_ESU","dd_IPM_pop","alpha","Rmax","R_pop_IPM","S_IPM"))
## @knitr
dev.off()


#--------------------------------------------------------------------------------
# Time series of observed and fitted total spawners or smolts for each pop
#--------------------------------------------------------------------------------

mod_name <- "SMS_Ricker"
life_stage <- "S"   # "S" = spawners, "M" = smolts

# dev.new(width=13,height=8)
png(filename=here("analysis", "results", paste0(life_stage, "_fit_", mod_name, ".png")),
    width=13*0.9, height=8*0.9, units="in", res=200, type="cairo-png")

## @knitr plot_spawner_smolt_ts
life_cycle <- unlist(strsplit(mod_name, "_"))[1]
dat <- switch(life_cycle, SS = fish_data_SS, SMS = fish_data_SMS)
SR_fun <- unlist(strsplit(mod_name, "_"))[2]

par(mfrow=c(3,4), mar=c(1,3,4.1,1), oma=c(4.1,3.1,0,0))

N_obs <- dat[,paste0(life_stage, "_obs")]
N_IPM <- do.call(extract1, list(as.name(mod_name), life_stage))
tau <- do.call(extract1, list(as.name(mod_name), 
                    switch(life_cycle, SS = "tau",
                           SMS = switch(life_stage, M = "tau_M", S = "tau_S"))))
N_obs_IPM <- N_IPM * rlnorm(length(N_IPM), 0, tau)

c1 <- "slategray4"
c1t <- transparent(c1, trans.val = 0.5)
c1tt <- transparent(c1, trans.val = 0.7)

for(i in levels(dat$pop))
{
  yi <- dat$year[dat$pop==i]
  plot(yi, N_obs[dat$pop==i], pch = "",
       xlim = range(dat$year),
       ylim = range(pmax(N_obs[dat$pop==i], 1),
                    colQuantiles(N_obs_IPM[,dat$pop==i], probs = c(0.025,0.975)), na.rm = T), 
       las = 1, xaxt = "n", yaxs = "i", yaxt = "n", xlab = "", ylab = "", log = "y")
  axis(side = 1, at = yi[yi %% 5 == 0], cex.axis = 1.2)
  rug(yi[yi %% 5 != 0], ticksize = -0.02)
  # tck <- maglab(10^par("usr")[3:4], log = TRUE)
  # axis(2, at = tck$labat, labels = tck$exp, cex.axis = 1.2, las = 1)
  magaxis(2, minorn = 0, crunch = FALSE, tcl = -0.4, las = 1, cex.axis = 1.2)
  mtext(i, side = 3, line = 0.5, cex = par("cex")*1.5)
  if(par("mfg")[2] == 1) 
    mtext(switch(life_stage, M = "Smolts", S = "Spawners"), 
          side = 2, line = 3.5, cex = par("cex")*1.5)
  if(par("mfg")[1] == par("mfg")[3]) 
    mtext("Year", side = 1, line = 3, cex = par("cex")*1.5)
  lines(yi, colMedians(N_IPM[,dat$pop==i]), col = c1, lwd = 3)
  polygon(c(yi, rev(yi)), 
          c(colQuantiles(N_IPM[,dat$pop==i], probs = 0.025), 
            rev(colQuantiles(N_IPM[,dat$pop==i], probs = 0.975))),
          col = c1t, border = NA)
  polygon(c(yi, rev(yi)), 
          c(colQuantiles(N_obs_IPM[,dat$pop==i], probs = 0.025), 
            rev(colQuantiles(N_obs_IPM[,dat$pop==i], probs = 0.975))),
          col = c1tt, border = NA)
  points(yi, N_obs[dat$pop==i], pch=16, cex = 1.5)
}

rm(list = c("mod_name","life_stage","life_cycle","dat","SR_fun",
            "N_IPM","N_obs_IPM","N_obs","c1","c1t","c1tt","yi","tau"))
## @knitr
dev.off()


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
     ylim = range(colQuantiles(phi, probs = c(0.025,0.975))), xaxs = "i", xaxt = "n",
     xlab = "Brood year", ylab = "Productivity anomaly")
abline(h = 0, col = "gray")
polygon(c(y, rev(y)), 
        c(colQuantiles(phi, probs = 0.025), rev(colQuantiles(phi, probs = 0.975))),
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

# dev.new(width=6,height=8)
png(filename=here("analysis","results",paste0("phi_", mod_name, ".png")),
    width=6, height=8, units="in", res=200, type="cairo-png")

## @knitr plot_phi
y <- sort(unique(fish_data_SMS$year))
phi_M <- do.call(extract1, list(as.name(mod_name), "phi_M"))
phi_MS <- do.call(extract1, list(as.name(mod_name), "phi_MS"))

c1 <- "slategray4"
c1t <- transparent(c1, trans.val = 0.5)

par(mfcol = c(2,1), mar = c(5.1,4.1,2,2),  oma = c(0,0.1,0,0))

# Smolt recruitment
plot(y, colMedians(phi_M), type = "n", las = 1, cex.axis = 1.2, cex.lab = 1.5,
     ylim = range(colQuantiles(phi_M, probs = c(0.025,0.975))), xaxs = "i", xaxt = "n",
     xlab = "Brood year", ylab = "Productivity anomaly", main = "Smolt recruitment")
abline(h = 0, col = "gray")
polygon(c(y, rev(y)), 
        c(colQuantiles(phi_M, probs = 0.025), rev(colQuantiles(phi_M, probs = 0.975))),
        col = c1t, border = NA)
lines(y, colMedians(phi_M), lwd = 3)
axis(side = 1, at = y[y %% 5 == 0], cex.axis = 1.2)
rug(y[y %% 5 != 0], ticksize = -0.02)

# SAR
plot(y, colMedians(phi_MS), type = "n", las = 1, cex.axis = 1.2, cex.lab = 1.5,
     ylim = range(colQuantiles(phi_MS, probs = c(0.025,0.975))), xaxs = "i", xaxt = "n",
     xlab = "Outmigration year", ylab = "Productivity anomaly", main = "SAR")
abline(h = 0, col = "gray")
polygon(c(y, rev(y)), 
        c(colQuantiles(phi_MS, probs = 0.025), rev(colQuantiles(phi_MS, probs = 0.975))),
        col = c1t, border = NA)
lines(y, colMedians(phi_MS), lwd = 3)
axis(side = 1, at = y[y %% 5 == 0], cex.axis = 1.2)
rug(y[y %% 5 != 0], ticksize = -0.02)

rm(list = c("mod_name","y","phi_M","phi_MS","c1","c1t"))
## @knitr
dev.off()


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



