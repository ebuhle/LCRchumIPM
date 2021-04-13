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
source(here("analysis","R","01_LCRchumIPM_data.R"))
# load plotting functions
source(here("analysis","R","03_LCRchumIPM_plots.R"))
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
# FIGURES SPECIFIC TO LCR CHUM IPM
# Figures for generic IPM_SS_pp models are in 04_LCRchumIPM_extra-plots.R
#===========================================================================

#--------------------------------------------------------------------
# Life-cycle multiplot
# S-R curves (spawners to smolts)
# Posterior distributions of fecundity, survival and Mmax parameters
# Time series of smolt productivity process errors and SAR
#--------------------------------------------------------------------

mod_name <- "LCRchum_Ricker"
save_plot <- FALSE

## @knitr multiplot
LCRchumIPM_multiplot(mod = get(mod_name), SR_fun = strsplit(mod_name, "_")[[1]][2], 
                     fish_data = fish_data_SMS, save_plot = save_plot,
                     filename = here("analysis","results",paste0("multiplot_",mod_name,".png")))
## @knitr

#--------------------------------------------------------------------------------
# Time series of observed and fitted total spawners or smolts for each pop
#--------------------------------------------------------------------------------

mod_name <- "LCRchum_Ricker"
life_stage <- "S"   # "S" = spawners, "M" = smolts
save_plot <- TRUE

## @knitr plot_spawner_smolt_ts
LCRchumIPM_MS_timeseries(mod = get(mod_name), life_stage = life_stage, 
                         fish_data = fish_data_SMS, show_plot = TRUE, save_plot = save_plot, 
                         filename = here("analysis","results",paste0(life_stage, "_fit_", mod_name, ".png")))
## @knitr

#--------------------------------------------------------------------------------
# Time series of observed and fitted spawner age structure for each pop
#--------------------------------------------------------------------------------

mod_name <- "LCRchum_Ricker"
save_plot <- TRUE

## @knitr plot_spawner_age_ts
LCRchumIPM_age_timeseries(mod = get(mod_name), fish_data = fish_data_SMS, 
                          show_plot = TRUE, save_plot = save_plot, 
                          filename = here("analysis","results",paste0("q_fit_", mod_name, ".png")))
## @knitr

#--------------------------------------------------------------------------------
# Time series of observed and fitted sex ratio for each pop
#--------------------------------------------------------------------------------

mod_name <- "LCRchum_Ricker"
save_plot <- TRUE

## @knitr plot_sex_ratio_ts
LCRchumIPM_sex_timeseries(mod = get(mod_name), fish_data = fish_data_SMS, 
                          show_plot = TRUE, save_plot = save_plot, 
                          filename = here("analysis","results",paste0("q_F_fit_", mod_name, ".png")))
## @knitr

#--------------------------------------------------------------------------------
# Time series of observed and fitted p_HOS for each pop
#--------------------------------------------------------------------------------

mod_name <- "LCRchum_Ricker"
save_plot <- FALSE

## @knitr plot_p_HOS_ts
LCRchumIPM_p_HOS_timeseries(mod = get(mod_name), fish_data = fish_data_SMS, 
                            show_plot = TRUE, save_plot = save_plot, 
                            filename = here("analysis","results",paste0("p_HOS_fit_", mod_name, ".png")))
## @knitr


#--------------------------------------------------------------------------------
# Distributions of observed and fitted fecundity by age
#--------------------------------------------------------------------------------

mod_name <- "LCRchum_Ricker"
save_plot <- FALSE

## @knitr plot_fecundity_fit
LCRchumIPM_fecundity_plot(get(mod_name), fish_data = fish_data_SMS, fecundity_data = fecundity_data,
                          filename = here("analysis","results",paste0("fecundity_fit_", mod_name, ".png")),
                          save_plot = save_plot)
## @knitr

#--------------------------------------------------------------------------------
# Lower Columbia chum spawner-egg-smolt-spawner #
# Observed and fitted distributions of "known" smolt and spawner 
# observation error SDs
#--------------------------------------------------------------------------------

mod_name <- "LCRchum_Ricker"
save_plot <- TRUE

## @knitr plot_obs_error_fit
LCRchumIPM_obs_error_plot(get(mod_name), fish_data = fish_data_SMS,
                          filename = here("analysis","results",paste0("tau_fit_", mod_name, ".png")),
                          save_plot = save_plot)
## @knitr


#===========================================================================
# TABLES
#===========================================================================

#--------------------------------------------------------------------------------
# Summary of 1-year-ahead forecasts by pop
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

