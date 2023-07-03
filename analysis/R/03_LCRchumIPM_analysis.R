#===========================================================================
# SETUP
#===========================================================================

## @knitr getting_started
if(.Platform$OS.type == "windows") options(device = "windows")
options(mc.cores = parallel::detectCores(logical = FALSE))

library(Hmisc)
library(dplyr)
library(tidyr)
library(salmonIPM)
library(rstan)
library(shinystan)
library(posterior)
library(distributional)
library(matrixStats)
library(yarrr)
library(zoo)
library(ggplot2)
library(ggdist)
theme_set(theme_bw(base_size = 16))
library(scales)
library(magicaxis)
library(viridis)
library(here)

# load data
source(here("analysis","R","01_LCRchumIPM_data.R"))
# load plotting functions
source(here("analysis","R","02_LCRchumIPM_plots.R"))
# load saved stanfit objects
if(file.exists(here("analysis","results","LCRchumIPM.RData")))
  load(here("analysis","results","LCRchumIPM.RData"))
## @knitr


#===========================================================================
# FIT RETROSPECTIVE MODELS
#===========================================================================

# # Density-independent
# ## @knitr fit_exp
# fit_exp <- salmonIPM(fish_data = fish_data,  fecundity_data = fecundity_data,
#                      ages = list(M = 1), stan_model = "IPM_LCRchum_pp", SR_fun = "exp",
#                      pars = c("mu_Emax","sigma_Emax","Emax"),
#                      include = FALSE, log_lik = TRUE,
#                      chains = 4, iter = 1500, warmup = 500,
#                      control = list(adapt_delta = 0.99, max_treedepth = 14))
# 
# ## @knitr print_fit_exp
# print(fit_exp, prob = c(0.05,0.5,0.95),
#       pars = c("eta_pop_EM","eta_year_EM","eta_year_MS","eta_pop_p","p",
#                "tau_M","tau_S","p_HOS","B_rate","E","S","M","s_EM","s_MS","q","LL"),
#       include = FALSE, use_cache = FALSE)
# ## @knitr

# # Beverton-Holt
# ## @knitr fit_BH
# fit_BH <- salmonIPM(stan_model = "IPM_LCRchum_pp", SR_fun = "BH", 
#                     par_models = list(s_MS ~ pop_type), 
#                     center = FALSE, scale = FALSE, ages = list(M = 1), 
#                     fish_data = fish_data, fecundity_data = fecundity_data,
#                     log_lik = TRUE, chains = 4, iter = 2000, warmup = 1000,
#                     control = list(adapt_delta = 0.95, max_treedepth = 14))
# 
# ## @knitr print_fit_BH
# print(fit_BH, prob = c(0.05,0.5,0.95),
#       pars = c("psi","Mmax","eta_year_M","eta_year_MS","eta_pop_p","mu_pop_alr_p","p","p_F",
#                "tau_M","tau_S","B_rate","E_hat","M","S","s_MS","q","q_F","q_origin","p_HOS","LL"), 
#       include = FALSE, use_cache = FALSE)
# ## @knitr

# Ricker
## @knitr fit_Ricker
fit_Ricker <- salmonIPM(stan_model = "IPM_LCRchum_pp", SR_fun = "Ricker", 
                        par_models = list(s_MS ~ pop_type), 
                        center = FALSE, scale = FALSE, ages = list(M = 1), 
                        fish_data = fish_data, fecundity_data = fecundity_data,
                        log_lik = TRUE, chains = 4, iter = 2000, warmup = 1000,
                        control = list(adapt_delta = 0.95, max_treedepth = 15))

## @knitr print_fit_Ricker
print(fit_Ricker, prob = c(0.05,0.5,0.95),
      pars = c("psi","Mmax","eta_year_M","eta_year_MS","eta_pop_p","mu_pop_alr_p","p","p_F",
               "tau_M","tau_S","B_rate","E_hat","M","S","s_MS","q","q_F","q_origin","p_HOS","LL"),
      include = FALSE, use_cache = FALSE)
## @knitr


# #--------------------------------------------------------------
# # Model selection using LOO
# #--------------------------------------------------------------
# 
# # Observationwise log-likelihood of each fitted model
# # Here an observation is a row of fish_data, and the total likelihood includes 
# # components for smolt abundance, spawner abundance, age-frequency, H/W-frequency,
# # and smolt and spawner observation error SDs
# ## @knitr loo
# LL <- lapply(list(exp = fit_exp, BH = fit_BH, Ricker = fit_Ricker),
#              loo::extract_log_lik, parameter_name = "LL", merge_chains = FALSE)
# 
# # Relative ESS of posterior draws of observationwise likelihood 
# r_eff <- lapply(LL_LCRchum, function(x) relative_eff(exp(x)))
# 
# # PSIS-LOO
# LOO <- lapply(1:length(LL), 
#               function(i) loo(LL[[i]], r_eff = r_eff[[i]]))
# names(LOO) <- names(LL)
# 
# ## Compare all three models
# loo_compare(LOO)
# 
# ## Exponential vs. Ricker
# loo_compare(LOO[c("exp","Ricker")])
# 
# ## Exponential vs. Beverton-Holt
# loo_compare(LOO[c("exp","BH")])
# 
# ## Beverton-Holt vs. Ricker
# loo_compare(LOO[c("BH","Ricker")])
# ## @knitr


#===========================================================================
# FIT PROSPECTIVE FORECASTING MODELS
#===========================================================================

# # Beverton_Holt
# # no broodstock removals or hatchery smolt releases
# # @knitr foreH0_BH
# foreH0_BH <- salmonIPM(stan_model = "IPM_LCRchum_pp", SR_fun = "BH", 
#                        par_models = list(s_MS ~ pop_type), 
#                        center = FALSE, scale = FALSE, ages = list(M = 1), 
#                        fish_data = fish_data_foreH0, 
#                        fecundity_data = fecundity_data,
#                        chains = 4, iter = 2000, warmup = 1000,
#                        control = list(adapt_delta = 0.95, max_treedepth = 14))
# 
# ## @knitr print_foreH0_BH
# print(foreH0_BH, prob = c(0.05,0.5,0.95),
#       pars = c("psi","Mmax","eta_year_M","eta_year_MS","eta_pop_p","mu_pop_alr_p","p","p_F",
#                "tau_M","tau_S","B_rate","E_hat","M","S","s_MS","q","q_F","q_origin","p_HOS","LL"), 
#       include = FALSE, use_cache = FALSE)
# ## @knitr
# 
# # Beverton_Holt
# # broodstock removal rates and hatchery smolt releases at maximum observed
# # @knitr foreHmax_BH
# foreHmax_BH <- salmonIPM(stan_model = "IPM_LCRchum_pp", SR_fun = "BH", 
#                          par_models = list(s_MS ~ pop_type), 
#                          center = FALSE, scale = FALSE, ages = list(M = 1), 
#                          fish_data = fish_data_foreHmax, 
#                          fecundity_data = fecundity_data,
#                          chains = 4, iter = 2000, warmup = 1000,
#                          control = list(adapt_delta = 0.95, max_treedepth = 14))
# 
# ## @knitr print_foreHmax_BH
# print(foreHmax_BH, prob = c(0.05,0.5,0.95),
#       pars = c("psi","Mmax","eta_year_M","eta_year_MS","eta_pop_p","mu_pop_alr_p","p","p_F",
#                "tau_M","tau_S","B_rate","E_hat","M","S","s_MS","q","q_F","q_origin","p_HOS","LL"), 
#       include = FALSE, use_cache = FALSE)
# ## @knitr

# Ricker
# no broodstock removals or hatchery smolt releases
# @knitr fit_foreH0_Ricker
foreH0_Ricker <- salmonIPM(stan_model = "IPM_LCRchum_pp", SR_fun = "Ricker", 
                           par_models = list(s_MS ~ pop_type), 
                           center = FALSE, scale = FALSE, ages = list(M = 1), 
                           fish_data = fish_data_foreH0, 
                           fecundity_data = fecundity_data,
                           chains = 4, iter = 2000, warmup = 1000,
                           control = list(adapt_delta = 0.95, max_treedepth = 15))

## @knitr print_foreH0_Ricker
print(foreH0_Ricker, prob = c(0.05,0.5,0.95),
      pars = c("psi","Mmax","eta_year_M","eta_year_MS","eta_pop_p","mu_pop_alr_p","p","p_F",
               "tau_M","tau_S","B_rate","E_hat","M","S","s_MS","q","q_F","q_origin","p_HOS","LL"), 
      include = FALSE, use_cache = FALSE)
## @knitr

# Ricker
# broodstock removal rates and hatchery smolt releases at maximum observed
## @knitr fit_foreHmax_Ricker
foreHmax_Ricker <- salmonIPM(stan_model = "IPM_LCRchum_pp", SR_fun = "Ricker", 
                             par_models = list(s_MS ~ pop_type), 
                             center = FALSE, scale = FALSE, ages = list(M = 1), 
                             fish_data = fish_data_foreHmax, 
                             fecundity_data = fecundity_data,
                             chains = 4, iter = 2000, warmup = 1000,
                             control = list(adapt_delta = 0.95, max_treedepth = 15))

## @knitr print_foreHmax_Ricker
print(foreHmax_Ricker, prob = c(0.05,0.5,0.95),
      pars = c("psi","Mmax","eta_year_M","eta_year_MS","eta_pop_p","mu_pop_alr_p","p","p_F",
               "tau_M","tau_S","B_rate_all","E_hat","M","S","s_MS","q","q_F","q_origin","p_HOS"), 
      include = FALSE, use_cache = FALSE)
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
# Life-cycle multiplot
# Posterior distributions of fecundity, survival and Mmax parameters
# Reconstructed S-R curve (spawners to smolts)
# Time series of smolt productivity process errors and SAR
#--------------------------------------------------------------------

mod_name <- "fit_Ricker"
save_plot <- TRUE

if(save_plot) {
  png(filename = here("analysis","results", 
                      paste0("multiplot_", strsplit(mod_name, "_")[[1]][2], ".png")), 
      width=12, height=5.5, units="in", res=300, type = "cairo-png")
} else dev.new(width = 12, height = 5.5)

## @knitr multiplot
multiplot(mod = get(mod_name), SR_fun = strsplit(mod_name, "_")[[1]][2], 
          fish_data = fish_data)
## @knitr
if(save_plot) dev.off()

#--------------------------------------------------------------------------------
# Distributions of observed and fitted fecundity by age
#--------------------------------------------------------------------------------

mod_name <- "fit_Ricker"
save_plot <- TRUE

## @knitr plot_fecundity_fit
gg <- fecundity_plot(get(mod_name), fish_data = fish_data, fecundity_data = fecundity_data)
## @knitr

if(save_plot) {
  ggsave(filename = here("analysis","results",
                         paste0("fecundity_fit_", strsplit(mod_name, "_")[[1]][2], ".png")), 
         width=7, height=7, units="in", dpi=200, type = "cairo-png")
} else {
  dev.new(width=7,height=7)
  plot(gg)
}

#--------------------------------------------------------------------------------
# S-R parameters at population and ESU level 
#--------------------------------------------------------------------------------

mod_name <- "fit_Ricker"
save_plot <- TRUE

## @knitr plot_psi_Mmax
gg <- psi_Mmax_plot(mod = get(mod_name), fish_data)
## @knitr

if(save_plot) {
  ggsave(filename = here("analysis","results",
                         paste0("psi_Mmax_", strsplit(mod_name, "_")[[1]][2], ".png")), 
         width=8, height=7, units="in", dpi=300, type = "cairo-png")
} else {
  dev.new(width=8,height=7)
  plot(gg)
}

# #--------------------------------------------------------------------
# # Spawner-recruit curves for each pop with data and states
# #--------------------------------------------------------------------
# 
# mod_name <- "fit_Ricker"
# life_stage <- "M"   # "M" = smolts, "R" = adult recruits
# save_plot <- FALSE
# 
# ## @knitr SR_plot
# gg <- SR_plot(mod = get(mod_name), SR_fun = strsplit(mod_name, "_")[[1]][2],
#                          life_stage = life_stage, fish_data = fish_data)
# ## @knitr
# 
# if(save_plot) {
#   ggsave(filename=here("analysis","results",
#                        paste0("SR_", strsplit(mod_name, "_")[[1]][2], ".png")),
#          width=11, height=7, units="in", dpi=300)
# } else {
#   dev.new(width=11,height=7)
#   plot(gg)
# }

#-------------------------------------------------------------------------
# Time series of smolt productivity anomaly and natural and hatchery SAR
#-------------------------------------------------------------------------

mod_name <- "fit_Ricker"
save_plot <- TRUE

## @knitr plot_M_anomaly_SAR
gg <- smolt_SAR_ts(mod = get(mod_name), fish_data = fish_data)
## @knitr

if(save_plot) {
  ggsave(filename = here("analysis","results",
                         paste0("M_anomaly_SAR_fit_", strsplit(mod_name, "_")[[1]][2], ".png")), 
         width=7, height=7, units="in", dpi=300)
} else {
  dev.new(width=7,height=7)
  plot(gg)
}

#--------------------------------------------------------------------------------
# Straying matrix: probability of dispersal from each origin to each population
#--------------------------------------------------------------------------------

mod_name <- "fit_Ricker"
save_plot <- TRUE

## @knitr plot_p_origin
gg <- p_origin_plot(mod = get(mod_name), fish_data = fish_data)
## @knitr

if(save_plot) {
  ggsave(filename = here("analysis","results",
                         paste0("p_origin_", strsplit(mod_name, "_")[[1]][2], ".png")), 
         width=11, height=7, units="in", dpi=300)
} else {
  dev.new(width=11, height=7)
  plot(gg)
}

#--------------------------------------------------------------------------------
# Time series of observed and fitted total spawners or smolts for each pop
#--------------------------------------------------------------------------------

mod_name <- "fit_Ricker"
life_stage <- "S"   # "S" = spawners, "M" = smolts
save_plot <- FALSE

## @knitr smolt_spawner_ts
gg <- smolt_spawner_ts(mod = get(mod_name), life_stage = life_stage, 
                       fish_data = switch(head(unlist(strsplit(mod_name, "_")), 1),
                                          foreH0 = fish_data_foreH0, 
                                          foreHmax = fish_data_foreHmax,
                                          fish_data))
## @knitr

if(save_plot) {
  ggsave(filename = here("analysis","results",
                         paste0(life_stage, "_fit_", strsplit(mod_name, "_")[[1]][2], ".png")), 
         width=12, height=7, units="in", dpi=300)
} else {
  dev.new(width=12,height=7)
  plot(gg)
}

#--------------------------------------------------------------------------------
# Observed and fitted distributions of "known" smolt and spawner 
# observation error SDs
#--------------------------------------------------------------------------------

mod_name <- "fit_Ricker"
save_plot <- TRUE

if(save_plot) {
  png(filename = here("analysis","results",
                      paste0("tau_fit_", strsplit(mod_name, "_")[[1]][2], ".png")), 
      width=10, height=5, units="in", res=200, type = "cairo-png")
} else dev.new(width=10,height=5)

## @knitr plot_obs_error_fit
obs_error_plot(get(mod_name), fish_data = fish_data)
## @knitr

if(save_plot) dev.off()

#--------------------------------------------------------------------------------
# Time series of observed and fitted spawner age structure for each pop
#--------------------------------------------------------------------------------

mod_name <- "fit_Ricker"
save_plot <- TRUE

## @knitr plot_spawner_age_ts
gg <- age_timeseries(mod = get(mod_name), fish_data = fish_data)
## @knitr

if(save_plot) {
  ggsave(filename = here("analysis","results",
                         paste0("q_fit_", strsplit(mod_name, "_")[[1]][2], ".png")), 
         width=12, height=7, units="in", dpi=300)
} else {
  dev.new(width=12,height=7)
  plot(gg)
}

#--------------------------------------------------------------------------------
# Time series of observed and fitted sex ratio for each pop
#--------------------------------------------------------------------------------

mod_name <- "fit_Ricker"
save_plot <- TRUE

## @knitr plot_sex_ratio
gg <- plot_sex_ratio(mod = get(mod_name), fish_data = fish_data)
## @knitr

if(save_plot) {
  ggsave(filename = here("analysis", "results", 
                         paste0("q_F_fit_", strsplit(mod_name, "_")[[1]][2], ".png")),
         width=11, height=7, units="in", dpi=300)
} else {
  dev.new(width=11,height=7)
  plot(gg)
}

#--------------------------------------------------------------------------------
# Time series of observed and fitted p_HOS for each pop
#--------------------------------------------------------------------------------

mod_name <- "fit_Ricker"
save_plot <- TRUE

## @knitr plot_p_HOS
gg <- p_HOS_timeseries(mod = get(mod_name), fish_data = fish_data)
## @knitr

if(save_plot) {
  ggsave(filename = here("analysis","results",
                         paste0("p_HOS_fit_", strsplit(mod_name, "_")[[1]][2], ".png")), 
         width=11, height=7, units="in", dpi=300)
} else {
  dev.new(width=11, height=7)
  plot(gg)
}

#--------------------------------------------------------------------------------
# Conditioning forecast trajectories on time-averaged SAR anomalies
#--------------------------------------------------------------------------------

mod_name <- "foreH0_Ricker"
save_plot <- TRUE

if(save_plot) {
  png(filename = here("analysis","results",
                      paste0("SAR_fore_", strsplit(mod_name, "_")[[1]][2], ".png")), 
      width=7, height=7, units="in", res=300, type = "cairo-png")
} else dev.new(width=7,height=7)

## @knitr plot_SAR_fore
SAR_fore_plot(mod = get(mod_name), fish_data_fore = fish_data_fore, example_pop = "Hardy Creek")
## @knitr

if(save_plot) dev.off()

#--------------------------------------------------------------------------------
# Distributions of forecast spawner abundance under alternative scenarios
#--------------------------------------------------------------------------------

modH0_name <- "foreH0_Ricker"
modHmax_name <- "foreHmax_Ricker"
save_plot <- FALSE

## @knitr plot_S_fore
gg <- S_fore_plot(modH0 = get(modH0_name), modHmax = get(modHmax_name), 
                  fish_data_foreH0 = fish_data_foreH0, 
                  fish_data_foreHmax = fish_data_foreHmax,
                  pop_names = pop_names)
## @knitr

if(save_plot) {
  ggsave(filename = here("analysis","results",
                         paste0("S_fore_", strsplit(modH0_name, "_")[[1]][2], ".png")), 
         width=12, height=4, units="in", dpi=300)
} else {
  dev.new(width=12, height=4)
  plot(gg)
}

#------------------------------------------------------------------------------------------
# Distributions of forecast final : initial spawner abundance under alternative scenarios
#------------------------------------------------------------------------------------------

modH0_name <- "foreH0_Ricker"
modHmax_name <- "foreHmax_Ricker"
save_plot <- FALSE

## @knitr plot_StS0_fore
gg <- StS0_fore_plot(modH0 = get(modH0_name), modHmax = get(modHmax_name),
                     fish_data_foreH0 = fish_data_foreH0,
                     fish_data_foreHmax = fish_data_foreHmax,
                     pop_names = pop_names)
## @knitr

if(save_plot) {
  ggsave(filename = here("analysis","results",
                         paste0("StS0_fore_", strsplit(modH0_name, "_")[[1]][2], ".png")),
         width=12, height=4, units="in", dpi=300)
} else {
  dev.new(width=12, height=4)
  plot(gg)
}

#--------------------------------------------------------------------------------
# Distributions of forecast p_HOS under alternative scenarios
#--------------------------------------------------------------------------------

modH0_name <- "foreH0_Ricker"
modHmax_name <- "foreHmax_Ricker"
save_plot <- FALSE

## @knitr plot_p_HOS_fore
gg <- p_HOS_fore_plot(modH0 = get(modH0_name), modHmax = get(modHmax_name), 
                      fish_data_foreH0 = fish_data_foreH0, 
                      fish_data_foreHmax = fish_data_foreHmax,
                      pop_names = pop_names)
## @knitr

if(save_plot) {
  ggsave(filename = here("analysis","results",
                         paste0("p_HOS_fore_", strsplit(modH0_name, "_")[[1]][2], ".png")), 
         width=12, height=4, units="in", dpi=300)
} else {
  dev.new(width=12, height=4)
  plot(gg)
}

#--------------------------------------------------------------------------------
# Probability of quasi-extinction under alternative scenarios
#--------------------------------------------------------------------------------

modH0_name <- "foreH0_Ricker"
modHmax_name <- "foreHmax_Ricker"
save_plot <- FALSE

## @knitr plot_PQE_fore
gg <- PQE_plot(modH0 = get(modH0_name), modHmax = get(modHmax_name), 
               fish_data_foreH0 = fish_data_foreH0, 
               fish_data_foreHmax = fish_data_foreHmax,
               pop_names = pop_names, QET = 50)
## @knitr

if(save_plot) {
  ggsave(filename = here("analysis","results",
                         paste0("PQE_", strsplit(modH0_name, "_")[[1]][2], ".png")), 
         width=12, height=4, units="in", dpi=300)
} else {
  dev.new(width=12, height=4)
  plot(gg)
}

#--------------------------------------------------------------------------------
# Probability of recovery under alternative scenarios
#--------------------------------------------------------------------------------

modH0_name <- "foreH0_Ricker"
modHmax_name <- "foreHmax_Ricker"
save_plot <- FALSE

## @knitr plot_Precovery_fore
gg <- Precovery_plot(modH0 = get(modH0_name), modHmax = get(modHmax_name), 
                     fish_data_foreH0 = fish_data_foreH0, 
                     fish_data_foreHmax = fish_data_foreHmax,
                     pop_names = pop_names, recovery_targets = recovery_targets)
## @knitr

if(save_plot) {
  ggsave(filename = here("analysis","results",
                         paste0("Precovery_", strsplit(modH0_name, "_")[[1]][2], ".png")), 
         width=12, height=4, units="in", dpi=300)
} else {
  dev.new(width=12, height=4)
  plot(gg)
}


#===========================================================================
# TABLES
#===========================================================================

#--------------------------------------------------------------------------------
# Summary of 1-year ahead escapement forecasts by pop
#   (use current calendar year, i.e. pre-season forecast)
#--------------------------------------------------------------------------------

## @knitr forecast_df
draws <- as_draws_rvars(as.matrix(foreH0_Ricker, c("psi","Mmax","S")))
forecast_df <- fish_data_fore %>% 
  mutate(A = round(A/1000, 1), # convert to km
         S = draws$S, S50 = round(median(S)), 
         S05 = round(as.vector(quantile(S, 0.05))), 
         S95 = round(as.vector(quantile(S, 0.95))),
         Forecast = paste0(S50, " (", S05, ", ", S95, ")")) %>% 
  filter(year == year(Sys.Date())) %>% arrange(pop) %>% 
  mutate(psi = draws$psi, psi50 = round(median(psi), 2), 
         psi05 = round(as.vector(quantile(psi, 0.05)), 2), 
         psi95 = round(as.vector(quantile(psi, 0.95)), 2),
         psi_mci = paste0(psi50, " (", psi05, ", ", psi95, ")"),
         Mmax = draws$Mmax/1000, Mmax50 = round(median(Mmax), 1), # smolts/m -> mil/km
         Mmax05 = round(as.vector(quantile(Mmax, 0.05)), 1), 
         Mmax95 = round(as.vector(quantile(Mmax, 0.95)), 1),
         Mmax_mci = paste0(Mmax50, " (", Mmax05, ", ", Mmax95, ")")) %>% 
  filter(pop_type == "natural") %>% 
  select(pop, A, psi_mci, Mmax_mci, Forecast) %>% 
  rename(Population = pop, `Habitat (km)` = A,
         `Max egg-smolt <br> survival ($\\psi$)` = psi_mci, 
         `Smolt capacity <br> ($M_\\text{max}$ 10^6^ km^-1^)` = Mmax_mci)
## @knitr

