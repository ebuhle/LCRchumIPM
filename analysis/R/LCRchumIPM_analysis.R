options(device = windows)

library(salmonIPM)
library(rstan)
library(shinystan)
library(matrixStats)
library(dplyr)
library(reshape2)
library(yarrr)
library(corrplot)
library(magicaxis)
library(zoo)
library(here)


#===========================================================================
# DATA
#===========================================================================

# Mapping of location to population
location_pop <- read.csv(here("data","Location.Reach_Population.csv"), header = TRUE) %>% 
  rename(strata = Strata, location = Location.Reach, pop = Population)

# Mapping of disposition to hatchery vs. wild (i.e., broodstock vs. natural spawner)
disposition_HW <- read.csv(here("data","Disposition_HW.csv"), header = TRUE) %>% 
  rename(disposition = Disposition) %>% arrange(HW)

# Start dates of hatcheries associated with populations
hatchery_start_dates <- read.csv(here("data","Hatchery_Start_Dates.csv"), header = TRUE)

# Spawner abundance data
spawner_data <- read.csv(here("data","Data_ChumSpawnerAbundance_2019-12-12.csv"), header = TRUE) %>% 
  rename(species = Species, stage = LifeStage, year = Return.Yr., strata = Strata,
         location = Location.Reach, disposition = Disposition, method = Method,
         S_obs = Abund.Mean, SD = Abund.SD) %>% 
  mutate(pop = location_pop$pop[match(location, location_pop$location)],
         disposition_HW = disposition_HW$HW[match(disposition, disposition_HW$disposition)]) %>% 
  select(species:location, pop, disposition, disposition_HW, method:SD) %>% 
  arrange(strata, location, year)

# special sum function for use in aggregating S_obs
sum.agg <- function(x) {
  if(length(x) == 0) 0 
  else if(all(is.na(x))) NA 
  else sum(x, na.rm = TRUE)
}

spawner_data_agg <- aggregate(S_obs ~ species + stage + year + strata + location + pop + disposition_HW,
                              data = spawner_data, FUN = sum.agg, na.action = na.pass) %>%
  dcast(species + stage + year + strata + pop ~ disposition_HW, value.var = "S_obs",
        fun.aggregate = sum.agg) %>% rename(S_obs = W, B_take_obs = H) %>% 
  mutate(B_take_obs = replace(B_take_obs, is.na(B_take_obs), 0)) %>% 
  arrange(strata, pop, year)

spawner_data_agg$S_obs <- spawner_data_agg$S_obs + spawner_data_agg$B_take_obs

# Spawner age-, sex-, and origin-frequency (aka BioData)
bio_data <- read.csv(here("data","Data_ChumSpawnerBioData_2019-12-12.csv"), header = TRUE) %>% 
  rename(species = Species, stage = LifeStage, year = Return.Yr., strata = Strata,
         location = Location.Reach, disposition = Disposition, origin = Origin,
         sex = Sex, age = Age, count = Count) %>% 
  mutate(pop = location_pop$pop[match(location, location_pop$location)],
         origin_HW = ifelse(origin == "Natural_spawner", "W", "H"),
         count = ifelse(is.na(count), 0, count)) %>% 
  select(species:location, pop, disposition, origin, origin_HW, sex:count) %>%
  arrange(strata, location, year, origin, age, sex)

bio_data_age <- bio_data %>% 
  dcast(species + stage + year + strata + pop ~ age, value.var = "count", 
        fun.aggregate = sum)

bio_data_origin <- bio_data %>% 
  dcast(species + stage + year + strata + pop ~ origin_HW, value.var = "count", 
        fun.aggregate = sum)

# Fish data formatted for salmonIPM
# Use A = 1 for now (so Rmax in units of spawners)
fish_data <- full_join(spawner_data_agg, bio_data_age, 
                        by = c("species","stage","year","strata","pop")) %>% 
  full_join(bio_data_origin, by = c("species","stage","year","strata","pop")) %>% 
  mutate(B_take_obs = replace(B_take_obs, is.na(B_take_obs), 0)) %>% 
  rename_at(vars(contains("Age-")), funs(paste0(sub("Age-","n_age",.), "_obs"))) %>% 
  rename(n_H_obs = H, n_W_obs = W) %>% mutate(A = 1, fit_p_HOS = NA, F_rate = 0) %>% 
  select(species, stage, strata, pop, year, A, S_obs, n_age2_obs:n_W_obs, 
         fit_p_HOS, B_take_obs, F_rate) %>% arrange(strata, pop, year) 

# fill in fit_p_HOS
for(i in 1:nrow(fish_data)) {
  indx <- match(fish_data$pop[i], hatchery_start_dates$pop)
  fish_data$fit_p_HOS[i] <- ifelse((!is.na(indx) & 
                                     fish_data$year[i] >= hatchery_start_dates$start_brood_year[indx] + 2) |
                                     fish_data$n_H_obs[i] > 0, 1, 0)
}


#----------------------
# Data Exploration
#----------------------




#===========================================================================
# FIT RETROSPECTIVE MODELS
#===========================================================================

#--------------------------------------------------------------
# Spawner-to-spawner IPM
#--------------------------------------------------------------

# Density-independent
SS_exp <- salmonIPM(fish_data = fish_data, stan_model = "IPM_SS_pp", SR_fun = "exp",
                   pars = c("mu_alpha","sigma_alpha","alpha",
                            "sigma_phi","rho_phi","phi",
                            "mu_p","sigma_gamma","R_gamma","gamma",
                            "sigma_p","R_p","p",
                            "p_HOS","sigma","tau","S","R","q","LL"),
                   chains = 3, iter = 1500, warmup = 500,
                   control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))

print(SS_exp, prob = c(0.025,0.5,0.975),
      pars = c("alpha","phi","p_HOS","q","gamma","p","S","R","LL"), 
      include = FALSE, use_cache = FALSE)

launch_shinystan(SS_exp)


# Beverton-Holt
SS_BH <- salmonIPM(fish_data = fish_data, stan_model = "IPM_SS_pp", SR_fun = "BH",
                    pars = c("mu_alpha","sigma_alpha","alpha",
                             "mu_Rmax","sigma_Rmax","Rmax","rho_alphaRmax",
                             "sigma_phi","rho_phi","phi",
                             "mu_p","sigma_gamma","R_gamma","gamma",
                             "sigma_p","R_p","p",
                             "p_HOS","sigma","tau","S","R","q","LL"),
                    chains = 3, iter = 1500, warmup = 500,
                    control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))

print(SS_BH, prob = c(0.025,0.5,0.975),
      pars = c("alpha","Rmax","phi","p_HOS","q","gamma","p","S","R","LL"), 
      include = FALSE, use_cache = FALSE)

launch_shinystan(SS_BH)


# Ricker
SS_Ricker <- salmonIPM(fish_data = fish_data, stan_model = "IPM_SS_pp", SR_fun = "Ricker",
                       pars = c("mu_alpha","sigma_alpha","alpha",
                                "mu_Rmax","sigma_Rmax","Rmax","rho_alphaRmax",
                                "sigma_phi","rho_phi","phi",
                                "mu_p","sigma_gamma","R_gamma","gamma",
                                "sigma_p","R_p","p",
                                "p_HOS","sigma","tau","S","R","q","LL"),
                       chains = 3, iter = 1500, warmup = 500,
                       control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 13))

print(SS_Ricker, prob = c(0.025,0.5,0.975),
      pars = c("alpha","Rmax","phi","p_HOS","q","gamma","p","S","R","LL"), 
      include = FALSE, use_cache = FALSE)

launch_shinystan(SS_Ricker)



#===========================================================================
# FIGURES
#===========================================================================

#--------------------------------------------------------------------------------
# Time series of observed and fitted total spawners for each pop
#--------------------------------------------------------------------------------

mod_name <- "SS_BH"
dat <- fish_data

# dev.new(width=12,height=6)
png(filename=here("analysis","results",paste0("S_fit_", mod_name, ".png")),
    width=12*0.9, height=6*0.9, units="in", res=200, type="cairo-png")

par(mfrow=c(2,4), mar=c(1,2.5,4.1,1), oma=c(4.1,3.1,0,0))

S_IPM <- do.call(extract1, list(as.name(mod_name), "S"))
tau <- do.call(extract1, list(as.name(mod_name), "tau"))
S_obs_IPM <- S_IPM * rlnorm(length(S_IPM), 0, tau)
init_NA <- dat$year
for(i in levels(dat$pop))
  init_NA[dat$pop==i] <- init_NA[dat$pop==i] - min(init_NA[dat$pop==i]) + 1
init_NA <- is.na(dat$S_obs) & init_NA < 5

c1 <- "slategray4"
c1t <- transparent(c1, trans.val = 0.5)
c1tt <- transparent(c1, trans.val = 0.7)

for(i in levels(dat$pop))
{
  y1 <- dat$year[dat$pop==i]
  plot(y1, dat$S_obs[dat$pop==i], pch = "",
       xlim = range(dat$year),
       ylim = range(pmax(dat$S_obs[dat$pop==i], 1),
                    apply(S_obs_IPM[,dat$pop==i & !init_NA], 2, quantile, c(0.025,0.975)), na.rm = T), 
       cex.axis = 1.5, las = 1, xaxt = "n", yaxs = "i", yaxt = "n", xlab = "", ylab = "", log = "y")
  axis(side = 1, at = y1[y1 %% 5 == 0], cex.axis = 1.2)
  rug(y1[y1 %% 5 != 0], ticksize = -0.02)
  at <- maglab(10^par("usr")[3:4], log = T)
  axis(2, at$labat, cex.axis = 1.2, las = 1,
       labels = sapply(log10(at$labat), function(i) as.expression(bquote(10^ .(i)))))
  mtext(i, side = 3, line = 0.5, cex = par("cex")*1.5)
  if(par("mfg")[2] == 1) 
    mtext("Spawners", side = 2, line = 3.5, cex = par("cex")*1.5)
  if(par("mfg")[1] == 2) 
    mtext("Year", side = 1, line = 3, cex = par("cex")*1.5)
  lines(y1, apply(S_IPM[,dat$pop==i], 2, median), col = c1, lwd = 2)
  polygon(c(y1, rev(y1)), 
          c(apply(S_IPM[,dat$pop==i], 2, quantile, 0.025), 
            rev(apply(S_IPM[,dat$pop==i], 2, quantile, 0.975))),
          col = c1t, border = NA)
  polygon(c(y1, rev(y1)), 
          c(apply(S_obs_IPM[,dat$pop==i], 2, quantile, 0.025), 
            rev(apply(S_obs_IPM[,dat$pop==i], 2, quantile, 0.975))),
          col = c1tt, border = NA)
  points(y1, dat$S_obs[dat$pop==i], pch=16, cex = 1)
}

dev.off()
rm(list = c("mod_name","dat","S_IPM","S_obs_IPM","at",
            "c1","c1t","c1tt","y1","init_NA","tau"))


