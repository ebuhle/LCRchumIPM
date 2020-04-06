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


#===========================================================================
# DATA
#===========================================================================

# Mapping of location to population
location_pop <- read.csv(here("data","Location.Reach_Population.csv"), header = TRUE) %>% 
  rename(strata = Strata, location = Location.Reach, pop1 = Population1, pop2 = Population2)

# Mapping of disposition to hatchery vs. wild (i.e., broodstock vs. natural spawner)
disposition_HW <- read.csv(here("data","Disposition_HW.csv"), header = TRUE) %>% 
  rename(disposition = Disposition) %>% arrange(HW)

# Start dates of hatcheries associated with populations
hatchery_start_dates <- read.csv(here("data","Hatchery_Start_Dates.csv"), header = TRUE)

# Spawner abundance data
# Assumptions:
# (1) NAs in hatchery dispositions (incl. Duncan Channel) are really zeros
# (2) NAs in Duncan Creek from 2004-present are really zeros
# (3) All other NAs are real missing observations
spawner_data <- read.csv(here("data","Data_ChumSpawnerAbundance_2019-12-12.csv"), header = TRUE) %>% 
  rename(species = Species, stage = LifeStage, year = Return.Yr., strata = Strata,
         location = Location.Reach, disposition = Disposition, method = Method,
         S_obs = Abund.Mean, SD = Abund.SD) %>% 
  mutate(pop = location_pop$pop2[match(location, location_pop$location)],
         disposition_HW = disposition_HW$HW[match(disposition, disposition_HW$disposition)],
         S_obs = replace(S_obs, is.na(S_obs) & disposition_HW == "H", 0),
         S_obs = replace(S_obs, is.na(S_obs) & pop == "Duncan_Creek" & year >= 2004, 0)) %>% 
  select(species:location, pop, disposition, disposition_HW, method:SD) %>% 
  arrange(strata, location, year)

names_S_obs <- disposition_HW$disposition
names_B_take_obs <- disposition_HW$disposition[disposition_HW$HW == "H"]

spawner_data_agg <- aggregate(S_obs ~ species + stage + year + strata + pop + disposition,
                              data = spawner_data, FUN = sum, na.action = na.pass) %>%
  dcast(species + stage + year + strata + pop ~ disposition, value.var = "S_obs", 
        fun.aggregate = identity, fill = 0) %>% 
  add_column(S_obs = rowSums(select(., all_of(names_S_obs))),
             B_take_obs = rowSums(select(., all_of(names_B_take_obs)))) %>% 
  select(-all_of(names_S_obs)) %>% arrange(strata, pop, year)

# Spawner age-, sex-, and origin-frequency (aka BioData)
bio_data <- read.csv(here("data","Data_ChumSpawnerBioData_2019-12-12.csv"), header = TRUE) %>% 
  rename(species = Species, stage = LifeStage, year = Return.Yr., strata = Strata,
         location = Location.Reach, disposition = Disposition, origin = Origin,
         sex = Sex, age = Age, count = Count) %>% 
  mutate(pop = location_pop$pop2[match(location, location_pop$location)],
         origin_HW = ifelse(origin == "Natural_spawner", "W", "H"),
         count = ifelse(is.na(count), 0, count)) %>% 
  select(species:location, pop, disposition, origin, origin_HW, sex:count) %>%
  arrange(strata, location, year, origin, age, sex)

# age of wild spawners only
bio_data_age <- bio_data %>% filter(origin_HW == "W") %>%  
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

# drop cases with initial NAs in S_obs, even if bio data is present
init_NA <- function(x) { cumsum(!is.na(x)) > 0 }
fish_data <- fish_data %>% group_by(pop) %>% filter(init_NA(S_obs)) %>% as.data.frame()
                                                                  
# fill in fit_p_HOS
for(i in 1:nrow(fish_data)) {
  indx <- match(fish_data$pop[i], hatchery_start_dates$pop)
  fish_data$fit_p_HOS[i] <- ifelse((!is.na(indx) & 
                                     fish_data$year[i] >= hatchery_start_dates$start_brood_year[indx] + 2) |
                                     fish_data$n_H_obs[i] > 0, 1, 0)
}

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
SS_exp <- salmonIPM(fish_data = fish_data, stan_model = "IPM_SS_pp", SR_fun = "exp",
                   pars = c("mu_alpha","sigma_alpha","alpha",
                            "sigma_phi","rho_phi","phi",
                            "mu_p","sigma_gamma","R_gamma","gamma",
                            "sigma_p","R_p","p",
                            "p_HOS","sigma","tau","S","R","q","LL"),
                   chains = 3, iter = 1500, warmup = 500,
                   control = list(adapt_delta = 0.99, max_treedepth = 13))

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
                    control = list(adapt_delta = 0.99, max_treedepth = 13))

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
                       control = list(adapt_delta = 0.99, max_treedepth = 13))

print(SS_Ricker, prob = c(0.025,0.5,0.975),
      pars = c("alpha","Rmax","phi","p_HOS","q","gamma","p","S","R","LL"), 
      include = FALSE, use_cache = FALSE)

launch_shinystan(SS_Ricker)


#--------------------------------------------------------------
# Model selection using LOO
#--------------------------------------------------------------

# Observationwise log-likelihood of each fitted model
# Here an observation is a row of fish_data, and the total likelihood includes 
# components for spawner abundance, age-frequency, and hatchery/wild-frequency
LL <- lapply(list(exp = SS_exp, BH = SS_BH, Ricker = SS_Ricker),
             loo::extract_log_lik, parameter_name = "LL", merge_chains = FALSE)

# Relative ESS of posterior draws of observationwise likelihood 
r_eff <- lapply(LL, function(x) relative_eff(exp(x)))

# PSIS-LOO
LOO <- lapply(1:length(LL), function(i) loo(LL[[i]], r_eff = r_eff[[i]]))
names(LOO) <- names(LL)

## Compare all three models
loo_compare(LOO)

## Exponential vs. Ricker
loo_compare(LOO[c("exp","Ricker")])

## Exponential vs. Beverton-Holt
loo_compare(LOO[c("exp","BH")])

## Beverton-Holt vs. Ricker
loo_compare(LOO[c("BH","Ricker")])



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

mod_name <- "SS_Ricker"

SR_eval <- function(alpha, Rmax = NULL, S, SR_fun = "BH") 
{
  switch(SR_fun,
         exp = alpha*S,
         BH = alpha*S/(1 + alpha*S/Rmax),
         Ricker = alpha*S*exp(-alpha*S/(exp(1)*Rmax)))
}

mu_alpha <- as.vector(do.call(extract1, list(as.name(mod_name), "mu_alpha")))
mu_Rmax <- as.vector(do.call(extract1, list(as.name(mod_name), "mu_Rmax")))
S <- matrix(seq(0, quantile(fish_data$S_obs/fish_data$A, 0.9, na.rm = TRUE), length = 100),
            nrow = length(mu_alpha), ncol = 100, byrow = TRUE)
R_ESU_IPM <- SR_eval(alpha = exp(mu_alpha), Rmax = exp(mu_Rmax), S = S, 
                     SR_fun = substring(mod_name, 4))
alpha <- do.call(extract1, list(as.name(mod_name), "alpha"))
Rmax <- do.call(extract1, list(as.name(mod_name), "Rmax"))
R_pop_IPM <- sapply(1:ncol(alpha), function(i) {
  colMedians(SR_eval(alpha = alpha[,i], Rmax = Rmax[,i], S = S,
                     SR_fun = substring(mod_name, 4)))
  })

c1 <- "slategray4"
c1t <- transparent(c1, trans.val = 0.5)
c1tt <- transparent(c1, trans.val = 0.7)

# dev.new(width = 7, height = 7)
png(filename=here("analysis","results",paste0("SR_",mod_name,".png")), 
    width=7, height=7, units="in", res=200, type="cairo-png")

par(mfrow = c(2,2), mar = c(5.1,5.1,1,1))

# Recruits vs. spawners
plot(S[1,], colMedians(R_ESU_IPM), type = "l", lwd=3, col = c1, las = 1,
     cex.axis = 1.2, xaxs = "i", yaxs = "i", ylim = range(R_pop_IPM),
     xlab = "" , ylab = "")
for(i in 1:ncol(R_pop_IPM))
  lines(S[1,], R_pop_IPM[,i], col = c1t)
polygon(c(S[1,], rev(S[1,])), 
        c(colQuantiles(R_ESU_IPM, probs = 0.025), rev(colQuantiles(R_ESU_IPM, probs = 0.975))), 
        col = c1tt, border = NA)
title(xlab="Spawners", ylab="Recruits", line = 3.5, cex.lab = 1.8)
# text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "A", cex = 2)

# Recruits/spawner vs. spawners
plot(S[1,], colMedians(R_ESU_IPM)/S[1,], type = "l", lwd=3, col = c1, las = 1,
     cex.axis = 1.2, xlim = range(S[S > 0]), xaxs = "i", xlab = "" , 
     ylim = range(colQuantiles(R_ESU_IPM[,-1]/S[1,-1], probs = c(0.025,0.975))), 
     yaxs = "i", ylab = "")
for(i in 1:ncol(R_pop_IPM))
  lines(S[1,], R_pop_IPM[,i]/S[1,], col = c1t)
polygon(c(S[1,-1], rev(S[1,-1])),
        c(colQuantiles(R_ESU_IPM[,-1]/S[1,-1], probs = 0.025),
          rev(colQuantiles(R_ESU_IPM[,-1]/S[1,-1], probs = 0.975))),
        col = c1tt, border = NA)
title(xlab="Spawners", ylab="Recruits / spawner", line = 3.5, cex.lab = 1.8)
# text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "B", cex = 2)

# Posterior densities of log(alpha)
dd_IPM_ESU <- density(mu_alpha)
dd_IPM_pop <- vector("list", length(levels(fish_data$pop)))
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
dd_IPM_pop <- vector("list", length(levels(fish_data$pop)))
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

rm(list=c("mu_alpha","mu_Rmax","S","R_ESU_IPM","c1","c1t","c1tt",
          "dd_IPM_ESU","dd_IPM_pop","alpha","Rmax","R_pop_IPM"))
dev.off()


#--------------------------------------------------------------------------------
# Time series of observed and fitted total spawners for each pop
#--------------------------------------------------------------------------------

mod_name <- "SS_Ricker"

# dev.new(width=12,height=8)
png(filename=here("analysis","results",paste0("S_fit_", mod_name, ".png")),
    width=12*0.9, height=8*0.9, units="in", res=200, type="cairo-png")

par(mfrow=c(3,4), mar=c(1,2.5,4.1,1), oma=c(4.1,3.1,0,0))

S_IPM <- do.call(extract1, list(as.name(mod_name), "S"))
tau <- do.call(extract1, list(as.name(mod_name), "tau"))
S_obs_IPM <- S_IPM * rlnorm(length(S_IPM), 0, tau)

c1 <- "slategray4"
c1t <- transparent(c1, trans.val = 0.5)
c1tt <- transparent(c1, trans.val = 0.7)

for(i in levels(fish_data$pop))
{
  yi <- fish_data$year[fish_data$pop==i]
  plot(yi, fish_data$S_obs[fish_data$pop==i], pch = "",
       xlim = range(fish_data$year),
       ylim = range(pmax(fish_data$S_obs[fish_data$pop==i], 1),
                    colQuantiles(S_obs_IPM[,fish_data$pop==i], probs = c(0.025,0.975)), na.rm = T), 
       cex.axis = 1.5, las = 1, xaxt = "n", yaxs = "i", yaxt = "n", xlab = "", ylab = "", log = "y")
  axis(side = 1, at = yi[yi %% 5 == 0], cex.axis = 1.2)
  rug(yi[yi %% 5 != 0], ticksize = -0.02)
  at <- maglab(10^par("usr")[3:4], log = T)
  axis(2, at$labat, cex.axis = 1.2, las = 1,
       labels = sapply(log10(at$labat), function(i) as.expression(bquote(10^ .(i)))))
  mtext(i, side = 3, line = 0.5, cex = par("cex")*1.5)
  if(par("mfg")[2] == 1) 
    mtext("Spawners", side = 2, line = 3.5, cex = par("cex")*1.5)
  if(par("mfg")[1] == par("mfg")[3]) 
    mtext("Year", side = 1, line = 3, cex = par("cex")*1.5)
  lines(yi, colMedians(S_IPM[,fish_data$pop==i]), col = c1, lwd = 3)
  polygon(c(yi, rev(yi)), 
          c(colQuantiles(S_IPM[,fish_data$pop==i], probs = 0.025), 
            rev(colQuantiles(S_IPM[,fish_data$pop==i], probs = 0.975))),
          col = c1t, border = NA)
  polygon(c(yi, rev(yi)), 
          c(colQuantiles(S_obs_IPM[,fish_data$pop==i], probs = 0.025), 
            rev(colQuantiles(S_obs_IPM[,fish_data$pop==i], probs = 0.975))),
          col = c1tt, border = NA)
  points(yi, fish_data$S_obs[fish_data$pop==i], pch=16, cex = 1.5)
}

dev.off()
rm(list = c("mod_name","S_IPM","S_obs_IPM","at","c1","c1t","c1tt","yi","tau"))


#--------------------------------------------------------------------------------
# Shared process errors (brood year productivity anomalies)
#--------------------------------------------------------------------------------

mod_name <- "SS_Ricker"

# dev.new(width=7,height=5)
png(filename=here("analysis","results",paste0("phi_", mod_name, ".png")),
    width=7, height=5, units="in", res=200, type="cairo-png")

y <- sort(unique(fish_data$year))
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



