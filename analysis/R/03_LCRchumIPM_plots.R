#==========================================================================#
# FUNCTIONS TO MAKE CUSTOM PLOTS FOR LCR CHUM IPM
#==========================================================================#

#--------------------------------------------------------------------
# Life-cycle multiplot
# S-R curves (spawners to smolts)
# Posterior distributions of fecundity, survival and Mmax parameters
# Time series of smolt productivity process errors and SAR
#--------------------------------------------------------------------

# Spawner-recruit function
SR_eval <- function(alpha, Rmax = NULL, S, SR_fun) 
{
  switch(SR_fun,
         exp = alpha*S,
         BH = alpha*S/(1 + alpha*S/(Rmax)),
         Ricker = alpha*S*exp(-alpha*S/(exp(1)*Rmax)))
}

# Plot function
LCRchumIPM_multiplot <- function(mod, SR_fun, fish_data, save_plot = FALSE, filename = NULL)
{
  # fecundity
  mu_E <- extract1(mod, "mu_E")
  ages <- substring(names(select(fish_data, starts_with("n_age"))), 6, 6)
  delta_NG <- extract1(mod, "delta_NG")
  # egg deposition and egg-to-smolt survival
  A <- fish_data$A
  q <- extract1(mod, "q")
  q_F <- extract1(mod, "q_F")
  mu_psi <- extract1(mod, "mu_psi")
  psi <- extract1(mod, "psi")
  alpha <- apply(sweep(q, c(1,3), mu_E, "*"), 1:2, sum) * q_F * psi[,as.numeric(fish_data$pop)]
  alpha <- t(as.matrix(aggregate(t(alpha), list(pop = fish_data$pop), median)[,-1]))
  mu_alpha <- rowMeans(log(alpha))
  mu_Mmax <- as.vector(extract1(mod, "mu_Mmax"))
  Mmax <- extract1(mod, "Mmax")
  S <- colMedians(extract1(mod, "S"))
  SA_grid <- matrix(seq(0, quantile(S/A, 0.9, na.rm = TRUE), length = 100),
                    nrow = length(mu_alpha), ncol = 100, byrow = TRUE)
  M_ESU <- SR_eval(alpha = exp(mu_alpha), Rmax = exp(mu_Mmax)*1e3, S = SA_grid, SR_fun = SR_fun)
  M_pop <- sapply(1:ncol(Mmax), function(i) {
    colMedians(SR_eval(alpha = alpha[,i], Rmax = Mmax[,i]*1e3, S = SA_grid, SR_fun = SR_fun))
  })
  # smolt recruitment process errors
  y <- sort(unique(fish_data$year))
  eta_year_M <- extract1(mod, "eta_year_M")
  sigma_M <- extract1(mod, "sigma_M")
  zeta_M <- stan_mean(mod, "zeta_M")
  epsilon_M <- outer(sigma_M, zeta_M, "*")
  error_M <- eta_year_M[,as.numeric(factor(fish_data$year))] + epsilon_M
  # SAR
  eta_year_MS <- extract1(mod, "eta_year_MS")
  mu_MS <- extract1(mod, "mu_MS")
  s_hat_MS <- plogis(sweep(eta_year_MS, 1, qlogis(mu_MS), "+"))
  s_MS <- extract1(mod, "s_MS")
  # colors
  c1 <- "slategray4"
  c1t <- transparent(c1, trans.val = 0.3)
  c1tt <- transparent(c1, trans.val = 0.5)
  ac <- viridis(length(ages), end = 0.8, direction = -1, alpha = 0.5) 

  if(save_plot) {
    png(filename=filename, width=12, height=5.5, units="in", res=300, type="cairo-png")
  } else dev.new(width = 12, height = 5.5)
  
  par(mfrow = c(2,4), mar = c(5.1,5.1,1,0.5), oma = c(0,0,0,1))
  
  # Posterior distributions of fecundity by age
  dd_age <- lapply(as.data.frame(mu_E), density)
  
  plot(dd_age[[1]]$x, dd_age[[1]]$y, pch = "", las = 1, cex.axis = 1.2, cex.lab = 1.5,
       xlab = bquote("Mean fecundity (" * mu[italic(E)] * ")"), ylab = "", 
       xlim = range(sapply(dd_age, function(m) m$x)),
       ylim = range(sapply(dd_age, function(m) m$y)),
       xaxs = "i", yaxt = "n")
  for(a in 1:length(ages))
    polygon(dd_age[[a]], col = ac[a], border = NA)
  title(ylab = "Probability density", line = 2, cex.lab = 1.5)
  text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "A", cex = 1.5)
  legend("topright", paste("age", ages, "  "), cex = 1.2, fill = ac, border = NA,
         bty = "n", inset = c(-0.05,0))
  
  # Posterior distribution of non-green female fecundity offset
  dd_NG <- density(delta_NG)
  
  plot(dd_NG, type = "n", las = 1, cex.axis = 1.2, cex.lab = 1.5,
       xlab = bquote("Non-green fecundity offset (" * delta[NG] * ")"), ylab = "", main = "",
       xlim = c(0,1), xaxs = "i", yaxt = "n", xpd = NA)
  polygon(dd_NG, col = c1tt, border = NA)
  text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "B", cex = 1.5)
  
  # Posterior densities of psi
  dd_ESU <- density(mu_psi)
  dd_pop <- lapply(as.data.frame(psi), density)
  
  plot(dd_ESU$x, dd_ESU$y, pch = "", lwd = 3, col = c1, las = 1, 
       xaxs = "i", yaxt = "n", cex.axis = 1.2, cex.lab = 1.5, 
       xlab = bquote("Maximum egg-to-smolt survival (" * psi * ")"), ylab = "",
       xlim = c(0,1), ylim = range(dd_ESU$y, sapply(dd_pop, function(m) m$y)), xpd = NA)
  polygon(dd_ESU, col = c1tt, border = NA)
  for(i in 1:length(dd_pop))
    lines(dd_pop[[i]]$x, dd_pop[[i]]$y, col = c1t)
  text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "C", cex = 1.5)
  
  # Posterior densities of log(Mmax)
  dd_ESU <- density(mu_Mmax * log10(exp(1)))  # convert to base 10
  dd_pop <- lapply(as.data.frame(log10(Mmax)), density)
  
  plot(dd_ESU$x, dd_ESU$y, pch = "", lwd = 3, col = c1, las = 1, 
       xaxt = "n", yaxt = "n", cex.axis = 1.2, cex.lab = 1.5, 
       xlab = bquote("Maximum smolts" ~ km^-1 ~ "(" * italic(M)[max] * ")"), ylab = "",
       xlim = range(dd_ESU$x[dd_ESU$y > 0.02], sapply(dd_pop, function(m) range(m$x[m$y > 0.02]))),
       ylim = range(dd_ESU$y, sapply(dd_pop, function(m) m$y)), xpd = NA)
  polygon(dd_ESU, col = c1tt, border = NA)
  for(i in 1:length(dd_pop))
    lines(dd_pop[[i]]$x, dd_pop[[i]]$y, col = c1t)
  tck <- maglab(10^par("usr")[1:2], log = TRUE)
  axis(1, at = log10(tck$labat), labels = tck$labat, cex.axis = 1.2, hadj = 0.5)
  text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "D", cex = 1.5)
  
  # Spawner-recruit function
  plot(SA_grid[1,], colMedians(M_ESU)*1e-6, type = "l", lwd=3, col = c1, las = 1,
       cex.axis = 1.2, cex.lab = 1.5, xaxs = "i", yaxs = "i",
       ylim = range(M_pop*1e-6), xlab = bquote("Spawner density (" * km^-1 * ")"), 
       ylab = bquote("Smolt density (" * 10^6 ~ km^-1 * ")"))
  for(i in 1:ncol(M_pop)) 
    lines(SA_grid[1,], M_pop[,i]*1e-6, col = c1t)
  polygon(c(SA_grid[1,], rev(SA_grid[1,])), 
          c(colQuantiles(M_ESU, probs = 0.05), rev(colQuantiles(M_ESU, probs = 0.95)))*1e-6, 
          col = c1tt, border = NA)
  rug(S, col = c1t)
  text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "E", cex = 1.5)
  
  # Smolts-per-spawner function
  plot(SA_grid[1,], colMedians(M_ESU)/SA_grid[1,], type = "l", lwd=3, col = c1, las = 1,
       cex.axis = 1.2, cex.lab = 1.5, xaxs = "i", yaxs = "i", 
       ylim = range(sweep(M_pop[-1,], 1, SA_grid[1,-1], "/")), 
       xlab = bquote("Spawner density (" * km^-1 * ")"), ylab = "Smolts per spawner")
  for(i in 1:ncol(M_pop)) 
    lines(SA_grid[1,], M_pop[,i]/SA_grid[1,], col = c1t)
  polygon(c(SA_grid[1,], rev(SA_grid[1,])), 
          c(colQuantiles(M_ESU, probs = 0.05)/SA_grid[1,], 
            rev(colQuantiles(M_ESU, probs = 0.95)/SA_grid[1,])), 
          col = c1tt, border = NA)
  rug(S, col = c1t)
  text(par("usr")[2], par("usr")[4], adj = c(1.5,1.5), "F", cex = 1.5)
  
  # Smolt recruitment process errors
  plot(y, colMedians(eta_year_M), type = "n", las = 1, cex.axis = 1.2, cex.lab = 1.5,
       ylim = range(colQuantiles(eta_year_M, probs = c(0.05, 0.95)), 
                    colQuantiles(error_M, probs = c(0.05, 0.95))), 
       xaxs = "i", xaxt = "n", xlab = "Brood year", 
       ylab = "Smolt productivity anomaly", xpd = NA)
  polygon(c(y, rev(y)), 
          c(colQuantiles(eta_year_M, probs = 0.05), 
            rev(colQuantiles(eta_year_M, probs = 0.95))),
          col = c1tt, border = NA)
  lines(y, colMedians(eta_year_M), col = c1t, lwd = 3)
  for(j in levels(fish_data$pop))
    lines(fish_data$year[fish_data$pop == j], colMedians(error_M[,fish_data$pop == j]), col = c1t)
  axis(side = 1, at = y[y %% 5 == 0], cex.axis = 1.2)
  rug(y[y %% 5 != 0], ticksize = -0.02)
  text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "G", cex = 1.5)
  
  # SAR
  plot(y, colMedians(s_hat_MS)*100, type = "n", las = 1, cex.axis = 1.2, cex.lab = 1.5,
       ylim = range(colQuantiles(s_hat_MS, probs = 0.95), colMedians(s_MS))*100, 
       xaxs = "i", xaxt = "n", xlab = "Outmigration year", ylab = "SAR (%)")
  # mtext("SAR (%)", side = 2, line = 3.7, cex = par("cex")*1.5)
  polygon(c(y, rev(y)), 
          c(colQuantiles(s_hat_MS, probs = 0.05), 
            rev(colQuantiles(s_hat_MS, probs = 0.95)))*100,
          col = c1tt, border = NA)
  lines(y, colMedians(s_hat_MS)*100, col = c1t, lwd = 3)
  for(j in levels(fish_data$pop))
    lines(fish_data$year[fish_data$pop == j], colMedians(s_MS[,fish_data$pop == j])*100, col = c1t)
  axis(side = 1, at = y[y %% 5 == 0], cex.axis = 1.2)
  rug(y[y %% 5 != 0], ticksize = -0.02)
  text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "H", cex = 1.5)
  
  if(save_plot) dev.off()
}

#--------------------------------------------------------------------
# Spawner-to-smolt S-R curve for each pop with data and states
#--------------------------------------------------------------------

### NOT READY FOR PRIME TIME 
### b/c of orders-of-magnitude differences in axis scales among pops
### standardizing by habitat area may help
### ALSO overplotting data for Grays MS doesn't really work b/c of upstream smolts

LCRchumIPM_SR_plot <- function(mod, SR_fun, fish_data, save_plot = FALSE, 
                               show_plot = !save_plot, filename = NULL)
{
  mu_E <- extract1(mod, "mu_E")
  q <- extract1(mod, "q")
  q_F <- extract1(mod, "q_F")
  psi <- extract1(mod, "psi")
  f <- apply(sweep(q, c(1,3), mu_E, "*"), 1:2, sum) * q_F  # eggs per spawner
  alpha <- f * psi[,as.numeric(fish_data$pop)]
  alpha <- t(as.matrix(aggregate(t(alpha), list(pop = fish_data$pop), median)[,-1]))
  Mmax <- extract1(mod, "Mmax")
  S <- extract1(mod, "S")
  M <- extract1(mod, "M")*1e-6  # convert smolts to millions
  dat <- fish_data %>% 
    mutate(S_obs_L = qlnorm(0.05, log(S_obs), tau_S_obs),
           S_obs_U = qlnorm(0.95, log(S_obs), tau_S_obs),
           M_obs = M_obs*1e-6,  # convert smolts to millions
           M_obs_L = qlnorm(0.05, log(M_obs), tau_M_obs),
           M_obs_U = qlnorm(0.95, log(M_obs), tau_M_obs),
           S_L = colQuantiles(S, probs = 0.05), S_U = colQuantiles(S, probs = 0.95), 
           S = colMedians(S), M_L = colQuantiles(M, probs = 0.05), 
           M_U = colQuantiles(M, probs = 0.95), M = colMedians(M), f = colMedians(f))
  
  gg <- expand_grid(pop = 1:ncol(alpha), iter = 1:nrow(alpha)) %>% 
    mutate(pop = factor(levels(fish_data$pop)[pop], levels = levels(fish_data$pop)),
           alpha = as.vector(alpha), Mmax = as.vector(Mmax), 
           Sadj = tapply(dat$S, dat$pop, function(x) max(x)*1.5)[pop],
           Smax = tapply(dat$S_U, dat$pop, max)[pop]) %>%
    group_by(pop, iter, alpha, Mmax) %>% summarize(S = seq(0, min(Sadj,Smax), length = 100)) %>% 
    mutate(M = SR_eval(alpha = alpha, Rmax = Mmax*1e3, S = S, SR_fun = SR_fun)*1e-6) %>% 
    group_by(pop,S) %>% summarize(L = quantile(M,0.05), U = quantile(M,0.95), M = median(M)) %>%
    ggplot(aes(x = S, y = M)) +
    geom_ribbon(aes(ymin = L, ymax = U), fill = "slategray4", alpha = 0.5) +
    geom_line(col = "slategray4", alpha = 0.7, lwd = 1) +
    geom_pointrange(aes(x = S, y = M, ymin = M_L, ymax = M_U), data = dat,
                    pch = 16, fatten = 3, col = "slategray4", alpha = 0.7) +
    # geom_errorbarh(aes(y = M, xmin = S_L, xmax = S_U), data = dat,
    #                col = "slategray4", alpha = 0.7) +
    geom_pointrange(aes(x = S_obs, y = M_obs, ymin = M_obs_L, ymax = M_obs_U), data = dat,
                    pch = 16, fatten = 3, alpha = 0.7) +
    geom_errorbarh(aes(y = M_obs, xmin = M_obs_L, xmax = M_obs_U), data = dat, alpha = 0.7) +
    scale_x_continuous(expand = expansion(0,0)) + labs(x = "Spawners", y = "Smolts") +
    facet_wrap(vars(pop), ncol = 4, scales = "free") + theme_bw(base_size = 16) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(),
          strip.background = element_rect(fill = NA))

  if(show_plot) {
    dev.new(width=11,height=7)
    show(gg)
  }
  if(save_plot)   
    ggsave(filename=filename, width=11, height=7, units="in", dpi=300, type="cairo-png")
}


#--------------------------------------------------------------------------------
# Time series of observed and fitted total spawners or smolts for each pop
#--------------------------------------------------------------------------------

LCRchumIPM_MS_timeseries <- function(mod, life_stage = c("M","S"), fish_data, 
                                     save_plot = FALSE, show_plot = !save_plot, 
                                     filename = NULL)
{
  life_cycle <- strsplit(mod@model_name, "_")[[1]][2]
  N <- extract1(mod, life_stage)
  tau <- extract1(mod, switch(life_cycle, SS = "tau", LCRchum = paste0("tau_", life_stage)))
  if(life_cycle == "LCRchum" & life_stage == "M")
    N[,na.omit(fish_data$downstream_trap)] <- N[,na.omit(fish_data$downstream_trap)] + 
    N[,which(!is.na(fish_data$downstream_trap))]
  N_ppd <- N * rlnorm(length(N), 0, tau)
  
  gg <- fish_data %>% 
    mutate(N_obs = !!sym(paste0(life_stage, "_obs")),
           tau_obs = !!sym(paste0("tau_", life_stage, "_obs")),
           N_obs_L = qlnorm(0.05, log(N_obs), tau_obs),
           N_obs_U = qlnorm(0.95, log(N_obs), tau_obs),
           N_L = colQuantiles(N, probs = 0.05),
           N_U = colQuantiles(N, probs = 0.95),
           N_ppd_L = colQuantiles(N_ppd, probs = 0.05),
           N_ppd_U = colQuantiles(N_ppd, probs = 0.95),
           pch = ifelse(is.na(tau_obs), 1, 16)) %>% 
    ggplot(aes(x = year, y = N_obs)) +
    geom_ribbon(aes(ymin = N_L, ymax = N_U), fill = "slategray4", alpha = 0.5) +
    geom_ribbon(aes(ymin = N_ppd_L, ymax = N_ppd_U), fill = "slategray4", alpha = 0.3) +
    geom_point(aes(shape = pch), size = 2.5) + scale_shape_identity() +
    geom_errorbar(aes(ymin = N_obs_L, ymax = N_obs_U), width = 0) +
    labs(x = "Year", y = switch(life_stage, M = "Smolts (thousands)", S = "Spawners")) + 
    scale_y_log10(labels = function(y) y*switch(life_stage, M = 1e-3, S = 1)) + 
    facet_wrap(vars(pop), ncol = 4) + theme_bw(base_size = 16) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(),
          strip.background = element_rect(fill = NA))
  
  if(show_plot) {
    dev.new(width=11,height=7)
    show(gg)
  }
  if(save_plot)   
    ggsave(filename=filename, width=11, height=7, units="in", dpi=300, type="cairo-png")
}

#--------------------------------------------------------------------------------
# Time series of observed and fitted spawner age structure for each pop
#--------------------------------------------------------------------------------

LCRchumIPM_age_timeseries <- function(mod, fish_data, save_plot = FALSE, 
                                      show_plot = !save_plot, filename = NULL)
{
  q <- extract1(mod, "q")
  
  gg <- fish_data %>% 
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
    labs(x = "Year", y = "Proportion at age") + 
    facet_wrap(vars(pop), ncol = 4) + theme_bw(base_size = 16) + 
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
          strip.background = element_rect(fill = NA), 
          legend.box.margin = margin(0,-10,0,-15))
  
  if(show_plot) {
    dev.new(width=12,height=7)
    show(gg)
  }
  if(save_plot)
    ggsave(filename=filename, width=12, height=7, units="in", dpi=300, type="cairo-png")
}

#--------------------------------------------------------------------------------
# Time series of observed and fitted sex ratio for each pop
#--------------------------------------------------------------------------------

LCRchumIPM_sex_timeseries <- function(mod, fish_data, save_plot = FALSE, 
                                      show_plot = !save_plot, filename = NULL)
{
  q_F <- extract1(mod, "q_F")
  
  gg <- cbind(fish_data, colQuantiles(q_F, probs = c(0.05, 0.5, 0.95))) %>%
    mutate(n_MF_obs = n_M_obs + n_F_obs) %>% 
    cbind(., with(., binconf(x = n_F_obs, n = n_MF_obs))) %>%
    ggplot(aes(x = year, y = PointEst, ymin = Lower, ymax = Upper)) + 
    geom_abline(intercept = 0.5, slope = 0, color = "gray") + 
    geom_ribbon(aes(ymin = `5%`, ymax = `95%`), fill = "slategray4", alpha = 0.5) +
    geom_line(aes(y = `50%`), col = "slategray4", lwd = 1) +
    geom_point(pch = 16, size = 2.5) + geom_errorbar(width = 0) +
    labs(x = "Year", y = "Proportion female") +
    facet_wrap(vars(pop), ncol = 4) + theme_bw(base_size = 16) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
          strip.background = element_rect(fill = NA))
  
  if(show_plot) {
    dev.new(width=11,height=7)
    show(gg)
  }
  if(save_plot) 
    ggsave(filename=here("analysis", "results", paste0("q_F_fit_", mod_name, ".png")),
           width=11, height=7, units="in", dpi=300, type="cairo-png")
}

#--------------------------------------------------------------------------------
# Time series of observed and fitted p_HOS for each pop
#--------------------------------------------------------------------------------

LCRchumIPM_p_HOS_timeseries <- function(mod, fish_data, save_plot = FALSE, 
                                        show_plot = !save_plot, filename = NULL)
{
  p_HOS <- extract1(mod, "p_HOS")
  
  gg <- fish_data %>% 
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
    facet_wrap(vars(pop), ncol = 4) + theme_bw(base_size = 16) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = NA))
  
  if(show_plot) {
    dev.new(width=11,height=7)
    show(gg)
  }
  if(save_plot) 
    ggsave(filename=filename, width=11, height=7, units="in", dpi=300, type="cairo-png")
}

#--------------------------------------------------------------------------------
# Distributions of observed and fitted fecundity by age
#--------------------------------------------------------------------------------

LCRchumIPM_fecundity_plot <- function(mod, fish_data, fecundity_data,
                                      save_plot = FALSE, filename = NULL)
{
  ages <- substring(names(select(fish_data, starts_with("n_age"))), 6, 6)
  E_obs <- fecundity_data$E_obs
  E_seq <- seq(min(E_obs, na.rm = TRUE), max(E_obs, na.rm = TRUE), length = 500)
  mu_E <- extract1(mod, "mu_E")
  sigma_E <- extract1(mod, "sigma_E")
  E_fit <- array(NA, c(nrow(mu_E), length(E_seq), ncol(mu_E)))
  for(a in 1:length(ages))
    E_fit[,,a] <- sapply(E_seq, function(x) dnorm(x, mu_E[,a], sigma_E[,a]))
  
  c1 <- viridis(length(ages), end = 0.8, direction = -1) 
  c1t <- transparent(c1, trans.val = 0.5)
  c1tt <- transparent(c1, trans.val = 0.7)
  
  if(save_plot) {
    png(filename=filename, width=7, height=7, units="in", res=200, type="cairo-png")
  } else dev.new(width=7,height=7)

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
  
  if(save_plot) dev.off()
}


#--------------------------------------------------------------------------------
# Observed and fitted distributions of "known" smolt and spawner 
# observation error SDs
#--------------------------------------------------------------------------------

LCRchumIPM_obs_error_plot <- function(mod, fish_data, save_plot = FALSE, filename = NULL)
{
  tau_M_obs <- fish_data$tau_M_obs
  tau_M_seq <- seq(min(tau_M_obs, na.rm = TRUE), max(tau_M_obs, na.rm = TRUE), length = 500)
  mu_tau_M <- extract1(mod, "mu_tau_M")
  sigma_tau_M <- extract1(mod, "sigma_tau_M")
  tau_M_fit <- sapply(tau_M_seq, function(x) dlnorm(x, log(mu_tau_M), sigma_tau_M))
  
  tau_S_obs <- fish_data$tau_S_obs
  tau_S_seq <- seq(min(tau_S_obs, na.rm = TRUE), max(tau_S_obs, na.rm = TRUE), length = 500)
  mu_tau_S <- extract1(mod, "mu_tau_S")
  sigma_tau_S <- extract1(mod, "sigma_tau_S")
  tau_S_fit <- sapply(tau_S_seq, function(x) dlnorm(x, log(mu_tau_S), sigma_tau_S))
  
  c1 <- "slategray4"
  c1t <- transparent(c1, trans.val = 0.6)
  
  if(save_plot) {
    png(filename=filename, width=6, height=8, units="in", res=200, type="cairo-png")
  } else dev.new(width=6,height=8)
  
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
  
  if(save_plot) dev.off()
}
