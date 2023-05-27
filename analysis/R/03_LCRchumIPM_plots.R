#==========================================================================#
# FUNCTIONS TO MAKE CUSTOM PLOTS FOR LCR CHUM IPM
#==========================================================================#

#--------------------------------------------------------------------
# Life-cycle multiplot
# S-R curves (spawners to smolts)
# Posterior distributions of fecundity, survival and Mmax parameters
# Time series of smolt productivity process errors and SAR
#--------------------------------------------------------------------

# Plot function
LCRchumIPM_multiplot <- function(mod, SR_fun, fish_data)
{
  # non-hatchery populations and cases
  which_W_pop <- which(!grepl("Hatchery", levels(fish_data$pop)))
  which_W_obs <- which(!grepl("Hatchery", fish_data$pop))
  W_pop <- droplevels(fish_data$pop[which_W_obs])
  W_year <- fish_data$year[which_W_obs]
  # fecundity
  mu_E <- extract1(mod, "mu_E")
  ages <- substring(names(select(fish_data, starts_with("n_age"))), 6, 6)
  delta_NG <- extract1(mod, "delta_NG")
  # egg deposition and egg-to-smolt survival
  A <- fish_data$A[which_W_obs]
  q <- extract1(mod, "q")[,which_W_obs,]
  q_F <- extract1(mod, "q_F")[,which_W_obs]
  mu_psi <- extract1(mod, "mu_psi")
  psi <- extract1(mod, "psi")[,which_W_pop]
  alpha <- apply(sweep(q, c(1,3), mu_E, "*"), 1:2, sum) * q_F * psi[,as.numeric(W_pop)]
  alpha <- t(as.matrix(aggregate(t(alpha), list(pop = W_pop), median)[,-1]))
  mu_alpha <- rowMeans(log(alpha))
  mu_Mmax <- as.vector(extract1(mod, "mu_Mmax"))
  Mmax <- extract1(mod, "Mmax")[,which_W_pop]
  S <- colMedians(extract1(mod, "S"))[which_W_obs]
  SA_grid <- matrix(seq(0, quantile(S/A, 0.9, na.rm = TRUE), length = 100),
                    nrow = length(mu_alpha), ncol = 100, byrow = TRUE)
  M_ESU <- SR(SR_fun, alpha = exp(mu_alpha), Rmax = exp(mu_Mmax), S = SA_grid)/1000 # mil/km
  M_pop <- sapply(1:ncol(Mmax), function(i) {
    colMedians(SR(SR_fun, alpha = alpha[,i], Rmax = Mmax[,i], S = SA_grid))/1000 # mil/km
  })
  # smolt recruitment process errors
  y <- sort(unique(fish_data$year))
  eta_year_M <- extract1(mod, "eta_year_M")
  sigma_M <- extract1(mod, "sigma_M")
  zeta_M <- stan_mean(mod, "zeta_M")
  epsilon_M <- outer(sigma_M, zeta_M, "*")
  error_M <- eta_year_M[,as.numeric(factor(fish_data$year))] + epsilon_M
  error_M <- error_M[,which_W_obs]
  # SAR
  eta_year_MS <- extract1(mod, "eta_year_MS")
  mu_MS <- extract1(mod, "mu_MS")
  s_hat_MS <- plogis(sweep(eta_year_MS, 1, qlogis(mu_MS), "+"))
  s_MS <- extract1(mod, "s_MS")[,which_W_obs]
  # colors
  c1 <- "slategray4"
  c1t <- transparent(c1, trans.val = 0.3)
  c1tt <- transparent(c1, trans.val = 0.5)
  ac <- viridis(length(ages), end = 0.8, direction = -1, alpha = 0.5) 
  
  par(mfrow = c(2,4), mar = c(5.1,5.1,1,0.5), oma = c(0,0,0,1))
  
  # Posterior distributions of fecundity by age
  dd_age <- lapply(as.data.frame(mu_E), density)
  
  plot(dd_age[[1]]$x, dd_age[[1]]$y, pch = "", las = 1, cex.axis = 1.2, cex.lab = 1.5,
       xlab = bquote("Mean fecundity (" * mu[italic(E)] * ")"), ylab = "", 
       xlim = range(sapply(dd_age, function(m) m$x)),
       ylim = range(sapply(dd_age, function(m) m$y))*1.02,
       xaxs = "i", yaxs = "i", yaxt = "n")
  for(a in 1:length(ages))
    polygon(dd_age[[a]], col = ac[a], border = NA)
  title(ylab = "Probability density", line = 2, cex.lab = 1.5)
  text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "A", cex = 1.5)
  legend("topright", paste("age", ages, "  "), cex = 1.2, fill = ac, border = NA,
         bty = "n", inset = c(-0.05,0))
  
  # Posterior distribution of non-green female fecundity discount
  dd_NG <- density(delta_NG)
  
  plot(dd_NG, type = "n", las = 1, cex.axis = 1.2, cex.lab = 1.5,
       xlab = bquote("Non-green fecundity discount (" * delta[NG] * ")"), ylab = "", main = "", 
       xlim = c(0,1), ylim = range(dd_NG$y)*1.02, xaxs = "i", yaxs = "i", yaxt = "n", xpd = NA)
  polygon(dd_NG, col = c1tt, border = NA)
  text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "B", cex = 1.5)
  
  # Posterior densities of psi
  dd_ESU <- density(mu_psi)
  dd_pop <- lapply(as.data.frame(psi), density)
  
  plot(dd_ESU$x, dd_ESU$y, pch = "", lwd = 3, col = c1, las = 1, 
       xaxs = "i", yaxs = "i", yaxt = "n", cex.axis = 1.2, cex.lab = 1.5, 
       xlab = bquote("Maximum egg-to-smolt survival (" * psi * ")"), ylab = "",
       xlim = c(0,1), ylim = range(dd_ESU$y, sapply(dd_pop, function(m) m$y))*1.02, xpd = NA)
  polygon(dd_ESU, col = c1tt, border = NA)
  for(i in 1:length(dd_pop))
    lines(dd_pop[[i]]$x, dd_pop[[i]]$y, col = c1t)
  text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "C", cex = 1.5)
  
  # Posterior densities of log(Mmax)
  dd_ESU <- density(mu_Mmax * log10(exp(1)) - 3)  # convert to base 10, units of mil/km 
  dd_pop <- lapply(as.data.frame(log10(Mmax) - 3), density)
  
  plot(dd_ESU$x, dd_ESU$y, pch = "", lwd = 3, col = c1, las = 1, 
       xaxt = "n", yaxt = "n", yaxs = "i", cex.axis = 1.2, cex.lab = 1.5, 
       xlab = bquote("Maximum smolt density" ~ "(" * italic(M)[max] * ")"), ylab = "",
       xlim = range(dd_ESU$x[dd_ESU$y > 0.02], sapply(dd_pop, function(m) range(m$x[m$y > 0.02]))),
       ylim = range(dd_ESU$y, sapply(dd_pop, function(m) m$y))* 1.02, xpd = NA)
  polygon(dd_ESU, col = c1tt, border = NA)
  for(i in 1:length(dd_pop))
    lines(dd_pop[[i]]$x, dd_pop[[i]]$y, col = c1t)
  tck <- maglab(10^par("usr")[1:2], log = TRUE)
  axis(1, at = log10(tck$labat), labels = tck$labat, cex.axis = 1.2, hadj = 0.5)
  text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "D", cex = 1.5)
  
  # Spawner-recruit function
  plot(SA_grid[1,], colMedians(M_ESU), type = "l", lwd=3, col = c1, las = 1,
       cex.axis = 1.2, cex.lab = 1.5, xaxs = "i", yaxs = "i",
       ylim = range(M_pop), xlab = bquote("Spawner density (" * km^-1 * ")"), 
       ylab = bquote("Smolt density (" * 10^6 ~ km^-1 * ")"))
  for(i in 1:ncol(M_pop)) 
    lines(SA_grid[1,], M_pop[,i], col = c1t)
  polygon(c(SA_grid[1,], rev(SA_grid[1,])), 
          c(colQuantiles(M_ESU, probs = 0.05), rev(colQuantiles(M_ESU, probs = 0.95))), 
          col = c1tt, border = NA)
  rug(S/A, col = c1t)
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
  rug(S/A, col = c1t)
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
  for(j in levels(W_pop))
    lines(W_year[W_pop == j], colMedians(error_M[,W_pop == j]), col = c1t)
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
  for(j in levels(W_pop))
    lines(W_year[W_pop == j], colMedians(s_MS[,W_pop == j])*100, col = c1t)
  axis(side = 1, at = y[y %% 5 == 0], cex.axis = 1.2)
  rug(y[y %% 5 != 0], ticksize = -0.02)
  text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "H", cex = 1.5)
}

#--------------------------------------------------------------------------------
# Time series of SAR for natural populations and hatcheries
#--------------------------------------------------------------------------------

LCRchumIPM_SAR_timeseries <- function(mod, fish_data)
{
  logit <- rfun(qlogis)
  ilogit <- rfun(plogis)
  
  draws <- as.matrix(mod, c("s_MS","mu_MS","beta_MS","eta_year_MS")) %>% as_draws_rvars() %>% 
    mutate_variables(s_hat_MS_N = ilogit(logit(mu_MS) + eta_year_MS),
                     s_hat_MS_H = ilogit(logit(mu_MS) + eta_year_MS + beta_MS[1]))
  
  hyper <- data.frame(year = sort(unique(fish_data$year)), 
                      s_hat_MS_N = draws$s_hat_MS_N, s_hat_MS_H = draws$s_hat_MS_H)
  
  gg <- fish_data %>% cbind(s_MS = draws$s_MS) %>% 
    ggplot(aes(x = year, y = median(s_MS), group = pop, color = pop_type)) +
    geom_line(size = 0.7) + 
    geom_line(aes(x = year, y = median(s_hat_MS_N)), inherit.aes = FALSE,
              data = hyper, size = 1.5, col = alpha("slategray4", 0.8)) +
    geom_ribbon(aes(x = year, ymin = as.vector(quantile(s_hat_MS_N, 0.05)), 
                    ymax = as.vector(quantile(s_hat_MS_N, 0.95))),
                inherit.aes = FALSE, data = hyper, fill = alpha("slategray4", 0.3)) +
    geom_line(aes(x = year, y = median(s_hat_MS_H)), inherit.aes = FALSE,
              data = hyper, size = 1.5, col = alpha("salmon", 0.8)) +
    geom_ribbon(aes(x = year, ymin = as.vector(quantile(s_hat_MS_H, 0.05)), 
                    ymax = as.vector(quantile(s_hat_MS_H, 0.95))),
                inherit.aes = FALSE, data = hyper, fill = alpha("salmon", 0.3)) +
    scale_y_continuous(labels = function(x) x*100) +
    scale_color_discrete(type = c(natural = alpha("slategray4", 0.5), 
                                  hatchery = alpha("salmon", 0.6))) +
    labs(x = "Year", y = "SAR (%)", color = "Origin") +
    theme_bw(base_size = 16) +
    theme(panel.grid.minor = element_blank(), legend.position = c(0.86, 0.91))
  
  return(gg)
}


#--------------------------------------------------------------------------------
# Spawner-to-smolt S-R plot with fit, states, and observations for each pop
#--------------------------------------------------------------------------------

LCRchumIPM_SR_plot <- function(mod, SR_fun, life_stage, fish_data)
{
  n_grid <- 100
  cl <- 0.8
  qnt <- c((1 - cl)/2, 1 - (1 - cl)/2)
  
  logit <- function(x) log(x) - log(1 - x)
  ilogit <- function(x) exp(x) / (1 + exp(x))
  
  # S-R parameters, states and observations including reconstructed recruits
  draws <- as.matrix(mod, c("mu_E", "q", "q_F", "psi", "S", "M", "s_MS")) %>% 
    as_draws_rvars() %>% 
    mutate_variables(alpha = (q %**% mu_E) * q_F * psi[fish_data$pop], 
                     R = M * s_MS,
                     N = switch(life_stage, M = M, R = R))
  
  states_obs <- run_recon(fish_data) %>% 
    mutate(N_obs = switch(life_stage, M = fish_data$M_obs, R = R_obs)) %>% 
    select(pop, year, A, S_obs, N_obs) %>% 
    data.frame(alpha = draws$alpha, S = draws$S, N = draws$N) %>% 
    group_by(pop) %>% 
    mutate(S_q2 = as.vector(quantile(S, qnt[2])), S_upper = pmin(S_q2, max(median(S))*1.1),
           N_q2 = as.vector(quantile(N, qnt[2])), N_upper = pmin(N_q2, max(median(N)*1.1)))
  
  # spawner densities at which to evaluate S-R function
  S_grid <- states_obs %>% group_by(pop) %>% 
    summarize(A = mean(A), alpha = rvar_mean(alpha), 
              S = seq(0, max(S_upper, S_obs, na.rm = TRUE), length = n_grid))
  
  # posteriors of S-R fit with total process and proc + obs error (PPD)
  ppd <- as.matrix(mod, c("Mmax", "sigma_year_M", "rho_M", "sigma_M", "tau_M", "M",
                          "mu_MS", "sigma_year_MS", "rho_MS", "sigma_MS", "tau_S")) %>% 
    as_draws_rvars() %>% 
    mutate_variables(A = as_rvar(S_grid$A), S = S_grid$S,
                     alpha = S_grid$alpha, Mmax = rep(Mmax, each = n_grid),
                     M_hat = SR(SR_fun = SR_fun, alpha = alpha, Rmax = Mmax, S = S, A = A),
                     sd_year_M = sigma_year_M / sqrt(1 - rho_M^2),
                     sd_proc_M = sqrt(sd_year_M^2 + sigma_M^2),
                     sd_ppd_M = sqrt(sd_proc_M^2 + tau_M^2),
                     M_proc = rvar_rng(rlnorm, length(M_hat), log(M_hat), sd_proc_M),
                     M_ppd = rvar_rng(rlnorm, length(M_hat), log(M_hat), sd_ppd_M),
                     R_hat = M_hat * mu_MS,
                     sd_year_MS = sigma_year_MS / sqrt(1 - rho_MS^2),
                     sd_proc_MS = sqrt(sd_year_MS^2 + sigma_MS^2),
                     MS_proc = ilogit(rvar_rng(rnorm, length(M_hat), logit(mu_MS), sd_proc_M)),
                     R_proc = M_proc * MS_proc,
                     R_ppd = rvar_rng(rlnorm, length(R_proc), log(R_proc), tau_S))
  
  ppd <- S_grid %>% select(pop, S) %>% 
    data.frame(N_hat = switch(life_stage, M = ppd$M_hat, R = ppd$R_hat),
               N_proc = switch(life_stage, M = ppd$M_proc, R = ppd$R_proc),
               N_ppd = switch(life_stage, M = ppd$M_ppd, R = ppd$R_ppd)) 
  
  gg <- ppd %>% ggplot(aes(x = S, y = median(N_hat))) +
    geom_ribbon(aes(ymin = quantile(N_hat, qnt[1]), ymax = quantile(N_hat, qnt[2])), 
                fill = "slategray4", alpha = 0.4) +
    geom_ribbon(aes(ymin = quantile(N_proc, qnt[1]), ymax = quantile(N_proc, qnt[2])),
                fill = "slategray4", alpha = 0.3) +
    # geom_ribbon(aes(ymin = quantile(N_ppd, qnt[1]), ymax = quantile(N_ppd, qnt[2])),
    #             fill = "slategray4", alpha = 0.2) +
    geom_line(lwd = 1, col = "slategray4") +
    geom_segment(aes(x = quantile(S, qnt[1]), xend = S_upper,
                     y = median(N), yend = median(N)),
                 data = states_obs, col = "slategray4", alpha = 0.8) +
    geom_segment(aes(x = median(S), xend = median(S),
                     y = quantile(N, qnt[1]), yend = N_upper),
                 data = states_obs, col = "slategray4", alpha = 0.8) +
    # geom_segment(aes(x = S_obs, xend = median(S), y = N_obs, yend = median(N)),
    #              data = states_obs, col = "slategray4", alpha = 0.4) +
    geom_point(aes(x = median(S), y = median(N)), data = states_obs,
               pch = 21, size = 2, col = "slategray4", fill = "white") +
    geom_point(aes(x = S_obs, y = N_obs), data = states_obs,
               pch = 16, size = 2, alpha = 0.6) +
    scale_x_continuous(labels = label_number(scale = 1e-3), expand = c(0,0)) +
    scale_y_continuous(labels = label_number(scale = switch(life_stage, M = 1e-6, R = 1e-3)),
                       expand = c(0,0)) +
    labs(x = "Spawners (thousands)", 
         y = switch(life_stage, M = "Smolts (millions)", R = "Recruits (thousands)")) +
    facet_wrap(vars(pop), ncol = 4, scales = "free") + theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11),
          panel.grid = element_blank(), strip.background = element_rect(fill = NA),
          strip.text = element_text(margin = margin(b = 3, t = 3)))
  
  return(gg)
}


#--------------------------------------------------------------------------------
# Time series of observed and fitted total spawners or smolts for each pop
#--------------------------------------------------------------------------------

LCRchumIPM_MS_timeseries <- function(mod, life_stage = c("M","S"), fish_data)
{
  life_cycle <- strsplit(mod@model_name, "_")[[1]][2]
  N <- extract1(mod, life_stage)
  tau <- extract1(mod, switch(life_cycle, SS = "tau", LCRchum = paste0("tau_", life_stage)))
  if(life_cycle == "LCRchum" & life_stage == "M")
    N[,na.omit(fish_data$downstream_trap)] <- N[,na.omit(fish_data$downstream_trap)] + 
                                              N[,which(!is.na(fish_data$downstream_trap))]
  N_ppd <- N * rlnorm(length(N), 0, tau)
  year <- fish_data$year
  
  gg <- fish_data %>% 
    mutate(N_obs = !!sym(paste0(life_stage, "_obs")),
           tau_obs = !!sym(paste0("tau_", life_stage, "_obs")),
           N_obs_L = qlnorm(0.05, log(N_obs), tau_obs),
           N_obs_U = qlnorm(0.95, log(N_obs), tau_obs),
           N_L = colQuantiles(N, probs = 0.05),
           N_m = colMedians(N),
           N_U = colQuantiles(N, probs = 0.95),
           N_ppd_L = colQuantiles(N_ppd, probs = 0.05),
           N_ppd_U = colQuantiles(N_ppd, probs = 0.95),
           pch = ifelse(is.na(tau_obs), 1, 16)) %>% 
    filter(!grepl("Hatchery", pop)) %>% 
    ggplot(aes(x = year, y = N_obs)) +
    geom_ribbon(aes(ymin = N_L, ymax = N_U), fill = "slategray4", alpha = 0.5) +
    geom_ribbon(aes(ymin = N_ppd_L, ymax = N_ppd_U), fill = "slategray4", alpha = 0.3) +
    geom_line(aes(y = N_m), lwd = 1, col = "slategray4") +
    geom_point(aes(shape = pch), size = 2.5) + scale_shape_identity() +
    geom_errorbar(aes(ymin = N_obs_L, ymax = N_obs_U), width = 0) +
    labs(x = "Year", y = switch(life_stage, M = "Smolts (thousands)", S = "Spawners")) + 
    scale_x_continuous(breaks = round(seq(min(year), max(year), by = 5)[-1]/5)*5) +
    scale_y_log10(labels = function(y) y*switch(life_stage, M = 1e-3, S = 1)) + 
    facet_wrap(vars(pop), ncol = 4) + theme_bw(base_size = 16) +
    theme(panel.grid.minor = element_blank(), 
          strip.background = element_rect(fill = NA),
          strip.text = element_text(margin = margin(b = 3, t = 3)))
  
  return(gg)
}

#--------------------------------------------------------------------------------
# Time series of observed and fitted spawner age structure for each pop
#--------------------------------------------------------------------------------

LCRchumIPM_age_timeseries <- function(mod, fish_data)
{
  q <- extract1(mod, "q")
  year <- fish_data$year
  
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
    filter(!grepl("Hatchery", pop)) %>% 
    ggplot(aes(x = year, group = age, color = age, fill = age)) +
    geom_line(aes(y = q_age_m), lwd = 1, alpha = 0.8) +
    geom_ribbon(aes(ymin = q_age_L, ymax = q_age_U), color = NA, alpha = 0.3) +
    geom_point(aes(y = PointEst), pch = 16, size = 2.5, alpha = 0.8) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0, alpha = 0.8) +
    scale_color_manual(values = viridis(3, end = 0.8, direction = -1)) +
    scale_fill_manual(values = viridis(3, end = 0.8, direction = -1)) +
    scale_x_continuous(breaks = round(seq(min(year), max(year), by = 5)[-1]/5)*5) +
    labs(x = "Year", y = "Proportion at age") + 
    facet_wrap(vars(pop), ncol = 4) + theme_bw(base_size = 16) + 
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
          strip.background = element_rect(fill = NA),
          strip.text = element_text(margin = margin(b = 3, t = 3)), 
          legend.box.margin = margin(0,-10,0,-15))
  
  return(gg)  
}

#--------------------------------------------------------------------------------
# Time series of observed and fitted sex ratio for each pop
#--------------------------------------------------------------------------------

LCRchumIPM_sex_timeseries <- function(mod, fish_data)
{
  q_F <- as_draws_rvars(as.matrix(mod, "q_F"))
  year <- fish_data$year
  
  gg <- cbind(fish_data, q_F = q_F) %>%
    mutate(n_MF_obs = n_M_obs + n_F_obs) %>% 
    cbind(., with(., binconf(x = n_F_obs, n = n_MF_obs))) %>%
    filter(!grepl("Hatchery", pop)) %>% 
    ggplot(aes(x = year, y = PointEst, ymin = Lower, ymax = Upper)) + 
    geom_abline(intercept = 0.5, slope = 0, color = "gray") + 
    geom_ribbon(aes(ymin = t(quantile(q_F, 0.05)), ymax = t(quantile(q_F, 0.95))), 
                fill = "slategray4", alpha = 0.5) +
    geom_line(aes(y = median(q_F)), col = "slategray4", lwd = 1) +
    geom_point(pch = 16, size = 2.5) + geom_errorbar(width = 0) +
    scale_x_continuous(breaks = round(seq(min(year), max(year), by = 5)[-1]/5)*5) +
    labs(x = "Year", y = "Proportion female") +
    facet_wrap(vars(pop), ncol = 4) + theme_bw(base_size = 16) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
          strip.background = element_rect(fill = NA),
          strip.text = element_text(margin = margin(b = 3, t = 3)))
  
  return(gg)
}

#--------------------------------------------------------------------------------
# Time series of observed and fitted p_HOS for each pop
#--------------------------------------------------------------------------------

LCRchumIPM_p_HOS_timeseries <- function(mod, fish_data)
{
  p_HOS <- as_draws_rvars(as.matrix(mod, "p_HOS"))
  year <- fish_data$year
  
  gg <- fish_data %>% cbind(p_HOS = p_HOS) %>% 
    mutate(p_HOS_obs = binconf(n_H_obs, n_H_obs + n_W_obs, alpha = 0.1)) %>% 
    do.call(data.frame, .) %>% # unpack col with nested data frame
    filter(!grepl("Hatchery", pop)) %>% 
    ggplot(aes(x = year)) +
    geom_ribbon(aes(ymin = t(quantile(p_HOS, 0.05)), ymax = t(quantile(p_HOS, 0.95))), 
                fill = "slategray4", alpha = 0.5) +
    geom_line(aes(y = median(p_HOS)), col = "slategray4", lwd = 1) +
    geom_point(aes(y = p_HOS_obs.PointEst), pch = 16, size = 2.5) +
    geom_errorbar(aes(ymin = p_HOS_obs.Lower, ymax = p_HOS_obs.Upper), width = 0) +
    scale_x_continuous(breaks = round(seq(min(year), max(year), by = 5)[-1]/5)*5) +
    labs(x = "Year", y = bquote(italic(p)[HOS])) +
    facet_wrap(vars(pop), ncol = 4) + theme_bw(base_size = 16) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = NA),
          strip.text = element_text(margin = margin(b = 3, t = 3)))
  
  return(gg)
}

#--------------------------------------------------------------------------------
# Straying matrix: probability of straying from each origin to each population
#--------------------------------------------------------------------------------

LCRchumIPM_p_origin <- function(mod, fish_data)
{
  p_origin <- as_draws_rvars(as.array(mod, "p_origin"))
  
  gg <- data.frame(p_origin = p_origin) %>% 
    setNames(unique(fish_data$pop[fish_data$pop_type == "natural"])) %>% 
    cbind(origin = unique(fish_data$pop[fish_data$pop_type == "hatchery"])) %>% 
    pivot_longer(cols = -origin, names_to = "pop", values_to = "p_origin") %>% 
    mutate(pop = factor(pop, levels = levels(fish_data$pop))) %>% 
    ggplot(aes(xdist = p_origin, y = pop)) +
    ggdist::stat_eye(.width = c(0.5, 0.9), normalize = "groups") + 
    scale_y_discrete(limits = rev) + labs(x = "Straying rate", y = "") + 
    facet_wrap(vars(origin)) + theme_bw(base_size = 16) +
    theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = NA),
          strip.text = element_text(margin = margin(b = 3, t = 3)))
  
  return(gg)
}

#--------------------------------------------------------------------------------
# Distributions of observed and fitted fecundity by age
#--------------------------------------------------------------------------------

LCRchumIPM_fecundity_plot <- function(mod, fish_data, fecundity_data)
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
         labels = paste("age", ages[a]), cex = 1.8, col = c1[a], adj = 1)
  }
  title(xlab = "Fecundity", ylab = "Probability density", cex.lab = 1.9, line = 0, outer = TRUE)
}


#--------------------------------------------------------------------------------
# Observed and fitted distributions of "known" smolt and spawner 
# observation error SDs
#--------------------------------------------------------------------------------

LCRchumIPM_obs_error_plot <- function(mod, fish_data)
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
  
  par(mfcol = c(2,1), mar = c(5,5,0,1),  oma = c(0,0,1,0))
  
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
}
