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
multiplot <- function(mod, SR_fun, fish_data)
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
       xlab = bquote("Mean fecundity (" * mu[italic(E)] * ")"), ylab = NA, 
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
       xlab = bquote("Non-green fecundity discount (" * delta[NG] * ")"), ylab = NA, main = NA, 
       xlim = c(0,1), ylim = range(dd_NG$y)*1.02, xaxs = "i", yaxs = "i", yaxt = "n", xpd = NA)
  polygon(dd_NG, col = c1tt, border = NA)
  text(par("usr")[1], par("usr")[4], adj = c(-1,1.5), "B", cex = 1.5)
  
  # Posterior densities of psi
  dd_ESU <- density(mu_psi)
  dd_pop <- lapply(as.data.frame(psi), density)
  
  plot(dd_ESU$x, dd_ESU$y, pch = "", lwd = 3, col = c1, las = 1, 
       xaxs = "i", yaxs = "i", yaxt = "n", cex.axis = 1.2, cex.lab = 1.5, 
       xlab = bquote("Maximum egg-to-smolt survival (" * psi * ")"), ylab = NA,
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
       xlab = bquote("Maximum smolt density" ~ "(" * italic(M)[max] * ")"), ylab = NA,
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
       ylim = range(M_pop), xlab = bquote("Spawner density (" * 10^3 ~ km^-1 * ")"), 
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
       xlab = bquote("Spawner density (" * 10^3 ~ km^-1 * ")"), 
       ylab = bquote("Smolts per spawner (" * 10^3 * ")"))
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

SAR_timeseries <- function(mod, fish_data)
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
    geom_line(linewidth = 0.7) + 
    geom_line(aes(x = year, y = median(s_hat_MS_N)), inherit.aes = FALSE,
              data = hyper, linewidth = 1.5, col = alpha("slategray4", 0.8)) +
    geom_ribbon(aes(x = year, ymin = as.vector(quantile(s_hat_MS_N, 0.05)), 
                    ymax = as.vector(quantile(s_hat_MS_N, 0.95))),
                inherit.aes = FALSE, data = hyper, fill = alpha("slategray4", 0.3)) +
    geom_line(aes(x = year, y = median(s_hat_MS_H)), inherit.aes = FALSE,
              data = hyper, linewidth = 1.5, col = alpha("salmon", 0.8)) +
    geom_ribbon(aes(x = year, ymin = as.vector(quantile(s_hat_MS_H, 0.05)), 
                    ymax = as.vector(quantile(s_hat_MS_H, 0.95))),
                inherit.aes = FALSE, data = hyper, fill = alpha("salmon", 0.3)) +
    scale_y_continuous(labels = function(x) x*100) +
    scale_color_discrete(type = c(natural = alpha("slategray4", 0.5), 
                                  hatchery = alpha("salmon", 0.6))) +
    labs(x = "Year", y = "SAR (%)", color = "Origin") +
    theme(panel.grid.minor = element_blank(), legend.position = c(0.86, 0.91))
  
  return(gg)
}

#--------------------------------------------------------------------------------
# Spawner-to-smolt S-R plot with fit, states, and observations for each pop
#--------------------------------------------------------------------------------

SR_plot <- function(mod, SR_fun, life_stage, fish_data)
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

MS_timeseries <- function(mod, life_stage = c("M","S"), fish_data)
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
    facet_wrap(vars(pop), ncol = 4) + 
    theme(panel.grid.minor = element_blank(), 
          strip.background = element_rect(fill = NA),
          strip.text = element_text(margin = margin(b = 3, t = 3)))
  
  return(gg)
}

#--------------------------------------------------------------------------------
# Time series of observed and fitted spawner age structure for each pop
#--------------------------------------------------------------------------------

age_timeseries <- function(mod, fish_data)
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
    facet_wrap(vars(pop), ncol = 4) + 
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
          strip.background = element_rect(fill = NA),
          strip.text = element_text(margin = margin(b = 3, t = 3)), 
          legend.box.margin = margin(0,-10,0,-15))
  
  return(gg)  
}

#--------------------------------------------------------------------------------
# Time series of observed and fitted sex ratio for each pop
#--------------------------------------------------------------------------------

sex_timeseries <- function(mod, fish_data)
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
    facet_wrap(vars(pop), ncol = 4) + 
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
          strip.background = element_rect(fill = NA),
          strip.text = element_text(margin = margin(b = 3, t = 3)))
  
  return(gg)
}

#--------------------------------------------------------------------------------
# Time series of observed and fitted p_HOS for each pop
#--------------------------------------------------------------------------------

p_HOS_timeseries <- function(mod, fish_data)
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
    facet_wrap(vars(pop), ncol = 4) + 
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = NA),
          strip.text = element_text(margin = margin(b = 3, t = 3)))
  
  return(gg)
}

#--------------------------------------------------------------------------------
# Straying matrix: probability of straying from each origin to each population
#--------------------------------------------------------------------------------

p_origin_plot <- function(mod, fish_data)
{
  p_origin <- as_draws_rvars(as.array(mod, "p_origin"))
  
  gg <- data.frame(p_origin = p_origin) %>% 
    setNames(unique(fish_data$pop[fish_data$pop_type == "natural"])) %>% 
    cbind(origin = unique(fish_data$pop[fish_data$pop_type == "hatchery"])) %>% 
    pivot_longer(cols = -origin, names_to = "pop", values_to = "p_origin") %>% 
    mutate(pop = factor(pop, levels = levels(fish_data$pop))) %>% 
    ggplot(aes(xdist = p_origin, y = pop)) +
    stat_eye(.width = c(0.5, 0.9), normalize = "groups", 
             color = "slategray4", fill = alpha("slategray4", 0.5)) + 
    scale_y_discrete(limits = rev) + labs(x = "Straying rate", y = "") + 
    facet_wrap(vars(origin)) + 
    theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = NA),
          strip.text = element_text(margin = margin(b = 3, t = 3)))
  
  return(gg)
}

#--------------------------------------------------------------------------------
# Distributions of observed and fitted fecundity by age
#--------------------------------------------------------------------------------

fecundity_plot <- function(mod, fish_data, fecundity_data)
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
         xlab = NA, ylab = NA, main = NA, xaxs = "i", yaxt = "n", bty = "n")
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

obs_error_plot <- function(mod, fish_data)
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
       ylab = "Probability density", main = NA)
  polygon(c(tau_M_seq, rev(tau_M_seq)),
          c(colQuantiles(tau_M_fit, probs = 0.05), rev(colQuantiles(tau_M_fit, probs = 0.95))),
          col = c1t, border = NA)
  lines(tau_M_seq, colMedians(tau_M_fit), col = c1, lwd = 3)
  
  # spawner observation error SD
  hist(tau_S_obs, 20, prob = TRUE, las = 1, cex.axis = 1.2, cex.lab = 1.5, 
       col = "lightgray", border = "white",
       ylim = c(0, max(colQuantiles(tau_S_fit, probs = 0.95))),
       xlab = bquote("Spawner observation error (" * tau[italic(S)] * ")"), 
       ylab = "Probability density", main = NA)
  polygon(c(tau_S_seq, rev(tau_S_seq)),
          c(colQuantiles(tau_S_fit, probs = 0.05), rev(colQuantiles(tau_S_fit, probs = 0.95))),
          col = c1t, border = NA)
  lines(tau_S_seq, colMedians(tau_S_fit), col = c1, lwd = 3)
}

#--------------------------------------------------------------------------------
# Conditioning forecast trajectories on time-averaged SAR anomalies
#--------------------------------------------------------------------------------

SAR_fore_plot <- function(mod, fish_data_fore, example_pop)
{
  dat <- fish_data_fore %>% 
    mutate(yr = as.numeric(factor(year))) %>% subset(pop == example_pop)
  
  rfindInterval <- rfun(findInterval)
  draws <- as_draws_rvars(as.matrix(mod, c("eta_year_MS", "S"))) %>% 
    mutate_variables(S = S[as.numeric(rownames(dat))], 
                     gmean_S = exp(rvar_mean(log(S[dat$forecast]))),
                     eta_year_MS = eta_year_MS[dat$yr], 
                     eta_mean_MS = rvar_mean(eta_year_MS[dat$forecast]),
                     q_MS = quantile(eta_mean_MS, (1:2)/3),
                     qnt_MS = rfindInterval(eta_mean_MS, q_MS) + 1)
  rdraws <- draws %>% resample_draws(ndraws = 100)
  
  layout(matrix(c(1:4), nrow = 2, byrow = TRUE), widths = c(7, 1, 7, 1))
  par(oma = c(2, 0, 0, 0))
  cols <- c(low = "firebrick1", med = "gold", high = "darkgreen")

  # eta_year_MS time series
  par(mar = c(2, 5, 1, 0.1))
  with(dat, {
    plot(min(year), 0, type = "n", 
         xlim = range(year), ylim = range(range(rdraws$eta_year_MS)),
         xlab = NA, ylab = quote("SAR anomaly" ~ (eta[italic(MS)]^year)), 
         xaxs = "i", las = 1, cex.axis = 1.2, cex.lab = 1.5)
    rug(year[year %% 10 != 0], ticksize = -0.01)
    for_each_draw(rdraws, 
                  segments(year[-nrow(dat)], eta_year_MS[-nrow(dat)], year[-1], eta_year_MS[-1],
                           col = alpha(ifelse(forecast, cols[qnt_MS], "slategray4"), 0.3)))
    legend("topleft", title = "Mean future SAR anomaly", legend = "", cex = 1.2, bty = "n")
    legend("topleft", title = "", legend = rev(names(cols)),
           cex = 1.2, fill = alpha(rev(cols), 0.5), border = NA, bty = "n")
  })
  
  # distribution of eta_mean_MS
  d <- density(draws_of(draws$eta_mean_MS))
  q_MS <- mean(draws$q_MS)
  xout <- sort(c(d$x, q_MS))
  d <- approx(d$x, d$y, xout = xout)
  xqnt <- findInterval(d$x, q_MS) + 1
  
  par(mar = c(2, 0, 1, 0.5))
  plot(d$y, d$x, type = "n", ylim = par("usr")[3:4],
       xaxt = "n", yaxt = "n", bty = "n", xlab = NA, ylab = NA)
  with(d, 
       for(q in 1:3)
         polygon(c(y[xqnt == q], rep(0, sum(xqnt == q))), c(x[xqnt == q], rev(x[xqnt == q])), 
                 border = NA, col = alpha(cols[q], 0.5)))
  
  # S time series for example_pop
  par(mar = c(3, 5, 1, 0.1))
  with(dat, {
    plot(year, log10(S_obs), type = "n", ylim = c(0, log10(max(range(rdraws$S)))),
         xaxs = "i", yaxt = "n", las = 1, cex.axis = 1.2, cex.lab = 1.5,
         xlab = NA, ylab = "Spawners")
    mtext("Year", side = 1, outer = TRUE, cex = par("cex")*1.5)
    tck <- maglab(10^par("usr")[3:4], log = TRUE)
    axis(2, at = log10(tck$labat), labels = tck$exp, las = 1, cex.axis = 1.2)
    rug(year[year %% 10 != 0], ticksize = -0.01)
    for_each_draw(rdraws,
                  segments(year[-nrow(dat)], log10(S[-nrow(dat)]), year[-1], log10(S[-1]),
                           col = alpha(ifelse(forecast, cols[qnt_MS], "slategray4"), 0.3)))
    points(year, log10(S_obs), pch = ifelse(is.na(tau_S_obs), 1, 16), cex = 1.5)
    segments(year, log10(qlnorm(0.05, log(S_obs), tau_S_obs)),
             year, log10(qlnorm(0.95, log(S_obs), tau_S_obs)))
  })
  
  # distributions of geometric mean S
  dens <- lapply(1:3, function(q) {
    log10_gmean_S <- draws %>% subset_draws(iter = which(draws_of(.$qnt_MS) == q)) %>% 
      .$gmean_S %>% log10()
    return(density(draws_of(log10_gmean_S)))
  })
  
  par(mar = c(3, 0, 1, 0.5))
  plot(dens[[1]]$y, dens[[1]]$x, type = "n", 
       xlim = c(0, max(sapply(dens, function(d) max(d$y)))), ylim = par("usr")[3:4], 
       xaxt = "n", yaxt = "n", bty = "n", xlab = NA, ylab = NA)
  for(q in 1:3)
    polygon(c(dens[[q]]$y, rep(0, length(dens[[q]]$y))), c(dens[[q]]$x, rev(dens[[q]]$x)), 
            border = NA, col = alpha(cols[q], 0.5))
}

#--------------------------------------------------------------------------------
# Distributions of forecast spawner abundance under alternative scenarios
#--------------------------------------------------------------------------------

S_fore_plot <- function(modH0, modHmax, fish_data_foreH0, fish_data_foreHmax, pop_names) 
{
  year <- as.numeric(factor(fish_data_foreH0$year))
  fore_years <- sort(unique(year[fish_data_foreH0$forecast]))
  n_draws <- (sum(modH0@sim$n_save - modH0@sim$warmup) %/% 3) * 3
  
  rfindInterval <- rfun(findInterval)
  drawsH0 <- as_draws_rvars(as.matrix(modH0, c("eta_year_MS","S"))) %>% 
    subset_draws(iter = 1:n_draws) %>% 
    mutate_variables(eta_mean_MS = rvar_mean(eta_year_MS[fore_years]),
                     q_MS = quantile(eta_mean_MS, (1:2)/3),
                     qnt_MS = rfindInterval(eta_mean_MS, q_MS) + 1)
  drawsHmax <- as_draws_rvars(as.matrix(modHmax, c("eta_year_MS","S"))) %>% 
    subset_draws(iter = 1:n_draws) %>% 
    mutate_variables(eta_mean_MS = rvar_mean(eta_year_MS[fore_years]),
                     q_MS = quantile(eta_mean_MS, (1:2)/3),
                     qnt_MS = rfindInterval(eta_mean_MS, q_MS) + 1)
  
  datH0 <- fish_data_foreH0 %>% 
    mutate(scenario = "No", qnt_MS = drawsH0$qnt_MS, S = drawsH0$S)
  datHmax <- fish_data_foreHmax %>% 
    mutate(scenario = "Yes", qnt_MS = drawsHmax$qnt_MS, S = drawsHmax$S)
  
  dat <- rbind(datH0, datHmax) %>% arrange(scenario, pop, year) %>% 
    filter(pop_type == "natural" & forecast) %>% 
    left_join(pop_names) %>% mutate(pop = recovery_pop) %>% # comment out to use model pops instead of recovery pops #
    select(scenario, pop, year, qnt_MS, S) %>% 
    group_by(scenario, pop, year) %>% 
    summarize(qnt_MS = unique(qnt_MS), S = rvar_sum(S)) %>% 
    mutate(S_low = subset_draws(S, iter = which(draws_of(qnt_MS[1]) == 1)),
           S_med = subset_draws(S, iter = which(draws_of(qnt_MS[1]) == 2)),
           S_high = subset_draws(S, iter = which(draws_of(qnt_MS[1]) == 3)),
           S = subset_draws(S, iter = 1:(n_draws/3))) %>%
    rename(S_all = S) %>% 
    pivot_longer(cols = contains("S_"), names_to = "SAR", names_prefix = "S_", values_to = "S") %>%
    mutate(SAR = factor(SAR, levels = c("all","low","med","high"))) %>%
    group_by(scenario, pop, SAR) %>% summarize(gmean_S = exp(rvar_mean(log(S))))
  
  cols <- c(all = "slategray4", low = "firebrick1", med = "gold", high = "darkgreen")
  
  gg <- dat %>% 
    ggplot(aes(x = scenario, ydist = gmean_S, color = SAR, fill = SAR)) +
    stat_eye(.width = c(0.5, 0.9), normalize = "groups", position = "dodge",
             point_size = 2, slab_alpha = 0.5, slab_color = NA) +
    scale_y_log10(labels = function(y) y) + 
    coord_cartesian(ylim = range(quantile(dat$gmean_S, c(0.05, 0.95)))) +
    scale_color_manual(values = cols) + scale_fill_manual(values = cols) +
    labs(x = "Hatchery scenario", y = "Geometric mean spawners") +
    facet_wrap(vars(pop), ncol = 3) + 
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), 
          strip.background = element_rect(fill = NA),
          strip.text = element_text(margin = margin(b = 3, t = 3)))
  
  return(gg)
}

#------------------------------------------------------------------------------------------
# Distributions of forecast final : initial spawner abundance under alternative scenarios
#------------------------------------------------------------------------------------------

StS0_fore_plot <- function(modH0, modHmax, fish_data_foreH0, fish_data_foreHmax, pop_names) 
{
  year <- as.numeric(factor(fish_data_foreH0$year))
  fore_years <- sort(unique(year[fish_data_foreH0$forecast]))
  n_draws <- (sum(modH0@sim$n_save - modH0@sim$warmup) %/% 3) * 3
  
  rfindInterval <- rfun(findInterval)
  drawsH0 <- as_draws_rvars(as.matrix(modH0, c("eta_year_MS","S"))) %>% 
    subset_draws(iter = 1:n_draws) %>% 
    mutate_variables(eta_mean_MS = rvar_mean(eta_year_MS[fore_years]),
                     q_MS = quantile(eta_mean_MS, (1:2)/3),
                     qnt_MS = rfindInterval(eta_mean_MS, q_MS) + 1)
  drawsHmax <- as_draws_rvars(as.matrix(modHmax, c("eta_year_MS","S"))) %>% 
    subset_draws(iter = 1:n_draws) %>% 
    mutate_variables(eta_mean_MS = rvar_mean(eta_year_MS[fore_years]),
                     q_MS = quantile(eta_mean_MS, (1:2)/3),
                     qnt_MS = rfindInterval(eta_mean_MS, q_MS) + 1)
  
  datH0 <- fish_data_foreH0 %>% 
    mutate(scenario = "No", qnt_MS = drawsH0$qnt_MS, S = drawsH0$S)
  datHmax <- fish_data_foreHmax %>% 
    mutate(scenario = "Yes", qnt_MS = drawsHmax$qnt_MS, S = drawsHmax$S)
  
  dat <- rbind(datH0, datHmax) %>% arrange(scenario, pop, year) %>% 
    filter(pop_type == "natural") %>% 
    left_join(pop_names) %>% mutate(pop = recovery_pop) %>% # comment out to use model pops instead of recovery pops #
    select(scenario, pop, year, forecast, qnt_MS, S) %>% 
    group_by(scenario, pop, year) %>% 
    summarize(forecast = any(forecast), qnt_MS = qnt_MS[1], S = rvar_sum(S)) %>% 
    mutate(S_low = subset_draws(S, iter = which(draws_of(qnt_MS[1]) == 1)),
           S_med = subset_draws(S, iter = which(draws_of(qnt_MS[1]) == 2)),
           S_high = subset_draws(S, iter = which(draws_of(qnt_MS[1]) == 3)),
           S = subset_draws(S, iter = 1:(n_draws/3))) %>%
    rename(S_all = S) %>% 
    pivot_longer(cols = contains("S_"), names_to = "SAR", names_prefix = "S_", values_to = "S") %>%
    mutate(SAR = factor(SAR, levels = c("all","low","med","high"))) %>%
    group_by(scenario, pop, SAR) %>% summarize(StS0 = tail(S,1) / tail(S[!forecast], 1))
  
  cols <- c(all = "slategray4", low = "firebrick1", med = "gold", high = "darkgreen")
  
  gg <- dat %>% 
    ggplot(aes(x = scenario, ydist = StS0, color = SAR, fill = SAR)) +
    geom_hline(yintercept = 1) +
    stat_eye(.width = c(0.5, 0.9), normalize = "groups", position = "dodge",
             point_size = 2, slab_alpha = 0.5, slab_color = NA) +
    scale_y_log10(breaks = 10^seq(-4, 4, by = 2), labels = function(y) y) + 
    coord_cartesian(ylim = range(quantile(dat$StS0, c(0.02, 0.98)))) +
    scale_color_manual(values = cols) + scale_fill_manual(values = cols) +
    labs(x = "Hatchery scenario", y = "Final / current spawners") +
    facet_wrap(vars(pop), ncol = 4) + 
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), 
          strip.background = element_rect(fill = NA),
          strip.text = element_text(margin = margin(b = 3, t = 3)))
  
  return(gg)
}

#--------------------------------------------------------------------------------
# Distributions of forecast p_HOS under alternative scenarios
#--------------------------------------------------------------------------------

p_HOS_fore_plot <- function(modH0, modHmax, fish_data_foreH0, fish_data_foreHmax, pop_names) 
{
  year <- as.numeric(factor(fish_data_foreH0$year))
  fore_years <- sort(unique(year[fish_data_foreH0$forecast]))
  n_draws <- (sum(modH0@sim$n_save - modH0@sim$warmup) %/% 3) * 3
  
  rfindInterval <- rfun(findInterval)
  drawsH0 <- as_draws_rvars(as.matrix(modH0, c("eta_year_MS","S","p_HOS"))) %>%
    subset_draws(iter = 1:n_draws) %>% 
    mutate_variables(eta_mean_MS = rvar_mean(eta_year_MS[fore_years]),
                     q_MS = quantile(eta_mean_MS, (1:2)/3),
                     qnt_MS = rfindInterval(eta_mean_MS, q_MS) + 1,
                     S_H = S*p_HOS)
  drawsHmax <- as_draws_rvars(as.matrix(modHmax, c("eta_year_MS","S","p_HOS"))) %>%
    subset_draws(iter = 1:n_draws) %>% 
    mutate_variables(eta_mean_MS = rvar_mean(eta_year_MS[fore_years]),
                     q_MS = quantile(eta_mean_MS, (1:2)/3),
                     qnt_MS = rfindInterval(eta_mean_MS, q_MS) + 1,
                     S_H = S*p_HOS)
  
  datH0 <- fish_data_foreH0 %>% 
    mutate(scenario = "No", qnt_MS = drawsH0$qnt_MS, S = drawsH0$S, S_H = drawsH0$S_H)
  datHmax <- fish_data_foreHmax %>% 
    mutate(scenario = "Yes", qnt_MS = drawsHmax$qnt_MS, S = drawsHmax$S, S_H = drawsHmax$S_H)
  
  dat <- rbind(datH0, datHmax) %>% arrange(scenario, pop, year) %>% 
    filter(pop_type == "natural" & forecast) %>% 
    left_join(pop_names) %>% mutate(pop = recovery_pop) %>% # comment out to use model pops instead of recovery pops #
    select(scenario, pop, year, qnt_MS, S, S_H) %>% 
    group_by(scenario, pop, year) %>% 
    summarize(qnt_MS = unique(qnt_MS), p_HOS = rvar_sum(S_H)/rvar_sum(S)) %>% 
    mutate(p_HOS_low = subset_draws(p_HOS, iter = which(draws_of(qnt_MS[1]) == 1)),
           p_HOS_med = subset_draws(p_HOS, iter = which(draws_of(qnt_MS[1]) == 2)),
           p_HOS_high = subset_draws(p_HOS, iter = which(draws_of(qnt_MS[1]) == 3)),
           p_HOS = subset_draws(p_HOS, iter = 1:(n_draws/3))) %>%
    rename(p_HOS_all = p_HOS) %>% 
    pivot_longer(cols = contains("p_HOS_"), names_to = "SAR", names_prefix = "p_HOS_", 
                 values_to = "p_HOS") %>%
    mutate(SAR = factor(SAR, levels = c("all","low","med","high"))) %>%
    group_by(scenario, pop, SAR) %>% summarize(p_HOS_mean = rvar_mean(p_HOS))
  
  cols <- c(all = "slategray4", low = "firebrick1", med = "gold", high = "darkgreen")
  
  gg <- dat %>% 
    ggplot(aes(x = scenario, ydist = p_HOS_mean, color = SAR, fill = SAR)) +
    stat_eye(.width = c(0.5, 0.9), normalize = "groups", position = "dodge",
             point_size = 2, slab_alpha = 0.5, slab_color = NA) +
    ylim(c(0,1)) + scale_color_manual(values = cols) + scale_fill_manual(values = cols) +
    labs(x = "Hatchery scenario", y = bquote("Mean" ~ italic(p)[HOS])) +
    facet_wrap(vars(pop), ncol = 3) + 
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), 
          strip.background = element_rect(fill = NA),
          strip.text = element_text(margin = margin(b = 3, t = 3)))
  
  return(gg)
}

#--------------------------------------------------------------------------------
# Probability of recovery under alternative scenarios
#--------------------------------------------------------------------------------

Precovery_plot <- function(modH0, modHmax, fish_data_foreH0, fish_data_foreHmax, 
                           pop_names, recovery_targets) 
{
  year <- as.numeric(factor(fish_data_foreH0$year))
  fore_years <- sort(unique(year[fish_data_foreH0$forecast]))
  n_draws <- (sum(modH0@sim$n_save - modH0@sim$warmup) %/% 3) * 3
  
  rfindInterval <- rfun(findInterval)
  drawsH0 <- as_draws_rvars(as.matrix(modH0, c("eta_year_MS","S"))) %>% 
    subset_draws(iter = 1:n_draws) %>% 
    mutate_variables(eta_mean_MS = rvar_mean(eta_year_MS[fore_years]),
                     q_MS = quantile(eta_mean_MS, (1:2)/3),
                     qnt_MS = rfindInterval(eta_mean_MS, q_MS) + 1)
  drawsHmax <- as_draws_rvars(as.matrix(modHmax,  c("eta_year_MS","S"))) %>% 
    subset_draws(iter = 1:n_draws) %>% 
    mutate_variables(eta_mean_MS = rvar_mean(eta_year_MS[fore_years]),
                     q_MS = quantile(eta_mean_MS, (1:2)/3),
                     qnt_MS = rfindInterval(eta_mean_MS, q_MS) + 1)
  
  datH0 <- fish_data_foreH0 %>% 
    mutate(scenario = "No", qnt_MS = drawsH0$qnt_MS, S = drawsH0$S)
  datHmax <- fish_data_foreHmax %>% 
    mutate(scenario = "Yes", qnt_MS = drawsHmax$qnt_MS, S = drawsHmax$S)
  
  dat <- rbind(datH0, datHmax) %>% arrange(scenario, pop, year) %>% 
    filter(pop_type == "natural" & forecast) %>% 
    left_join(pop_names) %>% 
    select(scenario, recovery_pop, year, qnt_MS, S) %>% 
    group_by(scenario, recovery_pop, year) %>% 
    summarize(qnt_MS = unique(qnt_MS), S = rvar_sum(S)) %>% 
    mutate(S_low = subset_draws(S, iter = which(draws_of(qnt_MS[1]) == 1)),
           S_med = subset_draws(S, iter = which(draws_of(qnt_MS[1]) == 2)),
           S_high = subset_draws(S, iter = which(draws_of(qnt_MS[1]) == 3)),
           S = subset_draws(S, iter = 1:(n_draws/3))) %>%
    rename(S_all = S) %>% 
    pivot_longer(cols = contains("S_"), names_to = "SAR", names_prefix = "S_", values_to = "S") %>%
    mutate(SAR = factor(SAR, levels = c("all","low","med","high"))) %>%
    left_join(recovery_targets) %>% 
    group_by(scenario, recovery_pop, SAR) %>% 
    mutate(gmean4_S = c(rep(as_rvar(NA), 3), 
                        rvar_apply(array(4:n()), .margin = 1, 
                                   .f = function(i) rvar_mean(S[(i-3):i]))),
           .after = S) %>% 
    summarize(recovered = rvar_any(gmean4_S > target[1], na.rm = TRUE)) %>% 
    ungroup() %>% 
    mutate(Precovery = binconf(sum(recovered), ndraws(recovered), alpha = 0.1)) %>% 
    do.call(data.frame, .) # unpack col with nested data frame
  
  cols <- c(all = "slategray4", low = "firebrick1", med = "gold", high = "darkgreen")
  
  gg <- dat %>% 
    ggplot(aes(x = scenario, y = Precovery.PointEst, color = SAR, fill = SAR)) +
    geom_col(position = "dodge", color = NA, alpha = 0.5) +
    geom_pointrange(aes(ymin = Precovery.Lower, ymax = Precovery.Upper), linewidth = 1, 
                    position = position_dodge(width = 0.9)) +
    scale_color_manual(values = cols) + scale_fill_manual(values = cols) +
    labs(x = "Hatchery scenario", y = "Probability of recovery") +
    facet_wrap(vars(recovery_pop), ncol = 3) + 
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), 
          strip.background = element_rect(fill = NA),
          strip.text = element_text(margin = margin(b = 3, t = 3)))
  
  return(gg)
}

#--------------------------------------------------------------------------------
# Probability of quasi-extinction under alternative scenarios
#--------------------------------------------------------------------------------

PQE_plot <- function(modH0, modHmax, fish_data_foreH0, fish_data_foreHmax, 
                     pop_names, QET) 
{
  year <- as.numeric(factor(fish_data_foreH0$year))
  fore_years <- sort(unique(year[fish_data_foreH0$forecast]))
  n_draws <- (sum(modH0@sim$n_save - modH0@sim$warmup) %/% 3) * 3
  
  rfindInterval <- rfun(findInterval)
  drawsH0 <- as_draws_rvars(as.matrix(modH0, c("eta_year_MS","S"))) %>% 
    subset_draws(iter = 1:n_draws) %>% 
    mutate_variables(eta_mean_MS = rvar_mean(eta_year_MS[fore_years]),
                     q_MS = quantile(eta_mean_MS, (1:2)/3),
                     qnt_MS = rfindInterval(eta_mean_MS, q_MS) + 1)
  drawsHmax <- as_draws_rvars(as.matrix(modHmax,  c("eta_year_MS","S"))) %>% 
    subset_draws(iter = 1:n_draws) %>% 
    mutate_variables(eta_mean_MS = rvar_mean(eta_year_MS[fore_years]),
                     q_MS = quantile(eta_mean_MS, (1:2)/3),
                     qnt_MS = rfindInterval(eta_mean_MS, q_MS) + 1)
  
  datH0 <- fish_data_foreH0 %>% 
    mutate(scenario = "No", qnt_MS = drawsH0$qnt_MS, S = drawsH0$S)
  datHmax <- fish_data_foreHmax %>% 
    mutate(scenario = "Yes", qnt_MS = drawsHmax$qnt_MS, S = drawsHmax$S)
  
  dat <- rbind(datH0, datHmax) %>% arrange(scenario, pop, year) %>% 
    filter(pop_type == "natural" & forecast) %>% 
    left_join(pop_names) %>% 
    select(scenario, recovery_pop, year, qnt_MS, S) %>% 
    group_by(scenario, recovery_pop, year) %>% 
    summarize(qnt_MS = unique(qnt_MS), S = rvar_sum(S)) %>% 
    mutate(S_low = subset_draws(S, iter = which(draws_of(qnt_MS[1]) == 1)),
           S_med = subset_draws(S, iter = which(draws_of(qnt_MS[1]) == 2)),
           S_high = subset_draws(S, iter = which(draws_of(qnt_MS[1]) == 3)),
           S = subset_draws(S, iter = 1:(n_draws/3))) %>%
    rename(S_all = S) %>% 
    pivot_longer(cols = contains("S_"), names_to = "SAR", names_prefix = "S_", values_to = "S") %>%
    mutate(SAR = factor(SAR, levels = c("all","low","med","high"))) %>%
    group_by(scenario, recovery_pop, SAR) %>% 
    mutate(gmean4_S = c(rep(as_rvar(NA), 3), 
                        rvar_apply(array(4:n()), .margin = 1, 
                                   .f = function(i) rvar_mean(S[(i-3):i]))),
           .after = S) %>% 
    summarize(quasi_extinct = rvar_any(gmean4_S < QET, na.rm = TRUE)) %>% 
    ungroup() %>% 
    mutate(PQE = binconf(sum(quasi_extinct), ndraws(quasi_extinct), alpha = 0.1)) %>% 
    do.call(data.frame, .) # unpack col with nested data frame
  
  cols <- c(all = "slategray4", low = "firebrick1", med = "gold", high = "darkgreen")
  
  gg <- dat %>% 
    ggplot(aes(x = scenario, y = PQE.PointEst, color = SAR, fill = SAR)) +
    geom_col(position = "dodge", color = NA, alpha = 0.5) +
    geom_pointrange(aes(ymin = PQE.Lower, ymax = PQE.Upper), linewidth = 1, 
                    position = position_dodge(width = 0.9)) +
    scale_color_manual(values = cols) + scale_fill_manual(values = cols) +
    labs(x = "Hatchery scenario", y = "Probability of quasi-extinction") +
    facet_wrap(vars(recovery_pop), ncol = 3) + 
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), 
          strip.background = element_rect(fill = NA),
          strip.text = element_text(margin = margin(b = 3, t = 3)))
  
  return(gg)
}

