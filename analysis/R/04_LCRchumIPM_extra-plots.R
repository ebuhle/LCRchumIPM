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


