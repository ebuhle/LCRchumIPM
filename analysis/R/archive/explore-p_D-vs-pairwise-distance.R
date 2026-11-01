d <- stan_data('IPM_LCRchum_pp', fish_data = fish_data, 
               fecundity_data = fecundity_data, ages = list(M = 1))
dr <- as.matrix(fit_Ricker, 'p_D') %>% as_draws_rvars()
dimnames(dr$p_D) <- list(levels(fish_data$pop)[d$which_O_pop],
                         levels(fish_data$pop)[d$which_D_pop])
dr <- dr %>% 
  mutate_variables(alr_p_D = sweep(log(p_D[,-ncol(p_D)]), 1, log(p_D[,ncol(p_D)]), "-"))
alr_p_D <- E(dr$alr_p_D)
dist_D <- as.matrix(pairwise_dist[rownames(dr$p_D), colnames(dr$p_D)])
A_D <- matrix(tapply(fish_data$A, fish_data$pop, mean)[colnames(dr$alr_p_D)],
              nrow = nrow(dr$alr_p_D), ncol = ncol(dr$alr_p_D), byrow = TRUE)
cols <- c("green","blue","brown","darkgray")

# p_D vs pairwise dist
windows()
plot(dist_D, dr$p_D, pch = "")
for(i in 1:nrow(dr$p_D))
  points(dist_D[i,], dr$p_D[i,], pch = 1, cex = 1.5, col = cols[i])
legend("topright", rownames(dr$p_D), pch = 1, pt.cex = 1.5, col = cols)

# alr(p_D) vs pairwise dist
lm1 <- lm(as.vector(alr_p_D) ~ as.vector(dist_D[,-ncol(dist_D)]))
summary(lm1)
# windows()
png(filename = here("analysis","results","archive","alr(p_D)-vs-distance.png"), 
    width=7, height=7, units="in", res=300, type = "cairo-png")
par(mar = c(5,5,1,1))
plot(dist_D[,-ncol(dist_D)], alr_p_D, pch = "", las = 1, cex.axis = 1.2, cex.lab = 1.5,
     xlab = "Distance (km)", ylab = bquote(alr(italic(p)[D])))
for(i in 1:nrow(alr_p_D))
  points(dist_D[i,-ncol(dist_D)], alr_p_D[i,], pch = 1, cex = 1.5, col = cols[i])
abline(coef(lm1))
legend("topright", rownames(dr$alr_p_D), pch = 1, cex = 1.2, pt.cex = 1.5, col = cols)
dev.off()

# alr(p_D) vs target habitat size
lm2 <- lm(as.vector(alr_p_D) ~ as.vector(A_D))
summary(lm2)
# windows()
png(filename = here("analysis","results","archive","alr(p_D)-vs-area.png"), 
    width=7, height=7, units="in", res=300, type = "cairo-png")
par(mar = c(5,5,1,1))
plot(A_D, alr_p_D, pch = "", las = 1, cex.axis = 1.2, cex.lab = 1.5,
     xlab = "Distance (km)", ylab = bquote(alr(italic(p)[D])))
for(i in 1:nrow(alr_p_D))
  points(A_D[i,], alr_p_D[i,], pch = 1, cex = 1.5, col = cols[i])
abline(coef(lm2))
legend("topright", rownames(dr$alr_p_D), pch = 1, cex = 1.2, pt.cex = 1.5, col = cols)
dev.off()


# posterior distributions of alr(p_D)
dat <- as.data.frame(t(dr$alr_p_D)) %>% 
  mutate(pop = factor(rownames(.), levels = levels(fish_data$pop))) %>% 
  pivot_longer(cols = -pop, names_to = "origin", values_to = "alr_p_D") %>% 
  mutate(origin = factor(origin, levels = levels(fish_data$pop)))

gg <- dat %>% 
  ggplot(aes(xdist = alr_p_D, y = pop)) +
  stat_eye(.width = c(0.5, 0.9), normalize = "groups", 
           color = "slategray4", fill = alpha("slategray4", 0.5)) + 
  scale_y_discrete(limits = rev) + labs(x = bquote(alr(italic(p)[D])), y = "") + 
  facet_wrap(vars(origin), nrow = 1) + 
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = NA),
        strip.text = element_text(margin = margin(b = 3, t = 3)))

windows(width = 11, height = 7)
show(gg)

# (a)symmetry in implied logistic normal prior on annual straying distribution
# mu_p_D ~ dirichlet(1)
# alr_p_D[year] = alr(mu_p_D) + zeta_D * sigma_D
# zeta_D ~ N(0,1)
mu_p_D <- gtools::rdirichlet(10000, rep(1,12))
alr_mu_p_D <- sweep(log(mu_p_D[,-ncol(mu_p_D)]), 1, log(mu_p_D[,ncol(mu_p_D)]), "-")
sigma_D <- abs(rnorm(nrow(mu_p_D), 0, 3)) # single SD for all entries
zeta_D <- array(rnorm(prod(dim(alr_mu_p_D)), 0, 1), dim = dim(alr_mu_p_D))
p_D <- alr_mu_p_D + sweep(zeta_D, 1, sigma_D, "*")
p_D <- exp(cbind(p_D, 0))
p_D <- sweep(p_D, 1, rowSums(p_D), "/")

windows()
par(mfrow = c(3,4), mar = c(3,3,1,1))
for(i in 1:ncol(mu_p_D))
  hist(mu_p_D[,i], prob = TRUE, xlim = c(0,1), xlab = "", ylab = "", main = "")

windows()
par(mfrow = c(3,4), mar = c(3,3,1,1))
for(i in 1:ncol(p_D))
  hist(p_D[,i], prob = TRUE, xlim = c(0,1), xlab = "", ylab = "", main = "")


