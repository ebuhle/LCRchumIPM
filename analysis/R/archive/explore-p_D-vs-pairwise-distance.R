d <- stan_data('IPM_LCRchum_pp', fish_data = fish_data, 
               fecundity_data = fecundity_data, ages = list(M = 1))
dr <- as.matrix(fit_Ricker, 'p_D') %>% as_draws_rvars()
p_D <- E(dr$p_D)
rownames(p_D) <- levels(fish_data$pop)[d$which_O_pop]
colnames(p_D) <- levels(fish_data$pop)[d$which_D_pop]
alr_p_D <- sweep(log(p_D[,-ncol(p_D)]), 1, log(p_D[,ncol(p_D)]), "-") 
dist_D <- as.matrix(pairwise_dist[rownames(p_D), colnames(p_D)])
cols <- c("green","blue","brown","darkgray")

# p_D vs pairwise dist
windows()
plot(dist_D, p_D, pch = "")
for(i in 1:nrow(p_D))
  points(dist_D[i,], p_D[i,], pch = 1, cex = 1.5, col = cols[i])
legend("topright", rownames(p_D), pch = 1, pt.cex = 1.5, col = cols)

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
legend("topright", rownames(p_D), pch = 1, cex = 1.2, pt.cex = 1.5, col = cols)
dev.off()
