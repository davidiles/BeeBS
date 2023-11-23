library(jagsUI)
library(tidyverse)

setwd("D:/Working_Files/1_Projects/Side_Projects/BeeBS/script/simulation")

rm(list=ls())

# ------------------------------
# Simulate data
# ------------------------------

n_sites <- 50
n_years <- 10

z_matrix <- matrix(NA, nrow=n_sites, ncol = n_years)

z_matrix

# first year occupancy
prob_occ_1 <- 0.3

z_matrix[,1] <-  rbinom(n_sites, 1, prob = prob_occ_1)

# 'persistence probability' of occupied sites
phi <- 0.5

# 'colonization probability' of unoccupied sites
gamma <- 0.1

# Loop through years and describe occupancy state
for (i in 1:n_sites){
  
  for (t in 2:n_years){
    z_matrix[i, t] <- rbinom(1,1,z_matrix[i, t-1]*phi + (1 - z_matrix[i, t-1])*gamma)
    
  }
}

# Fit model to data
sink("dynocc.jags")
cat("
    model {
    
    # ------------------------------
    # Priors
    # ------------------------------
    
    # Occupancy probability in first year
    psi[1] ~ dunif(0,1)  
    
    # Random effect priors for annual 'persistence probability' of occupied sites
    phi ~ dunif(0,1)

    # Random effect priors for annual 'colonization probability' of unoccupied sites
    gamma ~ dunif(0,1)

    # ------------------------------
    # Likelihood
    # ------------------------------
    
    # Loop through routes
    for (i in 1:n_route){
    
      # Occupancy state in first year
      z[i,1] ~ dbern(psi[1])
      
      # Loop through years and describe occupancy state
      for (t in 2:n_year){
        z[i, t] ~ dbern(z[i, t-1]*phi + (1 - z[i, t-1])*gamma)
      
      }
      
    }
  
    # Occupancy probability as a derived parameter
    for(t in 2:n_year) {
      psi[t] <- psi[t-1]*phi + (1-psi[t-1])*gamma
      
      # Annual growth rate defined as ratio of annual occupancy
      lambda[t-1] <- psi[t]/psi[t-1]
    }
    
    # Estimates of mean percent change per year based on annual psi
    trend <- 100 * ( pow( psi[n_year]/psi[1], 1/(n_year-1) )-1)
    

}
    ",fill = TRUE)
sink()

# Package data in JAGS
jags_data = list(n_year = ncol(z_matrix),
                 n_route = nrow(z_matrix),
                 z = z_matrix)

parameters.to.save = c("psi","phi","gamma","lambda","trend")
out <- jags(data = jags_data,
            parameters.to.save = parameters.to.save,
            inits = NULL,
            n.iter = 50000,
            n.burnin = 10000,
            n.thin = 1,
            model.file = "dynocc.jags",
            n.chains = 3,
            parallel = TRUE)

out

hist(out$sims.list$phi)
abline(v = phi, col = "blue", lwd = 2) # actual phi

hist(out$sims.list$gamma)
abline(v = gamma, col = "blue", lwd = 2) # actual phi

# --------------
# Estimates of annual occupancy (compared to observed)
# --------------
annual_psi = out$sims.list$psi %>%
  reshape2::melt() %>%
  rename(mc_sample = Var1, Year_Number = Var2, index = value) %>%
  group_by(Year_Number) %>%
  summarize(mean = mean(index),
            q025 = quantile(index,0.025),
            q05 = quantile(index,0.05),
            q50 = quantile(index,0.5),
            q95 = quantile(index,0.95),
            q975 = quantile(index,0.975))

obs_psi <- data.frame(Year_Number = 1:ncol(z_matrix),
                      psi_obs = colMeans(z_matrix,na.rm = TRUE))

annual_psi_plot = ggplot()+
  #geom_errorbar(data = annual_psi, aes(x = Year, ymin = q05, ymax = q95, y = q50), width = 0)+
  #geom_point(data = annual_psi, aes(x = Year,y = q50))+
  
  geom_ribbon(data = annual_psi, aes(x = Year_Number, ymin = q05, ymax = q95, y = q50), 
              fill = "dodgerblue", alpha = 0.2)+
  geom_line(data = annual_psi, aes(x = Year_Number,y = q50), col = "dodgerblue")+
  
  geom_point(data = obs_psi, aes(x = Year_Number, y = psi_obs))+
  ggtitle("Annual occurrence estimate")+
  ylab("Occurrence probability")+
  xlab("Year")+
  ylim(0,1)
print(annual_psi_plot)