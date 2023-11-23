library(jagsUI)
library(tidyverse)

setwd("D:/Working_Files/1_Projects/Side_Projects/BeeBS/script/simulation")

rm(list=ls())

sink("model_PObs.jags")
cat("
    model {
    
    # ------------------------------
    # Priors
    # ------------------------------
    
    # Occupancy probability in first year
    psi[1] ~ dunif(0,1)  
    
    # Random effect priors for annual 'persistence probability' of occupied sites
    phi_mean ~ dunif(0,1)
    logit_phi_mean <- logit(phi_mean)
    logit_phi_sigma ~ dunif(0,5)
    logit_phi_tau <- pow(logit_phi_sigma,-2)
    
    # Random effect priors for annual 'colonization probability' of unoccupied sites
    gamma_mean ~ dunif(0,1)
    logit_gamma_mean <- logit(gamma_mean)
    logit_gamma_sigma ~ dunif(0,5)
    logit_gamma_tau <- pow(logit_gamma_sigma,-2)
    
    # Draw year effects for phi and gamma
    for (t in 1:n_year){
    
      # Logit scale
      logit_phi[t]   ~ dnorm(logit_phi_mean,logit_phi_tau)       # Persistence probability at occupied sites
      logit_gamma[t] ~ dnorm(logit_gamma_mean,logit_gamma_tau)   # Colonization probability at unoccupied sites
    
      # Probability scale
      phi[t]   <- ilogit(logit_phi[t])
      gamma[t] <- ilogit(logit_gamma[t])
    }
    
    # ------------------------------
    # Likelihood
    # ------------------------------
    
    # Loop through routes
    for (i in 1:n_route){
    
      # Occupancy state in first year
      z[i,1] ~ dbern(psi[1])
      
      # Loop through years and describe occupancy state
      for (t in 2:n_year){
        z[i, t] ~ dbern(z[i, t-1]*phi[t] + (1 - z[i, t-1])*gamma[t])
      
      }
      
    }
  
    # Occupancy probability as a derived parameter
    for(t in 2:n_year) {
      psi[t] <- psi[t-1]*phi[t] + (1-psi[t-1])*gamma[t]
      
      # Annual growth rate defined as ratio of annual occupancy
      lambda[t-1] <- psi[t]/psi[t-1]
    }
    
    # Estimates of mean percent change per year based on annual psi
    trend <- 100 * ( pow( psi[n_year]/psi[1], 1/(n_year-1) )-1)
    

}
    ",fill = TRUE)
sink()

species_estimates <- read.csv(file = "../../output/species_estimates.csv")

logit_phi_sigma <- mean(species_estimates$logit_phi_sigma)
logit_gamma_sigma <- mean(species_estimates$logit_gamma_sigma)

precision_results <- expand.grid(psi_1 = c(0.1,0.9),
                                 logit_phi_mean = c(-2,2),
                                 logit_phi_sigma = mean(species_estimates$logit_phi_sigma),
                                 logit_gamma_mean = c(-2,2),
                                 logit_gamma_sigma = mean(species_estimates$logit_gamma_sigma),
                                 n_routes = c(10,80),
                                 n_years = c(5,40),
                                 simulation_rep = seq(1,10),
                                 trend_true = NA,
                                 trend_est = NA,
                                 trend_lcl = NA,
                                 trend_ucl = NA,
                                 trend_precision = NA)

for (i in 1:nrow(precision_results)){
  
  print(i)
  
  # ------------------------------
  # Simulate data
  # ------------------------------
  
  # first year occupancy
  psi_1 <- precision_results$psi_1[i]
  
  # 'persistence probability' of occupied sites
  phi_mean_logit <- precision_results$logit_phi_mean[i]
  phi_sd_logit <- precision_results$logit_phi_sigma[i] # Temporal variance of phi
  
  # 'colonization probability' of unoccupied sites
  gamma_mean_logit <- precision_results$logit_gamma_mean[i]
  gamma_sd_logit <- precision_results$logit_gamma_sigma[i]  # Temporal variance of gamma
  
  # Numerically calculate trend
  z_matrix <- matrix(NA, nrow=10000, ncol = 100)
  z_matrix[,1] <-  rbinom(nrow(z_matrix), 1, prob = psi_1)
  
  
  # Loop through years and describe occupancy state
  
  for (t in 2:ncol(z_matrix)){
    
    # Draw random values for phi and gamma this year
    phi_t <- plogis(rnorm(1,phi_mean_logit,phi_sd_logit))
    gamma_t <- plogis(rnorm(1,gamma_mean_logit,gamma_sd_logit))
    
    for (j in 1:nrow(z_matrix)){
      
      
      z_matrix[j, t] <- rbinom(1,1,z_matrix[j, t-1]*phi_t + (1 - z_matrix[j, t-1])*gamma_t)
      
    }
  }
  
  # Trends (true)
  psi_true <- trends_true <- c()
  trends_true <- c()
  for (t in 1:ncol(z_matrix)){
    psi_true[t] <- mean(z_matrix[,t])
    trends_true[t] <- 100 * ((psi_true[t]/psi_true[1])^(1/(t-1))-1)
  }
  
  z_matrix <- z_matrix[1:precision_results$n_routes[i],1:precision_results$n_years[i]]
  
  # Package data in JAGS
  jags_data = list(n_year = ncol(z_matrix),
                   n_route = nrow(z_matrix),
                   z = z_matrix)
  
  parameters.to.save = c("psi","phi","gamma","lambda","trend")
  out <- jags(data = jags_data,
              parameters.to.save = parameters.to.save,
              inits = NULL,
              n.iter = 5000,
              n.burnin = 1000,
              n.thin = 1,
              model.file = "model_PObs.jags",
              n.chains = 3,
              parallel = TRUE)
  
  # Results for this simulation run
  precision_results$trend_true[i] <- trends_true[precision_results$n_years[i]]
  precision_results$trend_est[i] <- out$mean$trend
  precision_results$trend_lcl[i] <- out$q2.5$trend
  precision_results$trend_ucl[i] <- out$q97.5$trend
  precision_results$trend_precision[i] <- out$q97.5$trend - out$q2.5$trend
  
}

precision_results_summary <- precision_results %>%
  group_by(psi_1,logit_phi_mean,logit_phi_sigma,logit_gamma_mean,logit_gamma_sigma,n_routes,n_years) %>%
  summarize_all(mean)


ggplot(data = precision_results_summary, aes(x = n_years, 
                                     y = trend_true,
                                     col = factor(psi_1)))+
  geom_line()+
  ylab("Trend")+
  xlab("Number of years")+
  scale_color_manual(values = viridis::viridis(length(unique(precision_results$n_routes))),
                     name = "# routes")+
  ggtitle("Precision analysis (trend)")+
  facet_grid(.~logit_phi_mean+logit_gamma_mean)

  