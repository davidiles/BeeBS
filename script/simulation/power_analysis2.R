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

precision_results <- expand.grid(psi_1 = c(0.1,0.9),
                                 change_scenario = c(0.7),
                                 n_routes = c(10,20,40,80),
                                 n_years = c(10),
                                 simulation_rep = seq(1,100),
                                 trend_true = NA,
                                 trend_est = NA,
                                 trend_lcl = NA,
                                 trend_ucl = NA,
                                 trend_precision = NA,
                                 prob_decline = NA)

for (i in 1:nrow(precision_results)){
  
  print(i)
  
  # ------------------------------
  # Simulate data
  # ------------------------------
  
  # first year occupancy
  psi_1 <- precision_results$psi_1[i]
  
  # Required logit-linear trend to achieve 30% decline in abundance in year 10
  psi_10 <- psi_1*precision_results$change_scenario[i]
  logit_slope <- (qlogis(psi_10) - qlogis(psi_1))/9
  
  trend_true <- 100 * ((psi_10/psi_1)^(1/(10-1))-1)
  
  # Numerically calculate trend
  z_matrix <- matrix(NA, nrow=precision_results$n_routes[i], ncol = precision_results$n_years[i])
  z_matrix[,1] <-  rbinom(nrow(z_matrix), 1, prob = psi_1)
  
  # Loop through years and describe occupancy state
  for (t in 2:ncol(z_matrix)){
    
    psi_t <- plogis(qlogis(psi_1) + logit_slope*(t-1))
    z_matrix[, t] <- rbinom(nrow(z_matrix),1,psi_t)
      
  }
  
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
  precision_results$trend_true[i] <- trend_true
  precision_results$trend_est[i] <- out$mean$trend
  precision_results$trend_lcl[i] <- out$q2.5$trend
  precision_results$trend_ucl[i] <- out$q97.5$trend
  precision_results$trend_precision[i] <- out$q97.5$trend - out$q2.5$trend
  precision_results$prob_decline[i] <- mean(out$sims.list$trend<0)
  precision_results$Rmax[i] <- max(unlist(out$Rhat))
}

# Omit results where models did not converge
precision_results[which(precision_results$Rmax > 1.1),c("trend_est","trend_lcl","trend_ucl")] <- NA
precision_results <- na.omit(precision_results)

# Does 95% credible interval overlap true trend?
precision_results$cov <- (precision_results$trend_lcl < precision_results$trend_true) & (precision_results$trend_ucl > precision_results$trend_true)

# Labels for figures
precision_results$route_label <- paste0("# routes = ",precision_results$n_routes)
precision_results$psi_label <- paste0("Initial occupancy = ",precision_results$psi_1)

# ----------------------------------------------------
# Plot results for each simulation run
# ----------------------------------------------------

trend_true <- mean(precision_results$trend_true)
ggplot()+
  geom_hline(yintercept = 0, col = "gray80", linetype = 2)+
  geom_hline(yintercept = trend_true)+
  
  geom_errorbar(data = precision_results, aes(x = simulation_rep, 
                                              ymin=trend_lcl, 
                                              ymax = trend_ucl,
                                              col = factor(cov)), width = 0)+
  geom_point(data = precision_results, aes(x = simulation_rep, 
                                           y=trend_est, 
                                           col = factor(cov)))+
  scale_colour_manual(values=c("orangered","dodgerblue"), guide = "none")+
  xlab("Simulation number")+
  ylab("Trend")+
  ggtitle("Precision analysis (trend)")+
  facet_grid(psi_label~route_label)

mean(precision_results$cov)

# ----------------------------------------------------
# Estimated probability the trend is less than zero
# ----------------------------------------------------

facet_labels <- precision_results %>%
  group_by(route_label,psi_label) %>%
  summarize(mean_p = round(mean(prob_decline),2),
            p_label = paste0("Mean P(trend<0) = ", round(mean(prob_decline),2)))
ggplot()+
  geom_hline(yintercept = 0, col = "gray80", linetype = 2)+
  geom_hline(yintercept = trend_true)+
  
  geom_errorbar(data = precision_results, aes(x = simulation_rep, 
                                              ymin=trend_lcl, 
                                              ymax = trend_ucl,
                                              col = prob_decline), width = 0)+
  geom_point(data = precision_results, aes(x = simulation_rep, 
                                           y=trend_est, 
                                           col = prob_decline))+
  
  geom_text(data = facet_labels, aes(x = 0, y = -40, label = p_label), hjust=0)+
  scale_colour_gradientn(colors = viridis(10), limits = c(0,1), name = "Prob trend < 0")+
  
  
  xlab("Simulation number")+
  ylab("Trend")+
  ggtitle("Precision analysis (trend)")+
  facet_grid(psi_label~route_label)

# ----------------------------------------------------
# Summarize mean precision (width of 95% credible interval) for each scenario
# ----------------------------------------------------

precision_mean <- precision_results %>%
  group_by(route_label,psi_label) %>%
  summarize(CRI_width = mean(trend_precision))

ggplot()+
  
  geom_bar(data = precision_mean, aes(x = route_label,y=CRI_width),
           stat = "identity")+
  
  xlab("# routes surveyed per year")+
  ylab("Width of 95% credible interval\n(uncertainty of trend estimate)")+
  ggtitle("Uncertainty in trend estimate")+
  facet_grid(psi_label~.)
