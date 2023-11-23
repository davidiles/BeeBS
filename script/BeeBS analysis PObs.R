# *****************************************************************
# *****************************************************************
# STEP 1: PREPARE SCRIPT
# *****************************************************************
# *****************************************************************

# ----------------------------------------
# Install/load necessary packages
# ----------------------------------------
my.packs = c('jagsUI',"ggplot2",'reshape2',
             'scales','tidyverse',
             'ggrepel','ggthemes','ggpubr')
if (any(!my.packs %in% installed.packages()[, 'Package']))install.packages(my.packs[which(!my.packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)
lapply(my.packs, require, character.only = TRUE)

setwd("D:/Working_Files/1_Projects/Side_Projects/BeeBS/script")

rm(list=ls())
theme_set(theme_bw()+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))

# ----------------------------------------
# Load dataset
# ----------------------------------------

dat = read.csv("../data/CWS Bombus count data_2017-2021.csv")
dat = subset(dat, !is.na(Species))

dat$Route_Number = dat$Route %>% as.factor() %>% as.numeric()
dat$Year_Number = dat$Year - min(dat$Year) + 1

# Links "year" to "year number"
year_table = unique(dat[,c("Year","Year_Number")])


# *****************************************************************
# *****************************************************************
# STEP 2: EXAMINE RAW DATA / SUMMARY STATISTICS
# *****************************************************************
# *****************************************************************

# Relative abundance of species
species_summary = dat %>% 
  group_by(Species) %>%
  summarize(total = sum(count, na.rm = TRUE),
            mean = mean(count, na.rm = TRUE),
            prop_obs = mean(count > 0, na.rm = TRUE)) %>%
  arrange(desc(total))

dat$Species <- factor(dat$Species, levels = species_summary$Species)
species_summary$Species <- factor(species_summary$Species, levels = species_summary$Species)

# Relative abundance in dataset
plot_1 <- ggplot(data = species_summary)+
  geom_bar(aes(x = Species, y = mean), fill = "dodgerblue", stat = "identity")+
  ggtitle("Relative abundance\n(across all years)")+
  ylab("Mean count per route")+
  xlab("Species")+
  theme()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
plot_1

pdf(paste0("../output/figures/summary_plots/summary_plot_1.pdf"),width=8,height=4)
print(plot_1)
dev.off()

annual_summary = na.omit(dat) %>%
  group_by(Species,Year) %>%
  summarize(mean = mean(count, na.rm = TRUE),
            prop_obs = mean(count>0, na.rm = TRUE),
            n_routes_surveyed = length(unique(Route)))

# Which species to plot
# Species_to_plot <- c("vancouverensis","mixtus","mckayi","frigidus","flavifrons","lapponicus")
Species_to_plot <- species_summary$Species[1:9] # plot 9 most abundant species

# Raw observations: separate lines for each route
plot_2 <- ggplot(data = subset(na.omit(dat), Species %in% Species_to_plot))+
  geom_point(aes(x = Year, y = count), col = "dodgerblue",alpha = 0.5)+
  geom_line(aes(x = Year, y = count, col = Route), alpha = 0.5, linetype = 2)+
  ggtitle("Observed counts")+
  ylab("Count")+
  xlab("Year")+
  facet_wrap(Species~., scales = "free_y")+
  scale_color_manual(values = rep("dodgerblue",length(unique(dat$Route))),
                     guide = "none")+
  theme()
plot_2

pdf(paste0("../output/figures/summary_plots/summary_plot_2.pdf"),width=8,height=6)
print(plot_2)
dev.off()

# Annual mean counts per route
plot_3 <- ggplot(data = subset(annual_summary, Species %in% Species_to_plot))+
  geom_point(aes(x = Year, y = mean, size = n_routes_surveyed),col = "dodgerblue")+
  geom_line(aes(x = Year, y = mean), col = "dodgerblue", linetype = 2, alpha = 0.5)+
  
  ggtitle("Mean counts per route")+
  ylab("Mean count")+
  xlab("Year")+
  facet_wrap(Species~., scales = "free_y")+
  scale_size_continuous(name = "# routes\nsurveyed", 
                        limits = c(5,20),
                        range = c(0.5,5))+
  theme()
plot_3

pdf(paste0("../output/figures/summary_plots/summary_plot_3.pdf"),width=8,height=6)
print(plot_3)
dev.off()

# Proportion of routes on which each species was observed
plot_4 <- ggplot(data = subset(annual_summary, Species %in% Species_to_plot))+
  geom_point(aes(x = Year, y = prop_obs, size = n_routes_surveyed),col = "dodgerblue")+
  geom_line(aes(x = Year, y = prop_obs), col = "dodgerblue", linetype = 2, alpha = 0.5)+
  
  ggtitle("Occurrence probability")+
  ylab("Occurence probability")+
  xlab("Year")+
  facet_wrap(Species~.)+
  scale_size_continuous(name = "# routes\nsurveyed", 
                        limits = c(5,20),
                        range = c(0.5,5))+
  theme()
plot_4

pdf(paste0("../output/figures/summary_plots/summary_plot_4.pdf"),width=8,height=6)
print(plot_4)
dev.off()

# *****************************************************************
# *****************************************************************
# STEP 3: DEFINE MODEL IN JAGS LANGUAGE
# *****************************************************************
# *****************************************************************

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


# *****************************************************************
# *****************************************************************
# STEP 4: LOOP THROUGH SPECIES, FIT MODEL, EXTRACT ESTIMATES
# *****************************************************************
# *****************************************************************
species_estimates <- data.frame()

for (sp in species_summary$Species){
  
  # ----------------------------------------
  # Create a 'data package' for JAGS
  # ----------------------------------------
  
  sp_dat = dat %>% subset(Species == sp)
  sp_dat$presence <- as.numeric(sp_dat$count > 0)
  
  # Convert to "wide" format.  Each row is a route, each column is a year
  # Occupancy status each year
  z = sp_dat[,c("Route","Year","presence")] %>% 
    spread(Year, presence)
  rownames(z) <- z$Route
  z = z[,-1]
  
  # Package data in JAGS
  jags_data = list(n_year = max(sp_dat$Year_Number),
                   n_route = max(sp_dat$Route_Number),
                   z = z)
  
  # ----------------------------------------
  # Fit model using JAGS
  # ----------------------------------------
  
  parameters.to.save = c("logit_phi_mean","logit_phi_sigma",
                         "logit_gamma_mean","logit_gamma_sigma",
                         "psi","phi","gamma","lambda","trend")
  out <- jags(data = jags_data,
              parameters.to.save = parameters.to.save,
              inits = NULL,
              n.iter = 110000,
              n.burnin = 10000,
              n.thin = 100,
              model.file = "model_PObs.jags",
              n.chains = 3,
              parallel = TRUE)

  out
  
  # ----------------------------------------
  # Extract (mean) parameter estimates to be used in power analysis
  # ----------------------------------------
  
  species_estimates <- rbind(species_estimates,
                             data.frame(Species = sp,
                                  trend = out$mean$trend,
                                  psi_1 = out$mean$psi[1],
                                  logit_phi_mean = out$mean$logit_phi_mean,
                                  logit_phi_sigma = out$mean$logit_phi_sigma,
                                  logit_gamma_mean = out$mean$logit_gamma_mean,
                                  logit_gamma_sigma = out$mean$logit_gamma_sigma))
  
  
  save(out, file = paste0("../output/",sp,"_out_PObs.Rdata"))
  
  # ----------------------------------------
  # Extract output and plot
  # ----------------------------------------
  
  # load(file = paste0("../output/",sp,"_out_PObs.Rdata"))
  
  # ----------------------------------------
  # Extract and plot annual indices
  # ----------------------------------------
  
  annual_psi = out$sims.list$psi %>%
    reshape2::melt() %>%
    rename(mc_sample = Var1, Year_Number = Var2, index = value) %>%
    group_by(Year_Number) %>%
    summarize(sp = sp,
              mean = mean(index),
              q025 = quantile(index,0.025),
              q05 = quantile(index,0.05),
              q50 = quantile(index,0.5),
              q95 = quantile(index,0.95),
              q975 = quantile(index,0.975)) %>%
    full_join(year_table)
  
  obs_psi <- sp_dat %>%
    group_by(Year) %>%
    summarize(psi_obs = mean(presence,na.rm = TRUE))
  
  annual_psi_plot = ggplot()+
    geom_ribbon(data = annual_psi, aes(x = Year, ymin = q05, ymax = q95, y = q50), 
                fill = "dodgerblue", alpha = 0.2)+
    geom_line(data = annual_psi, aes(x = Year,y = q50), col = "dodgerblue")+
    
    geom_point(data = obs_psi, aes(x = Year, y = psi_obs))+
    ggtitle("Annual occurrence estimate")+
    ylab("Occurrence probability")+
    xlab("Year")+
    ylim(0,1)
  print(annual_psi_plot)
  
  
  # ----------------------------------------
  # Estimate of trend
  # ----------------------------------------
  
  trend_df <- reshape2::melt(out$sims.list$trend)
  
  # rescale y axis
  trend_df$trend_rescale <- log(trend_df$value/100 + 1)
  axis_breaks = log(c(0,25,100,400)/100 + 1)
  axis_breaks = c(-axis_breaks,axis_breaks) %>% unique() %>% sort()
  axis_labels = paste0(100*(exp(axis_breaks)-1),"%")
  axis_labels[which(axis_breaks > 0)] <- paste0("+",axis_labels[which(axis_breaks > 0)])
  
  # Median trend estimate
  trend_q50 <- median(trend_df$value) %>% round(1)
  trend_q025 <- quantile(trend_df$value,0.025) %>% round(1)
  trend_q975 <- quantile(trend_df$value,0.975) %>% round(1)
  prob_decline <- mean(trend_df$value < 0) %>% round(2)
  
  trend_plot <- ggplot(trend_df, aes(y = trend_rescale, x = 1))+
    
    geom_hline(yintercept = 0, linetype = 2)+
    geom_violin(fill = "dodgerblue", alpha = 0.2, col = "dodgerblue",
                draw_quantiles = c(0.025,0.5,0.975))+
    scale_y_continuous(breaks = axis_breaks, labels = axis_labels, limits = range(axis_breaks))+
    scale_x_continuous(breaks = c(0,1,2), labels = rep("",3), limits = c(0,2))+
    
    ylab("% change per year")+
    xlab("")+
    ggtitle("5-year trend estimate")+
    geom_text(aes(x = 0, y = max(axis_breaks), 
                  label = paste0("Trend = ",trend_q50, "% per year\n(95% CI = ",trend_q025,"% to ", trend_q975,"%)\n\nProb. decline = ",prob_decline)),
              hjust = 0,
              vjust = 1,
              col = "dodgerblue",
              size = 3)+
    theme(axis.ticks.x = element_blank())
  
  
  # ----------------------------------------
  # Combine into single figure and save as png
  # ----------------------------------------
  
  species_plot <- ggarrange(annual_psi_plot, trend_plot,
            ncol = 2, nrow = 1, widths = c(0.7,0.3), align = "hv")
  
  
  species_plot_annotated <- annotate_figure(species_plot, 
                                            top =  text_grob(sp, color = "dodgerblue", face = "bold", size = 18))
  
  png(paste0("../output/figures/",sp,"_PObs.png"),width=10,height=6, units = "in", res = 500)
  print(species_plot_annotated)
  dev.off()
  
}

write.csv(species_estimates,file = "../output/species_estimates.csv",row.names = FALSE)
