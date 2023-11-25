# *****************************************************************
# *****************************************************************
# Analysis of stop-level BeeBS trends
# Stops are considered "sites"
# *****************************************************************
# *****************************************************************

# ----------------------------------------
# Install/load necessary packages
# ----------------------------------------
my.packs = c("ggplot2",'scales','tidyverse',
             'ggrepel','ggthemes','ggpubr',
             "viridis","sf",
             "terra",'stars',"tidyterra",
             "basemaps","ggtext",
             "jagsUI")
if (any(!my.packs %in% installed.packages()[, 'Package']))install.packages(my.packs[which(!my.packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)
lapply(my.packs, require, character.only = TRUE)

setwd("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/X_other_projects/BeeBS/script")

rm(list=ls())

AEA_proj <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-125 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "

# *****************************************************************
# STEP 1a: LOAD STOP-LEVEL DATA AND DO SOME PRELIMINARY PROCESSING
# *****************************************************************

dat = read.csv("../data/CWS Bombus survey data_Masterlist2017-2021.csv")

dat$Site <- paste0(dat$Latitude,"-",dat$Longitude) %>% factor() %>% as.numeric()
dat$Site_Year <- paste0(dat$Latitude,"-",dat$Longitude,"-",dat$Year) %>% factor() %>% as.numeric()

# How many stops per route?
n_stops <- dat %>% 
  group_by(Route) %>%
  summarize(n_stops = length(unique(Site))) %>%
  arrange(n_stops)

dat$Year_Number = dat$Year - min(dat$Year) + 1

# Links "Year" to "Year number"
Year_table = unique(dat[,c("Year","Year_Number")]) %>%
  arrange(Year)

# *****************************************************************
# STEP 1b: EVALUATE THE CLOSEST DISTANCES BETWEEN STOPS 
# (SHOULD ANY BE COMBINED?)
# *****************************************************************

# Unique survey locations
site_df <- dat[,c("Site","Longitude","Latitude")] %>%
  unique() %>%
  st_as_sf(coords=c("Longitude","Latitude"),crs=4326, remove = FALSE) %>%
  st_transform(AEA_proj)
site_df$Site

nrow(site_df) # 266 stops

# Distance between stops in metres
dists <- st_distance(site_df) 
diag(dists) <- NA
min_dists <- apply(dists,1,function(x) min(x,na.rm = TRUE)) %>% round(2)

range(min_dists)
sum(min_dists <= 1000) # 10 stops are less than 1000m from other stops

site_df$min_dist <- min_dists

# ------------------------------------
# Loop through stops, evaluate nearest site, and re-assign if necessary
# ------------------------------------

for (i in 1:nrow(site_df)){
  
  dists <- st_distance(site_df[i,],site_df) %>% as.numeric()
  dists[i] <- NA
  min_dist <- min(dists,na.rm = TRUE)
  if (min_dist >= 1000) next
  
  if (min_dist<1000){
    
    # Reassign Site
    old_Sites <- site_df$Site[which(dists <= 1000 | is.na(dists))]
    new_Site <- site_df$Site[which(dists <= 1000 | is.na(dists))][1]
    dat$Site[which(dat$Site %in% old_Sites )] <- new_Site
    
    # Reassign Lat/Long
    dat$Latitude[which(dat$Site %in% old_Sites )] <- site_df$Latitude[which(dists <= 1000 | is.na(dists))][1]
    dat$Longitude[which(dat$Site %in% old_Sites )] <- site_df$Longitude[which(dists <= 1000 | is.na(dists))][1]
    
  }
  
}

# ------------------------------------
# Create revised site_df dataframe; confirm 
# ------------------------------------
dat$Site <- factor(dat$Site) %>% as.numeric()
# Unique survey locations
site_df <- dat[,c("Site","Longitude","Latitude")] %>%
  unique() %>%
  st_as_sf(coords=c("Longitude","Latitude"),crs=4326, remove = FALSE) %>%
  st_transform(AEA_proj)

nrow(site_df) # 261 stops

# Distance between stops in metres
dists <- st_distance(site_df) 
diag(dists) <- NA
min_dists <- apply(dists,1,function(x) min(x,na.rm = TRUE)) %>% round(2)

range(min_dists)
sum(min_dists <= 1000) # Should be zero if code above worked

# *****************************************************************
# STEP 2: PREPARE STOP-LEVEL DATA FOR EACH SPECIES
# *****************************************************************

species_list <- table(dat$Species) %>% sort(decreasing = TRUE) %>% names()

Sites <- unique(dat$Site)
Site_dat <- expand.grid(Site = Sites,
                        Year = min(Year_table$Year):max(Year_table$Year)) %>%
  arrange(Site,Year)
Site_dat$Site_Year <- paste0(Site_dat$Site,"-",Site_dat$Year)

# Route membership of each site
Site_Route <- dat[,c("Site","Route")] %>% 
  unique()
duplicated_Routes <- Site_Route$Site[duplicated(Site_Route$Site)]


count_matrix <- matrix(NA, nrow = nrow(Site_dat), ncol = length(species_list))
rownames(count_matrix) <- Site_dat$Site_Year
colnames(count_matrix) <- species_list

dat$Site_Year <- paste0(dat$Site,"-",dat$Year)

for (i in 1:nrow(Site_dat)){
  
  # Was the site surveyed this year?
  idat <- subset(dat, Site_Year == Site_dat$Site_Year[i])
  if (nrow(idat)>0) count_matrix[i,] <- 0
  
  for (sp in species_list){
    
    sp_dat <- subset(dat, Site_Year == Site_dat$Site_Year[i] & Species == sp)
    if (nrow(sp_dat)==0) next
    
    count_matrix[i,sp] <- count_matrix[i,sp] + nrow(sp_dat)
  }
}

Site_dat     # location/year of each survey
count_matrix # Species counts associated with each survey

# *****************************************************************
# STEP 2: DEFINE TREND MODELS IN JAGS LANGUAGE
# *****************************************************************

sink("model_PObs_StopLevel.jags")
cat("
    model {
    
    # ------------------------------
    # Priors
    # ------------------------------
    
    # Occupancy probability in first Year
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
    
    # Draw Year effects for phi and gamma
    for (t in 1:n_Year){
    
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
    for (i in 1:n_Site){
    
      # Occupancy state in first Year
      z[i,1] ~ dbern(psi[1])
      
      # Loop through Years and describe occupancy state
      for (t in 2:n_Year){
        z[i, t] ~ dbern(z[i, t-1]*phi[t] + (1 - z[i, t-1])*gamma[t])
      
      }
      
    }
  
    # Occupancy probability as a derived parameter
    for(t in 2:n_Year) {
      psi[t] <- psi[t-1]*phi[t] + (1-psi[t-1])*gamma[t]
      
      # Annual growth rate defined as ratio of annual occupancy
      lambda[t-1] <- psi[t]/psi[t-1]
    }
    
    # Estimates of mean percent change per Year based on annual psi
    trend_GeometricMean <- 100 * ( pow( psi[n_Year]/psi[1], 1/(n_Year-1) )-1)

}
    ",fill = TRUE)
sink()

# *****************************************************************
# STEP 3: LOOP THROUGH SPECIES, FIT MODEL, EXTRACT ESTIMATES
# *****************************************************************

species_estimates <- data.frame()

for (sp in species_list){
  
  # ---------
  # Generate plots of abundance over time at each site
  # ---------
  
  sp_dat <- Site_dat %>% mutate(Species = sp, 
                                count = count_matrix[,sp],
                                presence = as.numeric(count_matrix[,sp]>0))
  
  mean_per_stop <- sp_dat %>%
    group_by(Year) %>%
    summarize(mean_count = mean(count,na.rm = TRUE),
              n_sites_surveyed = sum(!is.na(count)))
  
  ggplot()+
    geom_line(data = sp_dat, aes(x = Year, y = jitter(count,amount=0.1), col = factor(Site)))+
    geom_point(data = mean_per_stop,
               aes(x = Year, y = mean_count, size = n_sites_surveyed))+
    scale_color_manual(values = rep("dodgerblue",length(unique(sp_dat$Site))), guide = "none")+
    scale_size_continuous(limits = c(0,max(mean_per_stop$n_sites_surveyed)), name = "# sites surveyed")+
    
    theme_bw()+
    ylab("Count")+
    ggtitle(sp)
  
  # ---------
  # Prepare occupancy data
  # ---------
  
  # Convert to "wide" format.  Each row is a site, each column is a Year
  z = sp_dat[,c("Site","Year","presence")] %>% 
    spread(Year, presence)
  rownames(z) <- z$Route
  z = z[,-1]
  
  
  # Package data in JAGS
  jags_data = list(n_Year = nrow(Year_table),
                   n_Site = length(unique(sp_dat$Site)),
                   z = z)
  
  # ----------------------------------------
  # Fit model using JAGS
  # ----------------------------------------
  
  parameters.to.save = c("trend_GeometricMean",
                         "logit_phi_mean","logit_phi_sigma",
                         "logit_gamma_mean","logit_gamma_sigma",
                         "psi","phi","gamma","lambda")
  out <- jags(data = jags_data,
              parameters.to.save = parameters.to.save,
              inits = NULL,
              n.iter = 110000,
              n.burnin = 10000,
              n.thin = 100,
              model.file = "model_PObs_StopLevel.jags",
              n.chains = 3,
              parallel = TRUE)
  
  out$mcmc.info$elapsed.mins # about 1 minute
  
  # ----------------------------------------
  # Extract (mean) parameter estimates to be used in power analysis
  # ----------------------------------------
  
  species_estimates <- rbind(species_estimates,
                             data.frame(Species = sp,
                                        trend_GeometricMean = out$mean$trend_GeometricMean,
                                        psi_1 = out$mean$psi[1],
                                        logit_phi_mean = out$mean$logit_phi_mean,
                                        logit_phi_sigma = out$mean$logit_phi_sigma,
                                        logit_gamma_mean = out$mean$logit_gamma_mean,
                                        logit_gamma_sigma = out$mean$logit_gamma_sigma))
  
  
  #save(out, file = paste0("../output/",sp,"_out_PObs.Rdata"))
  
  # ----------------------------------------
  # Extract output and plot
  # ----------------------------------------
  
  # load(file = paste0("../output/",sp,"_out_PObs.Rdata"))
  
  # ----------------------------------------
  # Calculate logit-linear trend (trajectory)
  # ----------------------------------------
  logit_linear <- data.frame()
  logit_linear_pred <- data.frame()
  
  for (j in 1:out$mcmc.info$n.samples){
    fit <- lm(qlogis(out$sims.list$psi[j,])~Year_table$Year_Number)
    logit_linear <- rbind(logit_linear, data.frame(
      Intercept = as.numeric(coef(fit)[1]),
      Slope = as.numeric(coef(fit)[1])))
    
    logit_linear_pred <- rbind(logit_linear_pred, data.frame(
      samp = j,
      Year_Number = 1:jags_data$n_Year,
      pred = plogis(predict(fit))
      ))
    
  }
  
  logit_linear_pred_summary <- logit_linear_pred %>%
    group_by(Year_Number) %>%
    summarize(q025 = quantile(pred,0.025),
              q500 = quantile(pred,0.500),
              q975 = quantile(pred,0.975))  %>% left_join(Year_table)
  
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
    full_join(Year_table)
  
  obs_psi <- sp_dat %>%
    group_by(Year) %>%
    summarize(psi_obs = mean(presence,na.rm = TRUE))
  
  annual_psi_plot = ggplot()+
    
    geom_line(data = logit_linear_pred_summary, aes(x = Year,y = q500), col = "dodgerblue", size=1)+
    geom_ribbon(data = logit_linear_pred_summary, aes(x = Year,ymin = q025, ymax = q975), fill = "dodgerblue", col = "transparent", alpha = 0.2)+
    
    geom_errorbar(data = annual_psi, aes(x = Year, ymin = q05, ymax = q95, y = q50), 
                col = "darkblue", width=0)+
    geom_point(data = annual_psi, aes(x = Year,y = q50), col = "darkblue")+
    
    
    geom_point(data = obs_psi, aes(x = Year, y = psi_obs), shape = 1, size=3)+

    ggtitle("Annual occurrence estimate")+
    ylab("Occurrence probability")+
    xlab("Year")+
    ylim(0,1)+
    theme_bw() +
    theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank())
  print(annual_psi_plot)
  
  # ----------------------------------------
  # Estimate of log-linear trend
  # ----------------------------------------
  
  
  # Median trend_LogitLinear estimate
  trend_LogitLinear_q50 <- median(logit_linear$Slope) %>% round(1)
  trend_LogitLinear_q025 <- quantile(logit_linear$Slope,0.025) %>% round(1)
  trend_LogitLinear_q975 <- quantile(logit_linear$Slope,0.975) %>% round(1)
  prob_decline <- mean(logit_linear$Slope < 0) %>% round(2)
  
  trend_LogitLinear_plot <- ggplot(data = logit_linear, aes(y = Slope, x = 1))+
    
    geom_hline(yintercept = 0, linetype = 2)+
    geom_violin(fill = "dodgerblue", alpha = 0.2, col = "dodgerblue",
                draw_quantiles = c(0.025,0.5,0.975))+
    #scale_y_continuous(breaks = axis_breaks, labels = axis_labels, limits = range(axis_breaks))+
    scale_x_continuous(breaks = c(0,1,2), labels = rep("",3), limits = c(0,2))+
    
    #ylab("% change per Year")+
    xlab("")+
    ggtitle("5-Year Trend")+
    # geom_text(aes(x = 0, y = max(axis_breaks), 
    #               label = paste0("Trend = ",trend_LogitLinear_q50, "% per Year\n(95% CI = ",trend_LogitLinear_q025,"% to ", trend_LogitLinear_q975,"%)\n\nProb. decline = ",prob_decline)),
    #           hjust = 0,
    #           vjust = 1,
    #           col = "dodgerblue",
    #           size = 3)+
    theme(axis.ticks.x = element_blank())+
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # ----------------------------------------
  # Combine into single figure and save as png
  # ----------------------------------------
  
  species_plot <- ggarrange(annual_psi_plot, trend_GeometricMean_plot,
                            ncol = 2, nrow = 1, widths = c(0.7,0.3), align = "hv")
  
  
  species_plot_annotated <- annotate_figure(species_plot, 
                                            top =  text_grob(sp, color = "dodgerblue", face = "bold", size = 18))
  print(species_plot_annotated)
  
  png(paste0("../output/figures/",sp,"_PObs.png"),width=10,height=6, units = "in", res = 500)
  print(species_plot_annotated)
  dev.off()
  
}
