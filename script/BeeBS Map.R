# *****************************************************************
# *****************************************************************
# Map of species occurrence
# *****************************************************************
# *****************************************************************

# ----------------------------------------
# Install/load necessary packages
# ----------------------------------------
my.packs = c("ggplot2",'scales','tidyverse',
             'ggrepel','ggthemes','ggpubr',
             "viridis","sf",
             "terra",'stars',"tidyterra")
if (any(!my.packs %in% installed.packages()[, 'Package']))install.packages(my.packs[which(!my.packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)
lapply(my.packs, require, character.only = TRUE)

setwd("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/X_other_projects/BeeBS/script")

rm(list=ls())
theme_set(theme_bw()+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))

# ----------------------------------------
# Load spatial layers
# ----------------------------------------

AEA_proj <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-125 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "

BCR <- st_read("../data/Spatial/BCR/BCR_Terrestrial_master.shp")  %>%
  st_make_valid() %>%
  group_by(PROVINCE_S) %>%
  summarize(geometry = st_union(geometry))%>%
  st_transform(st_crs(AEA_proj))

# Land cover of Canada 2020 (not reprojected... faster to reproject grid)
#lcc2020 <- rast("../data/Spatial/LandCoverCanada2020/landcover-2020-classification.tif")
AMT <- rast("../data/Spatial/AnnualMeanTemperature/wc2.1_30s_bio_1.tif")
# ----------------------------------------
# Load raw data
# ----------------------------------------

dat = read.csv("../data/CWS Bombus survey data_Masterlist2017-2021.csv")
dat = subset(dat, !is.na(Species))

species_list <- unique(dat$Species)

site_list <- dat[,c("Latitude","Longitude")] %>%
  unique() %>%
  st_as_sf(coords=c("Longitude","Latitude"),crs=4326, remove = FALSE) %>%
  st_transform(AEA_proj)%>%
  st_intersection(BCR)

# Bounding box around study area (+/- 200 km in every direction)
study_box <- site_list %>%
  st_buffer(200000)%>%
  st_bbox() %>% 
  st_as_sfc()

# ----------------------------------------
# Crop lcc2020 to reasonable area
# ----------------------------------------

AMT <- AMT %>% crop(st_transform(study_box,crs(.))) 

# ----------------------------------------
# Plot example
# ----------------------------------------

sp <- "vancouverensis"

sp_dat <- subset(dat, Species == sp) %>%
  st_as_sf(coords=c("Longitude","Latitude"),crs=4326, remove = FALSE)%>%
  st_transform(AEA_proj)

# At each site, calculate number of years species was detected

sp_ndet <- site_list %>% mutate(n_Year = 0)

for (i in 1:nrow(site_list)){
  
  ndet <- st_intersection(sp_dat,site_list[i,] %>% 
                            st_buffer(100)) %>%
    as.data.frame() %>%
    summarize(n_Year = length(unique(Year)))
  sp_ndet$n_Year[i] <- ndet$n_Year
}


ggplot(BCR)+
  geom_spatraster(data = AMT)+
  geom_sf(data = BCR, fill = "transparent", col = "gray30") +
  
  geom_sf(data = subset(sp_ndet, n_Year > 0), aes(size = n_Year), col = "black")+
  geom_sf(data = subset(sp_ndet, n_Year == 0), col = "gray80", size = 0.5)+
  coord_sf(clip = "on",
           xlim = range(as.data.frame(st_coordinates(site_list))$X),
           ylim = range(as.data.frame(st_coordinates(site_list))$Y))+
  scale_fill_gradientn(colors = inferno(10)[3:9],guide="none", na.value = "transparent")+
  scale_size_continuous(name = "# years\nspecies\ndetected", range = c(1,6))+
  theme_bw()+
  ggtitle(sp)
