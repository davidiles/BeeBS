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
             "terra",'stars',"tidyterra",
             "basemaps","ggtext")
if (any(!my.packs %in% installed.packages()[, 'Package']))install.packages(my.packs[which(!my.packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)
lapply(my.packs, require, character.only = TRUE)

setwd("C:/Users/IlesD/OneDrive - EC-EC/Iles/Projects/X_other_projects/BeeBS/script")

rm(list=ls())
theme_set(theme_bw()+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))


AEA_proj <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-125 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "

# ----------------------------------------
# Load raw data
# ----------------------------------------

dat = read.csv("../data/CWS Bombus survey data_Masterlist2017-2021.csv")
dat = subset(dat, !is.na(Species))
species_list <- unique(dat$Species)

# ----------------------------------------
# Unique survey locations
# ----------------------------------------

site_list <- dat[,c("Latitude","Longitude")] %>%
  unique() %>%
  st_as_sf(coords=c("Longitude","Latitude"),crs=4326, remove = FALSE) %>%
  st_transform(AEA_proj)

# ----------------------------------------
# Bounding box for plotting area (in AEA_proj crs)
# ----------------------------------------
study_box <- xlim = c(-850000,30000),
ylim = c(2200000,3400000)

# ----------------------------------------
# Load spatial layers
# ----------------------------------------

BCR <- st_read("../data/Spatial/BCR/BCR_Terrestrial_master.shp")  %>%
  st_make_valid() %>%
  group_by(PROVINCE_S) %>%
  summarize(geometry = st_union(geometry))%>%
  st_transform(AEA_proj)

study_box <- bind_rows(site_list,
                       subset(BCR, PROVINCE_S == "YUKON"))%>%
  st_buffer(100000)%>%
  st_transform(AEA_proj) %>%
  st_bbox() %>% st_as_sfc()

# Elevation of north america (downloaded from http://www.cec.org/files/atlas_layers/0_reference/0_03_elevation/elevation_tif.zip)
elevation <- rast("../data/Spatial/NA_elevation/na_elevation.tif") %>%
  crop(st_transform(study_box,crs(.))) 

# Lakes and rivers shapefile downloaded from: https://www.hydrosheds.org/products/hydrolakes
lakes <- st_read("../data/Spatial/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.shp") %>%
  subset(Continent == "North America" & Lake_area > 10) %>%
  st_make_valid() %>%
  st_intersection(st_transform(study_box,crs(.))) %>%
  st_transform(AEA_proj)
  
rivers <- st_read("../data/Spatial/NA_rivers/hydrography_l_rivers_v2.shp") %>%
  st_intersection(st_transform(study_box,crs(.))) %>%
  st_transform(AEA_proj)



# ************************************************
# ************************************************
# Generate a plot for an example species
# ************************************************
# ************************************************

sp <- "vancouverensis"

sp_dat <- subset(dat, Species == sp) %>%
  st_as_sf(coords=c("Longitude","Latitude"),crs=4326, remove = FALSE)%>%
  st_transform(AEA_proj)

# ----------------------------------------
# At each site, calculate number of years species was detected
# ----------------------------------------

sp_ndet <- site_list %>% mutate(n_Year = 0)

for (i in 1:nrow(site_list)){
  
  ndet <- st_intersection(sp_dat,site_list[i,] %>% 
                            st_buffer(100)) %>%
    as.data.frame() %>%
    summarize(n_Year = length(unique(Year)))
  sp_ndet$n_Year[i] <- ndet$n_Year
}

# ----------------------------------------
# Map
# ----------------------------------------

colpal <- colorRampPalette(c("#FEFEFE", "#FBF7E2", "#FCF8D0", "#EEF7C2", "#CEF2B0", "#94E5A0", "#51C987", "#18A065", "#008C59", "#007F53", "#006344"))

sp_plot <- ggplot(data = BCR)+
  geom_spatraster(data = elevation)+
  geom_sf(data = subset(BCR,PROVINCE_S == "YUKON"), fill = "transparent", col = "white", linewidth = 1.5, alpha = 0.2) +
  geom_sf(data = rivers, col = "darkblue", alpha = 0.3)+
  #geom_sf(data = lakes2, col = "#59F3F3", fill = "#59F3F3")+
  
  geom_sf(data = sp_ndet, col = "gray90")+
  geom_sf(data = subset(sp_ndet, n_Year > 0), col = "black", size = 0.5)+
  coord_sf(clip = "off",
           #xlim = range(as.data.frame(st_coordinates(study_box))$X),
           #ylim = range(as.data.frame(st_coordinates(study_box))$Y))+
           xlim = c(-850000,30000),
           ylim = c(2200000,3400000))+
  scale_fill_gradientn(colors = rev(colpal(10))[4:9],guide="none", na.value = "#36454f")+
  scale_size_continuous(name = "# years\nspecies\ndetected", range = c(1,6))+
  theme_bw()+
  ggtitle(sp)+
  theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(5,10,5,-20),legend.title.align=0.5,
        legend.title = element_markdown(lineheight=.9,hjust = "left"))+
  annotate(geom="text",x=346000,y=1850000, label= paste0(sp),lineheight = .85,hjust = 0,size=6,fontface =2) +
  annotate(geom="text",x=346000,y=1400000, label= paste0("Prepared on ",Sys.Date()),size=3,lineheight = .75,hjust = 0,color="#3b3b3b")

print(sp_plot)

png(paste0("../output/Maps/",sp,".png"), width=6.5, height=8, units="in", res=300, type="cairo")
print(sp_plot)
dev.off()

