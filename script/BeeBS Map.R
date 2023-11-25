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

AEA_proj <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-125 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "

# ----------------------------------------
# Load raw data
# ----------------------------------------

dat = read.csv("../data/CWS Bombus survey data_Masterlist2017-2021.csv")

dat$stop_id <- paste0(dat$Latitude,"-",dat$Longitude) %>% factor() %>% as.numeric()
dat$stop_year <- paste0(dat$Latitude,"-",dat$Longitude,"-",dat$Year) %>% factor() %>% as.numeric()

# ----------------------------------------
# Unique survey locations
# ----------------------------------------

site_list <- dat[,c("Longitude","Latitude")] %>%
  unique() %>%
  st_as_sf(coords=c("Longitude","Latitude"),crs=4326, remove = FALSE) %>%
  st_transform(AEA_proj)

# ----------------------------------------
# Bounding box for plotting area (in AEA_proj crs)
# ----------------------------------------
study_box <- st_bbox(c(xmin = -850000, xmax = 30000, 
                       ymax = 2200000, ymin = 3400000), 
                     crs = st_crs(AEA_proj)) %>%
  st_as_sfc()

# ----------------------------------------
# Load spatial layers
# ----------------------------------------

BCR <- st_read("../data/Spatial/BCR/BCR_Terrestrial_master.shp")  %>%
  st_make_valid() %>%
  group_by(PROVINCE_S) %>%
  summarize(geometry = st_union(geometry))%>%
  st_transform(AEA_proj) %>%
  st_intersection(study_box)

# Elevation of north america (downloaded from http://www.cec.org/files/atlas_layers/0_reference/0_03_elevation/elevation_tif.zip)
elevation <- rast("../data/Spatial/NA_elevation/na_elevation.tif") %>%
  crop(st_transform(study_box,crs(.))) %>%
  terra::project(AEA_proj) %>%
  crop(study_box)

# # Lakes and rivers shapefile downloaded from: https://www.hydrosheds.org/products/hydrolakes
#  lakes <- st_read("../data/Spatial/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.shp") %>%
#    subset(Continent == "North America" & Lake_area > 1) %>%
#    st_make_valid() %>%
#    st_intersection(st_transform(study_box,crs(.))) %>%
#    st_transform(AEA_proj)

rivers <- st_read("../data/Spatial/NA_rivers/hydrography_l_rivers_v2.shp") %>%
  st_intersection(st_transform(study_box,crs(.))) %>%
  st_transform(AEA_proj) %>%
  st_intersection(study_box)



# ************************************************
# ************************************************
# Generate a plot for an example species
# ************************************************
# ************************************************
species_list <- unique(dat$Species) %>% na.omit()
table(dat$Species) %>% sort()
for (sp in species_list){
  
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
      summarize(n_Year = length(unique(Year)),
                n_count = sum())
    sp_ndet$n_Year[i] <- ndet$n_Year
    
  }
  
  stop_counts <- dat %>%
    as.data.frame() %>%
    group_by(stop_year) %>%
    summarize(count = sum(Species == sp,na.rm = TRUE))
  
  mean_count <- mean(stop_counts$count)
  if (mean_count > 0.1) mean_count <- signif(mean_count,2)
  if (mean_count < 0.1) mean_count <- signif(mean_count,1)
  
  mean_occ <- mean(stop_counts$count>0) %>% round(2)
  n_detections <- nrow(sp_dat)
  
  # ----------------------------------------
  # Map
  # ----------------------------------------
  
  colpal <- colorRampPalette(c("#FEFEFE", "#FBF7E2", "#FCF8D0", "#EEF7C2", "#CEF2B0", "#94E5A0", "#51C987", "#18A065", "#008C59", "#007F53", "#006344"))
  
  sp_plot <- ggplot(data = BCR)+
    geom_spatraster(data = elevation)+
    geom_sf(data = subset(BCR,PROVINCE_S == "YUKON"), fill = "gray90", col = "gray90", linewidth = 1.5, alpha = 0.2) +
    geom_sf(data = rivers, col = "darkblue", alpha = 0.3)+
    #geom_sf(data = lakes2, col = "#59F3F3", fill = "#59F3F3")+
    
    geom_sf(data = sp_ndet, col="gray50", size = 2)+
    geom_sf(data = sp_ndet, col = "white", size = 1.5)+
    geom_sf(data = subset(sp_ndet, n_Year > 0), col = "dodgerblue", size = 0.9)+
    geom_sf(data = study_box, col = "black", fill = "transparent",size = 1.5)+
    
    coord_sf(clip = "off",
             xlim = range(as.data.frame(st_coordinates(study_box))$X),
             ylim = range(as.data.frame(st_coordinates(study_box))$Y))+
    scale_fill_gradientn(colors = rev(colpal(10))[4:9],guide="none", na.value = "#36454f")+
    scale_size_continuous(name = "# years\nspecies\ndetected", range = c(1,6))+
    theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title.x=element_blank(), 
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(), 
          axis.text.y=element_blank(), 
          axis.ticks.y=element_blank())+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    theme(legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(5,10,5,-20),
          legend.title.align=0.5,
          legend.title = element_markdown(lineheight=.9,hjust = "left"))+
    annotate(geom="text",x=50000,y=3400000, label= paste0("Bombus\n",sp),lineheight = .85,hjust = 0,vjust=1,size=6,fontface =2) +
    
    annotate(geom="text",x=50000,y=2900000, label= paste0("Number of detections: ", n_detections),lineheight = .85,hjust = 0,size=4,fontface =1) +
    annotate(geom="text",x=50000,y=2850000, label= paste0("Mean PObs: ", mean_occ),lineheight = .85,hjust = 0,size=4,fontface =1) +
    annotate(geom="text",x=50000,y=2800000, label= paste0("Mean count per stop: ", mean_count),lineheight = .85,hjust = 0,size=4,fontface =1) +
    annotate(geom="text",x=50000,y=2750000, label= "Trend (occurrence):    ____",lineheight = .85,hjust = 0,size=4,fontface =1) +
    annotate(geom="text",x=50000,y=2220000, label= paste0("Prepared on ",Sys.Date()),size=3,lineheight = .75,hjust = 0,color="#3b3b3b")
  
  
  # ----------------------------------------
  # Add photo
  # ----------------------------------------
  photo_path <- paste0("../species_photos/",sp,".png")
  if (file.exists(photo_path)){
    img <-  readPNG(photo_path)
    img_dim <- dim(img)
    img_left <- 50000
    img_bottom <- 2950000
    img_top <- 3270000
    img_right <- (img_top-img_bottom)*(dim(img)[2]/dim(img)[1])+img_left
    
    sp_plot <- sp_plot + annotation_raster(img, 
                                           ymin = img_bottom,ymax= img_top,
                                           xmin = img_left,xmax = img_right)
    }
  sp_plot
  
  png(paste0("../output/Maps/",sp,".png"), width=12, height=8, units="in", res=300, type="cairo")
  print(sp_plot)
  dev.off()
}

