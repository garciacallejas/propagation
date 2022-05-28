library(rgbif)
library(tidyverse)
library("rnaturalearth")
library(sf)
library(ggspatial)
library(countrycode)
library(CoordinateCleaner)

# sources:
# https://rpubs.com/julianlavila/GbifCleanning

# -------------------------------------------------------------------------
# set grid size in meters
# grid.size <- 100000

# set time interval for observations - those older than this year will be discarded
year.first.obs <- 2010

# maximum number of observations gathered
# limit is 100000
max.obs <- 100000

# plot figure of species records?
plot.obs <- FALSE

# -------------------------------------------------------------------------
# NZ grid
# NZ_grid <- st_read(paste("data/NZ_grid_",grid.size/1e3,"km.shp",sep=""))
# 
# # shapefile for plotting
# NZ <- st_read('../datasets/spatial_data/NZ_main_islands.shp')
# NZ2 <- st_transform(NZ, crs= st_crs(27200))

# list of interactions with names of taxa
int.list <- read.csv("../datasets/plant-bird interactions and traits/plant_bird_interactions.csv")
sp.list <- unique(c(int.list$PLANTSPECIES,int.list$BIRDSPECIES))

sp.list <- str_replace(sp.list,"_"," ")

# taxa identified at the genus level often have other species of the same genus
# with observations (e.g. there are "Pseudopanax sp" and "Pseudopanax arboreus", etc)
# so it is not feasible to simply dump together all species in the "sp" observations.
# For now, remove them.
genus <- grep(" sp",sp.list)
sp.list <- sp.list[-genus]

# simplify subspecies/varieties
sp.list <- sort(unique(gsub("\\_.*","",sp.list)))

# -------------------------------------------------------------------------
# gather observations from every species

################### NOTE ######################
# Some species do not return valid observations in occ_search
# and I had to look for them manually, either using their taxonkey or name:
# Androstoma empetrifolium -> empetrifolia
# Callaeas cinerea -> cinereum
# Carduelis flammea -> not after 2010
# Citrus paradisi -> not after 2010

# Dacrydium cupressinum (taxonkey 5286163)
# Dysoxylum spectabile
# Tropaelum speciosum
# Neopanax colensoi (synonim Pseudopanax colensoi)
# Cornus capitata (7161042)
###############################################

# other errors
# 12  32  36  39  40  86  87  98 115 120 161 171 203 218
for(i.sp in 1:length(sp.list)){
  
  print(paste(Sys.time(),"- downloading data from sp:",sp.list[i.sp]))
  
  # account for some outlaws
  if(sp.list[i.sp] == "Dacrydium cupressinum"){
    dat <- occ_search(
      taxonKey = 5286163,
      # scientificName = sp.list[i.sp],
      limit = max.obs, # hard limit is 100k
      country = "NZ",
      hasCoordinate = T,
      fields=c("species", "decimalLongitude", "decimalLatitude", "countryCode",
               "gbifID",
               "coordinateUncertaintyInMeters", "year",
               "basisOfRecord", "institutionCode",
               "protocol")
    )
  }else if(sp.list[i.sp] == "Cornus capitata"){
    
  }else{
    #obtain data from GBIF via rgbif
    dat <- occ_search(
      # taxonKey = 2495905,
      scientificName = sp.list[i.sp],
      limit = max.obs, # hard limit is 100k
      country = "NZ",
      hasCoordinate = T,
      fields=c("species", "decimalLongitude", "decimalLatitude", "countryCode",
               "gbifID",
               "coordinateUncertaintyInMeters", "year",
               "basisOfRecord", "institutionCode",
               "protocol")
    )
  }
  
  # names(dat$data) #a lot of columns
 if(!is.null(dat$data)){
  #   
  #   #select columns of interest
  #   # weird behaviour in which sometimes columns differ among species
  #   if("coordinateUncertaintyInMeters" %in% names(dat$data)){
  # 
  #   dat <- dat$data %>%
  #     dplyr::select(species, decimalLongitude, decimalLatitude, countryCode,
  #                   gbifID,
  #                   coordinateUncertaintyInMeters, year,
  #                   basisOfRecord, 
  #                   institutionCode, 
  #                   protocol, 
  #                   samplingProtocol)
  #   }else{
  #     dat <- dat$data %>%
  #       dplyr::select(species, decimalLongitude, decimalLatitude, countryCode,
  #                     gbifID, year,
  #                     basisOfRecord, 
  #                     institutionCode, 
  #                     protocol, 
  #                     samplingProtocol)
  #   }
    # names(dat)
    
    # clean up records without coordinates
    # and by date
    dat2 <- dat$data %>%
      filter(!is.na(decimalLongitude))%>%
      filter(!is.na(decimalLatitude)) %>%
      filter(year >= year.first.obs)
    
    if(nrow(dat2)>0){
    
    #plot data to get an overview
    # wm <- borders("world", colour="gray70", fill="gray80")
    # ggplot()+ 
    #   coord_fixed()+ 
    #   wm +
    #   geom_point(data = dat, 
    #              aes(x = decimalLongitude, y = decimalLatitude),
    #              colour = "red", size = 0.9)+
    #   theme_bw()+ 
    #   labs(title = paste("Spatial distribution of records",dat$species))

    # -------------------------------------------------------------------------
    # more clean up of records with potential errors in different fields
    
    #convert country code from ISO2c to ISO3c
    dat2$countryCode <-  countrycode(dat2$countryCode, 
                                    origin =  'iso2c', 
                                    destination = 'iso3c')
    
    #flag problems
    dat2 <- data.frame(dat2)
    flags <- clean_coordinates(x = dat2, 
                               lon = "decimalLongitude", lat = "decimalLatitude",
                               countries = "countryCode", 
                               species = "species",
                               tests = c("capitals", "centroids", "equal",
                                         "gbif", "institutions",
                                         "zeros")) # most test are on by default
    
    # summary(flags)
    
    # exclude problematic records
    dat_cl <- dat2[flags$.summary,]
    
    # uncertainty of spatial coordinates
    # hist(dat_cl$coordinateUncertaintyInMeters / 1000, breaks = 20)
    
    # remove uncertain records
    if("coordinateUncertaintyInMeters" %in% names(dat_cl)){
      dat_cl <- dat_cl %>% 
        dplyr::filter(coordinateUncertaintyInMeters / 1000 <= 100 | is.na(coordinateUncertaintyInMeters))
    }    
    
    # remove suspicious data sources 
    # table(dat_cl$basisOfRecord)
    # dat_cl <- filter(dat_cl, basisOfRecord == "HUMAN_OBSERVATION")
    # dat_cl$year <- as.factor(dat_cl$year)
    # dat_cl %>% 
    #   group_by(year) %>% 
    #   summarise(records=n()) %>% 
    #   ggplot(aes(year, records, group = 1))+
    #   geom_line(color="steelblue")+
    #   labs(title=paste("Number of records of:",dat_cl$species, "by year"),
    #        y="",x="")+
    #   theme_bw()
    
    # -------------------------------------------------------------------------
    # write to disk
    my.sp <- str_replace(sp.list[i.sp]," ","_")
  
    write.csv2(dat_cl,paste("results/sp_observations/",
                               my.sp,".csv",sep=""),
               row.names = FALSE)
    
    }# if >0 observations
  }# if !is.null
}# for i.sp

# -------------------------------------------------------------------------
# which are not there?
my.sp.files <- list.files("results/sp_observations/")
my.sp.files <- str_replace(my.sp.files,"_"," ")
my.sp.files <- substr(my.sp.files,1,nchar(my.sp.files)-4)

remaining <- sp.list[which(!(sp.list %in% my.sp.files))]

# -------------------------------------------------------------------------
# assign a grid_id to each observation within the limits

# first, set the crs of the original data
# projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
# obs <- st_as_sf(x = dat_cl,                         
#                 coords = c("decimalLongitude", "decimalLatitude"),
#                 crs = projcrs)
# 
# # then, transform to NZMG
# obs.nz <- st_transform(obs, crs= st_crs(27200))
# observations_id <- st_intersection(obs.nz,NZ_grid)
# 
# # convert to standard dataframe
# obs_id_df <- observations_id[,c("cell_id","species","year","geometry","gbifID",
#                                 "protocol","institutionCode")] %>%
#   mutate(lat = sf::st_coordinates(.)[,2],
#          lon = sf::st_coordinates(.)[,1],
#          crs = "NZMG") %>%
#   sf::st_set_geometry(NULL)
# -------------------------------------------------------------------------
# if(plot.obs){
#   sp.plot <- ggplot()+
#     geom_sf(data = NZ2, fill = 'white', lwd = 0.05) +
#     # geom_sf(NZ_grid, fill = 'transparent', lwd = 0.3, mapping = aes()) +
#     coord_sf(datum = NA) +
#     geom_point(data = obs_id_df, aes(x = lon, y = lat),
#                colour = "cyan4", size = 1.2, alpha=0.6)+
#     annotation_scale(location = "br", width_hint = 0.1,line_width = 0.5) +
#     annotation_north_arrow(location = "tr", which_north = "true",
#                            pad_x = unit(0.5, "cm"),
#                            pad_y = unit(0.5, "cm"),
#                            height = unit(1, "cm"),
#                            width = unit(1, "cm"), # 0.2 # 0.3
#                            style = north_arrow_fancy_orienteering)+
#     labs(x = "", y = "",
#          title = obs_id_df$species[1],
#          caption = paste("Source: GBif\nAuthor: David Garcia-Callejas\nDate:",
#                          format(Sys.time(), '%Y-%m-%d')))+
#     theme(plot.title=element_text(size=14),
#           plot.caption = element_text(size = 7, color="grey60"),
#           panel.background= element_rect(fill = "grey96")) +
#     theme(axis.ticks = element_blank(), 
#           axis.text.x = element_blank(),
#           axis.text.y = element_blank()) +
#     NULL
#   
#   ggsave(paste("results/images/",my.sp,".pdf",
#                sep=""),
#          plot = sp.plot,
#          device = cairo_pdf,
#          width = 6, height = 6,dpi = 300)
#   
# }# if plot
