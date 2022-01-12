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
# NZ grid
NZ_grid <- st_read("data/NZ_grid.shp")

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
# Androstoma empetrifolia
# Dacrydium cupressinum (taxonkey 5286163)
# Dysoxylum spectabile
# Tropaelum speciosum
# Neopanax colensoi (synonim Pseudopanax colensoi)
# Cornus capitata (7161042)
###############################################

for(i.sp in 1:length(sp.list)){
  
  print(paste(Sys.time(),"- downloading data from sp:",sp.list[i.sp]))
  
  #obtain data from GBIF via rgbif
  # TODO: update limit for the real thing
  dat <- occ_search(
    # taxonKey = 7161042,
    scientificName = sp.list[i.sp],
    limit = 10000, # hard limit is 100k
    country = "NZ",
    hasCoordinate = T
  )
  
  # names(dat$data) #a lot of columns
  if(!is.null(dat$data)){
    
    #select columns of interest
    # weird behaviour in which sometimes columns differ among species
    if("coordinateUncertaintyInMeters" %in% names(dat$data)){

    dat <- dat$data %>%
      dplyr::select(species, decimalLongitude, decimalLatitude, countryCode,
                    gbifID,
                    coordinateUncertaintyInMeters, year,
                    basisOfRecord, institutionCode)
    }else{
      dat <- dat$data %>%
        dplyr::select(species, decimalLongitude, decimalLatitude, countryCode,
                      gbifID, year,
                      basisOfRecord, institutionCode)
    }
    # names(dat)
    
    # clean up records without coordinates
    dat <- dat %>%
      filter(!is.na(decimalLongitude))%>%
      filter(!is.na(decimalLatitude))
    
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
    
    # more clean up of records with potential errors in different fields
    
    #convert country code from ISO2c to ISO3c
    dat$countryCode <-  countrycode(dat$countryCode, 
                                    origin =  'iso2c', 
                                    destination = 'iso3c')
    
    #flag problems
    dat <- data.frame(dat)
    flags <- clean_coordinates(x = dat, 
                               lon = "decimalLongitude", lat = "decimalLatitude",
                               countries = "countryCode", 
                               species = "species",
                               tests = c("capitals", "centroids", "equal",
                                         "gbif", "institutions",
                                         "zeros")) # most test are on by default
    
    # summary(flags)
    
    # exclude problematic records
    dat_cl <- dat[flags$.summary,]
    
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
    # assign a grid_id to each observation within the limits
    
    projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
    obs <- st_as_sf(x = dat_cl,                         
                    coords = c("decimalLongitude", "decimalLatitude"),
                    crs = projcrs)
    observations_id <- st_intersection(obs,NZ_grid)
    
    # convert to standard dataframe
    obs_id_df <- observations_id[,c("grid_id","species","year","geometry","gbifID")] %>%
      mutate(lat = sf::st_coordinates(.)[,2],
             lon = sf::st_coordinates(.)[,1]) %>%
      sf::st_set_geometry(NULL)
    
    # -------------------------------------------------------------------------
    # write to disk
    my.sp <- str_replace(sp.list[i.sp]," ","_")
  
    write.csv2(obs_id_df,paste("results/sp_observations/",
                               my.sp,
                               ".csv",sep=""),row.names = FALSE)
    
  }# if !is.null
}# for i.sp

# -------------------------------------------------------------------------
# plot test

# ggplot() +
#   geom_sf(NZ_grid, fill = 'transparent', lwd = 0.3, mapping = aes()) +
#   # geom_text(data = grid_lab, aes(x = X, y = Y, label = grid_id), size = 3, color = "darkred") +
#   # coord_sf(datum = NA)  +
#   geom_sf(data = observations_id,
#              colour = "cyan4", size = 1.2, alpha=0.6)+
#   NULL

# -------------------------------------------------------------------------
# another way of plotting the data

# world <- ne_countries(scale = "medium", returnclass = "sf")
# 
# coords<-dat_cl %>% 
#   dplyr::select(species,decimalLongitude, decimalLatitude ) %>%
#   rename(Longitude=decimalLongitude, Latitude=decimalLatitude)
# 
# coords_sf<- coords %>% 
#   st_as_sf(coords=(2:3))
# 
# bbox_new <- st_bbox(coords_sf)
# 
# ggplot()+
#   geom_sf(data = world, fill= NA) +
#   coord_sf(xlim=c(bbox_new[1],bbox_new[3]), 
#            ylim=c(bbox_new[2],bbox_new[4]))+
#   geom_point(data = dat_cl, aes(x = decimalLongitude, y = decimalLatitude),
#              colour = "cyan4", size = 1.2, alpha=0.6)+
#   annotation_scale(location = "br", width_hint = 0.1,line_width = 0.5) +
#   annotation_north_arrow(location = "tr", which_north = "true", 
#                          pad_x = unit(1.5, "cm"), 
#                          pad_y = unit(0.5, "cm"),
#                          height = unit(1, "cm"),
#                          width = unit(1, "cm"), # 0.2 # 0.3
#                          style = north_arrow_fancy_orienteering)+
#   labs(x = "", y = "",
#        title =paste("Record distribution of:", dat_cl$species),
#        caption = paste("Source: GBif\nAuthor: yep\nDate:",
#                        format(Sys.time(), '%d %B, %Y')))+
#   theme(plot.title=element_text(size=19),
#         plot.caption = element_text(size = 7, color="grey60"),
#         panel.background= element_rect(fill = "grey96"))
