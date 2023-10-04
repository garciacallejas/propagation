
# tidy species and environment data

# INPUTS
# - species observations

# OUTPUTS
# - site x observations dataframe

# -------------------------------------------------------------------------
library(sf)
library(stars)
# library(raster)
library(terra)
library(tidyverse)
# remotes::install_github("hagc/rasterB")
# library(rasterB)
# -------------------------------------------------------------------------
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
# sp.list <- sort(unique(gsub("\\_.*","",sp.list)))

birds.tier1 <- read.csv2("data/birds_tier1.csv")

# to test concordance with the grid
# this should be already in WGS84
NZ <- st_read('../datasets/spatial_data/NZ_main_islands.shp')

# -------------------------------------------------------------------------
# cell size in meters
grid.size <- 100000

# where are species observations
sp.path <- "results/sp_observations"

# where are environmental rasters
env.path <- "../datasets/NZEnvDS_NZMG/"

# NZ grid
NZ_grid <- st_read(paste("data/NZ_grid_",grid.size/1e3,"km.shp",sep=""))

# -------------------------------------------------------------------------

sp.files <- list.files(path = sp.path,
                       pattern = ".csv",
                       full.names = T)

# -------------------------------------------------------------------------

tidy.obs <- list()

for(i.sp in 1:length(sp.files)){

  my.file <- read.csv2(sp.files[i.sp])

  # first, set the crs of the original data
  projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  obs <- st_as_sf(x = my.file,
                  coords = c("decimalLongitude", "decimalLatitude"),
                  crs = projcrs)

  # then, transform to NZTM2000
  obs.nz <- st_transform(obs, crs= st_crs(2193))
  observations_id <- st_intersection(obs.nz,NZ_grid)

  # convert to standard dataframe
  obs_id_df <- observations_id[,c("cell_id","species","year","geometry","gbifID",
                                  "protocol","institutionCode")] %>%
    mutate(lat = sf::st_coordinates(.)[,2],
           lon = sf::st_coordinates(.)[,1],
           crs = "NZTM2000") %>%
    sf::st_set_geometry(NULL)
  
  tidy.obs[[i.sp]] <- obs_id_df %>%
    group_by(species,cell_id) %>%
    summarise(observations = n()) %>%
    dplyr::select(cell_id, species, observations)

}

obs.df <- bind_rows(tidy.obs) %>%
  arrange(cell_id, species) 

# obs.df.wide <- obs.df %>% pivot_wider(names_from = species,
#                                       values_from = observations, values_fill = 0) %>%
#   arrange(cell_id)

sp.list <- sort(unique(obs.df$species))

# -------------------------------------------------------------------------
# add TIER1 data

# I manually checked that most bird sp from tier1 that have interactions
# (i.e. are in sp.list), are correctly identified. The birds that do not appear
# in sp.list are species that have no recorded frugivory interactions.

# these are the only exceptions
birds.tier1$species[birds.tier1$species == "Callaeas wilsoni"] <- "Callaeas cinereus"
birds.tier1$species[birds.tier1$species == "Carduelis chloris"] <- "Chloris chloris"
birds.tier1$species[birds.tier1$species == "Carduelis flammea"] <- "Acanthis flammea"
birds.tier1$species[birds.tier1$species == "Philesturnus rufusater"] <- "Philesturnus carunculatus"

# keep only frugivores
birds.tier1.2 <- subset(birds.tier1, species %in% sp.list)
tier1.sp <- sort(unique(birds.tier1.2$species))

tier1.obs <- list()

for(i.b in 1:length(tier1.sp)){
  
  my.obs <- subset(birds.tier1.2,species == tier1.sp[i.b])
  my.obs.clean <- subset(my.obs, !is.na(x) | !is.na(y))
  
  obs <- st_as_sf(x = my.obs.clean,
                  coords = c("x", "y"),
                  crs = 2193)
  
  observations_id <- st_intersection(obs,NZ_grid)
  observations_id$gbifID <- NA_character_
  observations_id$protocol <- "tier1_DOC"
  observations_id$institutionCode <- "tier1_DOC"
  
  # convert to standard dataframe
  obs_id_df <- observations_id[,c("cell_id","species","year","geometry","gbifID",
                                  "protocol","institutionCode")] %>%
    mutate(lat = sf::st_coordinates(.)[,2],
           lon = sf::st_coordinates(.)[,1],
           crs = "NZTM2000") %>%
    sf::st_set_geometry(NULL)
  
  tier1.obs[[i.b]] <- obs_id_df %>%
    group_by(species,cell_id) %>%
    summarise(observations = n()) %>%
    dplyr::select(cell_id, species, observations)
}

tier1.obs.df <- bind_rows(tier1.obs)

tier1.num.sp <- tier1.obs.df %>%
  group_by(cell_id) %>%
  summarise(num_sp_tier1DOC = length(unique(species)))

# -------------------------------------------------------------------------
# combine gbif and tier1 observations
all.obs <- bind_rows(obs.df,tier1.obs.df) %>% 
  group_by(cell_id,species) %>%
  summarise(obs = sum(observations)) 
names(all.obs)[3] <- "observations"
all.obs$species <- stringr::str_replace(all.obs$species," ","_")

all.obs.wide <- all.obs %>%
  pivot_wider(names_from = species,
              values_from = observations, values_fill = 0) %>%
  arrange(cell_id)

# -------------------------------------------------------------------------
# write to disk

write.csv2(all.obs, paste("results/sp_observations_long_",grid.size/1e3,"km.csv",sep=""),row.names = F)
# write.csv2(all.obs.wide, paste("results/sp_observations_",grid.size/1e3,"km.csv",sep=""),row.names = F)

# -------------------------------------------------------------------------

env.files <- list.files(path = env.path, 
                        pattern='.tif', 
                        full.names=TRUE)

# to change the resolution in the raster objects - much faster than transforming first
# original resolution is 100m
resolution.factor <- grid.size/100

# get the original grid
# transform to NZTM2000 projection (https://epsg.io/2193)
NZ2 <- st_transform(NZ, crs= st_crs(2193))
# st_crs(NZ2)$proj4string
# st_crs(NZ2)$units_gdal

# generate spatial grid
grid <- st_as_stars(st_bbox(NZ2), dx = grid.size, dy = grid.size)
grid <- st_as_sf(grid)
grid <- grid[NZ2,]
grid$cell_id <- 1:nrow(grid)

env.list <- list()
for(i.file in 1:length(env.files)){
  
  env.raster <- terra::rast(env.files[i.file])
  
  # sometimes crs is not set in the original file. I took this from the first raster
  if(crs(env.raster) == ""){
    terra::crs(env.raster) <- "+proj=nzmg +lat_0=-41 +lon_0=173 +x_0=2510000 +y_0=6023150 +ellps=intl +units=m +no_defs"
  }
  
  # and we need to reproject it to NZTM2000
  # this is currently the operation that takes the longest
  env.raster <- terra::project(env.raster,"EPSG:2193")
  # now rescale to the lower resolution. Terra works really fast,
  # no need for the previous "rasterB" # env.scaled2 <- rasterB::aggregateB(env.raster,fact = resolution.factor,fun = mean)
  env.scaled <- terra::aggregate(env.raster,fact = resolution.factor,fun = mean,na.rm = TRUE)
  
  # convert to dataframe
  env.df <- as.points(env.scaled) %>% as.data.frame()
  # the as.points function does not return coordinates, so retrieve them as well
  coords <- crds(env.scaled)
  env.df$x <- coords[,"x"]
  env.df$y <- coords[,"y"]
  
  # now, join to the grid
  env.sf <- st_as_sf(env.df,coords = c("x","y"), crs = st_crs(2193))
  env.grid <- grid %>% st_join(env.sf,join = st_intersects)
  # convert back to dataframe
  env.df.clean <- data.frame(cell_id = env.grid$cell_id)
  env.df.clean[,2] <- as.vector(env.grid[,3])[1]
  names(env.df.clean)[2] <- names(env.grid)[3]
  env.list[[i.file]] <- env.df.clean
}

env.df <- purrr::reduce(env.list,dplyr::left_join)
env.df.2 <- left_join(env.df,tier1.num.sp)
# ggplot() +
#   geom_sf(data = NZ2, fill = 'white', lwd = 0.05) +
#   # geom_sf(data = pts, color = 'red', size = 1.7) +
#   geom_sf(data = grid,fill = "transparent", lwd = 0.3) +
#   geom_sf(data = env.sf, aes(fill = elevation)) +
#   # geom_text(data = grid_with_labels,
#   #           aes(x = X, y = Y, label = cell_id),
#   #           size = 2) +
#   coord_sf(datum = NA)  +
#   labs(x = "") +
#   labs(y = "")

# -------------------------------------------------------------------------
env.df.2$num_sp_tier1DOC[which(is.na(env.df.2$num_sp_tier1DOC))] <- 0
write.csv2(env.df.2,paste("results/environmental_factors_",grid.size/1e3,"km.csv",sep=""),row.names = F)

# -------------------------------------------------------------------------
# trait data
trait.data <- read.csv2("results/trait_data.csv")
trait.data.sp <- unique(trait.data$species)
# dont need this here, I just need to widen the data
trait.data$n <- NULL
trait.data$sd.value <- NULL
trait.data.wide <- pivot_wider(trait.data, names_from = trait,values_from = mean.value)

obs.sp <- stringr::str_replace(sp.list," ","_")

trait.no.obs <- trait.data.sp[which(!trait.data.sp %in% obs.sp)]
obs.no.trait <- obs.sp[which(!obs.sp %in% trait.data.sp)]

trait.clean <- data.frame(species = obs.sp)
trait.clean$guild <- trait.data.wide$guild[match(trait.clean$species,trait.data.wide$species)]
trait.clean$status <- trait.data.wide$status[match(trait.clean$species,trait.data.wide$species)]

# three exotic plants have NA values
trait.clean$guild[which(is.na(trait.clean$guild))] <- "plants"
trait.clean$status[which(is.na(trait.clean$status))] <- "Exotic"

# lowercase
trait.clean$status <- tolower(trait.clean$status)

# -------------------------------------------------------------------------
write.csv2(trait.clean,"results/species_traits.csv",row.names = F)


# -------------------------------------------------------------------------
# test changed resolution

# 
# # add x-y coords of the centroid
# grid_with_labels <- st_centroid(grid) %>% cbind(st_coordinates(.))
# 
# # -------------------------------------------------------------------------
# ggplot() +
#   geom_sf(data = NZ2, fill = 'white', lwd = 0.05) +
#   # geom_sf(data = env.sf,aes(color = humidity_meanAnn)) +
#   geom_sf(data = grid, fill = 'transparent', lwd = 0.3) +
#   # geom_sf(data = env.grid, aes(fill = humidity_meanAnn), lwd = 0.3) +
# 
#   # geom_text(data = grid_with_labels,
#   #           aes(x = X, y = Y, label = cell_id),
#   #           size = 2) +
# 
#   coord_sf(datum = NA)  +
#   labs(x = "") +
#   labs(y = "")
