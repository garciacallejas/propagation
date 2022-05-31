
# tidy species and environment data

# INPUTS
# - species observations

# OUTPUTS
# - site x observations dataframe

# -------------------------------------------------------------------------
library(sf)
library(stars)
library(raster)
library(tidyverse)
# remotes::install_github("hagc/rasterB")
library(rasterB)
# -------------------------------------------------------------------------
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

  # then, transform to NZMG
  obs.nz <- st_transform(obs, crs= st_crs(27200))
  observations_id <- st_intersection(obs.nz,NZ_grid)

  # convert to standard dataframe
  obs_id_df <- observations_id[,c("cell_id","species","year","geometry","gbifID",
                                  "protocol","institutionCode")] %>%
    mutate(lat = sf::st_coordinates(.)[,2],
           lon = sf::st_coordinates(.)[,1],
           crs = "NZMG") %>%
    sf::st_set_geometry(NULL)
  
  tidy.obs[[i.sp]] <- obs_id_df %>%
    group_by(species,cell_id) %>%
    summarise(observations = n()) %>%
    select(cell_id, species, observations)

}

obs.df <- bind_rows(tidy.obs) %>%
  arrange(cell_id, species) 
obs.df.wide <- obs.df %>% pivot_wider(names_from = cell_id,
                                      values_from = observations, values_fill = 0) %>%
  arrange(species)

# -------------------------------------------------------------------------
# write to disk

write.csv2(obs.df.wide, paste("results/sp_observations_",grid.size/1e3,"km.csv",sep=""),row.names = F)

# -------------------------------------------------------------------------

env.files <- list.files(path = env.path, 
                        pattern='.tif', 
                        full.names=TRUE)

# to change the resolution in the raster objects - much faster than transforming first
# original resolution is 100m
resolution.factor <- grid.size/100

# get the original grid
# transform to NZMG projection (https://epsg.io/27200)
NZ2 <- st_transform(NZ, crs= st_crs(27200))
# st_crs(NZ2)$proj4string
# st_crs(NZ2)$units_gdal

# generate spatial grid
grid <- st_as_stars(st_bbox(NZ2), dx = grid.size, dy = grid.size)
grid <- st_as_sf(grid)
grid <- grid[NZ2,]
grid$cell_id <- 1:nrow(grid)

env.list <- list()
for(i.file in 1:length(env.files)){
  env.raster <- raster(env.files[i.file])
  env.scaled <- rasterB::aggregateB(env.raster,fact = resolution.factor)
  
  # convert to dataframe
  env.df <- rasterToPoints(env.scaled) %>% as.data.frame()
  # now, join to the grid
  env.sf <- st_as_sf(env.df,coords = c("x","y"), crs = st_crs(27200))
  env.grid <- grid %>% st_join(env.sf,join = st_intersects)
  # convert back to dataframe
  env.df.clean <- data.frame(cell_id = env.grid$cell_id)
  env.df.clean[,2] <- as.vector(env.grid[,3])[1]
  names(env.df.clean)[2] <- names(env.grid)[3]
  env.list[[i.file]] <- env.df.clean
}

env.df <- purrr::reduce(env.list,dplyr::left_join)

# -------------------------------------------------------------------------

write.csv2(env.df,paste("results/environmental_factors_",grid.size/1e3,"km.csv",sep=""),row.names = F)

# -------------------------------------------------------------------------
# test changed resolution

# 
# # add x-y coords of the centroid
# grid_with_labels <- st_centroid(grid) %>% cbind(st_coordinates(.))
# 
# # -------------------------------------------------------------------------
ggplot() +
  geom_sf(data = NZ2, fill = 'white', lwd = 0.05) +
  # geom_sf(data = env.sf,aes(color = humidity_meanAnn)) +
  geom_sf(data = grid, fill = 'transparent', lwd = 0.3) +
  # geom_sf(data = env.grid, aes(fill = humidity_meanAnn), lwd = 0.3) +

  # geom_text(data = grid_with_labels,
  #           aes(x = X, y = Y, label = cell_id),
  #           size = 2) +

  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "")
