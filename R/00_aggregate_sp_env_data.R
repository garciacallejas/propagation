
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
sp.path <- "results/sp_observations/"

# where are environmental rasters
env.path <- "../datasets/NZEnvDS_NZMG/"

# -------------------------------------------------------------------------

# sp.files <- list.files(path = sp.path, 
#                        pattern = paste(grid.size/1e3,"km",sep=""),
#                        full.names = T)
# 
# # -------------------------------------------------------------------------
# 
# tidy.obs <- list()
# 
# for(i.sp in 1:length(sp.files)){
#   
#   my.file <- read.csv2(sp.files[i.sp])
#   
#   tidy.obs[[i.sp]] <- my.file %>%
#     group_by(species,cell_id) %>%
#     summarise(observations = n())
#   
# }
# 
# obs.df <- bind_rows(tidy.obs)

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

write.csv2(env.df,"results/environmental_factors_100km.csv")



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
  # geom_sf(data = grid, fill = 'transparent', lwd = 0.3) +
  geom_sf(data = env.grid, aes(fill = humidity_meanAnn), lwd = 0.3) +

  # geom_text(data = grid_with_labels,
  #           aes(x = X, y = Y, label = cell_id),
  #           size = 2) +

  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "")


# -------------------------------------------------------------------------
# write to disk

# write.csv2(obs.df, paste("results/sp_observations_",grid.size/1e3,"km.csv",sep=""))


