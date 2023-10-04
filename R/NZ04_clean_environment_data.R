
library(tidyverse)
library(sf)
library(stars)
library(terra)

# -------------------------------------------------------------------------
env.path <- "../datasets/NZEnvDS_NZMG/"
NZ <- st_read('../datasets/spatial_data/NZ_main_islands.shp')

grid.size <- 10000

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
env.df.long <- pivot_longer(env.df,2:ncol(env.df),names_to = "environmental_factor",values_to = "value")
# env.df.2 <- left_join(env.df,tier1.num.sp)
# ggplot() +
#   geom_sf(data = NZ2, fill = 'white', lwd = 0.05) +
#   # geom_sf(data = pts, color = 'red', size = 1.7) +
#   geom_sf(data = grid,fill = "transparent", lwd = 0.3) +
#   # geom_sf(data = env.sf, aes(fill = elevation)) +
#   # geom_text(data = grid_with_labels,
#   #           aes(x = X, y = Y, label = cell_id),
#   #           size = 2) +
#   coord_sf(datum = NA)  +
#   labs(x = "") +
#   labs(y = "")

# -------------------------------------------------------------------------
# env.df.2$num_sp_tier1DOC[which(is.na(env.df.2$num_sp_tier1DOC))] <- 0
write.csv2(env.df.long,paste("data/environmental_factors_long_",grid.size/1e3,"km.csv",sep=""),row.names = F)
