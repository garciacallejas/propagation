
# add habitat cover to grid cells

# -------------------------------------------------------------------------

library(tidyverse)
library(terra)
library(sf)
library(exactextractr)

# -------------------------------------------------------------------------

# cell size in meters
grid.size <- 10000

# NZ <- st_read('../datasets/spatial_data/NZ_main_islands.shp')
# transform to NZTM2000 projection (https://epsg.io/2193)
# NZ2 <- st_transform(NZ, crs= st_crs(2193))

NZ_grid <- st_read(paste("data/NZ_grid_",grid.size/1e3,"km.shp",sep=""))

habitat.path <- "/home/david/Work/datasets/global_habitat_types/iucn_habitatclassification_composite_lvl1_ver004/"
habitat.raster <- terra::rast(paste(habitat.path,"iucn_habitatclassification_composite_lvl1_ver004.tif",sep=""))
codes <- read.csv2("/home/david/Work/datasets/global_habitat_types/habitat_codes.csv")

# first, reproject NZ grid to WGS84 in order to crop the habitat raster
NZ.wgs84 <- st_transform(NZ_grid,crs = st_crs(4326))
hab.nz <- crop(habitat.raster, ext(NZ.wgs84) + .01)
# if I try to reproject the full habitat raster to NZTM2000 it runs out of memory
# habitat.raster <- terra::project(habitat.raster,"EPSG:2193")

# now, reproject back to NZTM
# this takes a bit of time
# it is crucial to use method = "near", otherwise it will interpolate categories
hab.nz <- terra::project(hab.nz,"EPSG:2193", method = "near")
# add the labels of the land-use categories
levels(hab.nz) <- codes[,c(2,1)]
# test the frequencies in the full raster
nz.categories <- freq(hab.nz)

# and now, summarise habitat types per grid cell
habitat.freq <- exact_extract(hab.nz,NZ_grid,"frac")
names(habitat.freq) <- nz.categories$value
habitat.freq <- cbind(data.frame(cell.id = 1:nrow(habitat.freq)),habitat.freq)

habitat.long <- habitat.freq %>% pivot_longer(2:ncol(habitat.freq),
                                              names_to = "habitat_type",
                                              values_to = "frequency")

habitat.long$habitat_type <- recode(habitat.long$habitat_type,
                                    "1. Forest" = "forest",
                                    "11 Marine Deep Ocean Floor (Benthic and Demersal)" = "marine_benthic",
                                    "12 Marine Intertidal" = "marine_intertidal",
                                    "14 Artificial - Terrestrial" = "anthropic_habitats",
                                    "2. Savanna" = "savanna",
                                    "3. Shrubland" = "shrubland",
                                    "4. Grassland" = "grassland",
                                    "5. Wetlands (inland)" = "wetland",
                                    "6. Rocky Areas (e.g., inland cliffs, mountain peaks)" = "rocky_area",
                                    "8. Desert" = "desert",
                                    "9. Marine Neritic" = "marine_neritic")

# -------------------------------------------------------------------------

write.csv2(habitat.long,paste("data/land_use_frequencies_",grid.size/1e3,"km.csv",sep=""),row.names = F)

