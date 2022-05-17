

# create 1km grid over a NZ shapefile
library(sf)
library(tidyverse)

# this should be already in WGS84
NZ <- st_read('../datasets/spatial_data/NZ_main_islands.shp')

# -------------------------------------------------------------------------

# create 1km grid - in these units, 0.1 is 1km
# TODO check - that is not right, it should be .0001 I think
grid_1 <- st_make_grid(NZ, cellsize = c(.001, .001)) %>% 
  st_sf(grid_id = 1:length(.))

# crop grid, keep only cells with land
grid_1_land <- st_intersection(grid_1,NZ)

# create labels for each grid_id
grid_lab <- st_centroid(grid_1_land) %>% cbind(st_coordinates(.))

# view the polygons and grid
# ggplot() +
#   # geom_sf(data = NZ, fill = 'white', lwd = 0.05) +
#   # geom_sf(data = pts, color = 'red', size = 1.7) +
#   geom_sf(data = grid_1_land, fill = 'transparent', lwd = 0.3) +
#   # geom_text(data = grid_lab, aes(x = X, y = Y, label = grid_id), size = 2) +
#   coord_sf(datum = NA)  +
#   labs(x = "") +
#   labs(y = "")

# -------------------------------------------------------------------------

st_write(grid_lab,"data/NZ_grid_labs.csv")
st_write(grid_1_land,"data/NZ_grid.shp")

# -------------------------------------------------------------------------
# select a small set of cells for testing

# my.cells <- c(45,46)
# 
# grid.subset <- subset(grid_1_land, grid_id %in% my.cells)
# 
# # ggplot() +
# #   geom_sf(data = grid.subset, fill = 'transparent', lwd = 0.3) +
# #   coord_sf(datum = NA)
# 
# st_write(grid.subset,"data/NZ_grid_2cells.shp")





