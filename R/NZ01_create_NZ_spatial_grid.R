

# create square grid over a NZ shapefile
library(sf)
library(stars)
library(tidyverse)

# this should be already in WGS84
NZ <- st_read('data/NZ_main_islands.shp')

# -------------------------------------------------------------------------
# set cell size in meters

# grid.size <- 100000 # 100km
grid.size <- 10000 # 10km
# grid.size <- 1000 # 1km

# -------------------------------------------------------------------------

# transform to NZTM2000 projection (https://epsg.io/2193)
NZ2 <- st_transform(NZ, crs= st_crs(2193))
# st_crs(NZ2)$proj4string
# st_crs(NZ2)$units_gdal

# generate spatial grid
grid <- st_as_stars(st_bbox(NZ2), dx = grid.size, dy = grid.size)
grid <- st_as_sf(grid)
grid <- grid[NZ2,]
grid$cell_id <- 1:nrow(grid)

# add x-y coords of the centroid
grid_with_labels <- st_centroid(grid) %>% cbind(st_coordinates(.))

grid.coords <- data.frame(cell_id = grid_with_labels$cell_id,
                          lat_centroid = grid_with_labels$Y,
                          lon_centroid = grid_with_labels$X)

# Plot
# plot(st_geometry(grid), axes = TRUE, reset = FALSE)
# plot(st_geometry(NZ2), border = "darkred", add = TRUE)

# with ggplot
# grid.plot <- ggplot() +
#   geom_sf(data = NZ2, fill = 'white', lwd = 0.05) +
#   # geom_sf(data = pts, color = 'red', size = 1.7) +
#   geom_sf(data = grid, fill = 'transparent', lwd = 0.3) +
#   geom_text(data = grid_with_labels,
#             aes(x = X, y = Y, label = cell_id),
#             size = 2) +
#   coord_sf(datum = NA)  +
#   labs(x = "") +
#   labs(y = "")
# ggsave("results/images/NZ_grid_10km.png",
#        plot = grid.plot, width = 25, height = 25)

# -------------------------------------------------------------------------

st_write(grid_with_labels,paste("data/NZ_grid_",grid.size/1e3,"km.csv",sep=""),append = F)
st_write(grid,paste("data/NZ_grid_",grid.size/1e3,"km.shp",sep=""),append = F)
write.csv2(grid.coords,paste("data/NZ_grid_coords_",grid.size/1e3,"km.csv",sep=""))

# -------------------------------------------------------------------------

# older approach, takes much longer
# grid_1 <- st_make_grid(NZt, cellsize = c(grid.size, grid.size)) %>% 
#   st_sf(grid_id = 1:length(.))
# 
# # crop grid, keep only cells with land
# grid_1_land <- st_intersection(grid_1,NZt)
# 
# # create labels for each grid_id
# grid_lab <- st_centroid(grid_1_land) %>% cbind(st_coordinates(.))

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





