library(tidyverse)
library(terra)
library(sf)

# -------------------------------------------------------------------------
grid.size <- 10

# -------------------------------------------------------------------------
nz.ecoregions <- sf::read_sf("../datasets/spatial_data/ecoregions/Ecological_Districts/Ecological_Districts.shp")

NZ.grid <- st_read(paste("data/NZ_grid_",grid.size,"km.shp",sep=""))
cell.id <- sort(unique(NZ.grid$cell_id))

# -------------------------------------------------------------------------
nz.ecoregions.2 <- st_transform(nz.ecoregions, crs= st_crs(2193))
grid.centroids <- st_centroid(NZ.grid)

centroid.ecoregions <- st_join(grid.centroids,nz.ecoregions.2) %>% 
  st_drop_geometry()  %>% 
  select(cell_id,ECOLOGIC_2) %>% 
  rename("ecoregion" = "ECOLOGIC_2") %>%
  arrange(cell_id)

# -------------------------------------------------------------------------
# some centroids are outside any ecoregion - these are surely centroids that
# fall in water bodies
# assign to them the ecoregion most observed from their neighbours
max.dist <- 30000

NZ.distances <- st_distance(grid.centroids, grid.centroids)
NZ.distances <- units::drop_units(NZ.distances)
rownames(NZ.distances) <- cell.id
colnames(NZ.distances) <- cell.id

for(i.cell in 1:nrow(centroid.ecoregions)){
  if(is.na(centroid.ecoregions$ecoregion[i.cell])){
    my.cell.distances <- NZ.distances[cell.id[i.cell],]
    sublandscape.ids <- which(my.cell.distances <= max.dist)
    sublandscape.ids <- sublandscape.ids[-which(sublandscape.ids == i.cell)]
    my.neigh.ecoregions <- centroid.ecoregions$ecoregion[which(centroid.ecoregions$cell_id %in% sublandscape.ids)]
    my.neigh.ecoregions <- na.omit(my.neigh.ecoregions)
    imputed.ecoregion <- names(table(my.neigh.ecoregions)[1])
    centroid.ecoregions$ecoregion[i.cell] <- imputed.ecoregion
  }
}

# -------------------------------------------------------------------------
write.csv2(centroid.ecoregions, paste("data/NZ_grid_ecoregions_",grid.size,"km.csv",sep=""),row.names = F)



