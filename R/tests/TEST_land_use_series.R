library(tidyverse)
library(terra)
library(sf)
library(foreach)
library(doParallel)

# set number of cores -----------------------------------------------------

# workers <- 4
# cl <- makeCluster(workers)
# # register the cluster for using foreach
# registerDoParallel(cl)

# -------------------------------------------------------------------------
grid.size <- 10

# -------------------------------------------------------------------------

lucas <- sf::st_read("../datasets/spatial_data/land-use/lucas/lucas-nz-land-use-map-1990-2008-2012-2016-v011.shp")
lucas <- lucas[,c("LUCID_2016","LUCID_1990","AREA_HA","geometry")]
lucas.categories <- sort(unique(lucas$LUCID_2016))
# no SUBID in 1990, so compare general categories
# around 20% of change in these polygons
sum(lucas$LUCID_2016 != lucas$LUCID_1990)/nrow(lucas)

# to attach to the 1km grid cells
NZ.grid <- st_read(paste("data/NZ_grid_",grid.size,"km.shp",sep=""))

# test
# NZ.grid <- NZ.grid[1:10,]

# -------------------------------------------------------------------------
# reproject the lucas shp
lucas <- st_transform(lucas,"EPSG:2193")

# start_time <- Sys.time()
# t1 <- st_intersection(NZ.grid[i.cell,],lucas)
# end_time <- Sys.time()
# 
# end_time - start_time
# 
# s1 <- Sys.time()
# t2 <- st_crop(lucas,NZ.grid[i.cell,])
# e1 <- Sys.time()
# e1 - s1

# this is surely not the best praxis, to do it cell by cell, 
# but I have not found a quick and efficient way
# area.1990.list <- foreach(i.cell = 1:nrow(NZ.grid), 
#                           # .combine=comb.fun, 
#                           .packages = 'tidyverse') %dopar% {
#                             
#                             lucas.cell <- st_crop(lucas,NZ.grid[i.cell,])
#                             
#                             lucas.1990.freq <- lucas.cell %>% 
#                               select(LUCID_1990,AREA_HA) %>%
#                               group_by(LUCID_1990) %>% 
#                               summarise(area = sum(AREA_HA)) %>%
#                               mutate(rel.area.1990 = area/sum(area),cell_id = i.cell) %>%
#                               select(cell_id,LUCID_1990,rel.area.1990) %>% 
#                               st_drop_geometry()
#                             
#                             write.csv2(lucas.1990.freq,
#                                        paste("/home/david/Work/Projects/NZ/propagation/results/land_use_tests/1990_cells/cell_",
#                                              sprintf("%04d",i.cell),".csv",sep=""),
#                                        row.names = F)
#                             lucas.1990.freq
#                             
#                           }

area.1990.list <- list()
area.2016.list <- list()

for(i.cell in NZ.grid$cell_id){
  lucas.cell <- st_crop(lucas,NZ.grid[i.cell,])
  
  lucas.1990.freq <- lucas.cell %>%
    select(LUCID_1990,AREA_HA) %>%
    group_by(LUCID_1990) %>%
    summarise(area = sum(AREA_HA)) %>%
    mutate(rel.area.1990 = area/sum(area),cell_id = i.cell) %>%
    select(cell_id,LUCID_1990,rel.area.1990) %>%
    st_drop_geometry()
  
  lucas.2016.freq <- lucas.cell %>%
    select(LUCID_2016,AREA_HA) %>%
    group_by(LUCID_2016) %>%
    summarise(area = sum(AREA_HA)) %>%
    mutate(rel.area.2016 = area/sum(area),cell_id = i.cell) %>%
    select(cell_id,LUCID_2016,rel.area.2016) %>%
    st_drop_geometry()
  
  area.1990.list[[length(area.1990.list)+1]] <- lucas.1990.freq 
  area.2016.list[[length(area.2016.list)+1]] <- lucas.2016.freq 
  
  write.csv2(lucas.1990.freq,
             paste("/home/david/Work/Projects/NZ/propagation/results/land_use_tests/1990_cells/cell_",
                   sprintf("%04d",i.cell),".csv",sep=""),
             row.names = F)
  write.csv2(lucas.2016.freq,
             paste("/home/david/Work/Projects/NZ/propagation/results/land_use_tests/2016_cells/cell_",
                   sprintf("%04d",i.cell),".csv",sep=""),
             row.names = F)
}
# 
land.use.1990 <- bind_rows(area.1990.list)
land.use.2016 <- bind_rows(area.2016.list)
write.csv2(land.use.1990,"results/land_use_tests/land_use_1990.csv",row.names = F)
write.csv2(land.use.2016,"results/land_use_tests/land_use_2016.csv",row.names = F)

# 
# land.use.total <- expand.grid(cell_id = NZ.grid$cell_id,
#                               land_cover_category = lucas.categories)
# land.use.total <- left_join(land.use.total,land.use.1990,by = c("cell_id","land_cover_category" = "LUCID_1990"))
# land.use.total <- left_join(land.use.total,land.use.2016,by = c("cell_id","land_cover_category" = "LUCID_2016"))

# not sure this is the way
# land.use.change <- land.use.total %>% 
#   drop_na() %>%
#   group_by(land_cover_category) %>%
#   summarise(freq.change = sum(rel.area.2016) - sum(rel.area.1990))

# lucas2016.grid <- terra::intersect(NZ.grid,lucas2016)
