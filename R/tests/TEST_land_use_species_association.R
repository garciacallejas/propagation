
library(tidyverse)
library(sf)
library(terra)

# -------------------------------------------------------------------------
# NOTE:
# assume that a species present in a cell with a given land-use category
# in a given ecoregion is present in that category if it has >30% cover in the cell
min.category.cover <- 0.3

# -------------------------------------------------------------------------
grid.size <- 10

sp.obs <- read.csv2(paste("data/sp_observations_long_",grid.size,"km.csv",
                          sep=""))
nz.ecoregions <- read.csv2(paste("data/NZ_grid_ecoregions_",grid.size,"km.csv",
                                 sep=""))

land.use.1990 <- read.csv2("results/land_use_tests/land_use_1990.csv")
land.use.2016 <- read.csv2("results/land_use_tests/land_use_2016.csv")

categories <- data.frame(category = sort(unique(land.use.2016$LUCID_2016)),
                         category_char = c("Natural forest","Pre-1990 planted forest",
                                           "Post-1989 planted forest", "Grassland with woody biomass",
                                           "Grassland - high producing", "Grassland - low producing",
                                           "Cropland - perennial","Cropland - annual", "Wetland - open water",
                                           "Wetland - vegetated non forest", "Settlements", "Other"
                         ))
land.use.2016$category <- categories$category_char[match(land.use.2016$LUCID_2016,categories$category)]

# -------------------------------------------------------------------------
# i should do it more elegantly
cell.id <- sort(unique(nz.ecoregions$cell_id))
category.names <- sort(unique(categories$category_char))
sp <- sort(unique(sp.obs$species))
ecoregions <- sort(unique(nz.ecoregions$ecoregion))

sp.presences.list <- list()

cells.ecoregion.categories <- left_join(nz.ecoregions,land.use.2016)

for(i.cat in 1:length(category.names)){
  for(i.reg in 1:length(ecoregions)){
    
    # 1 - cells with this combination of category and ecoregion
    my.ecoregion.category <- subset(cells.ecoregion.categories,
                                    category == category.names[i.cat] &
                                      rel.area.2016 & min.category.cover &
                                      ecoregion == ecoregions[i.reg])
    # if there are cells, check every sp
    if(nrow(my.ecoregion.category)>0){
      
      my.eco.cat.sp <- left_join(my.ecoregion.category,sp.obs) %>% 
        drop_na() %>%
        arrange(cell_id, species)
      
      for(i.sp in 1:length(sp)){
        
        my.sp.presences <- sum(my.eco.cat.sp$species == sp[i.sp])
        my.sp.absences <- nrow(my.ecoregion.category) - my.sp.presences
        
        sp.presences.list[[length(sp.presences.list)+1]] <- 
          data.frame(species = sp[i.sp],
                     land.use.category = category.names[i.cat],
                     ecoregion = ecoregions[i.reg],
                     presences = my.sp.presences,
                     absences = my.sp.absences)
        
      }# species
    }# if ecoregion - category
    
  }# ecoregion
}# land-use category

species.presences <- bind_rows(sp.presences.list)
species.presences$presence.frequency <- species.presences$presences/(species.presences$presences + species.presences$absences)

# -------------------------------------------------------------------------
write.csv2(species.presences,"results/species_land_use_ecoregion.csv",row.names = F)

