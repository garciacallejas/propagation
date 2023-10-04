
library(tidyverse)

# -------------------------------------------------------------------------
grid.size <- 10

land.use.2016 <- read.csv2("results/land_use_tests/land_use_2016.csv")
land.use.1990 <- read.csv2("results/land_use_tests/land_use_1990.csv")
ecoregions <- read.csv2(paste("data/NZ_grid_ecoregions_",grid.size,"km.csv",sep=""))
species.presences.2016 <- read.csv2("results/species_land_use_ecoregion.csv")

categories <- data.frame(category = sort(unique(land.use.2016$LUCID_2016)),
                         category_char = c("Natural forest","Pre-1990 planted forest",
                                           "Post-1989 planted forest", "Grassland with woody biomass",
                                           "Grassland - high producing", "Grassland - low producing",
                                           "Cropland - perennial","Cropland - annual", "Wetland - open water",
                                           "Wetland - vegetated non forest", "Settlements", "Other"
                         ))

# -------------------------------------------------------------------------

territory.1990 <- left_join(ecoregions,land.use.1990)
territory.1990$land.use.category <- categories$category_char[match(territory.1990$LUCID_1990,
                                                                   categories$category)]

# should i consider only local introductions from 90 to 2016? 
# that is, species not present in 90 but present in 16
# or also local extinctions? perhaps the easiest way would be a SDM with
# land-use as covariate...

species.presences.1990 <- species.presences.2016

# do it sequentially because otherwise it would be too big a dataframe to handle
sp.obs.1990 <- list()

cell.id <- sort(unique(ecoregions$cell_id))

for(i.cell in cell.id){
  my.cell <- subset(territory.1990,cell_id == i.cell)
  my.ecoregions <- sort(unique(my.cell$ecoregion))
  for(i.reg in 1:length(my.ecoregions)){
    my.cell.reg <- subset(my.cell, ecoregion == my.ecoregions[i.reg])
    my.categories <- sort(unique(my.cell.reg$land.use.category))
    for(i.cat in 1:length(my.categories)){
      my.cat <- my.categories[i.cat]
      # my.sp <- species.presences.2016$species[which(species.presences.2016$ecoregion == my.ecoregions[i.reg] &
      #                                           species.presences.2016$land.use.category == my.cat)]
    }# i.cat
  }# for ecoregion

}# for i.cell









