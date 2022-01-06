
library(tidyverse)
library(sf)

# -------------------------------------------------------------------------

orig.data <- read.csv("../datasets/plant-bird interactions and traits/plant_bird_interactions.csv")

# small issues
# these names were automatically updated in the occ_search
orig.data$BIRDSPECIES[which(orig.data$BIRDSPECIES == "Callaeas_cinerea")] <- "Callaeas_cinereus"
# Callaeas_wilsoni is a subsp of Callaeas_cinereus
orig.data$BIRDSPECIES[which(orig.data$BIRDSPECIES == "Callaeas_wilsoni")] <- "Callaeas_cinereus"
# Carduelis chloris is Chloris chloris
orig.data$BIRDSPECIES[which(orig.data$BIRDSPECIES == "Carduelis_chloris")] <- "Chloris_chloris"
# Carduelis flammea is Acanthis flammea
orig.data$BIRDSPECIES[which(orig.data$BIRDSPECIES == "Carduelis_flammea")] <- "Acanthis_flammea"
# cyanoramphus forbesi only appears in Chatham Island
# Larus novaehollandiae is Chroicocephalus novaehollandiae
orig.data$BIRDSPECIES[which(orig.data$BIRDSPECIES == "Larus_novaehollandiae")] <- "Chroicocephalus_novaehollandiae"
# Mohoua novaeseelandiae is Finschia novaeseelandiae
orig.data$BIRDSPECIES[which(orig.data$BIRDSPECIES == "Mohoua_novaeseelandiae")] <- "Finschia_novaeseelandiae"
# Philesturnus rufusater is found only in offshore islands and sanctuaries
# Strygops habroptilus is Strygops habroptila
orig.data$BIRDSPECIES[which(orig.data$BIRDSPECIES == "Strigops_habroptilus")] <- "Strigops_habroptila"

bird.sp <- sort(unique(orig.data$BIRDSPECIES))
# simplify subspecies/varieties
bird.sp <- str_replace(bird.sp,"_"," ")
bird.sp <- sort(unique(gsub("\\_.*","",bird.sp)))
bird.sp <- str_replace(bird.sp," ","_")

# -------------------------------------------------------------------------
sp.files <- list.files(path = "results/sp_observations",
                       pattern = "*.csv",
                       full.names = T)
sp.obs <- sp.files %>% map_dfr(read_csv2)
sp.obs <- subset(sp.obs, !is.na(year))

NZ.grid <- st_read("data/NZ_grid.shp")

# taken directly from the grid. It may be the case that there are cells
# without observations, i.e. not represented in sp.obs
cell.id <- sort(unique(NZ.grid$grid_id))

# -------------------------------------------------------------------------

# issues:
# set a temporal threshold for observations?
# min(sp.obs$year,na.rm = T)

min.year <- 2010
sp.obs.2 <- subset(sp.obs, year >= min.year)

# assume that bird populations from adjacent cells are linked
# in next iterations, this can be improved by accounting for differences in dispersal ability
# possibly inferred by looking at morphological traits

# -------------------------------------------------------------------------
# keep the set of all species
sp.obs.2$species <- str_replace(sp.obs.2$species," ","_")
all.sp <- sort(unique(sp.obs.2$species))

# now, a few bird species are not present in the two main islands, which is 
# the territory I am considering. 
# so, remove them
bird.sp <- bird.sp[which(bird.sp %in% all.sp)]

# -------------------------------------------------------------------------
# I need to know which grid cells are adjacent

# check, plot the grid and the ids
# grid_lab <- st_centroid(NZ.grid) %>% cbind(st_coordinates(.))
# grid.plot <- ggplot() +
#   # geom_sf(data = NZ, fill = 'white', lwd = 0.05) +
#   # geom_sf(data = pts, color = 'red', size = 1.7) +
#   geom_sf(data = NZ.grid, fill = 'transparent', lwd = 0.3) +
#   geom_text(data = grid_lab, aes(x = X, y = Y, label = grid_id), size = 2) +
#   coord_sf(datum = NA)  +
#   labs(x = "") +
#   labs(y = "")
# grid.plot

# this function returns the positions of the touching cells, not their id
adjacent.list <- st_touches(NZ.grid)

adjacent.cells.matrix <- matrix(FALSE,
                                nrow = length(cell.id),
                                ncol = length(cell.id),
                                dimnames = list(cell.id,cell.id))

# populate a matrix that tells me whether cells in position i,j are adjacent.
# Id of the cells are in matrix names
for(i.row in 1:length(cell.id)){
  # for(i.col in 1:length(cell.id)){
    my.cell.adjacents <- adjacent.list[[i.row]]
    adjacent.cells.matrix[i.row,my.cell.adjacents] <- TRUE
  # }
}

# -------------------------------------------------------------------------
# build the full block matrix, considering only cells with observations

represented.cells <- cell.id[which(cell.id %in% unique(sp.obs$grid_id))]

matrix.n.rows <- length(all.sp)*length(represented.cells)
matrix.names <- rep(all.sp,length(represented.cells))

block.matrix <- matrix(0,
                       nrow = matrix.n.rows,
                       ncol = matrix.n.rows,
                       dimnames = list(matrix.names,matrix.names))

# -------------------------------------------------------------------------
# populate the block matrix

# conditions for linking:
# 1 - plant-bird interaction in the same cell
# 2 - only for birds: population of the same species in an adjacent cell

# there are a lot of cells to check
# I need some auxiliary data structures to do it efficiently



