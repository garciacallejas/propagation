
library(tidyverse)
library(sf)

# -------------------------------------------------------------------------

orig.data <- read.csv("../datasets/plant-bird interactions and traits/plant_bird_interactions.csv")

# small issues

# these names were updated in the gbif occ_search
# Actinidia deliciosa is a variety of Actinidia chinensis
orig.data$PLANTSPECIES[which(orig.data$PLANTSPECIES == "Actinidia_deliciosa")] <- "Actinidia_chinensis"
# Androstoma empetrifolium -> a
orig.data$PLANTSPECIES[which(orig.data$PLANTSPECIES == "Androstoma_empetrifolium")] <- "Androstoma_empetrifolia"
# citrus paradisi is not retrieved
# citrus sinensis is not retrieved
# dacrydium cupressinum... downloaded manually
# Dendrobenthamia_capitata is cornus capitata
orig.data$PLANTSPECIES[which(orig.data$PLANTSPECIES == "Dendrobenthamia_capitata")] <- "Cornus_capitata"
# Dysoxylum_spectabile... downloaded manually
# Eriobotrya japonica is thapiolepis loquata
orig.data$PLANTSPECIES[which(orig.data$PLANTSPECIES == "Eriobotrya_japonica")] <- "Rhaphiolepis_loquata"
# Leucopogon_fraseri is Styphelia_nesophila
orig.data$PLANTSPECIES[which(orig.data$PLANTSPECIES == "Leucopogon_fraseri")] <- "Styphelia_nesophila"
# Phormium_cookianum is Phormium_colensoi
orig.data$PLANTSPECIES[which(orig.data$PLANTSPECIES == "Phormium_cookianum")] <- "Phormium_colensoi"
# Phyllocladus_alpinus is Phyllocladus_trichomanoides
orig.data$PLANTSPECIES[which(orig.data$PLANTSPECIES == "Phyllocladus_alpinus")] <- "Phyllocladus_trichomanoides"
# Phytolacca_octandra is Phytolacca_icosandra
orig.data$PLANTSPECIES[which(orig.data$PLANTSPECIES == "Phytolacca_octandra")] <- "Phytolacca_icosandra"
# Piper_excelsum is Macropiper_excelsum
orig.data$PLANTSPECIES[which(orig.data$PLANTSPECIES == "Piper_excelsum")] <- "Macropiper_excelsum"
# Pseudopanax_colensoi is Neopanax_colensoi
orig.data$PLANTSPECIES[which(orig.data$PLANTSPECIES == "Pseudopanax_colensoi")] <- "Neopanax_colensoi"
orig.data$PLANTSPECIES[which(orig.data$PLANTSPECIES == "Pseudopanax_colensoi_var_colensoi")] <- "Neopanax_colensoi"
orig.data$PLANTSPECIES[which(orig.data$PLANTSPECIES == "Pseudopanax_colensoi_var_ternatus")] <- "Neopanax_colensoi"
# Pseudopanax chathamicus only appears in small islands, very rare
# Psidium_guajava is not so common, but should be there
# Solanum_nodiflorum is Solanum_americanum
orig.data$PLANTSPECIES[which(orig.data$PLANTSPECIES == "Solanum_nodiflorum")] <- "Solanum_americanum"
# Streblus_heterophyllus is Paratrophis_microphylla
orig.data$PLANTSPECIES[which(orig.data$PLANTSPECIES == "Streblus_heterophyllus")] <- "Paratrophis_microphylla"
# Tropaeolum_speciosum downloaded manually

plant.sp <- sort(unique(orig.data$PLANTSPECIES))

# simplify subspecies/varieties
plant.sp <- str_replace(plant.sp,"_"," ")
plant.sp <- sort(unique(gsub("\\_.*","",plant.sp)))
plant.sp <- str_replace(plant.sp," ","_")

# -------------------------------------------------------------------------

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

NZ.grid <- st_read("data/NZ_grid_2cells.shp")

# taken directly from the grid. It may be the case that there are cells
# without observations, i.e. not represented in sp.obs
cell.id <- sort(unique(NZ.grid$grid_id))

# -------------------------------------------------------------------------

# issues:
# set a temporal threshold for observations?
# min(sp.obs$year,na.rm = T)

min.year <- 2010

# duplicated, but for name consistency
represented.cells <- cell.id

# THIS IS A TEST FOR ONLY 2 CELLS
sp.obs.2 <- subset(sp.obs, year >= min.year & grid_id %in% represented.cells)

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

# for plants, I discard taxa identified at the genus level,
# and a couple of very rare species (Pseudopanax chathamicus,
# Jasminum officinale, citrus paradisi, citrus sinensis)
plant.sp <- plant.sp[which(plant.sp %in% all.sp)]
# discarded.plants <- plant.sp[which(!(plant.sp %in% all.sp))]

# -------------------------------------------------------------------------
# with this clean set of plant and bird species,
# I need to update the list of interactions to only account for these sp
clean.int.data <- subset(orig.data, BIRDSPECIES %in% bird.sp & PLANTSPECIES %in% plant.sp)

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
# a dataframe with observations of species (plant here, birds below) per cell

plant.sp.cells <- subset(sp.obs.2,species %in% plant.sp)
plant.sp.cells <- plant.sp.cells[,c("grid_id","species")]

plant.sp.cells.2 <- plant.sp.cells %>% group_by(grid_id,species) %>% summarise(obs = n())

# only in the small set of cells for this example
plant.sp.cells.subset <- subset(plant.sp.cells.2,grid_id %in% represented.cells)

plant.sp.cells.wide.subset <- pivot_wider(plant.sp.cells.subset,
                                          names_from = grid_id,
                                          values_from = obs,
                                          values_fill = 0)

# -------------------------------------------------------------------------
bird.sp.cells <- subset(sp.obs.2,species %in% bird.sp)
bird.sp.cells <- bird.sp.cells[,c("grid_id","species")]

bird.sp.cells.2 <- bird.sp.cells %>% group_by(grid_id,species) %>% summarise(obs = n())

# only in the small set of cells for this example
bird.sp.cells.subset <- subset(bird.sp.cells.2,grid_id %in% represented.cells)

bird.sp.cells.wide.subset <- pivot_wider(bird.sp.cells.subset,
                                         names_from = grid_id,
                                         values_from = obs,
                                         values_fill = 0)

# now, species selected are only those that appear in the two cells
all.sp.subset <- c(unique(bird.sp.cells.wide.subset$species),
             unique(plant.sp.cells.wide.subset$species))

# -------------------------------------------------------------------------
# build the full block matrix, considering only cells with observations

matrix.n.rows <- length(all.sp.subset)*length(represented.cells)
df.names <- expand.grid(all.sp.subset,represented.cells)
matrix.names <- paste(df.names[,1],df.names[,2],sep="-")

block.matrix <- matrix(0,
                       nrow = matrix.n.rows,
                       ncol = matrix.n.rows,
                       dimnames = list(matrix.names,matrix.names))

# -------------------------------------------------------------------------
# populate the block matrix
# I need some auxiliary data structures to do it efficiently
num.sp <- length(all.sp.subset)

# a metaweb adjacency matrix
adjacency.matrix <- matrix(0,length(all.sp.subset),length(all.sp.subset),
                           dimnames = list(all.sp.subset,all.sp.subset))
diag(adjacency.matrix) <- 1

for(i.obs in 1:nrow(clean.int.data)){
  adjacency.matrix[clean.int.data$PLANTSPECIES[i.obs],clean.int.data$BIRDSPECIES[i.obs]] <- 1
  adjacency.matrix[clean.int.data$BIRDSPECIES[i.obs],clean.int.data$PLANTSPECIES[i.obs]] <- 1
}# for i.obs

# -------------------------------------------------------------------------
# conditions for linking:
# 1 - plant-bird interaction in the same cell

# go through each cell (i.e. the diagonal of the block matrix)
for(i.cell in 1:length(represented.cells)){
  
  # the template is the metaweb adjacency matrix, with all interactions marked
  # so I need to "prune it" to keep only the interactions in each cell
  my.adj.matrix <- adjacency.matrix
  
  # go through each position and check whether both sp appear in that cell
  for(i.row in 1:nrow(my.adj.matrix)){
    for(i.col in 1:ncol(my.adj.matrix)){
      # only if they interact, of course
      if(my.adj.matrix[i.row,i.col] == 1){
        
        # check row species
        row.sp <- rownames(my.adj.matrix)[i.row]
        # check col species
        col.sp <- colnames(my.adj.matrix)[i.col]
        
        my.cell <- as.character(represented.cells[i.cell])
        
        # check presences in this cell
        # check again consistency of cells - there is at least one cell (86)
        # in which there are birds but no plants, and thus it is not recorded
        # in plant.sp.cells.wide, so any plant in that cell should have obs = 0
        if(row.sp %in% bird.sp){
          my.row <- which(bird.sp.cells.wide.subset$species == row.sp)
          if(my.cell %in% names(bird.sp.cells.wide.subset)){
            row.sp.obs <- as.numeric(bird.sp.cells.wide.subset[my.row,my.cell])
          }else{
            row.sp.obs <- 0
          }
        }else{
          my.row <- which(plant.sp.cells.wide.subset$species == row.sp)
          if(my.cell %in% names(plant.sp.cells.wide.subset)){
            row.sp.obs <- as.numeric(plant.sp.cells.wide.subset[my.row,my.cell])
          }else{
            row.sp.obs <- 0
          }
        }
        
        if(col.sp %in% bird.sp){
          my.col <- which(bird.sp.cells.wide.subset$species == col.sp)
          if(my.cell %in% names(bird.sp.cells.wide.subset)){
            col.sp.obs <- as.numeric(bird.sp.cells.wide.subset[my.col,my.cell])
          }else{
            col.sp.obs <- 0
          }
        }else{
          my.col <- which(plant.sp.cells.wide.subset$species == col.sp)
          if(my.cell %in% names(plant.sp.cells.wide.subset)){
            col.sp.obs <- as.numeric(plant.sp.cells.wide.subset[my.col,my.cell])
          }else{
            col.sp.obs <- 0
          }
        }
        
        # if any is missing, mark the matrix position as 0
        if(row.sp.obs == 0 | col.sp.obs == 0){
          my.adj.matrix[i.row,i.col] <- 0
        }
        
      }# if matrix position == 1
    }# for i.col
  }# for i.row
  
  # overwrite this diagonal of the block matrix with the updated adjacency matrix
  init.row <- 1 + (num.sp * (i.cell - 1))
  init.col <- init.row
  end.row <- num.sp + (num.sp * (i.cell - 1))
  end.col <- end.row
  
  block.matrix[init.row:end.row,init.col:end.col] <- my.adj.matrix
  
}# for i.cell

# 2 - only for birds: population of the same species in an adjacent cell
# there are a lot of cells to check

# go through each bird observation (bird.sp.cell.wide)
# check adjacent cells (adjacent.cells.matrix) 
# and mark (block.matrix) as linked the present cell and the adjacents

for(i.sp in 1:nrow(bird.sp.cells.wide.subset)){
  for(i.obs in 2:ncol(bird.sp.cells.wide.subset)){
    if(bird.sp.cells.wide.subset[i.sp,i.obs]>0){
      
      # this is the cell id
      my.cell.id <- as.numeric(names(bird.sp.cells.wide.subset)[i.obs])
      # and the number
      my.cell.num <- which(represented.cells == my.cell.id)
      # same for the adjacent cells: get id and number
      my.adjacent.ids <- as.numeric(names(which(adjacent.cells.matrix[my.cell.num,])))
      my.adjacent.num <- which(represented.cells %in% my.adjacent.ids)
      # likewise, which is the position of my species
      my.sp.num <- which(all.sp.subset == bird.sp.cells.wide.subset$species[i.sp])
      
      # diagonal (this cell, this species)
      diag.row <- num.sp * (my.cell.num - 1) + my.sp.num 
      block.matrix[diag.row,diag.row] <- 1
      
      # non-diagonals (this species linked in other cells)
      # I am only filling the upper diagonal, since the matrix is symmetric
      non.diag.row <- num.sp * (my.cell.num - 1) + my.sp.num
      non.diag.cols <- num.sp * (my.adjacent.num - 1) + my.sp.num
      
      block.matrix[non.diag.row,non.diag.cols] <- 1
      
    }# if >0 observations
  }# for each grid cell
}# for each bird sp

# -------------------------------------------------------------------------
# save the block matrix to disk
save(block.matrix,file = "results/community_block_matrix_2CELLS.Rdata")

# and the clean interaction data
write.csv2(clean.int.data,"data/plant_bird_clean_interaction_data_2CELLS.csv",row.names = F)
write.csv2(bird.sp.cells.wide.subset,"data/bird_cell_presences_2CELLS.csv",row.names = F)
write.csv2(plant.sp.cells.wide.subset,"data/plant_cell_presences_2CELLS.csv",row.names = F)

