

# TODO this is going to be deprecated, and the spatial matrix 
# will be created on the fly to avoid memory issues 
# see the new "get_subset_communicability" for details

library(tidyverse)
library(sf)

# -------------------------------------------------------------------------

sp.obs <- read.csv2("results/plant_bird_observations_clean.csv")
sp.int <- read.csv2("results/plant_bird_interactions_clean.csv")

# plant.sp <- sort(unique(sp.int$PLANTSPECIES))
# 
# plant.df <- data.frame(species = plant.sp)
# plant.df$species <- str_replace(plant.df$species,"_"," ")

# modify to select grid of different size
sp.cells <- read.csv2("results/sp_observations_long_100km.csv")
NZ.grid <- st_read("data/NZ_grid_100km.shp")

# -------------------------------------------------------------------------
all.sp <- sort(unique(sp.obs$species))
plant.sp <- sort(unique(sp.int$PLANTSPECIES))
bird.sp <- sort(unique(sp.int$BIRDSPECIES))

# taken directly from the grid. It may be the case that there are cells
# without observations, i.e. not represented in sp.obs
cell.id <- sort(unique(NZ.grid$cell_id))

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

represented.cells <- cell.id[which(cell.id %in% unique(sp.cells$cell_id))]

matrix.n.rows <- length(all.sp)*length(represented.cells)
df.names <- expand.grid(all.sp,represented.cells)
matrix.names <- paste(df.names[,1],df.names[,2],sep="-")

block.matrix <- matrix(0,
                       nrow = matrix.n.rows,
                       ncol = matrix.n.rows,
                       dimnames = list(matrix.names,matrix.names))

# -------------------------------------------------------------------------
# populate the block matrix
# I need some auxiliary data structures to do it efficiently
num.sp <- length(all.sp)

# a metaweb adjacency matrix
adjacency.matrix <- matrix(0,length(all.sp),length(all.sp),
                           dimnames = list(all.sp,all.sp))
diag(adjacency.matrix) <- 1

# -------------------------------------------------------------------------
# a dataframe with observations of species (plant here, birds below) per cell
for(i.obs in 1:nrow(sp.int)){
  adjacency.matrix[sp.int$PLANTSPECIES[i.obs],sp.int$BIRDSPECIES[i.obs]] <- 1
  adjacency.matrix[sp.int$BIRDSPECIES[i.obs],sp.int$PLANTSPECIES[i.obs]] <- 1
}# for i.obs

plant.sp.cells <- subset(sp.cells,species %in% plant.sp)
# plant.sp.cells <- plant.sp.cells[,c("cell_id","species")]

# plant.sp.cells.2 <- plant.sp.cells %>% group_by(cell_id,species) %>% summarise(obs = n())

plant.sp.cells.wide <- pivot_wider(plant.sp.cells,
                                  names_from = cell_id,
                                  values_from = observations,
                                  values_fill = 0)

# -------------------------------------------------------------------------
bird.sp.cells <- subset(sp.cells,species %in% bird.sp)
# bird.sp.cells <- bird.sp.cells[,c("cell_id","species")]
# 
# bird.sp.cells.2 <- bird.sp.cells %>% group_by(cell_id,species) %>% summarise(obs = n())

bird.sp.cells.wide <- pivot_wider(bird.sp.cells,
                                  names_from = cell_id,
                                  values_from = observations,
                                  values_fill = 0)

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
          my.row <- which(bird.sp.cells.wide$species == row.sp)
          if(my.cell %in% names(bird.sp.cells.wide)){
            row.sp.obs <- as.numeric(bird.sp.cells.wide[my.row,my.cell])
          }else{
            row.sp.obs <- 0
          }
        }else{
          my.row <- which(plant.sp.cells.wide$species == row.sp)
          if(my.cell %in% names(plant.sp.cells.wide)){
            row.sp.obs <- as.numeric(plant.sp.cells.wide[my.row,my.cell])
          }else{
            row.sp.obs <- 0
          }
        }
        
        if(col.sp %in% bird.sp){
          my.col <- which(bird.sp.cells.wide$species == col.sp)
          if(my.cell %in% names(bird.sp.cells.wide)){
            col.sp.obs <- as.numeric(bird.sp.cells.wide[my.col,my.cell])
          }else{
            col.sp.obs <- 0
          }
        }else{
          my.col <- which(plant.sp.cells.wide$species == col.sp)
          if(my.cell %in% names(plant.sp.cells.wide)){
            col.sp.obs <- as.numeric(plant.sp.cells.wide[my.col,my.cell])
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

for(i.sp in 1:nrow(bird.sp.cells.wide)){
  for(i.obs in 2:ncol(bird.sp.cells.wide)){
    if(bird.sp.cells.wide[i.sp,i.obs]>0){
      
      # this is the cell id
      my.cell.id <- as.numeric(names(bird.sp.cells.wide)[i.obs])
      # and the number
      my.cell.num <- which(represented.cells == my.cell.id)
      # same for the adjacent cells: get id and number
      my.adjacent.ids <- as.numeric(names(which(adjacent.cells.matrix[my.cell.num,])))
      my.adjacent.num <- which(represented.cells %in% my.adjacent.ids)
      # likewise, which is the position of my species
      my.sp.num <- which(all.sp == bird.sp.cells.wide$species[i.sp])
      
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
save(block.matrix,file = "results/community_block_matrix.Rdata")

# and the clean presence data
write.csv2(bird.sp.cells.wide,"data/bird_cell_presences.csv",row.names = F)
write.csv2(plant.sp.cells.wide,"data/plant_cell_presences.csv",row.names = F)

