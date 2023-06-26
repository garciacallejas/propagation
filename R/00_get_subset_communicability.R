
# -------------------------------------------------------------------------
library(foreach)
library(doParallel)
library(tidyverse)
library(sf)
source("R/communicability_network.R")
source("R/communicability.R")

# -------------------------------------------------------------------------
# this will likely be a standalone script soon

# dispersal rate based on distance and hand-wing index of the species
disp.fun <- function(dist,hwi){
  my.kern <- exp(-dist*(1/hwi))
  return(my.kern)
}

# set number of cores -----------------------------------------------------

workers <- 4
cl <- makeCluster(workers)
# register the cluster for using foreach
registerDoParallel(cl)

# -------------------------------------------------------------------------

# which source to use?
observation.source <- "occurrence_prob"
# observation.source <- "observed_records"

# grid size?
grid.size <- 10

# maximum distance for the submatrices (in km)
max.dist <- 3*grid.size

sp.int.orig <- read.csv2("results/plant_bird_interactions_clean.csv")
sp.disp <- read.csv2(paste("results/species_dispersal_rates_",grid.size,"km.csv",sep=""))

# what to assign to birds without HWI measure? 
hwi <- read.csv2("results/bird_hand_wing_index.csv")

# # note: convert dispersal to meters
# sp.disp$dispersal.distance <- sp.disp$dispersal.distance * 1e3

if(observation.source == "observed_records"){
  sp.cells <- read.csv2(paste("data/sp_observations_long_",grid.size,"km.csv",
                              sep=""))
}else if(observation.source == "occurrence_prob"){
  sp.cells <- read.csv2(paste("results/model_occurrences_",grid.size,"km.csv",
                              sep=""))
  # this is to harmonise datasets
  names(sp.cells)[which(names(sp.cells) == "presence")] <- "observations"
  # keep presences only
  sp.cells <- subset(sp.cells, observations == 1)
  
}

# careful with this, ensure that the ids are stored as character with a leading zero
# sp.cells$cell_id <- formatC(sp.cells$cell_id, width = 5, format = "d", flag = "0")

NZ.grid <- st_read(paste("data/NZ_grid_",grid.size,"km.shp",sep=""))
# NZ.grid$cell_id <- formatC(NZ.grid$cell_id, width = 5, format = "d", flag = "0")

# -------------------------------------------------------------------------
# sanity check - ensure the sets of species are consistent
sp.int <- subset(sp.int.orig, PLANTSPECIES %in% sp.cells$species & 
                         BIRDSPECIES %in% sp.cells$species)

# -------------------------------------------------------------------------
plant.sp <- sort(unique(sp.int$PLANTSPECIES))
bird.sp <- sort(unique(sp.int$BIRDSPECIES))
all.sp <- c(plant.sp,bird.sp)
# max.disp <- max(sp.disp$dispersal.distance,na.rm = T) # already in meters

# taken directly from the grid. It may be the case that there are cells
# without observations, i.e. not represented in sp.obs
cell.id <- sort(unique(NZ.grid$cell_id))

# -------------------------------------------------------------------------
# centroid of the cells
NZ.centroid <- st_centroid(NZ.grid) %>% cbind(st_coordinates(.))

# distance matrix between all centroids
NZ.distances <- st_distance(NZ.centroid, NZ.centroid)
NZ.distances <- units::drop_units(NZ.distances)
rownames(NZ.distances) <- cell.id
colnames(NZ.distances) <- cell.id

# convert to KM
NZ.distances <- NZ.distances/1e3

# -------------------------------------------------------------------------
# other data cleaning necessary for calculating the local community matrices

num.sp <- length(all.sp)

# a metaweb adjacency matrix
adjacency.matrix <- matrix(0,length(all.sp),length(all.sp),
                           dimnames = list(all.sp,all.sp))
diag(adjacency.matrix) <- 1

# a dataframe with observations of species (plant here, birds below) per cell
for(i.obs in 1:nrow(sp.int)){
  adjacency.matrix[sp.int$PLANTSPECIES[i.obs],sp.int$BIRDSPECIES[i.obs]] <- 1
  adjacency.matrix[sp.int$BIRDSPECIES[i.obs],sp.int$PLANTSPECIES[i.obs]] <- 1
}# for i.obs

plant.sp.cells <- subset(sp.cells,species %in% plant.sp)
plant.sp.cells.wide <- pivot_wider(plant.sp.cells,
                                   names_from = cell_id,
                                   values_from = observations,
                                   # values_from = presence,
                                   values_fill = 0)

bird.sp.cells <- subset(sp.cells,species %in% bird.sp)
bird.sp.cells.wide <- pivot_wider(bird.sp.cells,
                                  names_from = cell_id,
                                  values_from = observations,
                                  # values_from = presence,
                                  values_fill = 0)

# -------------------------------------------------------------------------
# this list will hold the aggregated communicability for each node (i.e.
# a local population of a given species)
# population.communicability.list <- vector(mode = "list", 
#                                           length = length(cell.id))

# -------------------------------------------------------------------------
# for each cell, obtain their neighbours within max.disp distance, 
# populate this subset of the landscape (including dispersal links),
# and obtain communicability of each local population

# tests
# i.cell <- 453
# cell.id <- 1:10

null.list <- foreach(i.cell = 1:length(cell.id), 
        # .combine=comb.fun, 
        .packages = 'tidyverse') %dopar% {
  
  source("/home/david/Work/Projects/NZ/propagation/R/communicability_network.R")
  source("/home/david/Work/Projects/NZ/propagation/R/communicability.R")
  
# for(i.cell in 1:length(cell.id)){
  
  # 1) create sub-landscape matrix and populate it
  # this matrix is "sublandscape.matrix"
  
  my.cell.distances <- NZ.distances[cell.id[i.cell],]
  sublandscape.ids <- which(my.cell.distances <= max.dist)
  
  # create sub-landscape
  matrix.n.rows <- length(all.sp)*length(sublandscape.ids)
  df.names <- expand.grid(all.sp,names(sublandscape.ids))
  matrix.names <- paste(df.names[,1],df.names[,2],sep="-")
  
  sublandscape.matrix <- matrix(0,
                         nrow = matrix.n.rows,
                         ncol = matrix.n.rows,
                         dimnames = list(matrix.names,matrix.names))
  
  # -------------------------------------------------------------------------
  # populate the sub-landscape
  
  # conditions for linking:
  # 1 - plant-bird interaction in the same cell
  
  # go through each cell of the sublandscape
  # (i.e. the diagonal of the block matrix)
  for(i.sub in 1:length(sublandscape.ids)){
    
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
          
          my.cell <- as.character(sublandscape.ids[i.sub])
          
          # check presences in this cell
          # check again consistency of cells - there is at least one cell
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
    init.row <- 1 + (num.sp * (i.sub - 1))
    init.col <- init.row
    end.row <- num.sp + (num.sp * (i.sub - 1))
    end.col <- end.row
    
    sublandscape.matrix[init.row:end.row,init.col:end.col] <- my.adj.matrix
    
  }# for i.sub
  
  # 2 - only for birds: population of the same species in a reachable cell
  # 2.1 set of bird species in this sublandscape
  birds.sublandscape <- sort(unique(bird.sp.cells$species[which(bird.sp.cells$cell_id %in% sublandscape.ids)]))
  
  # dispersal distances of these birds
  # birds.disp.sub <- subset(sp.disp,species %in% birds.sublandscape)
  sub.distances <- NZ.distances[sublandscape.ids,sublandscape.ids]
  
  # convert sublandscape distances to cell pairwise distances
  sub.dist.list <- sub.distances %>% 
    as_tibble(rownames = "cell1") %>%
    pivot_longer(-cell1,names_to = "cell2",values_to = "distance")
  sub.dist.list <- subset(sub.dist.list, distance > 0)
  
  # redo - simply apply the dispersal function to each pair of cells 
  # at a certain distance
  
  for(i.pair in 1:nrow(sub.dist.list)){
    for(i.bird in 1:length(birds.sublandscape)){
      
      my.sp <- birds.sublandscape[i.bird]
      my.cell.1 <- sub.dist.list$cell1[i.pair]
      my.cell.2 <- sub.dist.list$cell2[i.pair]
      bird.cells <- bird.sp.cells$cell_id[bird.sp.cells$species == my.sp]
      
      # 1) is this bird present in both cells?
      if(my.cell.1 %in% bird.cells & my.cell.2 %in% bird.cells){
        # cat(i.pair,"-",i.bird,"\n")
        my.distance <- sub.distances[as.character(my.cell.1),as.character(my.cell.2)]
        my.hwi <- hwi$hand.wing.index[which(hwi$species == my.sp)]
        
        my.dispersal <- disp.fun(my.distance,my.hwi)
          
        # assign dispersal links
        name.cell.1 <- paste(my.sp,"-",my.cell.1,sep="")
        name.cell.2 <- paste(my.sp,"-",my.cell.2,sep="")
        
        sublandscape.matrix[name.cell.1,name.cell.2] <- my.dispersal
        sublandscape.matrix[name.cell.2,name.cell.1] <- my.dispersal
        # cat(name.cell.1,"-",name.cell.2,"\n")
      }# if present
      
    }# for i.bird
  }# for i.pair
  
  # 2) communicability of each local population from the focal cell
  # i.e. i.cell
  # the thing is I cannot obtain these communicabilities alone, I need the whole
  # sublandscape
  sub.comm <- communicability_network(sublandscape.matrix,
                                      weighted = T,
                                      return.pairwise.comm = T)
  bin.comm <- sub.comm[[3]]
  w.comm <- sub.comm[[4]]
  
  # tidy functions only work with dataframes, but this still works, and is fast
  # dfb <- reshape2::melt(bin.comm,value.name = "binary.communicability")
  # dfw <- reshape2::melt(w.comm,value.name = "weighted.communicability")
  dfb <- as.data.frame(bin.comm) %>%
    mutate(Var1 = rownames(.)) %>%
    pivot_longer(cols = c(-Var1),names_to = "Var2",values_to = "binary.communicability")
  dfw <- as.data.frame(w.comm) %>%
    mutate(Var1 = rownames(.)) %>%
    pivot_longer(cols = c(-Var1),names_to = "Var2",values_to = "weighted.communicability")
  
  df1 <- left_join(dfb,dfw)
  
  df1$sp1 <- sub("\\-.*", "", df1$Var1)
  df1$cell.id.sp1 <- sub(".*-", "", df1$Var1)
  df1$sp2 <- sub("\\-.*", "", df1$Var2)
  df1$cell.id.sp2 <- sub(".*-", "", df1$Var2)  
  
  df1$scaled.binary.communicability <- scales::rescale(df1$binary.communicability)
  df1$scaled.weighted.communicability <- scales::rescale(df1$weighted.communicability)
  
  df1$diag <- ifelse(df1$sp1 == df1$sp2 & df1$cell.id.sp1 == df1$cell.id.sp2,TRUE,FALSE)
  
  df1$guild.sp1 <- ifelse(df1$sp1 %in% bird.sp, "birds", "plants")
  df1$guild.sp2 <- ifelse(df1$sp2 %in% bird.sp, "birds", "plants")
  
  df1 <- df1[,c("sp1","guild.sp1","cell.id.sp1",
                "sp2","guild.sp2","cell.id.sp2",
                "diag",
                "binary.communicability",
                "scaled.binary.communicability",
                "weighted.communicability",
                "scaled.weighted.communicability")]
  
  # I am only interested in the values for the focal cell - 
  # other cells are only partially represented, because they may be in the
  # border of this "sub-landscape".
  cell.df <- df1 %>% 
    filter(cell.id.sp1 == cell.id[i.cell]) %>%
    filter(diag == FALSE) %>%
    group_by(sp1,guild.sp1,cell.id.sp1) %>%
    summarise(population.bin.communicability = sum(binary.communicability),
              population.weighted.communicability = sum(weighted.communicability))
  
    # -------------------------------------------------------------------------

    write.csv2(cell.df,paste("results/communicability/NZ_networks/focal_",i.cell,"_",grid.size,"km.csv",sep=""))

    # -------------------------------------------------------------------------
    # there are two options: either writing the output of each iteration
    # independently, as above, or collecting all cell.df in a list as returned
    # by foreach. I don't want that second option, so I return NULL from each
    # iteration, simply to avoid saving the list elements.
    NULL

}# for i.cell

stopCluster(cl)


# population.communicability.df <- bind_rows(population.communicability.list)
# 
# # -------------------------------------------------------------------------
# 
# write.csv2(population.communicability.df,
#            paste("results/population_communicability_",grid.size,"km.csv",sep=""),
#            row.names = F)






