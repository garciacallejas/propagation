
# populates landscape matrices with dispersal links

# INPUTS
# - nested lists of presence/absence dataframes for every combination of
# landscape/network categories: sp.presence[[landscape.category]][[network.category]][[replicate]]
# "results/presence_dataframes.RData"
# - landscape matrices with local communities "wired"
# "results/sim_landscape_matrices/NETWORKCATEGORY_LANDSCAPECATEGORY_REPLICATE_RICHNESS_CELLS.RData
# - dispersal distances dataframe:
# "results/dispersal_kernels.csv"
# - cell coordinates and distances: "results/cell_coordinates.csv","results/cell_distances.csv"

# OUTPUTS
# - andscape matrices with dispersal, each stored individually
# results/sim_landscape_matrices/NETWORK_LANDSCAPE_DISPERSAL_REPLICATE_RICHNESS_CELLS.RData
# - dataframe with combinations of categories
# "results/landscapes_info.csv"

# -------------------------------------------------------------------------

library(tidyverse)
library(igraph)
library(NLMR)
library(sf)

# -------------------------------------------------------------------------
landscapes.path <- "results/sim_landscape_matrices/"

load("results/presence_dataframes.RData")
disp.df <- read.csv2("results/dispersal_kernels.csv")

cell.distances <- read.csv2("results/cell_distances.csv")
cell.coords <- read.csv2("results/cell_coordinates.csv")

# -------------------------------------------------------------------------
# get some parameters

landscape.categories <- sort(unique(names(sp.presence)))
network.categories <- sort(unique(names(sp.presence[[1]])))
dispersal.categories <- sort(unique(disp.df$dispersal.category))
replicates <- length(sp.presence[[1]][[1]])

landscape.rows <- max(sp.presence[[1]][[1]][[1]]$landscape.row)
landscape.cols <- max(sp.presence[[1]][[1]][[1]]$landscape.col)
cells <- landscape.rows * landscape.cols
sp.names <- sort(unique(disp.df$sp))
richness <- length(sp.names)

# -------------------------------------------------------------------------
# iterate through each generated landscape to add dispersal links

# i.net <- i.land <- i.rep <- i.disp <- 1

for(i.land in 1:length(landscape.categories)){
  for(i.net in 1:length(network.categories)){
    for(i.rep in 1:replicates){
      
      # -------------------------------------------------------------------
      # load landscape matrix
      
      landscape.name <- paste(network.categories[i.net],"_",
                              landscape.categories[i.land],"_",
                              "re",i.rep,"_",
                              richness,"sp_",
                              cells,"cells.RData",sep="")
      
      # my.landscape.matrix
      load(paste(landscapes.path,landscape.name,sep=""))
      
      # -------------------------------------------------------------------
      # add dispersal
      for(i.disp in 1:length(dispersal.categories)){
        
        landscape <- my.landscape.matrix
        
        # -----------------------------------------------------------------
        # obtain dispersal links
        
        my.disp <- subset(disp.df, 
                          dispersal.category == dispersal.categories[i.disp] &
                            replicate == i.rep)
        
        # expand the presence dataframe to see the cells to which
        # each sp can potentially disperse
        # by crossing presence information with cell distances
        # and with dispersal distances
        
        my.presence <- sp.presence[[i.land]][[i.net]][[i.rep]]
        my.presence <- subset(my.presence,presence == TRUE) %>%
          dplyr::select(sp,cell) %>%
          rename(cell_from = cell)
        my.presence.full <- expand_grid(my.presence,cell_to = 1:cells)
        my.presence.full.2 <- left_join(my.presence.full,cell.distances)
        my.presence.full.3 <- left_join(my.presence.full.2,my.disp[,c("sp","dispersal.distance")])
        my.presence.full.3 <- subset(my.presence.full.3,cell_from != cell_to)
        
        my.presence.full.3$dispersal.potential <- ifelse(my.presence.full.3$distance <= 
                                                           my.presence.full.3$dispersal.distance,
                                                         TRUE,FALSE)
        my.dispersal.potential <- subset(my.presence.full.3[,c("sp","cell_from",
                                                               "cell_to","dispersal.potential")],
                                         dispersal.potential == TRUE)
        
        # now, check which of the target cells have populations of the species
        # because only those will have realized dispersal
        my.presence$presence <- TRUE
        my.dispersal.realized <- left_join(my.dispersal.potential,my.presence,
                                           by = c("cell_to" = "cell_from",
                                                  "sp" = "sp")) %>%
          replace_na(list(presence = FALSE))
        my.dispersal.realized$dispersal <- as.logical(my.dispersal.realized$dispersal.potential * 
                                                        my.dispersal.realized$presence)
        my.dispersal.r2 <- subset(my.dispersal.realized[,c("sp","cell_from","cell_to","dispersal")],
                                  dispersal == TRUE)
        
        # remove symmetrical information
        my.dispersal.r2$duplicated <- FALSE
        for(i in 1:nrow(my.dispersal.r2)){
          if(!my.dispersal.r2$duplicated[i]){
            dup <- which(my.dispersal.r2$cell_from == my.dispersal.r2$cell_to[i] &
                           my.dispersal.r2$cell_to == my.dispersal.r2$cell_from[i] &
                           my.dispersal.r2$sp == my.dispersal.r2$sp[i])
            my.dispersal.r2$duplicated[dup] <- TRUE
          }# if
        }# for i
        
        # this is the distilled dataframe
        # including the dispersal coefficient to assign
        # i.e. 1/number of dispersing cells (irrespective of distance)
        
        # this may change in the future
        
        my.dispersal <- subset(my.dispersal.r2,duplicated == FALSE) %>%
          group_by(sp,cell_from) %>%
          mutate(dispersal.coef = 1/n()) %>%
          dplyr::select(sp,cell_from,cell_to,dispersal.coef)
        
        # ----------------------------------------------------------------
        # update landscape matrix
        
        # I only need to go through the upper triangle, since 
        # dispersal is symmetric
        # hence the weird nested for loops
        for(i.row in 1:(nrow(landscape)-1)){
          
          # which sp?
          my.sp <- i.row %% richness
          if(my.sp == 0){my.sp <- richness}
          
          for(i.col in (i.row+1):ncol(landscape)){
            
            # if diagonal element and not main diagonal,
            # it is a dispersal coefficient
            # if(i.row %% richness == i.col %% richness & i.row != i.col){
            if(i.row %% richness == i.col %% richness){
              
              # source and dest cell
              source.cell <- ceiling(i.col/richness)
              dest.cell <- ceiling(i.row/richness)
              
              # double check
              if(source.cell != dest.cell){
                
                valid.dispersal <- which(my.dispersal$sp == sp.names[my.sp] &
                                           (my.dispersal$cell_from == source.cell &
                                              my.dispersal$cell_to == dest.cell | 
                                              my.dispersal$cell_from == dest.cell &
                                              my.dispersal$cell_to == source.cell) )
                if(length(valid.dispersal) == 1){
                  
                  disp.coef <- my.dispersal$dispersal.coef[valid.dispersal]
                  
                  # fill the symmetric positions
                  landscape[i.row,i.col] <- disp.coef
                  landscape[i.col,i.row] <- disp.coef
                }# if realized dispersal
                
              }# if different cell
            }# if dispersal cell
            
          }# for i.col
        }# for i.row
        
        # -----------------------------------------------------------------
        # store landscape with dispersal
        
        # this is the name of the resulting landscape matrix
        landscape.dispersal.name <- paste(network.categories[i.net],"_",
                                          landscape.categories[i.land],"_",
                                          dispersal.categories[i.disp],"_",
                                          "re",i.rep,"_",
                                          richness,"sp_",
                                          cells,"cells.RData",sep="")
        
        save(landscape, file = paste(landscapes.path,landscape.dispersal.name,sep=""))
        
      }# for i.disp
    }# for i.rep
  }# for i.net
}# for i.land





