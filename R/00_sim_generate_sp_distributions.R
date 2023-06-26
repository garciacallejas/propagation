
# script to assign species presences to a simulated landscape

# INPUTS
# - list of simulated landscapes: "results/landscape_matrices.RData"
# - list of metawebs: "results/sim_degree_dist_networks.csv"
# - dataframe of species landscape suitability: "results/sp_suitability.csv"
# - dataframe of cell coordinates: "results/cell_coordinates.csv"

# OUTPUTS
# - list of dataframes with species presence per cell: 
# sp.presence[[landscape.category]][[generative.mdel]][[replicate]] 
# "results/sp_presence.RData"

# - list of supra-adjacency matrices without dispersal (interlayer links): 
# landscape.matrices[[landscape.category]][[generative.model]][[replicate]]
# "results/landscape_matrices.RData"
# -------------------------------------------------------------------------

# library(tidyverse)

# -------------------------------------------------------------------------
suit.df <- read.csv2("results/sp_suitability.csv")
load("results/landscape_matrices.RData")

# use the networks generated from a degree distribution gradient
sim.networks <- read.csv2("results/sim_degree_dist_networks.csv")

# load adjacency matrices
load("results/sim_degree_dist_network_matrices.RData")

cell.coords <- read.csv2("results/cell_coordinates.csv")
# -------------------------------------------------------------------------

landscape.categories <- names(landscape.list)
network.categories <- unique(sim.networks$generative.model)
network.replicates <- unique(sim.networks$replicate)
landscape.replicates <- 1:length(landscape.list[[1]])

landscape.rows <- nrow(landscape.list[[1]][[1]]) 
landscape.cols <- ncol(landscape.list[[1]][[1]])

cells <- landscape.rows * landscape.cols

sp.names <- unique(suit.df$sp)
num.sp <- length(sp.names)

# NOTE: I am assuming network.replicates == landscape.replicates, so that I generate
# that number of landscape matrices. Otherwise I would have to cross both, networks * landscapes,
# and that would mean even more computation and storage. 

if(identical(network.replicates,landscape.replicates)){
  replicates <- network.replicates
}else{
  replicates <- NA
}

# -------------------------------------------------------------------------

matrix.n.rows <- num.sp*cells
df.names <- expand.grid(sp.names,1:cells)
matrix.names <- paste(df.names[,1],df.names[,2],sep="-")

# create landscape matrix template
landscape.matrix.template <- matrix(0,
                                    nrow = matrix.n.rows,
                                    ncol = matrix.n.rows,
                                    dimnames = list(matrix.names,matrix.names))

# -------------------------------------------------------------------------

# dataframes of species presences per cell
sp.presence <- list()

# matrices
landscape.matrices <- list()

# go through each landscape category, network category, and replicate
for(i.land in 1:length(landscape.categories)){
  
  sp.presence[[i.land]] <- list()
  landscape.matrices[[i.land]] <- list()
  
  for(i.net in 1:length(network.categories)){
    
    sp.presence[[i.land]][[i.net]] <- list()
    landscape.matrices[[i.land]][[i.net]] <- list()
    
    for(i.rep in replicates){
      
      my.landscape.matrix <- landscape.matrix.template
      my.landscape <- landscape.list[[i.land]][[i.rep]]
      
      # list holding presences per cell
      presence.df <- list()
      
      for(i.row in 1:landscape.rows){
        for(i.col in 1:landscape.cols){
          
          # "reload" the metaweb for this replicate
          my.network <- sim.matrices[[i.net]][[i.rep]]
          
          # cells are referred to by its row and col numbers, but I need
          # a cell number for the block matrix (see below)
          cell.id <- (landscape.cols * (i.row - 1)) + i.col
          
          # suitability value of this cell
          landscape.value <- landscape.list[[i.land]][[i.rep]][i.row,i.col]
          
          # id of this cell
          cell.value <- cell.coords$cell[which(cell.coords$x == i.col &
                                                 cell.coords$y == i.row)]
          
          # presence dataframe for this cell
          cell.df <- tidyr::expand_grid(landscape.category = landscape.categories[i.land],
                                 network.category = network.categories[i.net],
                                 replicate = i.rep,
                                 landscape.row = i.row,
                                 landscape.col = i.col,
                                 cell = cell.value,
                                 sp = sp.names,
                                 presence = FALSE)
          
          # go through all sp, checking if it should be present in this cell
          # and update the df and network matrix of the cell
          
          # a species is present if the suitability value of the cell
          # is within its mean +- sd 
          
          for(i.sp in 1:num.sp){
            my.suit <- which(suit.df$sp == sp.names[i.sp])
            min.suit <- suit.df$optimum[my.suit] - suit.df$sd[my.suit]
            max.suit <- suit.df$optimum[my.suit] + suit.df$sd[my.suit]
            
            if(landscape.value >= min.suit & landscape.value <= max.suit){
              cell.df$presence[i.sp] <- TRUE
            }else{
              # if species not present, "prune" the matrix, setting the species'
              # elements to 0
              
              my.network[sp.names[i.sp],] <- 0
              my.network[,sp.names[i.sp]] <- 0
              
            }
          }# for i.sp
          presence.df[[length(presence.df)+1]] <- cell.df
          
          # now, add the cell network to the landscape matrix
          init.row <- 1 + (num.sp * (cell.id - 1))
          init.col <- init.row
          end.row <- num.sp + (num.sp * (cell.id - 1))
          end.col <- end.row

          my.landscape.matrix[init.row:end.row,init.col:end.col] <- my.network
          
        }# i.col
      }# i.row
      
      # append to matrices list
      # landscape.matrices[[i.land]][[i.net]][[i.rep]] <- my.landscape.matrix
      
      # save this matrix
      save(my.landscape.matrix,file = paste("results/sim_landscape_matrices/",network.categories[i.net],
                                            "_",landscape.categories[i.land],
                                            "_re",i.rep,
                                            "_",num.sp,"sp",
                                            "_",cells,"cells.RData",
                                            sep=""))

      # append to overall presences list
      sp.presence[[i.land]][[i.net]][[i.rep]] <- dplyr::bind_rows(presence.df)
    }# for i.rep
  }# for i.net
  names(sp.presence[[i.land]]) <- network.categories
}# for i.land
names(sp.presence) <- landscape.categories

# -------------------------------------------------------------------------

save(sp.presence,file = "results/presence_dataframes.RData")


