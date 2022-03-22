
# script to assign species presences to a simulated landscape

# INPUTS
# - list of simulated landscapes: "results/landscape_matrices.RData"
# - list of metawebs: "results/sim_degree_dist_networks.csv"
# - dataframe of species landscape suitability: "results/sp_suitability.csv"

# OUTPUTS
# - list of dataframes with species presence per cell: 
# sp.presence[[landscape.category]][[generative.mdel]][[replicate]] 
# "results/sp_presence.RData"

# - list of supra-adjacency matrices without dispersal (interlayer links): 
# landscape.matrices[[landscape.category]][[generative.model]][[replicate]]
# "results/landscape_matrices.RData"
# -------------------------------------------------------------------------

suit.df <- read.csv2("results/sp_suitability.csv")
load("results/landscape_matrices.RData")

# use the networks generated from a degree distribution gradient
sim.networks <- read.csv2("results/sim_degree_dist_networks.csv")

# -------------------------------------------------------------------------

landscape.categories <- names(landscape.list)
network.categories <- unique(sim.networks$generative.model)
network.replicates <- unique(sim.networks$replicate)
landscape.replicates <- 1:length(landscape.list[[1]])

landscape.rows <- nrow(landscape.list[[1]][[1]]) 
landscape.cols <- ncol(landscape.list[[1]][[1]])

cells <- landscape.rows * landscape.cols

sp.names <- unique(suit.df$sp)

# NOTE: I am assuming network.replicates == landscape.replicates, so that I generate
# that number of landscape matrices. Otherwise I would have to cross both, networks * landscapes,
# and that would mean even more computation and storage. 

if(identical(network.replicates,landscape.replicates)){
  replicates <- network.replicates
}else{
  replicates <- NA
}

# -------------------------------------------------------------------------

# TODO 50 sp * 50x50 cells is TOO BIG!!!

# create landscape matrix template
landscape.matrix.template <- matrix(0,nrow = cells * length(sp.names),
                                    ncol = cells * length(sp.names))

# -------------------------------------------------------------------------

# dataframes of species presences per cell
sp.presence <- list()

# matrices
landscape.matrices <- list()

# TODO for now, only sp.presence is filled

# go through each landscape category, network category, and replicate
for(i.land in 1:length(landscape.categories)){
  
  sp.presence[[i.land]] <- list()
  landscape.matrices[[i.land]] <- list()
  
  for(i.net in 1:length(network.categories)){
    
    sp.presence[[i.land]][[i.net]] <- list()
    landscape.matrices[[i.land]][[i.net]] <- list()
    
    for(i.rep in replicates){
      
      my.landscape <- landscape.list[[i.land]][[i.rep]]
      my.network <- subset(sim.networks, generative.model == network.categories[i.net] &
                             replicate == i.rep)
      
      # list holding presences per cell
      presence.df <- list()
      
      for(i.row in 1:landscape.rows){
        for(i.col in 1:landscape.cols){
          landscape.value <- landscape.list[[i.land]][[i.rep]][i.row,i.col]
          
          cell.df <- expand_grid(landscape.category = landscape.categories[i.land],
                                 network.category = network.categories[i.net],
                                 replicate = i.rep,
                                 landscape.row = i.row,
                                 landscape.col = i.col,
                                 sp = sp.names,
                                 presence = FALSE)
          for(i.sp in 1:length(sp.names)){
            my.suit <- which(suit.df$sp == sp.names[i.sp])
            min.suit <- suit.df$optimum[my.suit] - suit.df$sd[my.suit]
            max.suit <- suit.df$optimum[my.suit] + suit.df$sd[my.suit]
            
            if(landscape.value >= min.suit & landscape.value <= max.suit){
              cell.df$presence[i.sp] <- TRUE
            }
            
          }# for i.sp
          presence.df[[length(presence.df)+1]] <- cell.df
        }# i.col
      }# i.row
      
      # append to overall presences list
      sp.presence[[i.land]][[i.net]][[i.rep]] <- bind_rows(presence.df)
    }# for i.rep
    
    names(sp.presence[[i.land]]) <- network.categories
    names(landscape.matrices[[i.land]]) <- network.categories
  }# for i.net
  
  names(sp.presence) <- landscape.categories
  names(landscape.matrices) <- landscape.categories
}# for i.land





