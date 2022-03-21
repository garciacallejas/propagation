
# script to generate simulated networks according to different generative models

# INPUTS
# - richness
# - for increasing connectance: min, max values
# - categories in the gradient
# - replicates per category

# OUTPUTS
# - nested list of adjacency matrices: sim.matrices[[model]][[replicate]]
# "results/sim_network_matrices.Rdata"
# - dataframe with edge list of every network: "results/sim_networks.csv"

# -------------------------------------------------------------------------

library(tidyverse)
library(igraph)
source("R/horizontal_community_matrix.R")

# -------------------------------------------------------------------------
# start with random matrices of increasing connectance

richness <-  50
num.categories <- 10
num.category.replicates <- 10
min.connectance <- 0.1
max.connectance <- 0.9

# some constants
tau <- 1.5
diag.dom <- 0
# -------------------------------------------------------------------------

sp.names <- paste("sp",1:richness,sep="")

connectance.gradient <- seq(from = min.connectance,
                            to = max.connectance, 
                            length.out = num.categories)

generative.models <- paste("c",1:num.categories,sep="")

# this dataframe will hold the edge list of each matrix
sim.networks <- expand_grid(generative.model = generative.models,
                           replicate = 1:num.category.replicates,
                           node.from = sp.names,
                           node.to = sp.names,
                           value = 0) 

# this list will hold the actual matrices
sim.matrices <- list()
edge.lists <- list()
# -------------------------------------------------------------------------

for(i in 1:length(generative.models)){
  sim.matrices[[i]] <- list()
  for(j in 1:num.category.replicates){
    sim.matrices[[i]][[j]] <- horizontal_community_matrix(S = richness,
                                                          c = connectance.gradient[i],
                                                          tau = tau,
                                                          min.diag.dom = diag.dom,
                                                          restricted.positive = TRUE)
    dimnames(sim.matrices[[i]][[j]]) <- list(sp.names,sp.names)
    # fill up edge list
    my.edge.list <- expand_grid(generative.model = generative.models[i],
                                replicate = j,
                                node_from = sp.names,
                                node_to = sp.names,
                                value = 0)
    
    for(i.link in 1:nrow(my.edge.list)){
      my.edge.list$value[i.link] <- sim.matrices[[i]][[j]][my.edge.list$node_from[i.link],my.edge.list$node_to[i.link]]
    }
    edge.lists[[length(edge.lists)+1]] <- my.edge.list
  }
}
names(sim.matrices) <- generative.models

sim.networks <- bind_rows(edge.lists)

# -------------------------------------------------------------------------

save(sim.matrices,file = "results/sim_network_matrices.RData")
write.csv2(sim.networks,"results/sim_networks.csv",row.names = FALSE)




