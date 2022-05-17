# script to generate simulated networks from specified degree distributions - 
# sample them from Poisson distributions with increasing mean.

# INPUTS
# - richness
# - min, max lambda parameters of the Poisson distribution
# - categories in the gradient
# - replicates per category

# OUTPUTS
# - nested list of adjacency matrices: sim.degree.dist.matrices[[model]][[replicate]]
# "results/sim_degree_dist_network_matrices.Rdata"
# - dataframe with edge list of every network: "results/sim_degree_dist_networks.csv"
# - dataframe with network categories and the lambda parameter of each
# -------------------------------------------------------------------------

library(tidyverse)
library(extraDistr)
library(igraph)
library(gamlss.dist)

# -------------------------------------------------------------------------
param <- read.csv2("results/sim_landscape_matrices/parameters_v2.csv")

richness <-  param$richness
num.network.categories <- param$num.network.categories
num.category.replicates <- param$num.category.replicates

# poisson mean
min.lambda <- param$min.lambda
max.lambda <- param$max.lambda # this should vary with richness. for S = 50, 15 gives connectance = 0.3

# some constants for sampling interaction strengths
int.mean <- param$int.mean
int.sd <- param$int.sd
tau <- param$tau
diag.dom <- param$diag.dom
# -------------------------------------------------------------------------

sp.names <- paste("sp",sprintf("%02d", 1:richness),sep="")

degree.dist.gradient <- seq(from = min.lambda,
                            to = max.lambda, 
                            length.out = num.network.categories)

# generative.models <- paste("dd",1:num.network.categories,sep="")
generative.models <- paste("dd",sprintf("%02d", 1:num.network.categories),sep="")

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
connectances <- numeric(length(generative.models))

for(i in 1:length(generative.models)){
  sim.matrices[[i]] <- list()
  connectance.rep <- numeric(num.category.replicates)
  
  for(j in 1:num.category.replicates){
    
    my.dist <- extraDistr::rtpois(n = richness,lambda = degree.dist.gradient[i],a = 0)
    # make sum even
    if (sum(my.dist) %% 2 != 0) { my.dist[1] <- my.dist[1] + 1 }
    
    my.net <- igraph::sample_degseq(my.dist,method = "vl")
    my.matrix <- as.matrix(igraph::as_adjacency_matrix(my.net,type = "both"))
    # cat("lambda:",degree.dist.gradient[i],"- connectance:",(sum(my.matrix)/richness^2),"\n")
    connectance.rep[j] <- sum(my.matrix)/richness^2
    
    # assign interaction strengths according to an "extended" normal dist
    weights <- abs(gamlss.dist::rSHASHo(sum(my.dist), mu = int.mean, 
                                        sigma = int.sd, nu = 0, tau = tau))
    my.matrix[my.matrix == 1] <- weights
    diag(my.matrix) <- 1
    
    sim.matrices[[i]][[j]] <- my.matrix
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
  connectances[i] <- mean(connectance.rep)
}
names(sim.matrices) <- generative.models

sim.networks <- bind_rows(edge.lists)

# -------------------------------------------------------------------------

save(sim.matrices,file = "results/sim_degree_dist_network_matrices.RData")
write.csv2(sim.networks,"results/sim_degree_dist_networks.csv",row.names = FALSE)
write.csv2(data.frame(network.category = generative.models,
                      poisson.lambda = round(degree.dist.gradient,2),
                      connectance = round(connectances,3)),
           "results/network_gradient_categories.csv",
           row.names = FALSE)


