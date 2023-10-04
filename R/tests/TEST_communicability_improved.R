
# communicability improved - tests

library(tidyverse)
source("R/horizontal_community_matrix.R")
source("R/communicability.R")
source("R/auxiliary_functions/comm.R")
library(intsegration)
library(igraph)

# -------------------------------------------------------------------------

A <- matrix(runif(25),nrow = 5,dimnames = list(paste("sp",1:5,sep=""),
                                               paste("sp",1:5,sep="")))
# diag(A) <- 0
A[A < .7] <- 0

comm(A,normalised = T,weighted = F)

# A <- matrix(c(0,1.5,0,0,
#               1.5,0,4,1,
#               0,4,0,2,
#               0,1,2,0),nrow = 4,byrow = T,dimnames = list(c("u","q","r","v"),
#                                                           c("u","q","r","v")))

# -------------------------------------------------------------------------

# average communicability

# com.list <- communicability(A)
# bin.comm <- com.list[[1]]
# w.comm <- com.list[[2]]
# 
# # weird - this gives a lot of Inf values
# # g_all_shortest_paths <- intsegration::rcpp_floyd_flow(1/A)
# # Phi <- g_all_shortest_paths$F
# # diag(Phi) <- NA
# 
# # turns out to be because of the discontinuity of the Inf.
# inv.A <- 1/A
# # this works
# inv.A[is.infinite(inv.A)] <- .Machine$double.xmax
# all.shortest.pahts <- intsegration::rcpp_floyd_flow(inv.A)
# # flows
# phi.A <- all.shortest.pahts$F
# diag(Phi.A) <- 0
# 
# comm.phi <- communicability(phi.A)
# bin.comm.phi <- comm.phi[[1]]
# w.comm.phi <- comm.phi[[2]]
# 
# comm_ideal <- 1. / N / (N - 1) * sum(bin.comm.phi, na.rm = T)
# comm_obs <- 1. / N / (N - 1) * sum(bin.comm, na.rm = T)
# 
# comm_obs/comm_ideal

# follow the gce function from bertagnolli
# g <- igraph::graph_from_adjacency_matrix(A,weighted = T)
# igraph::E(g)$weight_inv <- 1. / igraph::E(g)$weight
# N <- igraph::gorder(g)
# A2 <- igraph::as_adjacency_matrix(g, attr = "weight")
# 
# x <- as.matrix(igraph::as_adjacency_matrix(g, attr = "weight_inv"))
# x[x == 0] <- .Machine$double.xmax         # almost as <- Inf
# g_all_shortest_paths2 <- intsegration::rcpp_floyd_flow(x)
# Phi2 <- g_all_shortest_paths2$F
# diag(Phi2) <- NA
# E_ideal <- 1. / N / (N - 1) * sum(g_all_shortest_paths2$F, na.rm = T)

