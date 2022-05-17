# example of global communication efficiency from Bertagnolli et al.
# https://doi.org/10.1038/s42005-021-00612-5


#devtools::install_github("gbertagnolli/intsegration")
library(intsegration)
library(gamlss.dist)
library(igraph)
# devtools::install_github("RadicalCommEcol/MultitrophicFun")
library(MultitrophicFun)
library(Rcpp)
# -------------------------------------------------------------------------

S <- 25
c <- .8

A <- horizontal_community_matrix(S = S,c = c,restricted.positive = T)
net.A <- igraph::graph_from_adjacency_matrix(A,weighted = T)
my.GCE <- my_GCE(g = net.A,directed = T,normalised = T)







