
library(tidyverse)
library(extraDistr)
library(igraph)
library(gamlss.dist)

# -------------------------------------------------------------------------

richness <-  30
num.categories <- 10
num.category.replicates <- 10

# exponetial decay rate
min.rate <- .75
max.rate <- .25 

# -------------------------------------------------------------------------

sp.names <- paste("sp",1:richness,sep="")

dispersal.rate.gradient <- seq(from = min.rate,
                            to = max.rate, 
                            length.out = num.categories)

# this list will hold the actual matrices
sim.matrices <- list()
edge.lists <- list()
# -------------------------------------------------------------------------

for(i in 1:nrow(dispersal.rate.gradient)){

my.dist <- rexp(n = richness,rate = dispersal.rate.gradient[i])

}