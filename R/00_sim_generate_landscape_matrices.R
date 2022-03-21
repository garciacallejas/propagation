
# script to generate simulated landscapes from a Null Landscape Model,
# according to a gradient in spatial autocorrelation

# INPUTS
# - number of columns and rows
# - categories in the gradient
# - replicates per category

# OUTPUTS
# - nested list: landscape.list[[category]][[replicate]]: "results/landscape_matrices.RData"

# -------------------------------------------------------------------------

library(tidyverse)

# remotes::install_github("ropensci/NLMR")
library(NLMR)

# remotes::install_github("ropensci/landscapetools")
library(landscapetools)

# for converting to matrix
library(raster)

# -------------------------------------------------------------------------
# generate a grid landscape with 0-1 habitat values in each cell
ncol = 50
nrow = 50

num.landscape.categories <- 10
num.category.replicates <- 10

# -------------------------------------------------------------------------

spatial.autocorr <- seq(0.1,0.9, length.out = num.landscape.categories)

# generate nested list
landscape.list <- list()

for(i in 1:length(spatial.autocorr)){
  landscape.list[[i]] <- list()
  
  for(i.rep in 1:num.category.replicates){
    my.landscape <- nlm_fbm(ncol = ncol,nrow = nrow,fract_dim = spatial.autocorr[i])
    my.matrix <- raster::as.matrix(my.landscape)
    landscape.list[[i]][[i.rep]] <- my.matrix
  }# for each replicate
  
}# for each landscape category
names(landscape.list) <- paste("sa_",spatial.autocorr,sep="")
# show_landscape(landscape.list,n_col = 3)

save(landscape.list,file = "results/landscape_matrices.RData")
