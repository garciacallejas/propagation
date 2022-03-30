
# script to generate simulated landscapes from a Null Landscape Model,
# according to a gradient in spatial autocorrelation

# INPUTS
# - number of columns and rows
# - categories in the gradient
# - replicates per category

# OUTPUTS
# - nested list: landscape.list[[category]][[replicate]]: "results/landscape_matrices.RData"
# - landscape categories and their autocorrelation values: "results/spatial_autocorrelation_categories.csv"
# - dataframe with cell coordinates and label: "results/cell_coordinates.csv"
# - distance dataframe between each pair of cells: "results/cell_distances.csv"

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
ncol = 20
nrow = 20

num.landscape.categories <- 10
num.category.replicates <- 10

# -------------------------------------------------------------------------

spatial.autocorr <- round(seq(0.1,0.9, length.out = num.landscape.categories),2)

# generate nested list
landscape.list <- list()

# check that mean suitability is conserved across categories
mean.suit <- expand_grid(spatial.autocorr = spatial.autocorr,replicate = 1:num.category.replicates,mean.suit = NA)

for(i in 1:length(spatial.autocorr)){
  landscape.list[[i]] <- list()
  
  for(i.rep in 1:num.category.replicates){
    my.landscape <- nlm_fbm(ncol = ncol,nrow = nrow,fract_dim = spatial.autocorr[i])
    # show_landscape(my.landscape)
    my.matrix <- raster::as.matrix(my.landscape)
    landscape.list[[i]][[i.rep]] <- my.matrix
    mean.suit$mean.suit[which(mean.suit$spatial.autocorr == spatial.autocorr[i] & mean.suit$replicate == i.rep)] <- mean(my.matrix)
  }# for each replicate
  
}# for each landscape category
names(landscape.list) <- paste("sa",1:num.landscape.categories,sep="")

# -------------------------------------------------------------------------
# is mean suitability maintained?
# yes it is

# -------------------------------------------------------------------------
# distance matrix calculation
land.df <- expand.grid(x = 1:ncol, y = 1:nrow)
distance.matrix <- as.matrix(dist(land.df))

land.df$cell <- 1:nrow(land.df)
distance.df <- reshape2::melt(distance.matrix)
names(distance.df) <- c("cell_from","cell_to","distance")

# -------------------------------------------------------------------------
write.csv2(data.frame(landscape.category = names(landscape.list),
                      spatial.autocorr.value = spatial.autocorr),
           "results/spatial_autocorrelation_categories.csv",
           row.names = FALSE)
save(landscape.list,file = "results/landscape_matrices.RData")

# cells and distances
write.csv2(land.df,"results/cell_coordinates.csv",row.names = FALSE)
write.csv2(distance.df,"results/cell_distances.csv",row.names = FALSE)
