
# calculate the communicability of simulated landscapes

# INPUTS
# - landscape matrices for each combination
# "results/sim_landscape_matrices/LANDSCAPE_NETWORK_DISPERSAL_REPLICATE_SP_CELLS.RData"
# - distances between cells

# OUTPUTS
# - one dataframe for each network, with pairwise communicability values
# "results/communicability/comm_LANDSCAPE_NETWORK_DISPERSAL_REPLICATE_SP_CELLS.csv"

# -------------------------------------------------------------------------

library(tidyverse)
library(igraph) # for path lenghts
source("R/communicability.R")
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# -------------------------------------------------------------------------
# read general data
param <- read.csv2("results/sim_landscape_matrices/parameters_v2.csv")

network.categories <- read.csv2("results/network_gradient_categories.csv")
landscape.categories <- read.csv2("results/spatial_autocorrelation_categories.csv")
dispersal.categories <-  read.csv2("results/dispersal_kernels.csv")
cell.distances <- read.csv2("results/cell_distances.csv")
cell.distances$cell_from <- as.character(cell.distances$cell_from)
cell.distances$cell_to <- as.character(cell.distances$cell_to)

# -------------------------------------------------------------------------
# recover factors

network.categories <- network.categories$network.category
landscape.categories <- landscape.categories$landscape.category
dispersal.categories <- unique(dispersal.categories$dispersal.category)
replicates <- param$num.category.replicates

richness <- param$richness
cells <- param$ncol * param$nrow

# -------------------------------------------------------------------------
# go through each replicate

i.land <- 1
i.net <- 1
i.disp <- 1
i.rep <- 1

for(i.land in 1:length(landscape.categories)){
  for(i.net in 1:length(network.categories)){
    for(i.disp in 1:length(dispersal.categories)){
      for(i.rep in 1:replicates){
        
        my.landscape.name <- paste("results/sim_landscape_matrices/",
                                   network.categories[i.net],"_",
                                   landscape.categories[i.land],"_",
                                   dispersal.categories[i.disp],"_",
                                   "re",i.rep,"_",
                                   richness,"sp_",cells,"cells.RData",sep="")
        my.df.name <- paste("results/communicability/comm_",
                            network.categories[i.net],"_",
                            landscape.categories[i.land],"_",
                            dispersal.categories[i.disp],"_",
                            "re",i.rep,"_",
                            richness,"sp_",cells,"cells.RData",sep="")
        load(my.landscape.name)
        
        communicability.matrices <- communicability(landscape) 
        
        # tidy functions only work with dataframes, but this still works, and is fast
        comm.df <- reshape2::melt(communicability.matrices[[1]],
                                  value.name = "binary.communicability")
        
        # extract cell of sp1 and cell of sp2
        comm.df$sp1 <- sub("\\-.*", "", comm.df$Var1)
        comm.df$cell1 <- sub(".*-", "", comm.df$Var1)
        comm.df$sp2 <- sub("\\-.*", "", comm.df$Var2)
        comm.df$cell2 <- sub(".*-", "", comm.df$Var2)
        
        comm.df$scaled.binary.communicability <- range01(comm.df$binary.communicability)
        
        # this should be valid because the two matrices have the same dimensions and names
        dfw <- reshape2::melt(communicability.matrices[[2]],
                              value.name = "weighted.communicability")
        comm.df$weighted.communicability <- dfw$weighted.communicability
        
        # -------------------------------------------------------------------------
        # get path lengths as well
        graph.D <- igraph::graph_from_adjacency_matrix(adjmatrix = communicability.matrices[[1]],
                                                       mode = "undirected",
                                                       weighted = "1",diag = FALSE)
        # this takes ~10min
        path.lengths <- igraph::distances(graph = graph.D,algorithm = "unweighted")
        
        # turning it to df is quick
        df.path.lengths <- reshape2::melt(path.lengths,value.name = "shortest.path.length")
        
        df.path.lengths$sp1 <- sub("\\-.*", "", df.path.lengths$Var1)
        df.path.lengths$cell1 <- sub(".*-", "", df.path.lengths$Var1)
        df.path.lengths$sp2 <- sub("\\-.*", "", df.path.lengths$Var2)
        df.path.lengths$cell2 <- sub(".*-", "", df.path.lengths$Var2) 
        
        comm.df$Var1 <- NULL
        comm.df$Var2 <- NULL
        
        df.path.lengths$Var1 <- NULL
        df.path.lengths$Var2 <- NULL
        
        # sp1,cell1,sp2,cell2 should be common
        comm.df <- left_join(comm.df,df.path.lengths)
        
        # -------------------------------------------------------------------------
        # get spatial distance between cells
        comm.df <- left_join(comm.df,cell.distances,by = c("cell1" = "cell_from",
                                                         "cell2" = "cell_to"))
        
        comm.df <- comm.df[,c("sp1","cell1",
                              "sp2","cell2","binary.communicability",
                              "scaled.binary.communicability",
                              "weighted.communicability",
                              "shortest.path.length",
                              "distance")]
        
        save(comm.df, file = my.df.name)
        
      }# for i.rep
    }# for i.disp
  }# for i.net
}# for i.land

# -------------------------------------------------------------------------

