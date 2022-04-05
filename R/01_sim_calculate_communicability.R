
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

# -------------------------------------------------------------------------
# read general data

# -------------------------------------------------------------------------
# recover factors

landscape.categories
network.categories
dispersal.categories
replicates

richness
cells

# -------------------------------------------------------------------------
# go through each replicate

for(i.land in 1:length(landscape.categories)){
  for(i.net in 1:length(network.categories)){
    for(i.disp in 1:length(dispersal.categories)){
      for(i.rep in 1:replicates){
        
        my.landscape.name <- paste("results/sim_landscape_matrices/",
                                   landscape.categories[i.land],"_",
                                   network.categories[i.net],"_",
                                   dispersal.categories[i.disp],"_",
                                   richness,"sp_",cells,"cells.RData",sep="")
        my.df.name <- paste("results/communicability/comm_",
                            landscape.categories[i.land],"_",
                            network.categories[i.net],"_",
                            dispersal.categories[i.disp],"_",
                            richness,"sp_",cells,"cells.RData",sep="")
        load(my.landscape.name)
        
        # TODO UPDATE WITH PROPER NAME
        communicability.matrices <- communicability(my.landscape) 
        
        # tidy functions only work with dataframes, but this still works, and is fast
        comm.df <- reshape2::melt(communicability.matrices[[1]],
                                  value.name = "binary.communicability")
        
        # TODO adapt to extract cell of sp1 and cell of sp2
        # df1$sp1 <- sub("\\-.*", "", df1$Var1)
        # df1$grid.id.sp1 <- sub(".*-", "", df1$Var1)
        # df1$sp2 <- sub("\\-.*", "", df1$Var2)
        # df1$grid.id.sp2 <- sub(".*-", "", df1$Var2)  
        
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
        
        # TODO adapt to extract cell of sp1 and cell of sp2
        # df.path.lengths$sp1 <- sub("\\-.*", "", df.path.lengths$Var1)
        # df.path.lengths$grid.id.sp1 <- sub(".*-", "", df.path.lengths$Var1)
        # df.path.lengths$sp2 <- sub("\\-.*", "", df.path.lengths$Var2)
        # df.path.lengths$grid.id.sp2 <- sub(".*-", "", df.path.lengths$Var2)  
        
        # sp1,cell1,sp2,cell2 should be common
        comm.df <- left_join(comm.df,df.path.lengths)
        
        # -------------------------------------------------------------------------
        # get spatial distance between cells
        comm.df <- left_join(comm.df,cell.distances,by(c("cell1" = "cell_from",
                                                         "cell2" = "cell_to")))
        
        comm.df <- comm.df[,c("sp1","cell1",
                              "sp2","cell2","binary.communicability",
                              "scaled.binary.communicability",
                              "weighted.communicability",
                              "shortest.path.length",
                              "spatial.distance")]
        
        save(comm.df, file = my.df.name)
        
      }# for i.rep
    }# for i.disp
  }# for i.net
}# for i.land

# -------------------------------------------------------------------------

