
# test how propagation metrics respond to removal of nodes, either random or targeted

# INPUTS

# OUTPUTS

# -------------------------------------------------------------------------

comm.path <- "results/extinction_sequences/communicability"
gce.path <- "results/extinction_sequences/"

# -------------------------------------------------------------------------
library(NLMR)
library(tidyverse)
library(igraph)
library(intsegration)
source("R/generate_landscape.R")
source("R/remove_sp_from_landscape.R")
source("R/communicability.R")
source("R/auxiliary_functions/GCE_weighted.R")
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# -------------------------------------------------------------------------
param <- read.csv2("results/sim_landscape_matrices/parameters_v2.csv")
suit.df <- read.csv2("results/sp_suitability.csv")

# -------------------------------------------------------------------------
# select one category of each factor for building the base network

landscape.cat <- "sa03"
network.cat <- "dd03"
dispersal.cat <- "dk03"

# -------------------------------------------------------------------------
# random, targeted_dispersal, targeted_degree, targeted_distribution

extinction.tests <- c("random","targeted_distribution","targeted_degree",
                      "targeted_dispersal")

# -------------------------------------------------------------------------

landscape.list <- generate_landscape(param = param,
                                     suit.df = suit.df,
                                     landscape.cat = landscape.cat,
                                     network.cat = network.cat,
                                     dispersal.cat = dispersal.cat)
landscape <- landscape.list[[1]]
sp.traits <- landscape.list[[2]]
cell.distances <- landscape.list[[3]]
cell.distances$cell_from <- as.character(cell.distances$cell_from)
cell.distances$cell_to <- as.character(cell.distances$cell_to)
num.sp <- nrow(sp.traits)

# -------------------------------------------------------------------------
# extinction sequences: for each sequence, calculate metrics 
# (gce and pairwise communicability/shortest path lengths)
# as well as some basic network properties

net.list <- list()
gce.list <- list()

for(i.test in 1:length(extinction.tests)){
  
  remaining.sp <- num.sp
  remaining.sp.traits <- sp.traits
  remaining.landscape <- landscape
  
  while(remaining.sp > 2){
    
    # index of the sp to remove
    if(extinction.tests[i.test] == "random"){
      to.remove <- sample(remaining.sp,1)
    }else if(extinction.tests[i.test] == "targeted_dispersal"){
      to.remove <- which(remaining.sp.traits$dispersal.distance ==
                           max(remaining.sp.traits$dispersal.distance))
    }else if(extinction.tests[i.test] == "targeted_degree"){
      to.remove <- which(remaining.sp.traits$degree ==
                           max(remaining.sp.traits$degree))
    }else if(extinction.tests[i.test] == "targeted_distribution"){
      to.remove <- which(remaining.sp.traits$presences ==
                           max(remaining.sp.traits$presences))
    }
    
    # if more than one sp
    if(length(to.remove)>1){
      to.remove <- to.remove[sample(length(to.remove),1)]
    }
    
    # name of the species to remove
    sp.to.remove <- remaining.sp.traits$sp[to.remove]
    
    # update landscape
    remaining.landscape <- remove_sp_from_landscape(remaining.landscape,
                                                    sp.name = sp.to.remove)
    # update traits and number
    remaining.sp <- remaining.sp - 1
    remaining.sp.traits <- subset(remaining.sp.traits, sp != sp.to.remove)
    
    # check that the resulting meta-adjacency matrix is valid, i.e. 
    # there are still interspecific links in the local communities
    # otherwise it does not make sense to calculate metrics
    num.interspecific.links <- sum(remaining.landscape != 0 & remaining.landscape != 1)
    
    # -------------------------------------------------------------------------
    # store some network metrics
    remaining.cells <- unique(substr(rownames(remaining.landscape),
                                     6,nchar(rownames(remaining.landscape))))
    
    net.list[[length(net.list)+1]] <- data.frame(extinction_sequence = extinction.tests[i.test],
                                                 species_removed = num.sp - remaining.sp,
                                                 num.intraspecific.links = length(diag(remaining.landscape)),
                                                 num.interspecific.links = num.interspecific.links,
                                                 num.dispersal.links = sum(remaining.landscape == 1) - length(diag(remaining.landscape)),
                                                 num.zeros = sum(remaining.landscape == 0),
                                                 num.cells = length(remaining.cells)) 
    
    # -------------------------------------------------------------------------
    # calculate metrics
    
    # if(num.interspecific.links > 0){
    #   
    #   communicability.matrices <- communicability(remaining.landscape) 
    #   
    #   # tidy functions only work with dataframes, but this still works, and is fast
    #   comm.df <- reshape2::melt(communicability.matrices[[1]],
    #                             value.name = "binary.communicability")
    #   
    #   # extract cell of sp1 and cell of sp2
    #   comm.df$sp1 <- sub("\\-.*", "", comm.df$Var1)
    #   comm.df$cell1 <- sub(".*-", "", comm.df$Var1)
    #   comm.df$sp2 <- sub("\\-.*", "", comm.df$Var2)
    #   comm.df$cell2 <- sub(".*-", "", comm.df$Var2)
    #   
    #   comm.df$scaled.binary.communicability <- range01(comm.df$binary.communicability)
    #   
    #   # this should be valid because the two matrices have the same dimensions and names
    #   dfw <- reshape2::melt(communicability.matrices[[2]],
    #                         value.name = "weighted.communicability")
    #   comm.df$weighted.communicability <- dfw$weighted.communicability
    #   
    #   # -------------------------------------------------------------------------
    #   # get path lengths as well
    #   graph.D <- igraph::graph_from_adjacency_matrix(adjmatrix = remaining.landscape,
    #                                                  mode = "undirected",
    #                                                  weighted = "1",diag = FALSE)
    #   # this takes ~10min
    #   path.lengths <- igraph::distances(graph = graph.D,algorithm = "unweighted")
    #   
    #   # turning it to df is quick
    #   df.path.lengths <- reshape2::melt(path.lengths,value.name = "shortest.path.length")
    #   
    #   df.path.lengths$sp1 <- sub("\\-.*", "", df.path.lengths$Var1)
    #   df.path.lengths$cell1 <- sub(".*-", "", df.path.lengths$Var1)
    #   df.path.lengths$sp2 <- sub("\\-.*", "", df.path.lengths$Var2)
    #   df.path.lengths$cell2 <- sub(".*-", "", df.path.lengths$Var2) 
    #   
    #   comm.df$Var1 <- NULL
    #   comm.df$Var2 <- NULL
    #   
    #   df.path.lengths$Var1 <- NULL
    #   df.path.lengths$Var2 <- NULL
    #   
    #   # sp1,cell1,sp2,cell2 should be common
    #   comm.df <- left_join(comm.df,df.path.lengths)
    #   
    #   # -------------------------------------------------------------------------
    #   # get spatial distance between cells
    #   comm.df <- left_join(comm.df,cell.distances,by = c("cell1" = "cell_from",
    #                                                      "cell2" = "cell_to"))
    #   
    #   comm.df <- comm.df[,c("sp1","cell1",
    #                         "sp2","cell2","binary.communicability",
    #                         "scaled.binary.communicability",
    #                         "weighted.communicability",
    #                         "shortest.path.length",
    #                         "distance")]
    #   comm.df$extinction_sequence <- extinction.tests[i.test]
    #   comm.df$species_removed <- num.sp - remaining.sp
    #   
    #   # -------------------------------------------------------------------------
    #   # global communication efficiency - a network level property
    #   
    #   landscape.graph <- igraph::graph_from_adjacency_matrix(remaining.landscape,
    #                                                          weighted = T)
    #   landscape.gce <- GCE_weighted(g = landscape.graph,normalised = T,
    #                                 directed = T)
    #   
    #   # -------------------------------------------------------------------------
    #   
    #   gce.list[[length(gce.list)+1]] <- data.frame(extinction_sequence = extinction.tests[i.test],
    #                                                species_removed = num.sp - remaining.sp,
    #                                                normalised.gce = landscape.gce$normalised) 
    #   
    #   # -------------------------------------------------------------------------
    #   
    #   save(comm.df, file = paste(comm.path,"/comm_",
    #                              extinction.tests[i.test],"_",
    #                              num.sp-remaining.sp,"removed.RData",sep=""))
    # }else{
    #   # artificially set the flag so that the while loop ends
    #   remaining.sp <- 1
    # }
  }# while remaining.sp
  
  cat(extinction.tests[i.test],"- completed\n")
  
}# for i.test

# gce.df <- bind_rows(gce.list)

# write.csv2(gce.df,file = paste(gce.path,"/gce_extinction_sequences.csv",sep=""))

net.df <- bind_rows(net.list)
write.csv2(net.df,file = paste(gce.path,"/network_properties_extinction_sequences.csv",sep=""))
