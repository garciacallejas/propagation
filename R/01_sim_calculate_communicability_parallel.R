
# calculate the communicability of simulated landscapes
# parallel version

# INPUTS
# - landscape matrices for each combination
# "results/sim_landscape_matrices/LANDSCAPE_NETWORK_DISPERSAL_REPLICATE_SP_CELLS.RData"
# - distances between cells

# OUTPUTS
# - one dataframe for each network, with pairwise communicability values
# "results/communicability/comm_LANDSCAPE_NETWORK_DISPERSAL_REPLICATE_SP_CELLS.csv"

# -------------------------------------------------------------------------
read.sim.landscape.path <- "/home/david/Work/datasets/NZ/results/sim_landscape_matrices/"
save.comm.path <- "/home/david/Work/datasets/NZ/results/communicability/"
save.net.comm.path <- "/home/david/Work/datasets/NZ/results/communicability/network_level/"

# -------------------------------------------------------------------------

library(foreach)
library(doParallel)

library(tidyverse)
library(igraph) # for path lenghts
source("R/auxiliary_functions/comm.R")
# range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# set number of cores -----------------------------------------------------

workers <- 3
cl <- makeCluster(workers)
# register the cluster for using foreach
registerDoParallel(cl)

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

# ID to loop over ---------------------------------------------------------

id <- expand.grid(network.categories,landscape.categories,
                  dispersal.categories,1:replicates)
id.char <- sort(paste(id[,1],"_",id[,2],"_",id[,3],"_",id[,4],sep=""))

# test
# i.id <- 256

# -------------------------------------------------------------------------
# iterate through each generated landscape 

results <- foreach(i.id = 1:length(id.char), 
                   # .combine=comb.fun, 
                   .packages = c('tidyverse',
                                 "igraph",
                                 # "intsegration",
                                 "reshape2")) %dopar% 
  {
    
    source("R/communicability.R")
    source("R/auxiliary_functions/comm.R")
    
    # recover landscape,network, and replicate from the ID
    # i.land <- sub(".*_", "", id.char[i.id])
    i.land <- as.numeric(substr(id.char[i.id],3,4))
    i.net <- as.numeric(substr(id.char[i.id],8,9))
    i.disp <- as.numeric(substr(id.char[i.id],13,14))
    i.rep <- as.numeric(substr(id.char[i.id],16,nchar(id.char[i.id])))
    
    my.landscape.name <- paste(read.sim.landscape.path,
                               network.categories[i.net],"_",
                               landscape.categories[i.land],"_",
                               dispersal.categories[i.disp],"_",
                               "re",i.rep,"_",
                               richness,"sp_",cells,"cells.RData",sep="")
    my.df.name <- paste(save.comm.path,"/comm_",
                        network.categories[i.net],"_",
                        landscape.categories[i.land],"_",
                        dispersal.categories[i.disp],"_",
                        "re",i.rep,"_",
                        richness,"sp_",cells,"cells.RData",sep="")
    load(my.landscape.name)
    
    communicability.metrics <- comm(landscape,
                                    normalised = FALSE,
                                    return.pairwise.comm = T) 
    
    # tidy functions only work with dataframes, but this still works, and is fast
    comm.df <- reshape2::melt(communicability.metrics[[2]],
                              value.name = "binary.communicability")
    
    # extract cell of sp1 and cell of sp2
    comm.df$sp1 <- sub("\\-.*", "", comm.df$Var1)
    comm.df$cell1 <- sub(".*-", "", comm.df$Var1)
    comm.df$sp2 <- sub("\\-.*", "", comm.df$Var2)
    comm.df$cell2 <- sub(".*-", "", comm.df$Var2)
    
    comm.df$scaled.binary.communicability <- scales::rescale(comm.df$binary.communicability,to = c(0,1))
    
    # this should be valid because the two matrices have the same dimensions and names
    dfw <- reshape2::melt(communicability.metrics[[3]],
                          value.name = "weighted.communicability")
    comm.df$weighted.communicability <- dfw$weighted.communicability
    
    # -------------------------------------------------------------------------
    # get path lengths as well
    graph.D <- igraph::graph_from_adjacency_matrix(adjmatrix = landscape,
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
    
    
    # -------------------------------------------------------------------------
    # global communication efficiency - a network level property
    
    # landscape.graph <- igraph::graph_from_adjacency_matrix(landscape,
    #                                                        weighted = T)
    # landscape.gce <- GCE_weighted(g = landscape.graph,normalised = T,
    #                               directed = T)

# -------------------------------------------------------------------------
    net.comm <- data.frame(landscape.category = landscape.categories[i.land],
                          network.category = network.categories[i.net],
                          dispersal.category = dispersal.categories[i.disp],
                          replicate = i.rep,
                          normalised.communicability = communicability.metrics[[1]])
    
# -------------------------------------------------------------------------

    write.csv2(net.comm,paste(save.net.comm.path,
                             network.categories[i.net],"_",
                             landscape.categories[i.land],"_",
                             dispersal.categories[i.disp],"_",
                             "re",i.rep,"_",
                             richness,"sp_",cells,"cells.csv",sep=""),row.names = F)
    save(comm.df, file = my.df.name)
    
  }# foreach id.char

stopCluster(cl)

# -------------------------------------------------------------------------

netcom.df <- list.files("/home/david/Work/datasets/NZ/results/communicability/network_level/",full.names = T) %>% map_dfr(read.csv2)
netcom.df$X <- NULL
write.csv2(netcom.df,"results/sim_network_level_communicability.csv",row.names = F)
