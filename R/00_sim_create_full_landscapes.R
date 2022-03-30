
# populates landscape matrices with dispersal links

# INPUTS
# - nested lists of presence/absence dataframes for every combination of
# landscape/network categories: sp.presence[[landscape.category]][[network.category]][[replicate]]
# "results/presence_dataframes.RData"
# - landscape matrices with local communities "wired"
# "results/sim_landscape_matrices/NETWORKCATEGORY_LANDSCAPECATEGORY_REPLICATE_RICHNESS_CELLS.RData
# - dispersal distances dataframe:
# "results/dispersal_kernels.csv"

# OUTPUTS
# - andscape matrices with dispersal, each stored individually
# results/sim_landscape_matrices/NETWORK_LANDSCAPE_DISPERSAL_REPLICATE_RICHNESS_CELLS.RData
# - dataframe with combinations of categories
# "results/landscapes_info.csv"

# -------------------------------------------------------------------------

library(tidyverse)
library(igraph)
library(NLMR)
library(sf)

# -------------------------------------------------------------------------
landscapes.path <- "results/sim_landscape_matrices/"

load("results/presence_dataframes.RData")
disp.df <- read.csv2("results/dispersal_kernels.csv")

# -------------------------------------------------------------------------
# get some parameters

landscape.categories <- sort(unique(names(sp.presence)))
network.categories <- sort(unique(names(sp.presence[[1]])))
dispersal.categories <- sort(unique(disp.df$dispersal.category))
replicates <- length(sp.presence[[1]][[1]])

landscape.rows <- max(sp.presence[[1]][[1]][[1]]$landscape.row)
landscape.cols <- max(sp.presence[[1]][[1]][[1]]$landscape.col)
cells <- landscape.rows * landscape.cols
sp.names <- sort(unique(disp.df$sp))
richness <- length(sp.names)

# -------------------------------------------------------------------------
# I need the distances between each cell of the landscapes
# since I know the number of cells, I don't need to load any created landscape
# rather, it's easier to build one and obtain distances on the fly

temp.land <- expand.grid(x = 1:landscape.cols, y = 1:landscape.rows)
distance.matrix <- dist(temp.land, diag=T, upper=T)

# -------------------------------------------------------------------------
# iterate through each generated landscape to add dispersal links

# i.net <- i.land <- i.rep <- i.disp <- 1

for(i.land in 1:length(landscape.categories)){
  for(i.net in 1:length(network.categories)){
    for(i.rep in 1:replicates){
      
      # -------------------------------------------------------------------
      # load landscape matrix
      
      landscape.name <- paste(network.categories[i.net],"_",
                              landscape.categories[i.land],"_",
                              "re",i.rep,"_",
                              richness,"sp_",
                              cells,"cells.RData",sep="")
      
      load(paste(landscapes.path,landscape.name,sep=""))
      
      # -------------------------------------------------------------------
      # add dispersal
      for(i.disp in 1:length(dispersal.categories)){
        
        landscape <- my.landscape.matrix
        
        # # 3 - generate matrix of distances between patches
        # # igraph structure from adjacency matrix
        # patch.graph <- igraph::graph_from_adjacency_matrix(adjmatrix = patch.connectivity,mode = "undirected",diag = FALSE)
        # # calculate distances between patches
        # patch.distances <- igraph::distances(graph = patch.graph,mode = "out")
        # 
        # # 4 - foraging and dispersal niches
        # 
        # # remap the niches to the range 0-d
        # foraging.niche <- c(0,niche,1)
        # foraging.niche <- foraging.niche - min(foraging.niche)
        # foraging.niche <- foraging.niche/max(foraging.niche)
        # foraging.niche <- foraging.niche * max.foraging.rate
        # foraging.niche <- foraging.niche[2:(length(foraging.niche)-1)]
        # 
        # dispersal.niche <- c(0,niche,1)
        # dispersal.niche <- dispersal.niche - min(dispersal.niche)
        # dispersal.niche <- dispersal.niche/max(dispersal.niche)
        # dispersal.niche <- dispersal.niche * max.dispersal.rate
        # dispersal.niche <- dispersal.niche[2:(length(dispersal.niche)-1)]
        # 
        # # create metacommunity matrix
        # D = matrix(0,nr=N*S,nc = N*S)
        # 
        # # 1. maximum distances reachable by any species
        # 
        # if(max.dispersal.rate > 0){
        #   
        #   sp.max.distance <- (max.dispersal.distance+0.1)/max(dispersal.niche)*dispersal.niche
        #   # matrix containing the dispersal coefficient of each sp for each distance
        #   dispersal.coef.data <- matrix(0,nrow=S,ncol=max(patch.distances))
        #   
        #   if(dispersal.coefs == "constant"){
        #     # every species disperses at least the nearest patch
        #     dispersal.coef.data[,1] <- max.dispersal.rate
        #     
        #     # the rest of the coefficients are constant within the reachable distances
        #     for(i.row in 1:nrow(dispersal.coef.data)){
        #       if(sp.max.distance[i.row]>2){
        #         dispersal.coef.data[i.row,2:floor(sp.max.distance[i.row])] <- max.dispersal.rate
        #       }# if reaches other patches
        #     }# for each row
        #     
        #     # divide dispersal across reachable patches
        #     if(ncol(dispersal.coef.data)>1){
        #       dispersal.coef.data <- t(apply(X = dispersal.coef.data,MARGIN = 1,FUN = function(x) x/sum(x>0)))
        #     }
        #     
        #   }else if(dispersal.coefs == "scaling"){
        #     
        #     # this formula is a line equation expanded to 1) find m and 2) find b, and these plugged into y=mx*b to get the values of dispersal(y) for a certain distance(x)
        #     for(i.sp in 1:nrow(dispersal.coef.data)){
        #       my.sp.distances <- ((max.dispersal.rate/(1-sp.max.distance[i.sp]))*1:floor(sp.max.distance[i.sp])) + (max.dispersal.rate*(-sp.max.distance[i.sp]))/(1-sp.max.distance[i.sp])
        #       if(length(my.sp.distances)<ncol(dispersal.coef.data)){
        #         my.sp.distances[(length(my.sp.distances)+1):ncol(dispersal.coef.data)] <- 0
        #       }
        #       dispersal.coef.data[i.sp,1:ncol(dispersal.coef.data)] <- my.sp.distances[1:ncol(dispersal.coef.data)] 
        #     }# for i.sp
        #     dispersal.coef.data[dispersal.coef.data<0] <- 0
        #     
        #     # divide dispersal effort across reachable patches
        #     for(i.row in 1:nrow(dispersal.coef.data)){
        #       if(sum(dispersal.coef.data[i.row,]>0)>1){
        #         dispersal.coef.data[i.row,] <- dispersal.coef.data[i.row,] * (max.dispersal.rate/sum(dispersal.coef.data[i.row,]))
        #       }# if more than one patch
        #     }# for i.row
        #     
        #   }# if-else dispersal is constant or scaling
        #   
        #   # 5 - assign dispersal coefficients
        #   # maintaining mass balance, i.e. the sum of the dispersal to all connected patches i.e. sum(d/N-1)
        #   # is equivalent to the loss of the source patch, i.e. -d
        #   # -d is the value of the diagonal patch, then.
        #   
        #   # update matrix D
        #   for(i.row in 1:nrow(D)){
        #     for(i.col in 1:ncol(D)){
        #       
        #       # if diagonal element and not main diagonal,
        #       # it is a dispersal coefficient
        #       if(i.row %% S == i.col %% S & i.row != i.col){
        #         
        #         # which sp?
        #         my.sp <- i.row %% S
        #         if(my.sp == 0){my.sp <- S}
        #         
        #         # source and dest patch
        #         source.patch <- ceiling(i.col/S) 
        #         dest.patch <- ceiling(i.row/S)
        #         
        #         # double check
        #         if(source.patch != dest.patch){
        #           
        #           # distance among them
        #           my.distance <- patch.distances[source.patch,dest.patch]
        #           
        #           # the value is the one from dispersal.coef.data (the net effort for a given distance)
        #           # divided by the number of patches at that distance
        #           
        #           # how many patches at this distance?
        #           patches.at.dist <- table(patch.distances[source.patch,])
        #           patches.at.dist <- as.integer(patches.at.dist[which(names(patches.at.dist) == my.distance)])
        #           
        #           # divide the effort from foraging.coef.data among all patches within a given distance
        #           dispersal.coef.weight <- dispersal.coef.data[my.sp,my.distance]/patches.at.dist
        #           
        #           D[i.row,i.col] <- dispersal.coef.weight
        #           
        #         }# if different patch
        #       }# if dispersal cell
        #       
        #     }# for i.col
        #   }# for i.row

        # -----------------------------------------------------------------
        # store landscape with dispersal
        
        # this is the name of the resulting landscape matrix
        landscape.dispersal.name <- paste(network.categories[i.net],"_",
                                          landscape.categories[i.land],"_",
                                          dispersal.categories[i.disp],"_",
                                          "re",i.rep,"_",
                                          richness,"sp_",
                                          cells,"cells.RData",sep="")
        
        save(landscape, paste(landscapes.path,landscape.dispersal.name,sep=""))
        
      }# for i.disp
    }# for i.rep
  }# for i.net
}# for i.land





