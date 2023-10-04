
# populates landscape matrices with dispersal links
# in this version, dispersal is correlated with species degree, 
# so dispersal rates are generated "on the fly", thus 
# the differences with the "standard" script
# the added code is basically copied from 01_sim_generate_dispersal_kernel.R

# INPUTS
# - nested lists of presence/absence dataframes for every combination of
# landscape/network categories: sp.presence[[landscape.category]][[network.category]][[replicate]]
# "results/presence_dataframes.RData"
# - landscape matrices with local communities "wired"
# "results/sim_landscape_matrices/NETWORKCATEGORY_LANDSCAPECATEGORY_REPLICATE_RICHNESS_CELLS.RData
# - dispersal distances dataframe:
# "results/dispersal_kernels_correlated.csv"
# - cell coordinates and distances: "results/cell_coordinates.csv","results/cell_distances.csv"

# OUTPUTS
# - landscape matrices with dispersal, each stored individually
# results/sim_landscape_matrices/NETWORK_LANDSCAPE_DISPERSAL_REPLICATE_RICHNESS_CELLS.RData
# - dataframe with combinations of categories
# "results/landscapes_info.csv"

# -------------------------------------------------------------------------

library(foreach)
library(doParallel)
library(tidyverse)
library(igraph)
library(NLMR)
library(sf)

# set number of cores -----------------------------------------------------

workers <- 6
cl <- makeCluster(workers)
# register the cluster for using foreach
registerDoParallel(cl)

# -------------------------------------------------------------------------
# if correlated, save to a different folder
landscapes.path <- "results/sim_landscape_matrices_correlated_deg_disp/"

load("results/presence_dataframes.RData")
cell.distances <- read.csv2("results/cell_distances.csv")
cell.coords <- read.csv2("results/cell_coordinates.csv")

# for dispersal
param <- read.csv2("results/sim_landscape_matrices/parameters_v3.csv")
richness <-  param$richness
num.dispersal.categories <- param$num.dispersal.categories
num.category.replicates <- param$num.category.replicates
deg.dist <- read.csv2("results/sim_degree_dist_networks.csv")

# exponetial rate
min.rate <- param$min.rate
max.rate <- param$max.rate
dispersal.rate.gradient <- seq(from = min.rate,
                               to = max.rate, 
                               length.out = num.dispersal.categories)

dispersal.categories <- paste("dk",sprintf("%02d", 1:num.dispersal.categories),sep="")
# -------------------------------------------------------------------------
# get some parameters

landscape.categories <- sort(unique(names(sp.presence)))
network.categories <- sort(unique(names(sp.presence[[1]])))
# dispersal.categories <- sort(unique(disp.df$dispersal.category))
replicates <- length(sp.presence[[1]][[1]])

landscape.rows <- max(sp.presence[[1]][[1]][[1]]$landscape.row)
landscape.cols <- max(sp.presence[[1]][[1]][[1]]$landscape.col)
cells <- landscape.rows * landscape.cols
sp.names <- sort(unique(deg.dist$node_from))
richness <- length(sp.names)

# ID to loop over ---------------------------------------------------------

id <- expand.grid(network.categories,landscape.categories,1:replicates)
id.char <- sort(paste(id[,1],"_",id[,2],"_",id[,3],sep=""))

# -------------------------------------------------------------------------
# iterate through each generated landscape to add dispersal links

results <- foreach(i.id = 1:length(id.char), 
                   # .combine=comb.fun, 
                   .packages = 'tidyverse') %dopar% {
                     
                     # recover landscape,network, and replicate from the ID
                     # i.land <- sub(".*_", "", id.char[i.id])
                     i.land <- as.numeric(substr(id.char[i.id],3,4))
                     i.net <- as.numeric(substr(id.char[i.id],8,9))
                     i.rep <- as.numeric(substr(id.char[i.id],11,nchar(id.char[i.id])))
                     
                     # -------------------------------------------------------------------
                     # load landscape matrix
                     
                     landscape.name <- paste(network.categories[i.net],"_",
                                             landscape.categories[i.land],"_",
                                             "re",i.rep,"_",
                                             richness,"sp_",
                                             cells,"cells.RData",sep="")
                     
                     # my.landscape.matrix
                       load(paste("results/sim_landscape_matrices/",landscape.name,sep=""))
                     
                     # -------------------------------------------------------------------
                     # add dispersal
                     for(i.disp in 1:length(dispersal.categories)){
                       
                       landscape <- my.landscape.matrix
                       
                       # -----------------------------------------------------------------
                       # obtain dispersal links
                       
                       # in this test, dispersal is correlated with species degree,
                       # so it is dynamically generated here
                       # the code is essentially the same from 01_sim_generate_dispersal_kernel
                       
                       my.dist <- rexp(n = richness,rate = dispersal.rate.gradient[i.disp])
                       my.deg <- subset(deg.dist, generative.model == network.categories[i.net] & replicate == i.rep) %>% 
                         group_by(node_from) %>%
                         summarise(degree = sum(value != 0)) %>%
                         arrange(desc(degree))
                       my.disp <- data.frame(sp = my.deg$node_from, dispersal.distance = sort(my.dist,decreasing = T))
                       
                       # expand the presence dataframe to see the cells to which
                       # each sp can potentially disperse
                       # by crossing presence information with cell distances
                       # and with dispersal distances
                       
                       my.presence <- sp.presence[[i.land]][[i.net]][[i.rep]]
                       my.presence <- subset(my.presence,presence == TRUE) %>%
                         dplyr::select(sp,cell) %>%
                         rename(cell_from = cell)
                       my.presence.full <- expand_grid(my.presence,cell_to = 1:cells)
                       my.presence.full.2 <- left_join(my.presence.full,cell.distances)
                       my.presence.full.3 <- left_join(my.presence.full.2,my.disp[,c("sp","dispersal.distance")])
                       my.presence.full.3 <- subset(my.presence.full.3,cell_from != cell_to)
                       
                       my.presence.full.3$dispersal.potential <- ifelse(my.presence.full.3$distance <= 
                                                                          my.presence.full.3$dispersal.distance,
                                                                        TRUE,FALSE)
                       my.dispersal.potential <- subset(my.presence.full.3[,c("sp","cell_from",
                                                                              "cell_to","dispersal.potential")],
                                                        dispersal.potential == TRUE)
                       
                       # now, check which of the target cells have populations of the species
                       # because only those will have realized dispersal
                       my.presence$presence <- TRUE
                       my.dispersal.realized <- left_join(my.dispersal.potential,my.presence,
                                                          by = c("cell_to" = "cell_from",
                                                                 "sp" = "sp")) %>%
                         replace_na(list(presence = FALSE))
                       my.dispersal.realized$dispersal <- as.logical(my.dispersal.realized$dispersal.potential * 
                                                                       my.dispersal.realized$presence)
                       my.dispersal.r2 <- subset(my.dispersal.realized[,c("sp","cell_from","cell_to","dispersal")],
                                                 dispersal == TRUE)
                       
                       # remove symmetrical information
                       my.dispersal.r2$duplicated <- FALSE
                       for(i in 1:nrow(my.dispersal.r2)){
                         if(!my.dispersal.r2$duplicated[i]){
                           dup <- which(my.dispersal.r2$cell_from == my.dispersal.r2$cell_to[i] &
                                          my.dispersal.r2$cell_to == my.dispersal.r2$cell_from[i] &
                                          my.dispersal.r2$sp == my.dispersal.r2$sp[i])
                           my.dispersal.r2$duplicated[dup] <- TRUE
                         }# if
                       }# for i
                       
                       # this is the distilled dataframe
                       # including the dispersal coefficient to assign
                       # i.e. 1/number of dispersing cells (irrespective of distance)
                       
                       # this may change in the future
                       
                       my.dispersal <- subset(my.dispersal.r2,duplicated == FALSE) %>%
                         group_by(sp,cell_from) %>%
                         mutate(dispersal.coef = 1/n()) %>%
                         dplyr::select(sp,cell_from,cell_to,dispersal.coef)
                       
                       # ----------------------------------------------------------------
                       # update landscape matrix
                       
                       # I only need to go through the upper triangle, since 
                       # dispersal is symmetric
                       # hence the weird nested for loops
                       for(i.row in 1:(nrow(landscape)-1)){
                         
                         # which sp?
                         my.sp <- i.row %% richness
                         if(my.sp == 0){my.sp <- richness}
                         
                         for(i.col in (i.row+1):ncol(landscape)){
                           
                           # if diagonal element and not main diagonal,
                           # it is a dispersal coefficient
                           # if(i.row %% richness == i.col %% richness & i.row != i.col){
                           if(i.row %% richness == i.col %% richness){
                             
                             # source and dest cell
                             source.cell <- ceiling(i.col/richness)
                             dest.cell <- ceiling(i.row/richness)
                             
                             # double check
                             if(source.cell != dest.cell){
                               
                               valid.dispersal <- which(my.dispersal$sp == sp.names[my.sp] &
                                                          (my.dispersal$cell_from == source.cell &
                                                             my.dispersal$cell_to == dest.cell | 
                                                             my.dispersal$cell_from == dest.cell &
                                                             my.dispersal$cell_to == source.cell) )
                               if(length(valid.dispersal) == 1){
                                 
                                 disp.coef <- my.dispersal$dispersal.coef[valid.dispersal]
                                 
                                 # fill the symmetric positions
                                 landscape[i.row,i.col] <- disp.coef
                                 landscape[i.col,i.row] <- disp.coef
                               }# if realized dispersal
                               
                             }# if different cell
                           }# if dispersal cell
                           
                         }# for i.col
                       }# for i.row
                       
                       
                       # -----------------------------------------------------------------
                       # remove absent species from every cell - this can be done simply
                       # by keeping only columns and rows with 1 entries in their diagonal,
                       # since each present species has 1 in the corresponding diagonal.
                       
                       present.sp <- diag(landscape) == 1
                       landscape <- landscape[present.sp,present.sp] 
                       
                       # -----------------------------------------------------------------
                       # store landscape with dispersal
                       
                       # this is the name of the resulting landscape matrix
                       landscape.dispersal.name <- paste(network.categories[i.net],"_",
                                                         landscape.categories[i.land],"_",
                                                         dispersal.categories[i.disp],"_",
                                                         "re",i.rep,"_",
                                                         richness,"sp_",
                                                         cells,"cells.RData",sep="")
                       
                       save(landscape, file = paste(landscapes.path,landscape.dispersal.name,sep=""))
                       
                     }# for i.disp
                   }# foreach id.char

stopCluster(cl)
