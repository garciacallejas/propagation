
# test how propagation metrics respond to removal of nodes, either random or targeted

# INPUTS

# OUTPUTS

# -------------------------------------------------------------------------
library(NLMR)
library(tidyverse)

# -------------------------------------------------------------------------
param <- read.csv2("results/sim_landscape_matrices/parameters_v2.csv")
suit.df <- read.csv2("results/sp_suitability.csv")

# -------------------------------------------------------------------------
# select one category of each factor for building the base network

landscape.cat <- "sa03"
network.cat <- "dd03"
dispersal.cat <- "dk03"

# -------------------------------------------------------------------------
sp.names <- unique(suit.df$sp)

richness <-  param$richness

# to fully replicate the workflow of building networks, generate the set of
# categories and take the exact value of the one chosen above.
num.network.categories <- param$num.network.categories
num.dispersal.categories <- param$num.dispersal.categories
num.landscape.categories <- param$num.landscape.categories

# landscape dimensions
ncol <- param$ncol
nrow <- param$nrow
cells <- nrow*ncol

# poisson mean
min.lambda <- param$min.lambda
max.lambda <- param$max.lambda # this should vary with richness. for S = 50, 15 gives connectance = 0.3

# exponetial rate for dispersal
min.rate <- param$min.rate
max.rate <- param$max.rate

# some constants for sampling interaction strengths
int.mean <- param$int.mean
int.sd <- param$int.sd
tau <- param$tau
diag.dom <- param$diag.dom
# -------------------------------------------------------------------------


# network degree distribution
degree.dist.gradient <- seq(from = min.lambda,
                            to = max.lambda, 
                            length.out = num.network.categories)
# generative.models <- paste("dd",1:num.network.categories,sep="")
generative.models <- paste("dd",sprintf("%02d", 1:num.network.categories),sep="")

landscape.lambda <- degree.dist.gradient[which(generative.models == network.cat)]

# landscape configuration
landscape.gradient <- round(seq(0.1,0.9, length.out = num.landscape.categories),2)
landscape.models <- paste("sa",sprintf("%02d", 1:num.landscape.categories),sep="")

landscape.sa <- landscape.gradient[which(landscape.models == landscape.cat)]

# dispersal kernel
dispersal.rate.gradient <- seq(from = min.rate,
                               to = max.rate, 
                               length.out = num.dispersal.categories)

dispersal.models <- paste("dk",sprintf("%02d", 1:num.dispersal.categories),sep="")

landscape.dispersal <- dispersal.rate.gradient[which(dispersal.models == dispersal.cat)]

# -------------------------------------------------------------------------
# first, create landscape matrix

my.landscape <- NLMR::nlm_fbm(ncol = ncol,nrow = nrow,fract_dim = landscape.sa)
# show_landscape(my.landscape)
my.matrix <- raster::as.matrix(my.landscape)


# -------------------------------------------------------------------------
# second, generate metaweb

my.dist <- extraDistr::rtpois(n = richness,lambda = landscape.lambda,a = 0)
# make sum even
if (sum(my.dist) %% 2 != 0) { my.dist[1] <- my.dist[1] + 1 }
my.net <- igraph::sample_degseq(my.dist,method = "vl")
my.matrix <- as.matrix(igraph::as_adjacency_matrix(my.net,type = "both"))
# assign interaction strengths according to an "extended" normal dist
weights <- abs(gamlss.dist::rSHASHo(sum(my.dist), mu = int.mean, 
                                    sigma = int.sd, nu = 0, tau = tau))
my.matrix[my.matrix == 1] <- weights
diag(my.matrix) <- 1

# -------------------------------------------------------------------------
# third, assign species presences and absences based on landscape suitability
landscape.n.rows <- richness*cells
df.names <- expand.grid(sp.names,1:cells)
landscape.names <- paste(df.names[,1],df.names[,2],sep="-")

# create landscape matrix template
landscape.template <- matrix(0,
                             nrow = landscape.n.rows,
                             ncol = landscape.n.rows,
                             dimnames = list(landscape.names,
                                             landscape.names))

my.landscape.matrix <- landscape.template
my.landscape <- landscape.list[[i.land]][[i.rep]]

# list holding presences per cell
presence.df <- list()

for(i.row in 1:landscape.rows){
  for(i.col in 1:landscape.cols){
    
    # "reload" the metaweb for this replicate
    my.network <- sim.matrices[[i.net]][[i.rep]]
    
    # cells are referred to by its row and col numbers, but I need
    # a cell number for the block matrix (see below)
    cell.id <- (landscape.cols * (i.row - 1)) + i.col
    
    # suitability value of this cell
    landscape.value <- landscape.list[[i.land]][[i.rep]][i.row,i.col]
    
    # id of this cell
    cell.value <- cell.coords$cell[which(cell.coords$x == i.col &
                                           cell.coords$y == i.row)]
    
    # presence dataframe for this cell
    cell.df <- tidyr::expand_grid(landscape.category = landscape.categories[i.land],
                                  network.category = network.categories[i.net],
                                  replicate = i.rep,
                                  landscape.row = i.row,
                                  landscape.col = i.col,
                                  cell = cell.value,
                                  sp = sp.names,
                                  presence = FALSE)
    
    # go through all sp, checking if it should be present in this cell
    # and update the df and network matrix of the cell
    
    # a species is present if the suitability value of the cell
    # is within its mean +- sd 
    
    for(i.sp in 1:num.sp){
      my.suit <- which(suit.df$sp == sp.names[i.sp])
      min.suit <- suit.df$optimum[my.suit] - suit.df$sd[my.suit]
      max.suit <- suit.df$optimum[my.suit] + suit.df$sd[my.suit]
      
      if(landscape.value >= min.suit & landscape.value <= max.suit){
        cell.df$presence[i.sp] <- TRUE
      }else{
        # if species not present, "prune" the matrix, setting the species'
        # elements to 0
        
        my.network[sp.names[i.sp],] <- 0
        my.network[,sp.names[i.sp]] <- 0
        
      }
    }# for i.sp
    presence.df[[length(presence.df)+1]] <- cell.df
    
    # now, add the cell network to the landscape matrix
    init.row <- 1 + (num.sp * (cell.id - 1))
    init.col <- init.row
    end.row <- num.sp + (num.sp * (cell.id - 1))
    end.col <- end.row
    
    my.landscape.matrix[init.row:end.row,init.col:end.col] <- my.network
    
  }# i.col
}# i.row





