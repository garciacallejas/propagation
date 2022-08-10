
# simulate small community and calculate communicability

# -------------------------------------------------------------------------
library(tidyverse)
source("R/auxiliary_functions/comm.R")
source("R/communicability.R")
source("R/horizontal_community_matrix.R")
# -------------------------------------------------------------------------
# three patches, three configurations. This needs to be handcrafted

S <- 4
n.blocks <- 3

connected.community <- matrix(1, nrow = S, ncol = S)
diag(connected.community) <- 0

full.block <- matrix(0,nrow = S*n.blocks, ncol = S*n.blocks)
full.block[1:S,1:S] <- connected.community
full.block[(S+1):(2*S),(S+1):(2*S)] <- connected.community
full.block[(2*S+1):(3*S),(2*S+1):(3*S)] <- connected.community

# set dispersal links
for(i in 1:nrow(full.block)){
  for(j in 1:ncol(full.block)){
    if(i%%S == j%%S & i!=j){
      full.block[i,j] <- 1
      full.block[j,i] <- 1
    }
  }
}

# -------------------------------------------------------------------------
# modified local structure

local.block <- full.block

local.block[2,1] <- local.block[1,2] <- 0
local.block[3,2] <- local.block[2,3] <- 0
local.block[4,1] <- local.block[1,4] <- 0

local.block[5,6] <- local.block[6,5] <- 0
local.block[8,5] <- local.block[5,8] <- 0
local.block[7,8] <- local.block[8,7] <- 0

local.block[9,11] <- local.block[11,9] <- 0
local.block[11,10] <- local.block[10,11] <- 0
local.block[10,12] <- local.block[12,10] <- 0


# -------------------------------------------------------------------------
# modified dispersal
dispersal.block <- matrix(0,nrow = S*n.blocks, ncol = S*n.blocks)
dispersal.block[1:S,1:S] <- connected.community
dispersal.block[(S+1):(2*S),(S+1):(2*S)] <- connected.community
dispersal.block[(2*S+1):(3*S),(2*S+1):(3*S)] <- connected.community

# set dispersal links
for(i in 1:nrow(dispersal.block)){
  for(j in 1:ncol(dispersal.block)){
    if(i%%S == 2 & j%%S == 2 & i!=j){
      dispersal.block[i,j] <- 1
      dispersal.block[j,i] <- 1
    }
  }
}

# -------------------------------------------------------------------------
# modified sp composition

composition.block <- full.block

composition.block[3,] <- composition.block[,3] <- 0
composition.block[4,] <- composition.block[,4] <- 0

composition.block[5,] <- composition.block[,5] <- 0
composition.block[7,] <- composition.block[,7] <- 0

composition.block[9,] <- composition.block[,9] <- 0
composition.block[10,] <- composition.block[,10] <- 0

# -------------------------------------------------------------------------

comm(full.block)
comm(dispersal.block)
comm(local.block)
comm(composition.block)
