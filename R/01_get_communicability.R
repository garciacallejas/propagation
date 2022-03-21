
# Script to obtain communicability matrix from interactions

# INPUTS: 
# - supra-adjacency matrix: "results/community_block_matrix.Rdata"
# - optional - interaction data: "data/plant_bird_clean_interaction_data.csv"
# - bird presences in cells: "data/bird_cell_presences.csv"
# - plant presences in cells: "data/plant_cell_presences.csv"

# OUTPUTS:
# - binary communicability matrix: "{external_path}/results/binary_communicability_matrix.Rdata"
# - weighted communicability matrix: "{external_path}/results/weighted_communicability_matrix.Rdata"
# - pairwise communicability dataframe: "{external_path}/results/communicability_pairwise.Rdata"

# NOTE:
# these outputs are too big for git/github. Save them externally and keep 
# the path {external_path} always the same

external_path <- "/home/david/Work/datasets/NZ/"

# -------------------------------------------------------------------------

library(data.table)
# library(tidyverse)
library(expm) # for matrix exponentiation
# library(igraph)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# -------------------------------------------------------------------------
load("results/community_block_matrix.Rdata")

# these are for building the list of interactions
# clean.int.data <- read.csv2("data/plant_bird_clean_interaction_data.csv")
bird.sp.cells.wide <- read.csv2("data/bird_cell_presences.csv")
plant.sp.cells.wide <- read.csv2("data/plant_cell_presences.csv")

# clean.int.data and bird/plant.sp.cells.wide should be consistent, CHECK
# all.sp <- sort(unique(c(clean.int.data$PLANTSPECIES,clean.int.data$BIRDSPECIES)))
bird.sp <- sort(unique(bird.sp.cells.wide$species))
plant.sp <- sort(unique(plant.sp.cells.wide$species))
all.sp <- sort(unique(c(bird.sp,plant.sp)))
num.sp <- length(all.sp)

represented.cells <- unique(c(names(bird.sp.cells.wide),names(plant.sp.cells.wide)))
represented.cells <- as.numeric(substr(represented.cells[-1],
                                       start = 2,
                                       stop = nchar(represented.cells[-1])))
# -------------------------------------------------------------------------

# is my matrix binary?
is.binary <- sum(block.matrix > 1) == 0
if(is.binary){
  # binary
  # beware, it takes a while for a 200M elements matrix - RAM is a limiting factor
  binary.comm.matrix <- expm(block.matrix)
  weighted.comm.matrix <- 0
}else{
  # binary and weighted
  binary.matrix <- ifelse(block.matrix != 0, 1, 0)
  binary.comm.matrix <- expm(binary.matrix)
  # TODO consider weighting these elements by community size

  # weighted
  
  scaled.comm.matrix <- range01(block.matrix)
  weighted.comm.matrix <- expm(scaled.comm.matrix)
  
  # abs.int.matrix <- abs(block.matrix)
  # degree.diag.matrix <- matrix(0,nrow = nrow(binary.matrix),
  #                              ncol = ncol(binary.matrix))
  # diag(degree.diag.matrix) <- rowSums(abs.int.matrix) # strength - suma de los pesos
  # 
  # d.inverse <- solve(degree.diag.matrix)
  # d.pow <- sqrt(d.inverse)
  # 
  # multiplied.matrix <- d.pow %*% abs.int.matrix %*% d.pow # eq. 2.2 of Crofts and Higham 2009
  # weighted.comm.matrix <- expm(multiplied.matrix)
}

save(binary.comm.matrix,file = paste(external_path,"results/binary_communicability_matrix.Rdata"),sep="")
save(weighted.comm.matrix,file = paste(external_path,"results/weighted_communicability_matrix.Rdata"),sep="")

rm(block.matrix)
gc(verbose = F)

# -------------------------------------------------------------------------
# convert matrix to dataframe

# load("results/binary_communicability_matrix.Rdata")
# load("results/weighted_communicability_matrix.Rdata")

# tidy functions only work with dataframes, but this still works, and is fast
df1 <- reshape2::melt(binary.comm.matrix,value.name = "binary.communicability")

df1$sp1 <- sub("\\-.*", "", df1$Var1)
df1$grid.id.sp1 <- sub(".*-", "", df1$Var1)
df1$sp2 <- sub("\\-.*", "", df1$Var2)
df1$grid.id.sp2 <- sub(".*-", "", df1$Var2)  

df1$scaled.binary.communicability <- range01(df1$binary.communicability)

df1$diag <- ifelse(df1$sp1 == df1$sp2 & df1$grid.id.sp1 == df1$grid.id.sp2,TRUE,FALSE)

df1$guild.sp1 <- ifelse(df1$sp1 %in% bird.sp, "birds", "plants")
df1$guild.sp2 <- ifelse(df1$sp2 %in% bird.sp, "birds", "plants")

df1 <- df1[,c("sp1","guild.sp1","grid.id.sp1",
              "sp2","guild.sp2","grid.id.sp2",
              "diag",
              "binary.communicability",
              "scaled.binary.communicability")]

if(!is.binary){
  # this should be valid because the two matrices have the same dimensions and names
  dfw <- reshape2::melt(weighted.comm.matrix,value.name = "weighted.communicability")
  df1$weighted.communicability <- dfw$weighted.communicability
}

# -------------------------------------------------------------------------

# it is unfeasible to save a csv file of this size
save(df1,file = paste(external_path,"results/communicability_pairwise.Rdata",sep=""))
# write_csv2(df1,"results/communicability_pairwise.csv")
