
# communicability_ij = exp(A)ij

# library(tidyverse)
library(expm) # for matrix exponentiation
# library(igraph)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# -------------------------------------------------------------------------
load("results/community_block_matrix.Rdata")

# these are for building the list of interactions
clean.int.data <- read.csv2("data/plant_bird_clean_interaction_data.csv")
bird.sp.cells.wide <- read.csv2("data/bird_cell_presences.csv")
plant.sp.cells.wide <- read.csv2("data/plant_cell_presences.csv")

# clean.int.data and bird/plant.sp.cells.wide should be consistent, CHECK
# all.sp <- sort(unique(c(clean.int.data$PLANTSPECIES,clean.int.data$BIRDSPECIES)))
all.sp <- sort(unique(c(bird.sp.cells.wide$species,plant.sp.cells.wide$species)))
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
  # beware, it takes at least > 20 min for a 200M elements matrix
  binary.comm.matrix <- expm(block.matrix)
  weighted.comm.matrix <- 0
  
}else{
  # binary and weighted
  binary.matrix <- ifelse(block.matrix != 0, 1, 0)
  binary.comm.matrix <- expm(binary.matrix)
  # TODO consider weighting these elements by community size
  
  # weighted
  abs.int.matrix <- abs(block.matrix)
  degree.diag.matrix <- matrix(0,nrow = nrow(binary.matrix),
                               ncol = ncol(binary.matrix))
  diag(degree.diag.matrix) <- rowSums(abs.int.matrix) # strength - suma de los pesos
  
  d.inverse <- solve(degree.diag.matrix)
  d.pow <- sqrt(d.inverse)
  
  multiplied.matrix <- d.pow %*% abs.int.matrix %*% d.pow # esto es eq. 2.2 de Crofts and Higham 2009
  weighted.comm.matrix <- expm(multiplied.matrix)
}

# -------------------------------------------------------------------------

if(is.binary){
  df1 <- expand.grid(sp1 = 1:length(all.sp),
                     grid.id.sp1 = represented.cells,
                     sp2 = 1:length(all.sp),
                     grid.id.sp2 = represented.cells,
                     binary.communicability = 0)
  
  for(i.pos in 1:nrow(df1)){
    
    my.sp1 <- which(all.sp == df1$sp1[i.pos])
    my.grid.id.sp1 <- df1$grid.id.sp1[i.pos]
    my.sp2 <- df1$sp2[i.pos]
    my.grid.id.sp2 <- df1$grid.id.sp2[i.pos]
    
    my.matrix.row <-  num.sp * (my.grid.id.sp1 - 1) + my.sp.num
    
    df1$binary.communicability[i.pos] <- binary.comm.matrix[df1$x[i.pos],df1$y[i.pos]]
    # df1$weighted.communicability[i.pos] <- weighted.comm.matrix[df1$x[i.pos],df1$y[i.pos]]
  }
  
  df1$scaled.binary.communicability <- range01(df1$binary.communicability)
  df1$diag <- ifelse(df1$x == df1$y,TRUE,FALSE)

}else{
  df1 <- expand.grid(sp1 = 1:length(all.sp),
                     grid.id.sp1 = represented.cells,
                     sp2 = 1:length(all.sp),
                     grid.id.sp2 = represented.cells,
                     binary.communicability = 0,
                     weighted.communicability = 0)
  
  for(i.pos in 1:nrow(df1)){
    df1$binary.communicability[i.pos] <- binary.comm.matrix[df1$x[i.pos],df1$y[i.pos]]
    df1$weighted.communicability[i.pos] <- weighted.comm.matrix[df1$x[i.pos],df1$y[i.pos]]
  }
  
  df1$scaled.binary.communicability <- range01(df1$binary.communicability)
  df1$diag <- ifelse(df1$x == df1$y,TRUE,FALSE)
}

# -------------------------------------------------------------------------

write.csv2(df1,"results/communicability_pairwise.csv",row.names = FALSE)

# -------------------------------------------------------------------------
# load("data/interaction_matrix.Rdata")
# num <- 100
# interaction.matrix[sample(length(interaction.matrix),num,FALSE)] <- abs(rnorm(num,0,5))
# # diag(interaction.matrix) <- 1
# richness <- nrow(interaction.matrix)
# 
# # binary
# lc.bin <- ifelse(interaction.matrix != 0, 1, 0)
# lc.comm.matrix <- expm(interaction.matrix)
# # TODO consider weighting these elements by community size
# 
# # weighted
# lc.diag.matrix <- matrix(0,nrow = nrow(interaction.matrix), 
#                              ncol = ncol(interaction.matrix))
# diag(lc.diag.matrix) <- rowSums(interaction.matrix) # strength - suma de los pesos
# 
# lc.d.inverse <- solve(lc.diag.matrix)
# lc.d.pow <- sqrt(lc.d.inverse)
# 
# lc.multiplied.matrix <- lc.d.pow %*% abs(interaction.matrix) %*% lc.d.pow # esto es eq. 2.2 de Crofts and Higham 2009
# lc.weighted.comm.matrix <- expm(lc.multiplied.matrix)
# 
# df1 <- expand.grid(x = 1:richness, y = 1:richness, 
#                    binary.communicability = 0,
#                    weighted.communicability = 0)
# 
# for(i.pos in 1:nrow(df1)){
#   df1$binary.communicability[i.pos] <- lc.comm.matrix[df1$x[i.pos],df1$y[i.pos]]
#   df1$weighted.communicability[i.pos] <- lc.weighted.comm.matrix[df1$x[i.pos],df1$y[i.pos]]
# }
# 
# df1$scaled.binary.communicability <- range01(df1$binary.communicability)
# df1$diag <- ifelse(df1$x == df1$y,TRUE,FALSE)
# 
# # -------------------------------------------------------------------------
# 
# write.csv2(df1,"results/communicability_pairwise.csv",row.names = FALSE)

