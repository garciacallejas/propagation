
# communicability_ij = exp(A)ij

library(tidyverse)
library(expm) # for matrix exponentiation
# library(igraph)

# -------------------------------------------------------------------------
load("data/interaction_matrix.Rdata")
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

richness <- nrow(interaction.matrix)
# -------------------------------------------------------------------------

# binary
binary.matrix <- ifelse(interaction.matrix != 0, 1, 0)
binary.comm.matrix <- expm(binary.matrix)
# TODO consider weighting these elements by community size

# weighted
abs.int.matrix <- abs(interaction.matrix)
degree.diag.matrix <- matrix(0,nrow = nrow(binary.matrix),
                             ncol = ncol(binary.matrix))
diag(degree.diag.matrix) <- rowSums(abs.int.matrix) # strength - suma de los pesos

d.inverse <- solve(degree.diag.matrix)
d.pow <- sqrt(d.inverse)

multiplied.matrix <- d.pow %*% abs.int.matrix %*% d.pow # esto es eq. 2.2 de Crofts and Higham 2009
weighted.comm.matrix <- expm(multiplied.matrix)

df1 <- expand.grid(x = 1:richness, y = 1:richness, 
                   binary.communicability = 0,
                   weighted.communicability = 0)

for(i.pos in 1:nrow(df1)){
  df1$binary.communicability[i.pos] <- binary.comm.matrix[df1$x[i.pos],df1$y[i.pos]]
  df1$weighted.communicability[i.pos] <- weighted.comm.matrix[df1$x[i.pos],df1$y[i.pos]]
}

df1$scaled.binary.communicability <- range01(df1$binary.communicability)
df1$diag <- ifelse(df1$x == df1$y,TRUE,FALSE)

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

