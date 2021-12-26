# generate network matrices to analyze propagation patterns

library(tidyverse)

# -------------------------------------------------------------------------
# create interaction matrix, first take

# 1 - richness

community.size <- 50

# 2 - type of structure
# specify connectance for now and assign links randomly

connectance <- .2

# 1 - type of interactions
# try negative interactions (competition)
# and negative diagonals with stronger intraspecific competition

num.interactions <- connectance * community.size^2
interaction.strengths <- (abs(rnorm(num.interactions, mean = 0, sd = 1)))*-1

# -------------------------------------------------------------------------

interaction.matrix <- matrix(0,nrow = community.size, ncol = community.size)
interaction.positions <- sample(length(interaction.matrix),
                                num.interactions,
                                replace = FALSE)
interaction.matrix[interaction.positions] <- interaction.strengths

diag(interaction.matrix) <- -1

# -------------------------------------------------------------------------

save(interaction.matrix, file = "data/interaction_matrix.Rdata")

# -------------------------------------------------------------------------

linear.chain.matrix <- matrix(0,nrow = community.size, ncol = community.size)

for(i.pos in 1:(nrow(linear.chain.matrix)-1)){
  linear.chain.matrix[i.pos,i.pos+1] <- 1
}

save(linear.chain.matrix,file = "data/linear_chain_matrix.Rdata")
