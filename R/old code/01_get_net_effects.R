
library(tidyverse)

# -------------------------------------------------------------------------

load("data/interaction_matrix.Rdata")
source("R/GLV_functions.R")
# -------------------------------------------------------------------------
# signs are the "real" ones in these functions, e.g. competition is negative
# interaction.matrix <- -interaction.matrix
# diag(interaction.matrix) <- -1
# -------------------------------------------------------------------------

# solve system
# need initial abundances, growth rates and the interaction matrix
# (positive values are positive interactions)
S <- nrow(interaction.matrix)

N0 <- rep(1,S) #initial conditions for species abundances
# N0 <- N0 / sum(N0)

r0 <- rep(1,S) # assume al growth rates are positive

# print(solve(interaction.matrix, -r0))
parameters <- list(r = r0, A = interaction.matrix)

eq.abund <- rootSolve::stode(y = N0,
                             func=GLV, 
                             parms = parameters,
                             positive = TRUE)[[1]]

# if there is a non-negative steady-state solution, it will be eq.abund
# if there is no positive solution, all eq.abund elements will be zero
if(sum(eq.abund)>0){
  
  # calculate jacobian
  J = rootSolve::jacobian.full(y=eq.abund,func=GLV,parms = parameters)
  
  # negative of the inverse jacobian: net effects matrix
  invJ <- -MASS::ginv(J)
  
}else{
  
  J <- matrix(0,nrow = nrow(interaction.matrix),
              ncol = ncol(interaction.matrix))
  
  invJ <- matrix(0,nrow = nrow(interaction.matrix),
                 ncol = ncol(interaction.matrix))
}

# -------------------------------------------------------------------------

df1 <- expand.grid(x = 1:S, y = 1:S, 
                   direct.effect = 0,
                   net.effect = 0)

for(i.pos in 1:nrow(df1)){
  df1$direct.effect[i.pos] <- J[df1$x[i.pos],df1$y[i.pos]]
  df1$net.effect[i.pos] <- invJ[df1$x[i.pos],df1$y[i.pos]]
}

df1$diag <- ifelse(df1$x == df1$y,TRUE,FALSE)

# -------------------------------------------------------------------------

write.csv2(df1,"results/interaction_effects_pairwise.csv",row.names = FALSE)

