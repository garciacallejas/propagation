
# compare net effects (dynamics) and communicability (structure)

library(tidyverse)
source("R/niche_model.R")
source("R/GLV_functions.R")
source("R/communicability.R")
source("R/horizontal_community_matrix.R")

# -------------------------------------------------------------------------

# how many simulated networks
num.nets <- 100
# richness
S <- seq(10,100,by = 10)
# connectance
C <- .25

net.type <- "horizontal"
# net.type <- "food web"

# -------------------------------------------------------------------------
# signs are the "real" ones in these functions, e.g. competition is negative
# interaction.matrix <- -interaction.matrix
# diag(interaction.matrix) <- -1
# -------------------------------------------------------------------------

net.metrics.list <- list()

for(i.s in 1:length(S)){
  my.s <- S[i.s]
  for(i.net in 1:num.nets){
    if(net.type == "horizontal"){
      A <- horizontal_community_matrix(S = my.s,c = C)
    }else{
      A <- niche_model(my.s,C)
    }
    
    # ensure zero diagonal?
    # diag(A) <- 0
    
    # parameters of the dynamical system
    N0 <- rep(1,my.s) #initial conditions for species abundances
    r0 <- rep(1,my.s) # assume al growth rates are positive
    parameters <- list(r = r0, A = A)
    
    # solve the dynamical system
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
      
      J <- matrix(0,nrow = nrow(A),
                  ncol = ncol(A))
      
      invJ <- matrix(0,nrow = nrow(A),
                     ncol = ncol(A))
    }# if-else feasible solution
    
    # calculate communicability
    comm <- communicability(A)
    bin.comm <- comm[[1]]
    w.comm <- comm[[2]]
    
    # generate dataframe
    df1 <- expand.grid(net = i.net,
                       structure = "niche",
                       sp1 = 1:my.s, sp2 = 1:my.s, 
                       direct.effect = 0,
                       net.effect = 0,
                       binary.communicability = 0,
                       weighted.communicability = 0)
    
    for(i.pos in 1:nrow(df1)){
      df1$richness <- my.s
      df1$direct.effect[i.pos] <- J[df1$sp1[i.pos],df1$sp2[i.pos]]
      df1$net.effect[i.pos] <- invJ[df1$sp1[i.pos],df1$sp2[i.pos]]
      df1$binary.communicability[i.pos] <- bin.comm[df1$sp1[i.pos],df1$sp2[i.pos]]
      df1$weighted.communicability[i.pos] <- w.comm[df1$sp1[i.pos],df1$sp2[i.pos]]
      
    }
    
    net.metrics.list[[length(net.metrics.list)+1]] <- df1
    
  }# for i.net
}# for i.s
net.metrics.df <- bind_rows(net.metrics.list)

# -------------------------------------------------------------------------
# cor(df1$net.effect,df1$weighted.communicability)

metrics.grouped <- net.metrics.df %>%
  group_by(richness) %>%
  summarise(cor.net.w = cor(abs(net.effect),abs(weighted.communicability)),
            cor.net.bin = cor(abs(net.effect),abs(binary.communicability)),
            cor.dir.w = cor(abs(direct.effect),abs(weighted.communicability)),
            cor.dir.bin = cor(abs(direct.effect),abs(binary.communicability)),)

