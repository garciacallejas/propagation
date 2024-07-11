# adapted from https://stefanoallesina.github.io/Sao_Paulo_School/multi.html#how-many-species-will-coexist
# and my own SME repo

library(tidyverse)
library(deSolve)

# model_stode = function(t,y,parms=NULL,A,D) {
#     dy = y*(1+A%*%y) + D%*%y
#     return(list(dy,1))
# }

# Generalized Lotka-Volterra model
GLV <- function(t, y, parameters){
    with(as.list(c(y, parameters)), {
        y[y < 10^-8] <- 0 # prevent numerical problems
        dydt <- y * (r + A %*% y)
        list(dydt)
    })
}

# include intraspecific density dependence explicitly (see Hale & Valdovinos 2021)
GLV_dd <- function(t, y, parameters){
  with(as.list(c(y, parameters)), {
    y[y < 10^-8] <- 0 # prevent numerical problems
    dydt <- y * (r + A %*% y - s * y) 
    list(dydt)
  })
}

# mutualistic bipartite dynamics
# besides all coefficients, include "num.plants" and "num.pols" in parameters
mutualistic_model <- function(t,y,parameters){
  with(as.list(c(y,parameters)), {
    
    num.plants <- length(rplants)
    num.pols <- length(rpols)
    
    y[y < 10^-10] <- 0
    
    dplants <- rep(0,num.plants)
    dpols <- rep(0,num.pols)
    
    for(i.plant in 1:num.plants){
      
      # sum over all pollinators
      sum.interactions <- 0
      for(i.pol in 1:num.pols) {
        pol.n.index <- i.pol+num.plants
        sum.interactions <- sum.interactions + B_plpo[i.plant,i.pol] * y[pol.n.index]
      }

      dplants[i.plant] <- (rplants[i.plant] + sum.interactions - (alphaplants[i.plant] + cplants[i.plant]*sum.interactions)*y[i.plant]) * y[i.plant]
    }
    
    for(i.pol in 1:num.pols){
      pol.n.index <- i.pol+num.plants
      # sum over all plants
      sum.interactions <- 0
      for(i.plant in 1:num.plants) {
        sum.interactions <- sum.interactions + B_popl[i.pol,i.plant] * y[i.plant]
      }
      
      dpols[i.pol] <- (rpols[i.pol] + sum.interactions - (alphapols[i.pol] + cpols[i.pol]*sum.interactions)*y[pol.n.index]) * y[pol.n.index]
    }
    
    list(c(dplants,dpols))
    
  })
}

# function to plot output
tidy_ODE_output <- function(out){
    out <- as.data.frame(out)
    # colnames(out) <- c("time", paste("sp", 1:(ncol(out) -1), sep = "_"))
    out <- as_tibble(out) %>% gather(species, density, -time)
    # pl <- ggplot(data = out) + 
    #     aes(x = time, y = density, colour = species) + 
    #     geom_line()
    # show(pl)
    return(out)
}
# general function to integrate GLV
integrate_GLV <- function(r, A, x0, maxtime = 100, steptime = 0.5){
    times <- seq(0, maxtime, by = steptime)
    parameters <- list(r = r, A = A)
    # solve numerically
    out <- ode(y = x0, times = times, 
               func = GLV, parms = parameters, 
               method = "ode45")
    # make into tidy form
    out2 <- tidy_ODE_output(out)
    return(out)
}

build_LDstable <- function(n){
    A <- matrix(0, n, n)
    A[upper.tri(A)] <- rnorm(n * (n - 1) / 2)
    # make symmetric
    A <- A + t(A)
    # now find the largest eigenvalue
    l1A <- max(eigen(A, only.values = TRUE, symmetric = TRUE)$values)
    if (l1A > 0){
        # set the diagonal to make it stable
        diag(A) <- diag(A) - l1A - 0.01
    }
    return(A)
}
