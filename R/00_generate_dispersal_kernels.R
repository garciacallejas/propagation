
# generate dispersal kernels for a set of species by sampling from 
# the exponential distribution

# INPUTS
# none

# OUTPUTS
# dataframe with dispersal distance estimated for every sp, 
out.file <- "results/dispersal_kernels.csv"

# -------------------------------------------------------------------------

library(tidyverse)

# -------------------------------------------------------------------------
param <- read.csv2("results/sim_landscape_matrices/parameters_v2.csv")

richness <-  param$richness
num.dispersal.categories <- param$num.dispersal.categories
num.category.replicates <- param$num.category.replicates

# exponetial rate
min.rate <- param$min.rate
max.rate <- param$max.rate

# -------------------------------------------------------------------------

sp.names <- paste("sp",sprintf("%02d", 1:richness),sep="")

dispersal.rate.gradient <- seq(from = min.rate,
                            to = max.rate, 
                            length.out = num.dispersal.categories)

dispersal.categories <- paste("dk",sprintf("%02d", 1:num.dispersal.categories),sep="")

# disp.df <- tidyr::expand_grid(sp = sp.names, dispersal.category = dispersal.categories,
#                               replicate = 1:num.category.replicates,
#                               dispersal.distance = NA)
# disp.df$exponential.rate <- dispersal.rate.gradient[as.numeric(sub(".*_", "", 
#                                                                    disp.df$dispersal.category))]
# disp.df <- arrange(disp.df[,c("sp","dispersal.category","exponential.rate","replicate","dispersal.distance")],
#                    sp,dispersal.category,replicate)

disp.df <- NULL

# -------------------------------------------------------------------------

for(i in 1:length(dispersal.rate.gradient)){
  for(i.rep in 1:num.category.replicates){
    
    my.dist <- rexp(n = richness,rate = dispersal.rate.gradient[i])
    
    my.data <- data.frame(sp = sp.names,
                          dispersal.category = dispersal.categories[i],
                          exponential.rate = dispersal.rate.gradient[i],
                          replicate = i.rep,
                          dispersal.distance = my.dist)
    disp.df <- bind_rows(disp.df,my.data)
  }
}

# -------------------------------------------------------------------------

write.csv2(disp.df,out.file,row.names = FALSE)

