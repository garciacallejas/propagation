
# generate dispersal kernels for a set of species by sampling from 
# the exponential distribution.
# It is possible to correlate the dispersal distance of species 
# with their degree, by setting the flag "degree.correlated" to TRUE.

# INPUTS
# simulated network structures: "results/sim_degree_dist_networks.csv"

# OUTPUTS
# dataframe with dispersal distance estimated for every sp
# -------------------------------------------------------------------------

library(tidyverse)

# -------------------------------------------------------------------------
# is dispersal correlated with species degree?
degree.correlated <- TRUE

out.file <- "results/dispersal_kernels.csv"
if(degree.correlated){
  out.file <- "results/dispersal_kernels_correlated.csv"
}

deg.dist <- read.csv2("results/sim_degree_dist_networks.csv")

param <- read.csv2("results/sim_landscape_matrices/parameters_v3.csv")
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

degree.cat <- sort(unique(deg.dist$generative.model))

for(i in 1:length(dispersal.rate.gradient)){
  for(i.rep in 1:num.category.replicates){
    
    my.dist <- rexp(n = richness,rate = dispersal.rate.gradient[i])
    
    # if dispersal rates are correlated with species degree, 
    # get the degrees, sort species, and assign dispersal rates
    # note that dispersal categories are equated to degree categories (dk01 ~ dd01, etc)
    if(degree.correlated){
      my.net <- subset(deg.dist,generative.model == degree.cat[i] & replicate == i.rep)
      my.deg <- my.net %>% group_by(node_from) %>%
        summarise(degree = sum(value != 0)) %>%
        arrange(desc(degree))
      
      my.data <- data.frame(sp = my.deg$node_from,
                            dispersal.category = dispersal.categories[i],
                            exponential.rate = dispersal.rate.gradient[i],
                            replicate = i.rep,
                            dispersal.distance = sort(my.dist,decreasing = T))
    }else{
      my.data <- data.frame(sp = sp.names,
                            dispersal.category = dispersal.categories[i],
                            exponential.rate = dispersal.rate.gradient[i],
                            replicate = i.rep,
                            dispersal.distance = my.dist)
    }
    
    disp.df <- bind_rows(disp.df,my.data)
  }
}

# -------------------------------------------------------------------------

write.csv2(disp.df,out.file,row.names = FALSE)

