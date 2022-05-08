
# statistical analyses on communicability of simulated landscapes

# INPUTS
# - individual communicability dataframes: "results/communicability/.."
# - categories of the different factors


# OUTPUTS


# -------------------------------------------------------------------------
library(tidyverse)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# -------------------------------------------------------------------------
# read general data
param <- read.csv2("results/sim_landscape_matrices/parameters_v2.csv")

network.categories <- read.csv2("results/network_gradient_categories.csv")
landscape.categories <- read.csv2("results/spatial_autocorrelation_categories.csv")
dispersal.categories <-  read.csv2("results/dispersal_kernels.csv")
cell.distances <- read.csv2("results/cell_distances.csv")
cell.distances$cell_from <- as.character(cell.distances$cell_from)
cell.distances$cell_to <- as.character(cell.distances$cell_to)

# -------------------------------------------------------------------------
# recover factors

network.categories <- network.categories$network.category
landscape.categories <- landscape.categories$landscape.category
dispersal.categories <- unique(dispersal.categories$dispersal.category)
replicates <- param$num.category.replicates

richness <- param$richness
cells <- param$ncol * param$nrow

# -------------------------------------------------------------------------
# idea: take a small sample of each dataframe (say 5-10%), create a df,
# and run the model on it. Bootstrap the process 100 or 1000 times.
# This should be theoretically valid, as these are indpendent instances
# of the same model, thus the interpretation of the coefficients does not vary.

i.model <- 1

sample.perc <- .05
n.models <- 100
my.models <- list()

my.files <- list.files("results/communicability",full.names = T)

# for(i.model in 1:n.models){
  
  pairs.data <- list()
  node.data <- list()
  my.files <- my.files[sample(1:length(my.files),20)]
  
  for(i.file in 1:length(my.files)){
    
    # get factor info
    my.net.cat <- substr(my.files[i.file],30,33)
    my.landscape.cat <- substr(my.files[i.file],35,38)
    my.dispersal.cat <- substr(my.files[i.file],40,43)
    
    # get communicability
    load(my.files[i.file])
    
    # subset
    valid.pairs <- subset(comm.df, sp1 != sp2 & weighted.communicability != 0)
    
    # aggregate per node
  my.node.data <- comm.df %>% group_by(sp1,cell1) %>%
    summarise(wc = sum(weighted.communicability),bc = sum(binary.communicability)) %>%
    mutate(network.category = my.net.cat,
           landscape.category = my.landscape.cat,
           dispersal.category = my.dispersal.cat) %>%
    arrange(cell1,sp1)
    
    # extract sample
    my.sample <- valid.pairs[sample(nrow(valid.pairs),round(nrow(valid.pairs)*sample.perc),F),]
    clean.sample <- data.frame(network.category = my.net.cat,
                               landscape.category = my.landscape.cat,
                               dispersal.category = my.dispersal.cat,
                               wc = my.sample[,c("weighted.communicability")])
    
    pairs.data[[length(pairs.data)+1]] <- clean.sample
    node.data[[length(node.data)+1]] <- my.node.data 
  }# for i.file
  
  pairs.df <- bind_rows(pairs.data)
  nodes.df <- bind_rows(node.data)
  
# }# for i.model

  pairs.means <- pairs.df %>%
    group_by(network.category,landscape.category,dispersal.category) %>%
    summarise(mean.wc = mean(wc),sd.wc = sd(wc)) 
  
  m1 <- lm(wc~network.category + landscape.category + dispersal.category,data = pairs.df)
  mm1 <- lm(mean.wc ~ network.category + landscape.category + dispersal.category,data = pairs.means)
  mn1 <- lm(wc ~ network.category + landscape.category + dispersal.category, data = nodes.df)
  
  plot.wc.bc <- ggplot(nodes.df) + 
    geom_point(aes(x = log(bc), y = log(wc))) +
    theme_bw() +
    NULL
  
  plot.node.land <- ggplot(nodes.df,aes(x = network.category,y = wc)) + 
    geom_boxplot(aes(fill = landscape.category)) +
    # facet_grid(dispersal.category~.)+
    theme_bw() +
    # ylim(c(0,1)) +
    NULL
  
  # plot.means.land <- ggplot(result.means) +
  #   geom_point(aes(x = landscape.category, y = mean.wc, fill = network.category), shape = 21) + 
  #   geom_errorbar(aes(x = landscape.category, ymin = mean.wc - sd.wc, ymax = mean.wc + sd.wc)) + 
  #   facet_grid(network.category~dispersal.category)+
  #   theme_bw() + 
  #   NULL
  # 
  # plot.means.net <- ggplot(result.means) +
  #   geom_point(aes(x = network.category, y = mean.wc, fill = landscape.category), shape = 21) + 
  #   geom_errorbar(aes(x = network.category, ymin = mean.wc - sd.wc, ymax = mean.wc + sd.wc)) + 
  #   facet_grid(landscape.category~dispersal.category)+
  #   theme_bw() + 
  #   NULL
  
# plot.land <- ggplot(result.df,aes(x = landscape.category,y = wc)) + 
#   geom_boxplot(aes(fill = network.category)) + 
#   theme_bw() +
#   ylim(c(0,1)) +
#   NULL
# 
# plot.net <- ggplot(result.df,aes(x = network.category,y = wc)) + 
#   geom_boxplot(aes(fill = landscape.category)) + 
#   theme_bw() +
#   ylim(c(0,1)) +
#   NULL


