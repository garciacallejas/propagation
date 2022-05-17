
# statistical analyses on communicability of simulated landscapes

# INPUTS
# - individual communicability dataframes: "results/communicability/.."
# - categories of the different factors


# OUTPUTS


# -------------------------------------------------------------------------
library(tidyverse)
library(patchwork)

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

# list of files with communicability of each landscape
my.files <- list.files("results/communicability",full.names = T,pattern = "*.RData")

# set output names
output.nodes <- "results/communicability_nodes.csv"
output.net <- "results/communicability_net.csv"

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

# another approach: take node-level aggregates, and in parallel, obtain 
# landscape-level aggregates: mean + sd communicability, total, and potentially
# network attributes.

# i.model <- 1

# sample.perc <- .05
# n.models <- 100
# my.models <- list()


# for(i.model in 1:n.models){
  
  # pairs.data <- list()
  node.data <- vector(mode = "list", length = length(my.files))
  net.data <- vector(mode = "list", length = length(my.files))
  # my.files <- my.files[sample(1:length(my.files),20)]
  
  for(i.file in 1:length(my.files)){
    
    # get factor info
    my.net.cat <- substr(my.files[i.file],30,33)
    my.landscape.cat <- substr(my.files[i.file],35,38)
    my.dispersal.cat <- substr(my.files[i.file],40,43)
    my.rep <- substr(my.files[i.file],47,48)
    my.rep <- as.numeric(stringr::str_remove(my.rep,"_"))
    
    # get communicability
    load(my.files[i.file])
    
    # subset
    valid.pairs <- subset(comm.df, sp1 != sp2 & weighted.communicability != 0)
    
    # aggregate per node
  my.node.data <- comm.df %>% 
    filter(sp1 != sp2) %>% 
    group_by(sp1,cell1) %>%
    summarise(wc = sum(weighted.communicability),bc = sum(binary.communicability)) %>%
    mutate(network.category = my.net.cat,
           landscape.category = my.landscape.cat,
           dispersal.category = my.dispersal.cat,
           replicate = my.rep) %>%
    arrange(cell1,sp1)
    
  my.net.data <- comm.df %>%
    filter(sp1 != sp2) %>%
    summarise(mean.wc = mean(weighted.communicability),
              sd.wc = sd(weighted.communicability),
              total.wc = sum(weighted.communicability)) %>%
    mutate(network.category = my.net.cat,
           landscape.category = my.landscape.cat,
           dispersal.category = my.dispersal.cat,
           replicate = my.rep)
  
    # extract sample
    # my.sample <- valid.pairs[sample(nrow(valid.pairs),round(nrow(valid.pairs)*sample.perc),F),]
    # clean.sample <- data.frame(network.category = my.net.cat,
    #                            landscape.category = my.landscape.cat,
    #                            dispersal.category = my.dispersal.cat,
    #                            wc = my.sample[,c("weighted.communicability")])
    
    # pairs.data[[length(pairs.data)+1]] <- clean.sample
  
    node.data[[i.file]] <- my.node.data 
    net.data[[i.file]] <- my.net.data
  }# for i.file
  
  # pairs.df <- bind_rows(pairs.data)
  nodes.df <- bind_rows(node.data)
  net.df <- bind_rows(net.data)
  
# -------------------------------------------------------------------------
  write.csv2(nodes.df,output.nodes, row.names = F)
  write.csv2(net.df,output.net, row.names = F)
# }# for i.model

  # nodes.df <- read.csv2("results/communicability_nodes.csv")
  # net.df <- read.csv2("results/communicability_net.csv")
  
  # pairs.means <- pairs.df %>%
  #   group_by(network.category,landscape.category,dispersal.category) %>%
  #   summarise(mean.wc = mean(wc),sd.wc = sd(wc)) 
  
  # m1 <- lm(wc~network.category + landscape.category + dispersal.category,data = pairs.df)
  # mm1 <- lm(mean.wc ~ network.category + landscape.category + dispersal.category,data = pairs.means)
  # mn1 <- lm(wc ~ network.category + landscape.category + dispersal.category, data = nodes.df)
  # 
  # plot.wc.bc <- ggplot(nodes.df) + 
  #   geom_point(aes(x = log(bc), y = log(wc))) +
  #   theme_bw() +
  #   NULL
  # 
  plot.node.land <- ggplot(nodes.df,aes(x = network.category,y = wc)) +
    geom_boxplot(aes(fill = landscape.category)) +
    facet_grid(dispersal.category~.)+
    theme_bw() +
    # ylim(c(0,1)) +
    NULL
  
  node.land.comm <- ggplot(nodes.df,aes(x = wc, y = landscape.category)) + 
    ggridges::geom_density_ridges(
      jittered_points = FALSE,
      position = ggridges::position_points_jitter(width = 0.05, height = 0),
      point_shape = '|', point_size = 1, point_alpha = 1, alpha = 0.1,
    ) + 
    xlim(c(0,100)) +
    xlab("node weighted communicability") +
    ggtitle("Landscape categories")
  
  node.net.comm <- ggplot(nodes.df,aes(x = wc, y = network.category)) + 
    ggridges::geom_density_ridges(
      jittered_points = FALSE,
      position = ggridges::position_points_jitter(width = 0.05, height = 0),
      point_shape = '|', point_size = 1, point_alpha = 1, alpha = 0.1,
    ) +
    xlim(c(0,100)) +
    
    xlab("node weighted communicability") +
    ggtitle("metaweb structure categories")
  
  node.disp.comm <- ggplot(nodes.df,aes(x = wc, y = dispersal.category)) + 
    ggridges::geom_density_ridges(
      jittered_points = FALSE,
      position = ggridges::position_points_jitter(width = 0.05, height = 0),
      point_shape = '|', point_size = 1, point_alpha = 1, alpha = 0.1,
    ) +
    xlim(c(0,100)) +
    xlab("node weighted communicability") +
    ggtitle("Dispersal kernel categories")
  
  category.densities <- node.land.comm/node.net.comm/node.disp.comm
  
  ggsave("results/images/node_communicability_densities2.pdf",
         plot = category.densities,
         device = cairo_pdf,
         width = 10, height = 10,dpi = 300)
  
  
# -------------------------------------------------------------------------
  
  net.df.agg <- net.df %>%
    group_by(network.category,landscape.category,dispersal.category) %>%
    summarise(mean_weighted_communicability = mean(mean.wc),
              sd_weighted_communicability = mean(sd.wc),
              total_weighted_communicability = mean(total.wc))

  mean.net.land <- ggplot(net.df.agg) +
    geom_boxplot(aes(x = landscape.category, 
                     y = mean_weighted_communicability)) + 
    theme_bw() +
    NULL
  
  mean.net.net <- ggplot(net.df.agg) +
    geom_boxplot(aes(x = network.category, 
                     y = mean_weighted_communicability)) + 
    theme_bw() +
    NULL
  
  mean.net.disp <- ggplot(net.df.agg) +
    geom_boxplot(aes(x = dispersal.category, 
                     y = mean_weighted_communicability)) + 
    theme_bw() +
    NULL
  
  # ggplot(iris, aes(x = Sepal.Length, y = Species)) +
  #   ggridges::geom_density_ridges(
  #     jittered_points = TRUE,
  #     position = ggridges::position_points_jitter(width = 0.05, height = 0),
  #     point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7,
  #   )
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


