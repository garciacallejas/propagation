

# -------------------------------------------------------------------------

library(tidyverse)
library(patchwork)
library(igraph)
library(ggraph)
library(tidygraph)

# -------------------------------------------------------------------------
load("results/sim_degree_dist_network_matrices.RData")
nets <- read.csv2("results/sim_degree_dist_networks.csv")
degree.values <- read.csv2("results/network_gradient_categories.csv")

# -------------------------------------------------------------------------

rep <- 1
lambda.values <- degree.values$poisson.lambda
categories <- degree.values$network.category

plot.list <- list()

for(i in 1:length(lambda.values)){
  my.net <- subset(nets, generative.model == categories[i] &
                     replicate == rep)
  my.matrix <- sim.matrices[[i]][[rep]]
  
  my.connectance <- round(sum(my.matrix != 0)/nrow(my.matrix)^2,2)
  
  # 1 - nodes/edges dataframes
  g  <- graph.adjacency(my.matrix,weighted=TRUE)
  my.edge.list <- get.data.frame(g)
  id.list <- unique(my.edge.list$from)
  
  my.net <- tbl_graph(edges = my.edge.list,
                      directed = TRUE)
  
  p1 <- ggraph(my.net) + 
    
    geom_edge_link0(aes(edge_width = weight),
                    edge_colour = "#A8A8A8",
                    # edge_width = 0.1,
                    edge_alpha = .5) +
    # geom_edge_link(aes(edge_width = weight)) +
    
    scale_edge_width(range = c(0.25, 1.5)) +
    
    # geom_node_point() +
    geom_node_point(fill = "darkorange",
                    colour = "#000000",
                    size = 3,
                    shape = 21,
                    stroke = 0.3) +
    ggtitle(paste(length(id.list), "sp\ndeg. dist with lambda:",lambda.values[i],"\nconnectance:",my.connectance)) +
    theme_graph() +
    # ggtitle(paste("plot",i.plot,sep = " ")) +
    theme(legend.position = "none") +
    theme(
      plot.title = element_text(size = 10),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA)
    )+
    NULL
  plot.list[[i]] <- p1
}

all.plots <- patchwork::wrap_plots(plot.list,ncol = 3)

# -------------------------------------------------------------------------
# ggsave("results/images/network_plots.pdf",
#        plot = all.plots,
#        device = cairo_pdf,
#        width = 8, height = 6,dpi = 300)

# -------------------------------------------------------------------------


