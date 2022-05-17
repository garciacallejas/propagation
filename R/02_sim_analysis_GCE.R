
# analyze GCE patterns in the simulations


# -------------------------------------------------------------------------

library(tidyverse)
library(ggridges)
library(patchwork)

# -------------------------------------------------------------------------

# param <- read.csv2("results/sim_landscape_matrices/parameters_v2.csv")
# 
# network.categories <- read.csv2("results/network_gradient_categories.csv")
# landscape.categories <- read.csv2("results/spatial_autocorrelation_categories.csv")
# dispersal.categories <-  read.csv2("results/dispersal_kernels.csv")
# cell.distances <- read.csv2("results/cell_distances.csv")
# cell.distances$cell_from <- as.character(cell.distances$cell_from)
# cell.distances$cell_to <- as.character(cell.distances$cell_to)

gce.df <- read.csv2("results/sim_GCE.csv")

# -------------------------------------------------------------------------

land.dist <- ggplot(gce.df, aes(x = normalised.gce, y = landscape.category)) + 
  geom_density_ridges()
net.dist <- ggplot(gce.df, aes(x = normalised.gce, y = network.category)) + 
  geom_density_ridges()
disp.dist <- ggplot(gce.df, aes(x = normalised.gce, y = dispersal.category)) + 
  geom_density_ridges()

# land.dist/net.dist/disp.dist

m1 <- lm(normalised.gce ~ landscape.category + network.category + dispersal.category, 
         data = gce.df)
