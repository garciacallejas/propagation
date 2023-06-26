
# NOTE: this is going to be deprecated, as everything is included in 03_individual...

# Script to plot communicability values in space, grid cells, and more

# INPUTS: 
# - pairwise communicability df: "{external_path}/results/pairwise_communicability.Rdata"
# - pairwise shortest path lengths df: "{external_path}/results/shortest_path_lengths.Rdata"
# - richness per cell df: "results/cell_richness.csv"
# - average communicability df: "results/average_communicability.csv"
# - NZ spatial grid shapefile: "data/NZ_grid.shp"

# OUTPUTS:
# - several plots: "results/images/*.pdf-png"

# NOTE:
# some inputs are too big for git/github. load them externally and keep 
# the path {external_path} always the same

external_path <- "/home/david/Work/datasets/NZ/"

# -------------------------------------------------------------------------

library(tidyverse)
library(sf)
library(patchwork)

# -------------------------------------------------------------------------

load(paste(external_path,"results/communicability_pairwise.Rdata",sep=""))
load(paste(external_path,"results/shortest_path_lengths.Rdata",sep=""))

cell.richnes <- read.csv2("results/cell_richness.csv")
avg.comm <- read.csv2("results/average_communicability.csv")

NZ.grid <- st_read("data/NZ_grid.shp")

# taken directly from the grid. It may be the case that there are cells
# without observations, i.e. not represented in sp.obs
cell.id <- sort(unique(NZ.grid$grid_id))

# grid_lab <- st_centroid(NZ.grid) %>% cbind(st_coordinates(.))
# 
# ggplot() +
#   # geom_sf(data = NZ, fill = 'white', lwd = 0.05) +
#   # geom_sf(data = pts, color = 'red', size = 1.7) +
#   geom_sf(data = NZ.grid, fill = 'transparent', lwd = 0.3) +
#   geom_text(data = grid_lab, aes(x = X, y = Y, label = grid_id), size = 2) +
#   coord_sf(datum = NA)  +
#   labs(x = "") +
#   labs(y = "")

# ggsave("results/images/NZ_grid.png",width = 6,height = 6)

# cell 46 includes christchurch
# propagation from plants in cell 46 to 1) birds everywhere, 2) plants everywhere

my.source.cell <- 46
my.source.guild <- "plants"

my.propagation.data <- subset(avg.comm,source.guild == my.source.guild & 
                                source.cell == my.source.cell)

my.propagation.data <- my.propagation.data[,c("recipient.cell","recipient.guild","scaled.avg.communicability")]
names(my.propagation.data) <- c("grid_id","guild","comm")
my.propagation.wide <- pivot_wider(my.propagation.data,names_from = guild,values_from = comm)

# -------------------------------------------------------------------------
# richness

cell.richness.wide <- pivot_wider(cell.richnes,names_from = guild,values_from = richness)
names(cell.richness.wide)[1] <- "grid_id"

# -------------------------------------------------------------------------
# join propagation/metrics data with grid

grid.plot.data <- left_join(NZ.grid, my.propagation.wide, by = "grid_id")

grid.metrics.data <- left_join(NZ.grid,cell.richness.wide, by = "grid_id")
grid.metrics.data[is.na(grid.metrics.data)] <- 0

# -------------------------------------------------------------------------
# plots

# richness
bird.richness.plot <- ggplot(data = grid.metrics.data) + 
  geom_sf(aes(fill = birds)) + 
  scale_fill_viridis_c() +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "") +
  theme_bw() +
  NULL

plant.richness.plot <- ggplot(data = grid.metrics.data) + 
  geom_sf(aes(fill = plants)) + 
  scale_fill_viridis_c() +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "") +
  theme_bw() +
  NULL

# propagation
bird.plot <- ggplot(data = grid.plot.data) + 
  geom_sf(aes(fill = birds)) + 
  scale_fill_viridis_c() +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "") +
  theme_bw() +
  NULL

plant.plot <- ggplot(data = grid.plot.data) + 
  geom_sf(aes(fill = plants)) + 
  scale_fill_viridis_c() +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "") +
  theme_bw() +
  NULL

metrics.plot <- plant.richness.plot + bird.richness.plot + plot_annotation(tag_levels = 'A')
prop.plot <- plant.plot + bird.plot + plot_annotation(tag_levels = 'A')

# -------------------------------------------------------------------------

ggsave(filename = "results/images/richness_plot.pdf",
       plot = metrics.plot,
       device = cairo_pdf,
       width = 10, height = 6,dpi = 300)
ggsave(filename = "results/images/richness_plot.png",
       plot = metrics.plot,
       # device = cairo_pdf,
       width = 10, height = 6,dpi = 300)

ggsave(filename = paste("results/images/propagation_plot_",my.source.cell,"_",my.source.guild,".pdf",sep=""),
       plot = prop.plot,
       device = cairo_pdf,
       width = 10, height = 6,dpi = 300)
ggsave(filename = paste("results/images/propagation_plot_",my.source.cell,"_",my.source.guild,".png",sep=""),
       plot = prop.plot,
       # device = cairo_pdf,
       width = 10, height = 6,dpi = 300)


