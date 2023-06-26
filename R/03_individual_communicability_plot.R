# Script to plot communicability values in space, grid cells, and more

# INPUTS: 
# - species/cell level communicability: "results/*level_communicability.csv"
# - species observations across cells: "results/sp_observations_long*"
# - NZ spatial grid shapefile: "data/NZ_grid_XXkm.shp"

# OUTPUTS:
# - several plots: "results/images/*.pdf-png"

# -------------------------------------------------------------------------

library(tidyverse)
# library(sf)
library(patchwork)
library(ggrepel)
library(colorblindr)
library(sf)

# -------------------------------------------------------------------------
grid.size <- 10 # km

# which source to use?
observation.source <- "occurrence_prob"

if(observation.source == "observed_records"){
  sp.obs <- read.csv2(paste("data/sp_observations_long_",grid.size,"km.csv",
                              sep=""))
}else if(observation.source == "occurrence_prob"){
  sp.obs <- read.csv2(paste("results/model_occurrences_",grid.size,"km.csv",
                              sep=""))
  # this is to harmonise datasets
  names(sp.obs)[which(names(sp.obs) == "presence")] <- "observations"
  # keep presences only
  sp.obs <- subset(sp.obs, observations == 1)
  
}

sp.comm <- read.csv2(paste("results/species_level_communicability_",grid.size,"km.csv",sep=""))
sp.comm$species <- factor(sp.comm$species, levels = sort(unique(sp.comm$species)))
sp.comm <- dplyr::arrange(sp.comm,desc(scaled.weighted.communicability)) %>%
  mutate(species.rank = rank(-scaled.weighted.communicability,ties.method = "first"))

cell.comm <- read.csv2(paste("results/cell_level_communicability_",grid.size,"km.csv",sep=""))
cell.comm <- dplyr::arrange(cell.comm,desc(scaled.weighted.communicability))
cell.comm$cell_id <- factor(cell.comm$cell_id, levels = unique(cell.comm$cell_id))

NZ.grid <- st_read(paste("data/NZ_grid_",grid.size,"km.shp",sep=""))

# -------------------------------------------------------------------------

sp.plot <- ggplot(data = sp.comm, aes(x = species.rank, 
                                     y = scaled.weighted.communicability)) + 
  geom_point(aes(fill = guild), shape = 21, size = 2) + 
  scale_fill_OkabeIto() +
  theme_bw() +
  scale_x_discrete(breaks=NULL) +
  geom_label_repel(data = subset(sp.comm, scaled.weighted.communicability > 0.25),
            aes(label=species),
            nudge_x = 50) +
  xlab("") + ylab("scaled weighted communicability") +
  NULL
# sp.plot

sp.hist <- ggplot(data = sp.comm, aes(x = scaled.weighted.communicability)) + 
  geom_density(aes(fill = guild), alpha = 0.5) + 
  theme_bw() +
  scale_fill_OkabeIto() +
  xlab("scaled weighted communicability") + ylab("density estimate") +
  NULL
# sp.hist

# -------------------------------------------------------------------------

cell.plot <- ggplot(data = cell.comm, aes(x = cell_id, 
                                       y = scaled.weighted.communicability)) + 
  geom_point() + 
  theme_bw() +
  # scale_x_discrete(breaks=NULL) +
  NULL
# cell.plot

cell.hist <- ggplot(data = cell.comm, aes(x = scaled.weighted.communicability)) + 
  geom_density() + 
  theme_bw() +
  # scale_x_discrete(breaks=NULL) +
  NULL
# cell.hist

# -------------------------------------------------------------------------

cell.richness <- sp.obs %>%
  group_by(cell_id) %>%
  summarise(richness = n())

grid_lab <- st_centroid(NZ.grid) %>% cbind(st_coordinates(.))

# ggplot() +
#   # geom_sf(data = NZ, fill = 'white', lwd = 0.05) +
#   # geom_sf(data = pts, color = 'red', size = 1.7) +
#   geom_sf(data = NZ.grid, fill = 'transparent', lwd = 0.3) +
#   geom_text(data = grid_lab, aes(x = X, y = Y, label = cell_id), size = 2) +
#   coord_sf(datum = NA)  +
#   labs(x = "") +
#   labs(y = "")

# ggsave("results/images/NZ_grid.png",width = 6,height = 6)

# -------------------------------------------------------------------------
# join propagation/metrics data with grid
cell.comm$cell_id <- as.numeric(as.character(cell.comm$cell_id))
grid.communicability.data <- left_join(NZ.grid, cell.comm)
grid.communicability.data.2 <- subset(grid.communicability.data,
                                      !is.na(scaled.weighted.communicability))
grid.lab.subset <- subset(grid_lab,cell_id %in% unique(grid.communicability.data.2$cell_id))

# -------------------------------------------------------------------------
# plots

# communicability
cell.comm.plot <- ggplot(data = grid.communicability.data.2) + 
  geom_sf(aes(fill = scaled.weighted.communicability)) + 
  # geom_text(data = grid.lab.subset, 
  #           aes(x = X, y = Y, label = cell_id), size = 2, color = "lightgrey") +
  scale_fill_viridis_c() +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "") +
  theme_bw() +
  NULL
cell.comm.plot

# -------------------------------------------------------------------------

ggsave("results/images/sp_communicability_distribution.png",
       plot = sp.plot, width = 8, height = 5)

ggsave("results/images/sp_communicability_density.png",
       plot = sp.hist, width = 8, height = 5)

# -------------------------------------------------------------------------

ggsave(paste("results/images/cell_communicability_",grid.size,"km.png",sep=""),
       plot = cell.comm.plot, width = 5, height = 5)

# ggsave(paste("results/images/cell_richness_",grid.size,"km.png",sep=""),
#        plot = richness.plot, width = 8, height = 8)
