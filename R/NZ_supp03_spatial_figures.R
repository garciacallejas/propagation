

# -------------------------------------------------------------------------

library(tidyverse)
# library(sf)
library(patchwork)
library(ggrepel)
library(colorblindr)
library(sf)

# -------------------------------------------------------------------------

grid.size <- 10 # km

sp.obs <- read.csv2(paste("results/sp_observations_long_",grid.size,"km.csv",sep=""))
sp.obs.record <- read.csv2(paste("results/sp_observations_record_type_long_",grid.size,"km.csv",sep=""))
sp.obs.status <- read.csv2(paste("results/sp_observations_status_long_",grid.size,"km.csv",sep=""))

NZ.grid <- st_read(paste("data/NZ_grid_",grid.size,"km.shp",sep=""))

grid.record.type <- left_join(NZ.grid,sp.obs.record) %>%
  filter(!is.na(guild))

grid.status <- left_join(NZ.grid,sp.obs.status) %>%
  filter(!is.na(guild))

# -------------------------------------------------------------------------

# it's best one plot per combination, otherwise plots are difficult to see

# record.bird.gbif <- subset(grid.record.type,guild == "birds" & record.type == "GBIF")
record.plant.nvs <- subset(grid.record.type,guild == "plants" & record.type == "NVS")

record.bird.gbif <- subset(grid.record.type,guild == "birds" & record.type == "GBIF")
record.plant.gbif <- subset(grid.record.type,guild == "plants" & record.type == "GBIF")

record.bird.tier1 <- subset(grid.record.type,guild == "birds" & record.type == "TIER1")
record.plant.tier1 <- subset(grid.record.type,guild == "plants" & record.type == "TIER1")

record.bird.gbif.plot <- ggplot(data = record.bird.gbif) +
  geom_sf(aes(fill = richness)) +
  scale_fill_viridis_c() +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "") +
  ggtitle("birds - GBIF") +
  theme_bw() +
  NULL
# record.bird.gbif.plot

record.plant.gbif.plot <- ggplot(data = record.plant.gbif) +
  geom_sf(aes(fill = richness)) +
  scale_fill_viridis_c() +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "") +
  ggtitle("plants - GBIF") +
  theme_bw() +
  NULL
# record.plant.gbif.plot

record.bird.tier1.plot <- ggplot(data = record.bird.tier1) +
  geom_sf(aes(fill = richness)) +
  scale_fill_viridis_c() +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "") +
  ggtitle("birds - Tier 1") +
  theme_bw() +
  NULL
# record.bird.tier1.plot

# plants NVS
record.plant.nvs.plot <- ggplot(data = record.plant.nvs) +
  geom_sf(aes(fill = richness)) +
  scale_fill_viridis_c() +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "") +
  ggtitle("plants - NVS") +
  theme_bw() +
  NULL
# record.plant.gbif.plot

record.type.plots <- record.bird.gbif.plot + record.bird.tier1.plot + record.plant.gbif.plot + record.plant.nvs.plot

# -------------------------------------------------------------------------

exotic.bird <- subset(grid.status,guild == "birds" & status == "Exotic")
exotic.plant <- subset(grid.status,guild == "plants" & status == "Exotic")

native.bird <- subset(grid.status,guild == "birds" & status == "Native")
native.plant <- subset(grid.status,guild == "plants" & status == "Native")

exotic.bird.plot <- ggplot(data = exotic.bird) +
  geom_sf(aes(fill = richness)) +
  scale_fill_viridis_c() +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "") +
  ggtitle("birds - exotic") +
  theme_bw() +
  NULL
# exotic.bird.plot

exotic.plant.plot <- ggplot(data = exotic.plant) +
  geom_sf(aes(fill = richness)) +
  scale_fill_viridis_c() +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "") +
  ggtitle("plants - exotic") +
  theme_bw() +
  NULL
# exotic.plant.plot

native.bird.plot <- ggplot(data = native.bird) +
  geom_sf(aes(fill = richness)) +
  scale_fill_viridis_c() +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "") +
  ggtitle("birds - native") +
  theme_bw() +
  NULL
# native.bird.plot

native.plant.plot <- ggplot(data = native.plant) +
  geom_sf(aes(fill = richness)) +
  scale_fill_viridis_c() +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "") +
  ggtitle("plants - native") +
  theme_bw() +
  NULL
# native.plant.plot

status.plots <- (native.bird.plot + exotic.bird.plot)/(native.plant.plot + exotic.plant.plot)
# status.plots

# -------------------------------------------------------------------------

ggsave("results/images/richness_record_type.png",
       plot = record.type.plots, width = 15, height = 7)

ggsave("results/images/richness_status.png",
       plot = status.plots, width = 15, height = 15)

