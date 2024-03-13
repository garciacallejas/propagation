# Script to plot communicability values in space, grid cells, and more

# INPUTS: 
# - population/species/cell level communicability: "results/*level_communicability.csv"
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
library(gghalves)
library(ggbeeswarm)

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

sp.traits <- read.csv2("data/trait_data.csv")
sp.status <- unique(sp.traits[,c("species","guild","status","IUCN_status","NZ_doc_status")])

pop.comm <- read.csv2(paste("results/population_level_communicability_",grid.size,"km.csv",sep=""))
sp.comm <- read.csv2(paste("results/species_level_communicability_",grid.size,"km.csv",sep=""))
sp.comm$species <- factor(sp.comm$species, levels = sort(unique(sp.comm$species)))
sp.comm <- dplyr::arrange(sp.comm,desc(scaled.weighted.communicability)) %>%
  mutate(species.rank = rank(-scaled.weighted.communicability,ties.method = "first")) %>%
  left_join(sp.status) %>%
  mutate(NZ_threatened = ifelse(NZ_doc_status %in% c("at risk","threatened","endangered","relict"),T,F),
         IUCN_threatened = ifelse(IUCN_status %in% c("CR","EN","LR/cd","LR/nt","NT","VU"),T,F))

cell.land.uses <- read.csv2(paste("data/land_use_frequencies_",grid.size,"km.csv",sep=""))
cell.comm <- read.csv2(paste("results/cell_level_communicability_",grid.size,"km.csv",sep=""))
cell.comm <- dplyr::arrange(cell.comm,desc(scaled.weighted.communicability))
cell.comm$cell_id <- factor(cell.comm$cell_id, levels = unique(cell.comm$cell_id))

NZ.grid <- st_read(paste("data/NZ_grid_",grid.size,"km.shp",sep=""))

all.sp <- sort(unique(pop.comm$species))
# -------------------------------------------------------------------------
# species presences
sp.presence.absences <- expand.grid(cell_id = unique(sp.obs$cell_id),
                                    species = all.sp) %>%
  left_join(sp.obs) %>% 
  replace_na(list(observations = 0)) %>%
  rename(presence = observations) %>%
  arrange(cell_id)

# -------------------------------------------------------------------------
# population averages

pop.averages <- pop.comm %>% 
  left_join(sp.presence.absences) %>%
  filter(presence == 1) %>%
  group_by(species,guild) %>%
  summarise(avg.comm = mean(population.weighted.communicability),
            sd.comm = sd(population.weighted.communicability)) %>%
  left_join(sp.status) %>%
  dplyr::arrange(desc(avg.comm)) %>%
  ungroup() %>%
  mutate(species.rank = rank(-avg.comm,ties.method = "first"),
         sp.labels = paste("italic('",gsub("_"," ",species),"')",sep=""),
         scaled.avg.comm = scales::rescale(avg.comm),
         NZ_threatened = ifelse(NZ_doc_status %in% c("at risk","threatened","endangered","relict"),T,F),
         IUCN_threatened = ifelse(IUCN_status %in% c("CR","EN","LR/cd","LR/nt","NT","VU"),T,F))

pop.plot.means <- ggplot(pop.averages, aes(x = guild, y = scaled.avg.comm)) + 
  geom_half_point(aes(color = status), 
                  transformation = position_quasirandom(width = 0.1),
                  side = "l", size = 0.5, alpha = 0.5) +
  geom_half_boxplot(aes(fill = status), side = "r", outlier.shape = NA) +
  scale_fill_OkabeIto(darken = .2) + 
  scale_color_OkabeIto(darken = .2) +
  theme_bw() +
  ylab("average weighted communicability") +
  NULL
# pop.plot.means

pop.avg.rank.plot <- ggplot(data = subset(pop.averages, guild == "birds"), aes(x = species.rank, 
                                      y = scaled.avg.comm)) + 
  # geom_point(aes(fill = guild), shape = 21, size = 4) +
  geom_point(aes(fill = NZ_threatened, shape = NZ_threatened, size = NZ_threatened)) +
  scale_shape_manual(values = c(21,25)) +
  scale_size_manual(values = c(4,6)) +
  scale_fill_OkabeIto(order = c(1,3)) +
  theme_bw() +
  # scale_x_discrete(breaks=NULL) +
  geom_label_repel(data = subset(pop.averages, NZ_threatened == T),
  # geom_label_repel(data = subset(pop.averages, avg.comm > 110000),
                   aes(label=sp.labels),
                   nudge_x = 50,
                   alpha = .8,
                   parse = T) +
  xlab("species rank - population level") + ylab("average population communicability") +
  scale_x_discrete(breaks=NULL) +
  theme(legend.position = "none") +
  NULL
# pop.avg.rank.plot

# -------------------------------------------------------------------------
# species-level ranking

sp.comm$sp.labels <- paste("italic('",gsub("_"," ",sp.comm$species),"')",sep="")

aggr.rank.plot <- ggplot(data = subset(sp.comm, guild == "birds"), aes(x = species.rank, 
                                     y = scaled.weighted.communicability)) + 
  # geom_point(aes(fill = guild), shape = 21, size = 4) +
  geom_point(aes(fill = NZ_threatened, shape = NZ_threatened, size = NZ_threatened)) +
  scale_shape_manual(values = c(21,25)) +
  scale_size_manual(values = c(4,6)) +
  scale_fill_OkabeIto(order = c(1,3)) +
  theme_bw() +
  scale_x_discrete(breaks=NULL) +
  geom_label_repel(data = subset(sp.comm, NZ_threatened == T),
  # geom_label_repel(data = subset(sp.comm, scaled.weighted.communicability > 0.187),
            aes(label=sp.labels),
            alpha = .8,
            nudge_x = 50,
            parse = T) +
  xlab("species rank - full territory") + ylab("aggregate species communicability") +
  # theme(legend.position = "none") +
  NULL
# aggr.rank.plot

sp.hist <- ggplot(data = sp.comm, aes(x = scaled.weighted.communicability)) + 
  geom_density(aes(fill = guild), alpha = 0.5) + 
  theme_bw() +
  scale_fill_OkabeIto() +
  xlab("aggregate species communicability") + ylab("density estimate") +
  NULL
# sp.hist

# -------------------------------------------------------------------------
# combined population and species-level communicability ranking
pop.avg.rank.plot.nolegend <- pop.avg.rank.plot + 
  theme(legend.position = "none") + ggtitle("A)")
sp.rank.plot.title <- aggr.rank.plot + ggtitle("B)") + 
  theme(legend.title = element_blank(), 
        legend.justification=c(.99,.99), legend.position=c(.99,.99))

combined.rank.plot <- pop.avg.rank.plot.nolegend + sp.rank.plot.title

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

# land use to compare with communicability
names(cell.land.uses)[1] <- "cell_id"
cell.prevalent.land.use <- cell.land.uses %>% 
  filter(!habitat_type %in% c("marine_benthic","marine_intertidal","marine_neritic")) %>%
  group_by(cell_id) %>%
  filter(frequency == max(frequency)) %>%
  select(-frequency)
cell.prevalent.land.use$habitat_type[which(cell.prevalent.land.use$habitat_type %in% 
                                             c("anthropic_habitats","grassland"))] <- "crops and pasture"
cell.prevalent.land.use$habitat_type[which(cell.prevalent.land.use$habitat_type %in% 
                                             c("rocky_area"))] <- "rocky areas"
cell.prevalent.land.use$habitat_type[which(cell.prevalent.land.use$habitat_type %in% 
                                             c("alpine_areas"))] <- "alpine areas"

grid.land.use.data <- left_join(NZ.grid,cell.prevalent.land.use)
# -------------------------------------------------------------------------
# plots

# communicability
cell.comm.plot <- ggplot(data = grid.communicability.data.2) + 
  geom_sf(aes(fill = scaled.weighted.communicability),color = NA) + 
  # geom_text(data = grid.lab.subset, 
  #           aes(x = X, y = Y, label = cell_id), size = 2, color = "lightgrey") +
  scale_fill_viridis_c(name = "communicability") +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "") +
  theme_bw() +
  theme(legend.justification=c(1,0), legend.position=c(.99,0.01)) +
  NULL
# cell.comm.plot

# land use
cell.land.use.plot <- ggplot(data = grid.land.use.data) + 
  geom_sf(aes(fill = habitat_type)) + 
  # geom_text(data = grid.lab.subset, 
  #           aes(x = X, y = Y, label = cell_id), size = 2, color = "lightgrey") +
  scale_fill_manual(name = "Land use",values = c("grey80","darkkhaki",
                                                 "darkolivegreen",
                                                 "grey20",
                                                 "darkolivegreen3",
                                                 "deepskyblue4")) +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "") +
  theme_bw() +
  theme(legend.justification=c(1,0), legend.position=c(.99,0.01)) +
  NULL
cell.land.use.plot

cell.plots <- cell.comm.plot + cell.land.use.plot

# -------------------------------------------------------------------------
# ggsave("results/images/population_communicability_distribution.png",
#        plot = pop.avg.rank.plot, width = 8, height = 5)
# 
# ggsave("results/images/population_communicability_by_guild_status.png",
#        plot = pop.plot.means, width = 8, height = 5)
# 
# 
# # -------------------------------------------------------------------------
# 
# ggsave("results/images/sp_communicability_distribution.png",
#        plot = aggr.rank.plot, width = 8, height = 5)
# 
# ggsave("results/images/combined_communicability_distribution.png",
#        plot = combined.rank.plot, width = 13, height = 6)
# 
# ggsave("results/images/sp_communicability_density.png",
#        plot = sp.hist, width = 8, height = 5)
# 
# # -------------------------------------------------------------------------
# 
# ggsave(paste("results/images/cell_communicability_",grid.size,"km.png",sep=""),
#        plot = cell.comm.plot, width = 5, height = 5)
# 
# ggsave(paste("results/images/cell_communicability_land_use_",grid.size,"km.png",sep=""),
#        plot = cell.plots, width = 9, height = 6)

# ggsave(paste("results/images/cell_richness_",grid.size,"km.png",sep=""),
#        plot = richness.plot, width = 8, height = 8)
