
# plot variability of communicability metric 
# depending on buffer from the focal cell

# -------------------------------------------------------------------------

library(tidyverse)
library(colorblindr)
library(ggrepel)
# -------------------------------------------------------------------------

grid.size <- 10

d.list <- list()

d.list[[1]] <- read.csv2("results/buffer_0_focal_1_10km.csv")
d.list[[2]] <- read.csv2("results/buffer_10_focal_1_10km.csv")
d.list[[3]] <- read.csv2("results/buffer_20_focal_1_10km.csv")
d.list[[4]] <- read.csv2("results/buffer_30_focal_1_10km.csv")
d.list[[5]] <- read.csv2("results/buffer_40_focal_1_10km.csv")

dist.comm <- bind_rows(d.list)

sp.traits <- read.csv2("data/trait_data.csv")
sp.status <- unique(sp.traits[,c("species","guild","status")])

# -------------------------------------------------------------------------

dist.comm$species <- factor(dist.comm$sp1, levels = sort(unique(dist.comm$sp1)))
dist.comm.rank <- dist.comm %>%
  group_by(buffer.distance) %>%
  mutate(scaled.weighted.communicability = scales::rescale(population.weighted.communicability)) %>%
  group_by(buffer.distance) %>%
  dplyr::arrange(desc(scaled.weighted.communicability)) %>%
  group_by(buffer.distance) %>%
  mutate(species.rank = rank(-scaled.weighted.communicability,ties.method = "first")) %>%
  filter(scaled.weighted.communicability > 0) %>%
  left_join(sp.status)

dist.comm.rank$sp.labels <- paste("italic('",gsub("_"," ",dist.comm.rank$species),"')",sep="")

dist.comm.rank <- dist.comm.rank %>%
  select(species,guild,buffer.distance,scaled.weighted.communicability,species.rank,sp.labels)

# -------------------------------------------------------------------------

rank.distances.plot <- ggplot(data = dist.comm.rank, 
                              aes(x = species.rank, 
                                  y = scaled.weighted.communicability), 
                              group = buffer.distance) + 
  geom_point(aes(fill = guild), shape = 21, size = 2) + 
  geom_line() +
  scale_fill_OkabeIto() +
  theme_bw() +
  scale_x_discrete(breaks=NULL) +
  geom_label_repel(data = subset(dist.comm.rank, species.rank < 6),
                   aes(label=sp.labels),
                   alpha = .8,
                   nudge_x = 50,
                   parse = T) +
  facet_grid(.~buffer.distance) +
  xlab("species rank") + ylab("population communicability") +
  NULL
# rank.distances.plot

# -------------------------------------------------------------------------

rank.sp.subset <- subset(dist.comm.rank, species.rank < 9)

rank.variation.plot <- ggplot(rank.sp.subset, aes(x = buffer.distance, y = species.rank, group = species)) +
  geom_line() +
  geom_point(aes(fill = species), shape = 21, size = 3) +
  scale_y_reverse() +
  labs(x = "maximum distance from focal cell (km)", y = "species communicability rank") +
  NULL
# rank.variation.plot

# -------------------------------------------------------------------------
ggsave("results/images/buffer_communicability_distributions.png",
       plot = rank.distances.plot, width = 16, height = 5)

ggsave("results/images/buffer_communicability_variations.png",
       plot = rank.variation.plot, width = 7, height = 5)
