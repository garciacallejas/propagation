
library(tidyverse)
# library(sf)
library(patchwork)
library(ggrepel)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# -------------------------------------------------------------------------

ind.com <- read.csv2("results/individual_aggregated_communicability.csv")

ind.ag <- ind.com %>%
  group_by(species) %>%
  summarise(agg.com = sum(aggregated.communicability)) %>%
  mutate(scaled.aggregated.communicability = range01(agg.com)) 
ind.ag$guild <- ind.com$guild[match(ind.ag$species,ind.com$species)]
ind.ag <- dplyr::arrange(ind.ag,desc(agg.com))
ind.ag$species <- factor(ind.ag$species, levels = unique(ind.ag$species))

cell.ag <- ind.com %>%
  group_by(cell) %>%
  summarise(agg.com = sum(aggregated.communicability)) %>%
  mutate(scaled.aggregated.communicability = range01(agg.com)) 
cell.ag <- dplyr::arrange(cell.ag,desc(agg.com))
cell.ag$cell <- factor(cell.ag$cell, levels = unique(cell.ag$cell))

# -------------------------------------------------------------------------

sp.plot <- ggplot(data = ind.ag, aes(x = species, 
                                     y = scaled.aggregated.communicability)) + 
  geom_point(aes(fill = guild), shape = 21) + 
  theme_bw() +
  scale_x_discrete(breaks=NULL) +
  geom_label_repel(data = subset(ind.ag, scaled.aggregated.communicability > 0.4),
            aes(label=species),
            nudge_x = 50) +
  xlab("") +
  NULL

sp.hist <- ggplot(data = ind.ag, aes(x = scaled.aggregated.communicability)) + 
  geom_density(aes(fill = guild), alpha = 0.5) + 
  theme_bw() +
  # scale_x_discrete(breaks=NULL) +
  NULL

# -------------------------------------------------------------------------

cell.plot <- ggplot(data = cell.ag, aes(x = cell, 
                                       y = scaled.aggregated.communicability)) + 
  geom_point() + 
  theme_bw() +
  # scale_x_discrete(breaks=NULL) +
  NULL

cell.hist <- ggplot(data = cell.ag, aes(x = scaled.aggregated.communicability)) + 
  geom_density() + 
  theme_bw() +
  # scale_x_discrete(breaks=NULL) +
  NULL

# -------------------------------------------------------------------------

ggsave("results/images/sp_communicability_distribution.png",
       plot = sp.plot, width = 8, height = 5)

ggsave("results/images/sp_communicability_density.png",
       plot = sp.hist, width = 8, height = 5)


ggsave("results/images/cell_communicability_distribution.png",
       plot = cell.plot, width = 12, height = 5)

ggsave("results/images/cell_communicability_density.png",
       plot = cell.hist, width = 8, height = 5)
