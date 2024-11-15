# tidy modeled co-occurrences from Otso Ovaskainen JSDM 

# INPUTS
# raw input from the JSDM
# - "results/plant_bird_predictions.RData"

# OUTPUTS
# - "results/model_occurrences_"

library(tidyverse)
library(sf)
library(patchwork)

# -------------------------------------------------------------------------

grid.size <- 10
occurrence.threshold <- 0.5

NZ.grid <- st_read(paste("data/NZ_grid_",grid.size,"km.shp",sep=""))

load("results/plant_bird_predictions.RData")

occurrence.prob.df <- as.data.frame(EpredY.xy)
occurrence.prob.df$x <- NULL
occurrence.prob.df$y <- NULL

occurrence.prob.df <- occurrence.prob.df %>%
  mutate(cell_id = rownames(occurrence.prob.df)) %>%
  pivot_longer(cols = c(-cell_id),names_to = "species",values_to = "prob")

presence.df <- occurrence.prob.df
presence.df$presence <- ifelse(presence.df$prob > occurrence.threshold,1,0)
presence.df$prob <- NULL
presence.df$cell_id <- as.numeric(presence.df$cell_id)

# test species
# occurrence.all.sp <- sort(unique(occurrence.prob.df$species))
# presences.all.sp <- sort(unique(presence.df$species[which(presence.df$presence == 1)]))

write.csv2(presence.df,paste("results/model_occurrences_",grid.size,"km.csv",sep=""),row.names = F)

# -------------------------------------------------------------------------
# more info: how many species per cell?

sp.per.cell <- presence.df %>% 
  filter(presence == 1) %>% 
  group_by(cell_id) %>%
  summarise(richness = n()) 

avg.values <- sp.per.cell %>% summarise(average.richness = mean(richness),
                                        median.richness = median(richness),
                                        sd.richness = sd(richness))

richness.plot <- ggplot(sp.per.cell) + 
  geom_histogram(aes(x = richness)) + 
  geom_vline(xintercept = avg.values$average.richness, linetype = "dashed", color = "darkred") +
  geom_vline(xintercept = avg.values$median.richness, linetype = "dashed", color = "darkorange") +
  theme_bw() +
  xlab("species richness in 10km grid cell") +
  ylab("number of grid cells") +
  NULL

ggsave("results/images/richness_per_cell.png",
       plot = richness.plot, width = 5, height = 5)
