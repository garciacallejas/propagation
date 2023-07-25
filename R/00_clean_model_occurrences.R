# clean model occurrences from otso 

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



