
# plot observations alongside model occurrences for every sp

library(tidyverse)
library(sf)
library(patchwork)

# -------------------------------------------------------------------------
grid.size <- 10
obs.traits <- read.csv2(paste("data/sp_observations_full_long_",grid.size,"km.csv",
                              sep=""))
obs.traits$cell_id <- formatC(obs.traits$cell_id, width = 5, format = "d", flag = "0")

sp.obs <- read.csv2(paste("data/sp_observations_long_",grid.size,"km.csv",
                            sep=""))
NZ.grid <- st_read(paste("data/NZ_grid_",grid.size,"km.shp",sep=""))

load("results/plant_bird_predictions.RData")
# dim(EpredY.xy)

occurrence.prob.df <- as.data.frame(EpredY.xy)
occurrence.prob.df$x <- NULL
occurrence.prob.df$y <- NULL

occurrence.prob.df <- occurrence.prob.df %>%
  mutate(cell_id = rownames(occurrence.prob.df)) %>%
  pivot_longer(cols = c(-cell_id),names_to = "species",values_to = "prob")

NZ.grid$cell_id <- formatC(NZ.grid$cell_id, width = 5, format = "d", flag = "0")

# -------------------------------------------------------------------------
# plot observations alongside occurrence probabilities for every sp

all.sp <- sort(unique(occurrence.prob.df$species))

for(i.sp in 1:length(all.sp)){
  my.sp <- subset(occurrence.prob.df,species == all.sp[i.sp])
  
  my.grid <- left_join(NZ.grid,my.sp[,c("cell_id","prob")])
  
  my.plot <- ggplot(data = my.grid) + 
    geom_sf(aes(fill = prob)) + 
    scale_fill_viridis_c() +
    coord_sf(datum = NA)  +
    labs(x = "") +
    labs(y = "") +
    theme_bw() +
    ggtitle(paste(all.sp[i.sp]," - occurrence probability")) +
    NULL
  
  my.traits <- subset(obs.traits,species == all.sp[i.sp])
  my.record.types <- sort(unique(my.traits$record.type))
  my.guild <- sort(unique(my.traits$guild))
  my.status <- sort(unique(my.traits$status))
  
  my.record.type.plots <- list()
  for(i.record.type in 1:length(my.record.types)){
    my.type.data <- subset(my.traits, record.type == my.record.types[i.record.type])
    my.type.obs <- left_join(NZ.grid,my.type.data[,c("cell_id","observations")])
    
    my.record.type.plots[[i.record.type]] <- ggplot(data = my.type.obs) + 
      geom_sf(aes(fill = observations)) + 
      scale_fill_viridis_c() +
      coord_sf(datum = NA)  +
      labs(x = "") +
      labs(y = "") +
      theme_bw() +
      ggtitle(paste(all.sp[i.sp],"-",my.record.types[i.record.type])) +
      NULL
  }# for record type
  
  my.record.type.plots[[length(my.record.type.plots)+1]] <- my.plot
  
  my.full.plot <- wrap_plots(my.record.type.plots)
  
  ggsave(paste("results/images/species_model_occurrences/",all.sp[i.sp],".pdf",sep=""),
         plot = my.full.plot,
         device = cairo_pdf,
         width = 15, height = 10,dpi = 300)
}





