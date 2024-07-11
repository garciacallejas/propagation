library(tidyverse)
library(sf)
library("rnaturalearth") 
library("rnaturalearthdata")
library(colorblindr)

# -------------------------------------------------------------------------
# this script prepares data for the validation analyses

# -------------------------------------------------------------------------

# cell size in km
grid.size <- 10

# NZ grid
NZ_grid <- st_read(paste("data/NZ_grid_",grid.size,"km.shp",sep=""))
NZ_contour <- st_read('data/NZ_main_islands.shp')

load("data/local_nz_networks/NZ_networks.Rdata")
nz.nets.metadata <- read.csv2("data/local_nz_networks/NZ_networks_metadata.csv") 
nz.nets.metadata <- nz.nets.metadata %>%
  mutate(net_id = paste("NZ",sprintf("%02d",1:nrow(nz.nets.metadata)),sep=""))

# metaweb
sp.int.orig <- read.csv2("results/plant_bird_interactions_clean.csv")
metaweb <- sp.int.orig %>% mutate(metaweb_interaction = 1) %>%
  dplyr::select(PLANTSPECIES,BIRDSPECIES,metaweb_interaction) %>% unique()
names(metaweb) <- c("plant.sp","bird.sp","metaweb_interaction")

# modelled occurrences
load("results/plant_bird_predictions.RData")

occurrence.prob.df <- as.data.frame(EpredY.xy)
occurrence.prob.df$x <- NULL
occurrence.prob.df$y <- NULL

occurrence.prob.df <- occurrence.prob.df %>%
  mutate(cell_id = rownames(occurrence.prob.df)) %>%
  pivot_longer(cols = c(-cell_id),names_to = "species",values_to = "prob")

# species guild/status
sp.traits <- read.csv2("data/trait_data.csv")
sp.status <- unique(sp.traits[,c("species","guild","status")])

# -------------------------------------------------------------------------
# network location in the cell grid - empirical networks

projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
nz.nets.sf <- st_as_sf(x = nz.nets.metadata,
                       coords = c("lon","lat"),
                       crs = projcrs)
nz.nets.sf <- st_transform(nz.nets.sf, crs= st_crs(2193))

nz.nets.cells <- st_intersection(nz.nets.sf,NZ_grid) %>% 
  dplyr::select(net_id,cell_id) %>%
  sf::st_set_geometry(NULL) %>%
  mutate(cell_id = as.character(cell_id))

# -------------------------------------------------------------------------
# rename a couple of sp, as in script NZ03
# and transform the matrices to lists of links
nz.empirical.link.lists <- list()

for(i.net in 1:length(nz.nets)){
  my.net.df <- nz.nets[[i.net]]
  
  names(my.net.df) <- str_replace(names(my.net.df)," ","_")
  rownames(my.net.df) <- str_replace(rownames(my.net.df)," ","_")
  
  # birds
  names(my.net.df)[which(names(my.net.df) == "Mohoua_novaeseelandiae")] <- "Finschia_novaeseelandiae"
  names(my.net.df)[which(names(my.net.df) == "Carduelis_flammea")] <- "Acanthis_flammea"
  
  # plants
  rownames(my.net.df)[which(names(my.net.df) == "Piper_excelsum")] <- "Macropiper_excelsum"
  rownames(my.net.df)[which(names(my.net.df) == "Streblus_heterophyllus")] <- "Paratrophis_microphylla"
  
  nz.nets[[i.net]] <- my.net.df
  
  my.birds <- colnames(my.net.df)
  my.plants <- rownames(my.net.df)
  
  my.interaction.list <- list()
  for(i.bird in 1:length(my.birds)){
    my.interaction.list[[length(my.interaction.list)+1]] <- data.frame(net_id = nz.nets.metadata$net_id[i.net],
                                                                       bird.sp = my.birds[i.bird],
                                                                       plant.sp = my.plants,
                                                                       interaction = my.net.df[,i.bird])
  }
  my.link.list <- bind_rows(my.interaction.list) 
  
  nz.empirical.link.lists[[length(nz.empirical.link.lists)+1]] <- my.link.list
}

nz.empirical.links <- bind_rows(nz.empirical.link.lists)
nz.empirical.links$interaction[nz.empirical.links$interaction != 0] <- 1
nz.empirical.links.wide <- pivot_wider(nz.empirical.links, names_from = net_id, 
                                       values_from = interaction,names_glue = "{net_id}_empirical_interaction")
# names(nz.empirical.links.wide) <- paste(names(nz.empirical.links.wide),"_empirical",sep="")
# -------------------------------------------------------------------------
# list of modelled species 
# small hack to remove leading zero
occurrence.prob.df$cell_id <- as.character(as.numeric(occurrence.prob.df$cell_id))

# modelled species pairs in each cell in which there is an empirical network
occurrence.net.cells <- occurrence.prob.df %>% 
  filter(cell_id %in% nz.nets.cells$cell_id) %>%
  left_join(sp.status[,c("species","guild")])

occurrence.plants <- occurrence.net.cells %>%
  filter(guild == "plants") %>% dplyr::select(cell_id,species,prob)
names(occurrence.plants) <- c("cell_id","plant.sp","plant.sp.prob")

occurrence.birds <- occurrence.net.cells %>%
  filter(guild == "birds") %>% dplyr::select(cell_id,species,prob)
names(occurrence.birds) <- c("cell_id","bird.sp","bird.sp.prob")

occurrence.interactions.list <- list()
  for(i.bird in 1:nrow(occurrence.birds)){
    my.cell.plants <- subset(occurrence.plants, cell_id == occurrence.birds$cell_id[i.bird])
    
    my.bird.interactions <- my.cell.plants
    my.bird.interactions$bird.sp <- occurrence.birds$bird.sp[i.bird]
    my.bird.interactions$bird.sp.prob <- occurrence.birds$bird.sp.prob[i.bird]
    
    my.bird.interactions <- my.bird.interactions %>% dplyr::select(cell_id,
                                                                   bird.sp,
                                                                   bird.sp.prob,
                                                                   plant.sp,
                                                                   plant.sp.prob)
    occurrence.interactions.list[[length(occurrence.interactions.list)+1]] <- my.bird.interactions
  }

occurrence.interactions.df <- bind_rows(occurrence.interactions.list)

# add metaweb value (1/0) to these pairs
occurrence.interactions.df <- occurrence.interactions.df %>% 
  left_join(metaweb) %>% replace_na(list(metaweb_interaction = 0))

# -------------------------------------------------------------------------
# join empirical and modelled networks
# in a grid cell there may be more than one empirical networks,
# hence the relationship many-to-many
occurrence.interactions.wide <- nz.nets.cells %>% 
  left_join(occurrence.interactions.df, relationship = "many-to-many") %>%
  dplyr::select(-cell_id) %>%
  pivot_wider(names_from = net_id,
              values_from = c(bird.sp.prob,plant.sp.prob),
              names_glue = "{net_id}_{.value}",
              names_sort = T,
              names_vary = "slowest") 

names.order <- data.frame(netname = sort(rep(nz.nets.cells$net_id,3)),
                          intname = rep(c("bird.sp.prob","plant.sp.prob","empirical_interaction"),length(nz.nets.cells))) %>%
  mutate(colname = paste(netname,intname,sep="_"))

full_interaction_df <- full_join(occurrence.interactions.wide,nz.empirical.links.wide) 
full_interaction_df <- full_interaction_df %>% dplyr::select(bird.sp,plant.sp,metaweb_interaction,names.order$colname)

# -------------------------------------------------------------------------

# write.csv2(nz.empirical.links,"data/local_nz_networks/empirical_networks_link_list.csv",row.names = F)
# write.csv2(occurrence.interactions.df,"data/local_nz_networks/modelled_networks_link.list.csv",row.names = F)
# write.csv2(nz.nets.cells,"data/net_id_grid_cells.csv",row.names = F)
# save(full_interaction_df,file = "data/local_nz_networks/local_networks_empirical_modelled_NZ.Rdata")

# -------------------------------------------------------------------------

# list of species per cell - empirical networks

sp.cell.list <- list()
for(i.net in 1:length(nz.nets)){
  my.net.df <- nz.nets[[i.net]]
  my.plants <- data.frame(species = rownames(my.net.df),guild = "plants")
  my.birds <- data.frame(species = colnames(my.net.df),guild = "birds")
  my.df <- data.frame(cell_id = nz.nets.cells$cell_id[i.net],
                      net_id = nz.nets.cells$net_id[i.net],
                      species = c(my.plants$species,my.birds$species),
                      guild = c(my.plants$guild,my.birds$guild))
  sp.cell.list[[length(sp.cell.list)+1]] <- my.df
}
sp.cell.df <- bind_rows(sp.cell.list)

# -------------------------------------------------------------------------
# plot the spatial location of the 16 empirical networks

nz.nets.metadata$location <- NA
nz.nets.metadata$location[c(1:5,10)] <- "Wellington area"
nz.nets.metadata$location[6:9] <- "Puhi-puhi river area"
nz.nets.metadata$location[12:14] <- "Nelson area"
nz.nets.metadata$location[11] <- "Windbag Valley"
nz.nets.metadata$location[15] <- "Banks Peninsula"
nz.nets.metadata$location[16] <- "Kowai Bush"

net.location.plot <- ggplot() + 
  geom_sf(data = NZ_contour, color = "grey80",fill = 'grey80', lwd = 0.1) +
  # geom_sf(color = "grey80", fill = "grey80") +
  geom_point(data = nz.nets.metadata, aes(x = lon, y = lat, fill = location), shape = 21, size = 2) +
  labs(x="",y="") +
  scale_color_OkabeIto() +
  theme_bw() +
  coord_sf(datum = NA) +
  NULL
# net.location.plot

# ggsave("results/images/local_nz_networks.png",
#        plot = net.location.plot, width = 6, height = 8)





