
library(pROC)
library(tidyverse)
library(sf)
library(gghalves)
library(ggbeeswarm)

# -------------------------------------------------------------------------
# this script performs the validation of the JSDM

# -------------------------------------------------------------------------
# species list
load("results/plant_bird_predictions.RData")

occurrence.prob.df <- as.data.frame(EpredY.xy)
occurrence.prob.df$x <- NULL
occurrence.prob.df$y <- NULL

occurrence.prob.df <- occurrence.prob.df %>%
  mutate(cell_id = rownames(occurrence.prob.df)) %>%
  pivot_longer(cols = c(-cell_id),names_to = "species",values_to = "prob")
occurrence.prob.df$cell_id <- as.numeric(occurrence.prob.df$cell_id)

presence.df <- occurrence.prob.df
presence.df$presence <- ifelse(presence.df$prob > 0.5,1,0)
presence.df$prob <- NULL
presence.df$cell_id <- as.numeric(presence.df$cell_id)

# metaweb
sp.int.orig <- read.csv2("results/plant_bird_interactions_clean.csv")
metaweb <- sp.int.orig %>% mutate(metaweb_interaction = 1) %>%
  dplyr::select(PLANTSPECIES,BIRDSPECIES,metaweb_interaction) %>% unique()
names(metaweb) <- c("plant.sp","bird.sp","metaweb_interaction")

# plant.sp <- sort(unique(metaweb$plant.sp))
# bird.sp <- sort(unique(metaweb$bird.sp))
# sp.list <- c(plant.sp,bird.sp)

sp.int <- subset(sp.int.orig, PLANTSPECIES %in% presence.df$species & 
                   BIRDSPECIES %in% presence.df$species)

# -------------------------------------------------------------------------
plant.sp <- sort(unique(sp.int$PLANTSPECIES))
bird.sp <- sort(unique(sp.int$BIRDSPECIES))
sp.list <- c(plant.sp,bird.sp)

# -------------------------------------------------------------------------
# read observed occurrences

grid.size <- 10
occurrence.threshold <- 0.5

NZ_grid <- st_read(paste("data/NZ_grid_",grid.size,"km.shp",sep=""))
projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

all.cells <- sort(unique(NZ_grid$cell_id))

nvs.obs <- read.csv("../datasets/DGCSpeciesDistribution_Approved29NOV22.csv")
nvs.otago <- read.csv("../datasets/DGCSpeciesDistribution_OtagoPeninsula30_NOV22.csv")

birds.tier1 <- read.csv2("data/birds_tier1.csv")
birds.tier1$species <- str_replace(birds.tier1$species," ","_")
birds.tier1.clean <- subset(birds.tier1,species %in% bird.sp)
# -------------------------------------------------------------------------
# limit year for NVS observations
nvs.year <- 2000

nvs.obs.full <- rbind(nvs.obs,nvs.otago)

names(nvs.obs.full)[1] <- "species"
names(nvs.obs.full)[which(names(nvs.obs.full) == "PlotEastingWG")] <- "decimalLongitude"
names(nvs.obs.full)[which(names(nvs.obs.full) == "PlotNorthingWG")] <- "decimalLatitude"
names(nvs.obs.full)[which(names(nvs.obs.full) == "ProjectStartDate")] <- "startyear"
names(nvs.obs.full)[which(names(nvs.obs.full) == "ProjectStopDate")] <- "endyear"

# species
# test
nvs.sp <- sort(unique(nvs.obs.full$species))
# nvs.nolist <- nvs.sp[which(!(nvs.sp %in% sp.list))]

nvs.obs.full$species <- str_replace(nvs.obs.full$species," ","_")

# year of nvs observations is the end year of the project unless it is empty
# in which case the start year of the project is assigned
nvs.obs.full$endyear[which(nvs.obs.full$endyear == "")] <- nvs.obs.full$startyear[which(nvs.obs.full$endyear == "")]
nvs.obs.full$year <- substr(nvs.obs.full$endyear,1,4)
nvs.obs.full$year <- as.numeric(nvs.obs.full$year)

# keep the useful fields and subset for the selected time interval
nvs.obs.clean <- nvs.obs.full[,c("species","decimalLongitude","decimalLatitude","year")]
nvs.obs.clean <- subset(nvs.obs.clean, year >= nvs.year & species %in% sp.list)

# add crs to the data
nvs.obs.grid <- st_as_sf(x = nvs.obs.clean,
                         coords = c("decimalLongitude", "decimalLatitude"),
                         crs = projcrs)

# then, transform to NZTM2000
nvs.obs.nz <- st_transform(nvs.obs.grid, crs = st_crs(2193))
nvs.obs.df <- st_intersection(nvs.obs.nz,NZ_grid) %>% 
  sf::st_set_geometry(NULL) %>%
  group_by(species,cell_id) %>%
  # group_by(species,cell_id,year) %>%
  summarise(observations = n()) %>%
  dplyr::select(cell_id, species, observations) %>%
  mutate(guild = "plants", dataset = "NVS")

tier1.clean <- subset(birds.tier1.clean, !is.na(x) | !is.na(y))

tier1.obs.grid <- st_as_sf(x = tier1.clean,
                           coords = c("x", "y"),
                           crs = 2193)
tier1.observations.id <- st_intersection(tier1.obs.grid,NZ_grid)

tier1.obs.df <- tier1.observations.id[,c("cell_id","species","year")] %>%
  mutate(lat = sf::st_coordinates(.)[,2],
         lon = sf::st_coordinates(.)[,1],
         crs = "NZTM2000") %>%
  sf::st_set_geometry(NULL) %>%
  group_by(species,cell_id) %>%
  # group_by(species,cell_id,year) %>%
  summarise(observations = n()) %>%
  dplyr::select(cell_id, species, observations) %>%
  mutate(guild = "birds", dataset = "TIER1")

# -------------------------------------------------------------------------
# modelled probabilities of occurrence

tier1.sites <- sort(unique(tier1.obs.df$cell_id))
nvs.sites <- sort(unique(nvs.obs.df$cell_id))

# generate species-specific modelled and observed dataframes
# i.sp <- 1

sp.included <- character()
sp.data.list <- list()
bird.data.wide.list <- list()
plant.data.wide.list <- list()
for(i.sp in 1:length(sp.list)){
  
  if(sp.list[i.sp] %in% bird.sp){
    my.obs <- subset(tier1.obs.df, species == sp.list[i.sp]) 
    my.sites <- tier1.sites
    my.guild <- "birds"
  
  }else{
    my.obs <- subset(nvs.obs.df, species == sp.list[i.sp]) 
    my.sites <- nvs.sites
    my.guild <- "plants"
  }
  
  # this is for storing observed and modelled data in a wide format
  my.wide.data <- data.frame(cell_id = my.sites) %>% 
    left_join(my.obs) %>%
    replace_na(list(observations = 0)) %>%
    mutate(obs = ifelse(observations > 0,1,0)) %>%
    dplyr::select(cell_id, obs)
  names(my.wide.data)[2] <- paste0(sp.list[i.sp],"__observed")
  
  # this in long format
  my.obs <- my.obs %>% 
    mutate(presence = 1) %>%
    dplyr::select(cell_id, species, presence)
  
  no.obs <- data.frame(cell_id = my.sites[which(!(my.sites %in% my.obs$cell_id))],
                       species = sp.list[i.sp],
                       presence = 0)
  
  my.sp.obs <- bind_rows(my.obs,no.obs) %>% arrange(cell_id)
  names(my.sp.obs) <- c("cell_id","species","observed")
  
  my.mod <- subset(occurrence.prob.df, species == sp.list[i.sp] & cell_id %in% my.sites)
  if(nrow(my.mod)>0){
    
    my.wide.mod <- data.frame(cell_id = my.sites) %>% 
      left_join(my.mod) %>%
      replace_na(list(prob = 0)) %>%
      dplyr::select(cell_id, prob)
    names(my.wide.mod)[2] <- paste0(sp.list[i.sp],"__modelled")
    
    if(sp.list[i.sp] %in% bird.sp){
      bird.data.wide.list[[length(bird.data.wide.list)+1]] <- left_join(my.wide.data,my.wide.mod)
    }else{
      plant.data.wide.list[[length(plant.data.wide.list)+1]] <- left_join(my.wide.data,my.wide.mod)
    }
    
    if(nrow(my.mod)<length(my.sites)){
    no.mod <- data.frame(cell_id = my.sites[which(!(my.sites %in% my.mod$cell_id))],
                         species = sp.list[i.sp],
                         prob = NA)
    my.sp.mod <- bind_rows(my.mod,no.mod) %>% arrange(cell_id)
    names(my.sp.mod) <- c("cell_id","species","modelled")
    my.sp.data <- left_join(my.sp.obs,my.sp.mod) %>% drop_na() %>%
      mutate(guild = my.guild)
    
    sp.data.list[[length(sp.data.list)+1]] <- my.sp.data
    sp.included[length(sp.included)+1] <- sp.list[i.sp]
    }# if needed to add modelled sites with absence
  }# if nrow modelled
}# for i.sp
names(sp.data.list) <- sp.included

bird.wide.data <- bird.data.wide.list %>% reduce(left_join)
plant.wide.data <- plant.data.wide.list %>% reduce(left_join)
# bird.wide.data <- bind_cols(bird.data.wide.list)
# plant.wide.data <- bind_cols(plant.data.wide.list)

# -------------------------------------------------------------------------
# auc analyses

auc.sp <- data.frame(species = sp.included, guild = NA, AUC = NA)

for(i.sp in 1:length(sp.included)){
  my.response <- sp.data.list[[i.sp]]$observed
  my.pred <- sp.data.list[[i.sp]]$modelled
  if(sum(my.response != 0)>0){
    auc.sp$AUC[i.sp] <- auc(predictor= my.pred,
                            response=my.response)
                            # direction="<")
    if(sp.included[i.sp] %in% bird.sp){
      auc.sp$guild[i.sp] <- "birds"
    }else{
      auc.sp$guild[i.sp] <- "plants"
    }
  }# if presences and absences
}# for i.sp
auc.sp.clean <- drop_na(auc.sp)

# double check
auc.birds <- list()
for(i.bird in 1:length(bird.data.wide.list)){
  my.data <- bird.data.wide.list[[i.bird]]
  my.sp <- substr(names(my.data)[2],1,nchar(names(my.data)[2]) - 10)
  my.response <- my.data[,2]
  my.pred <- my.data[,3]
  if(sum(my.response != 0)>0){
    auc.birds[[length(auc.birds)+1]] <- data.frame(guild = "birds", 
                                                   species = my.sp,
                                                   AUC = as.numeric(auc(predictor= my.pred,
                                                             response=my.response,direction="<")))
  }
}
auc.birds.df <- bind_rows(auc.birds)

# double check
auc.plants <- list()
for(i.plant in 1:length(plant.data.wide.list)){
  my.data <- plant.data.wide.list[[i.plant]]
  my.sp <- substr(names(my.data)[2],1,nchar(names(my.data)[2]) - 10)
  my.response <- my.data[,2]
  my.pred <- my.data[,3]
  if(sum(my.response != 0)>0){
    auc.plants[[length(auc.plants)+1]] <- data.frame(guild = "plants", 
                                                   species = my.sp,
                                                   AUC = as.numeric(auc(predictor= my.pred,
                                                                        response=my.response,direction="<")))
  }
}
auc.plants.df <- bind_rows(auc.plants)
auc.double.check <- bind_rows(auc.birds.df,auc.plants.df)

# -------------------------------------------------------------------------
# hist.auc <- ggplot(auc.sp.clean) + 
  geom_histogram(aes(x = AUC), fill = "grey40", color = "grey20",bins = 50) + 
  facet_grid(rows = vars(guild)) + 
  xlim(c(0,1)) +
  theme_bw() +
  ylab("number of species") +
  NULL
 # hist.auc
# -------------------------------------------------------------------------

sp.data.df <- bind_rows(sp.data.list) %>%
  mutate(obs = ifelse(observed == 0, "absence","presence")) 

guild.mod.plot <- ggplot(sp.data.df, aes(x = guild, y = modelled)) + 
  # geom_boxplot(aes(x = guild, y = modelled, fill = obs)) +
  geom_half_point(aes(color = obs),
                  transformation = position_quasirandom(width = 0.1),
                  side = "l", size = 0.5, alpha = 0.5) +
  geom_half_boxplot(aes(fill = obs), side = "r",outlier.size = 0.8) +
  theme_bw() +
  ylab("modelled probability of occurrence") +
  scale_fill_discrete(name = "observation") +
  scale_color_discrete(guide = "none") +
  NULL
# guild.mod.plot

# -------------------------------------------------------------------------
# ggsave("results/images/species_level_AUC.png",
#        plot = hist.auc, width = 5, height = 3)
# 
# ggsave("results/images/observed_modelled_probs.png",
#        plot = guild.mod.plot, width = 6, height = 4)

