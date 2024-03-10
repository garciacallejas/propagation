
# tidy species and interactions data
# IMPORTANT NOTE:
# this script will not work. It needs data from NVS and Tier1 which we cannot 
# share in its entirety, as we are not the owners. Therefore, it is here
# for reference, and the processed results are available in the "results" folder.

# NOTE 2: this script also integrates GBIF data that we did not use in the
# final analyses. It is left here for reference.

# INPUTS
# - species observations from different sources (GBIF, NVS, Tier1)

# OUTPUTS
# - site x observations dataframe

# -------------------------------------------------------------------------
library(sf)
library(stars)
# library(raster)
library(terra)
library(tidyverse)
# remotes::install_github("hagc/rasterB")
# library(rasterB)

# -------------------------------------------------------------------------
# interaction list
int.data <- read.csv2("results/plant_bird_interactions_clean.csv")
sp.list <- read.csv2("data/species_list.csv")

# -------------------------------------------------------------------------
# actual data
# NOTE: this is the data that we cannot share
nvs.obs <- read.csv("../datasets/DGCSpeciesDistribution_Approved29NOV22/DGCSpeciesDistribution_Approved29NOV22.csv")
nvs.otago <- read.csv("../datasets/DGCSpeciesDistribution_Approved29NOV22/DGCSpeciesDistribution_OtagoPeninsula30_NOV22.csv")
# nvs.all <- rbind(nvs.obs,nvs.otago)

gbif.obs <- read.csv2("results/species_observations_GBIF_2000.csv")

# NOTE: also this is not available
birds.tier1 <- read.csv2("data/birds_tier1.csv")

# -------------------------------------------------------------------------
# to test concordance with the grid
# this should be already in WGS84
NZ <- st_read('../datasets/spatial_data/NZ_main_islands.shp')

# -------------------------------------------------------------------------
# cell size in meters
grid.size <- 10000

# limit year for NVS observations
nvs.year <- 2000

# where are environmental rasters
env.path <- "../datasets/NZEnvDS_NZMG/"

# NZ grid
NZ_grid <- st_read(paste("data/NZ_grid_",grid.size/1e3,"km.shp",sep=""))

# -------------------------------------------------------------------------
# add cell ids to species observations

# first, set the crs of the original data
projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# -------------------------------------------------------------------------
# NVS data
# -------------------------------------------------------------------------
# in nvs, the colum PreferredSpeciesName can be discarded afaik, 
# gbif accepts the names I had
# td <- subset(nvs.obs, DGCSpeciesName != PreferredSpeciesName)

nvs.obs.full <- rbind(nvs.obs,nvs.otago)

names(nvs.obs.full)[1] <- "species"
names(nvs.obs.full)[which(names(nvs.obs.full) == "PlotEastingWG")] <- "decimalLongitude"
names(nvs.obs.full)[which(names(nvs.obs.full) == "PlotNorthingWG")] <- "decimalLatitude"
names(nvs.obs.full)[which(names(nvs.obs.full) == "ProjectStartDate")] <- "startyear"
names(nvs.obs.full)[which(names(nvs.obs.full) == "ProjectStopDate")] <- "endyear"

# species
# test
nvs.sp <- sort(unique(nvs.obs.full$species))
nvs.nolist <- nvs.sp[which(!(nvs.sp %in% sp.list))]

nvs.obs.full$species <- str_replace(nvs.obs.full$species," ","_")

# year of nvs observations is the end year of the project unless it is empty
# in which case the start year of the project is assigned
nvs.obs.full$endyear[which(nvs.obs.full$endyear == "")] <- nvs.obs.full$startyear[which(nvs.obs.full$endyear == "")]
nvs.obs.full$year <- substr(nvs.obs.full$endyear,1,4)
nvs.obs.full$year <- as.numeric(nvs.obs.full$year)

# keep the useful fields
nvs.obs.clean <- nvs.obs.full[,c("species","decimalLongitude","decimalLatitude","year")]

# -------------------------------------------------------------------------
# Filter by species. And by year?
nvs.obs.clean <- subset(nvs.obs.clean, year >= nvs.year & species %in% sp.list$species)

# add crs to the data
nvs.obs.grid <- st_as_sf(x = nvs.obs.clean,
                        coords = c("decimalLongitude", "decimalLatitude"),
                        crs = projcrs)

# then, transform to NZTM2000
nvs.obs.nz <- st_transform(nvs.obs.grid, crs = st_crs(2193))
nvs.observations.id <- st_intersection(nvs.obs.nz,NZ_grid)

# convert to standard dataframe
# and aggregate observations by cell
# obs_id_df <- observations_id[,c("cell_id","species","year","geometry",
#                                 "gbifID","institutionCode")] %>%
# nvs.obs.df <- nvs.observations.id[,c("cell_id","species","year")] %>%
#   mutate(lat = sf::st_coordinates(.)[,2],
#          lon = sf::st_coordinates(.)[,1],
#          crs = "NZTM2000") %>%
#   sf::st_set_geometry(NULL) %>%
#   group_by(species,cell_id) %>%
#   summarise(observations = n()) %>%
#   dplyr::select(cell_id, species, observations) %>%
#   mutate(dataset = "NVS")

nvs.cell.years <- nvs.observations.id %>%
  sf::st_set_geometry(NULL) %>%
  group_by(cell_id) %>%
  summarise(years_sampled = length(unique(year))) %>%
  mutate(dataset = "NVS")

nvs.obs.df <- nvs.observations.id[,c("cell_id","species","year")] %>%
  mutate(lat = sf::st_coordinates(.)[,2],
         lon = sf::st_coordinates(.)[,1],
         crs = "NZTM2000") %>%
  sf::st_set_geometry(NULL) %>%
  group_by(species,cell_id) %>%
  # group_by(species,cell_id,year) %>%
  summarise(observations = n()) %>%
  dplyr::select(cell_id, species, observations) %>%
  # dplyr::select(cell_id, species, year, observations) %>%
  # group_by(cell_id,species) %>%
  # summarise(observations = sum(observations),years_sampled = length(unique(year))) %>%
  mutate(dataset = "NVS")

nvs.num.sp <- nvs.obs.df %>%
  group_by(cell_id) %>%
  summarise(num_sp_NVS = length(unique(species)))

# -------------------------------------------------------------------------
# GBIF DATA
# -------------------------------------------------------------------------

# filter by institution code
gbif.institutions.to.remove <- c("","naturgucker","iNaturalist",
                                 "MO","K","F","E","S","W","Anymals.org")
gbif.obs.inst <- subset(gbif.obs, !(institutionCode %in% gbif.institutions.to.remove))

# filter by species
gbif.obs.clean <- subset(gbif.obs.inst, species %in% sp.list$species)

# add crs to the data
gbif.obs.grid <- st_as_sf(x = gbif.obs.clean,
                        coords = c("decimalLongitude", "decimalLatitude"),
                        crs = projcrs)

# then, transform to NZTM2000
gbif.obs.nz <- st_transform(gbif.obs.grid, crs = st_crs(2193))
gbif.observations.id <- st_intersection(gbif.obs.nz,NZ_grid)

gbif.cell.years <- gbif.observations.id %>%
  sf::st_set_geometry(NULL) %>%
  group_by(cell_id) %>%
  summarise(years_sampled = length(unique(year))) %>%
  mutate(dataset = "GBIF")

# convert to standard dataframe
# and aggregate observations by cell
# obs_id_df <- observations_id[,c("cell_id","species","year","geometry",
#                                 "gbifID","institutionCode")] %>%
gbif.obs.df <- gbif.observations.id[,c("cell_id","species","year")] %>%
  mutate(lat = sf::st_coordinates(.)[,2],
         lon = sf::st_coordinates(.)[,1],
         crs = "NZTM2000") %>%
  sf::st_set_geometry(NULL) %>%
  group_by(species,cell_id) %>%
  # group_by(species,cell_id,year) %>%
  summarise(observations = n()) %>%
  dplyr::select(cell_id, species, observations) %>%
  # dplyr::select(cell_id, species, year, observations) %>%
  # group_by(cell_id,species) %>%
  # summarise(observations = sum(observations),years_sampled = length(unique(year))) %>%
  # group_by(species,cell_id) %>%
  # summarise(observations = n()) %>%
  # dplyr::select(cell_id, species, observations) %>%
  mutate(dataset = "GBIF")

gbif.num.sp <- gbif.obs.df %>%
  group_by(cell_id) %>%
  summarise(num_sp_GBIF = length(unique(species)))

# -------------------------------------------------------------------------
# TIER1 data
# -------------------------------------------------------------------------

# I manually checked that most bird sp from tier1 that have interactions
# (i.e. are in sp.list), are correctly identified. The birds that do not appear
# in sp.list are species that have no recorded frugivory interactions.

# these are the only exceptions
birds.tier1$species[birds.tier1$species == "Callaeas wilsoni"] <- "Callaeas cinereus"
birds.tier1$species[birds.tier1$species == "Carduelis chloris"] <- "Chloris chloris"
birds.tier1$species[birds.tier1$species == "Carduelis flammea"] <- "Acanthis flammea"
birds.tier1$species[birds.tier1$species == "Philesturnus rufusater"] <- "Philesturnus carunculatus"

birds.tier1$species <- str_replace(birds.tier1$species," ","_")

# keep only frugivores
birds.tier1.2 <- subset(birds.tier1, species %in% sp.list$species)
# tier1.sp <- sort(unique(birds.tier1.2$species))
tier1.clean <- subset(birds.tier1.2, !is.na(x) | !is.na(y))

tier1.obs.grid <- st_as_sf(x = tier1.clean,
                        coords = c("x", "y"),
                        crs = 2193)
tier1.observations.id <- st_intersection(tier1.obs.grid,NZ_grid)

tier1.cell.years <- tier1.observations.id %>%
  sf::st_set_geometry(NULL) %>%
  group_by(cell_id) %>%
  summarise(years_sampled = length(unique(year))) %>%
  mutate(dataset = "TIER1")

tier1.obs.df <- tier1.observations.id[,c("cell_id","species","year")] %>%
  mutate(lat = sf::st_coordinates(.)[,2],
         lon = sf::st_coordinates(.)[,1],
         crs = "NZTM2000") %>%
  sf::st_set_geometry(NULL) %>%
  group_by(species,cell_id) %>%
  # group_by(species,cell_id,year) %>%
  summarise(observations = n()) %>%
  dplyr::select(cell_id, species, observations) %>%
  # dplyr::select(cell_id, species, year, observations) %>%
  # group_by(cell_id,species) %>%
  # summarise(observations = sum(observations),years_sampled = length(unique(year))) %>%
  # group_by(species,cell_id) %>%
  # summarise(observations = n()) %>%
  # dplyr::select(cell_id, species, observations) %>%
  mutate(dataset = "TIER1")

tier1.num.sp <- tier1.obs.df %>%
  group_by(cell_id) %>%
  summarise(num_sp_tier1DOC = length(unique(species)))

# -------------------------------------------------------------------------
# combine nvs, gbif, and tier1 observations

years.sampled.obs <- bind_rows(nvs.cell.years,gbif.cell.years,tier1.cell.years)

dataset.obs <- bind_rows(nvs.obs.df,gbif.obs.df,tier1.obs.df)

all.obs <- bind_rows(nvs.obs.df, gbif.obs.df,tier1.obs.df) %>% 
  group_by(cell_id,species) %>%
  summarise(obs = sum(observations)) 
names(all.obs)[3] <- "observations"

all.obs.wide <- all.obs %>%
  pivot_wider(id_cols = cell_id,
              names_from = species,
              values_from = observations, values_fill = 0) %>%
  arrange(cell_id)

# -------------------------------------------------------------------------
# write to disk

write.csv2(years.sampled.obs, paste("data/years_sampled_dataset_long_",grid.size/1e3,"km.csv",sep=""),row.names = F)
write.csv2(dataset.obs, paste("data/sp_observations_dataset_long_",grid.size/1e3,"km.csv",sep=""),row.names = F)
write.csv2(all.obs, paste("data/sp_observations_long_",grid.size/1e3,"km.csv",sep=""),row.names = F)
write.csv2(all.obs.wide, paste("data/sp_observations_",grid.size/1e3,"km.csv",sep=""),row.names = F)

# -------------------------------------------------------------------------

# put together richness values per cell according to different factors

tier1.obs.df.traits <- left_join(tier1.obs.df,sp.list)
tier1.obs.df.traits$record.type <- "TIER1"

gbif.obs.df.traits <- left_join(gbif.obs.df,sp.list)
gbif.obs.df.traits$record.type <- "GBIF"

nvs.obs.df.traits <- left_join(nvs.obs.df,sp.list)
nvs.obs.df.traits$record.type <- "NVS"

obs.traits <- rbind(tier1.obs.df.traits,gbif.obs.df.traits,nvs.obs.df.traits)
obs.traits$dataset <- NULL

# 1 - observation type (GBIF/Tier1/NVS) and guild
obs.record.type <- obs.traits %>%
  group_by(cell_id,guild,record.type) %>%
  summarise(richness = length(unique(species)))

# 2 - status (native/exotic) and guild
obs.status <- obs.traits %>%
  group_by(cell_id,guild,status) %>%
  summarise(richness = length(unique(species)))

write.csv2(obs.traits,paste("data/sp_observations_full_long_",grid.size/1e3,"km.csv",sep=""),row.names = F)
write.csv2(obs.record.type,paste("data/sp_observations_record_type_long_",grid.size/1e3,"km.csv",sep=""),row.names = F)
write.csv2(obs.status,paste("data/sp_observations_status_long_",grid.size/1e3,"km.csv",sep=""),row.names = F)

# -------------------------------------------------------------------------
# update metaweb

updated.sp.list <- sort(unique(all.obs$species))
updated.interactions <- subset(int.data,(PLANTSPECIES %in% updated.sp.list | 
                                 BIRDSPECIES %in% updated.sp.list))

write.csv2(updated.interactions, 
           "data/interactions_occurring_species.csv",
           row.names = F)

# -------------------------------------------------------------------------
# test changed resolution

# 
# # add x-y coords of the centroid
# grid_with_labels <- st_centroid(grid) %>% cbind(st_coordinates(.))
# 
# # -------------------------------------------------------------------------
# ggplot() +
#   geom_sf(data = NZ2, fill = 'white', lwd = 0.05) +
#   # geom_sf(data = env.sf,aes(color = humidity_meanAnn)) +
#   geom_sf(data = grid, fill = 'transparent', lwd = 0.3) +
#   # geom_sf(data = env.grid, aes(fill = humidity_meanAnn), lwd = 0.3) +
# 
#   # geom_text(data = grid_with_labels,
#   #           aes(x = X, y = Y, label = cell_id),
#   #           size = 2) +
# 
#   coord_sf(datum = NA)  +
#   labs(x = "") +
#   labs(y = "")
