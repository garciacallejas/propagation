
# test species presences in otso's occurrences,
# this is because there seem to be too few birds in "bird.data" 
# from the script 04_statistical_analyses

library(tidyverse)

# -------------------------------------------------------------------------
grid.size <- 10
occurrence.threshold <- 0.5

# -------------------------------------------------------------------------
# read and fill trait data

sp.traits <- read.csv2("data/trait_data.csv")
sp.traits.2 <- sp.traits
# fill by hand some body mass values that are not present in the gathered dataset
# source: https://en.wikipedia.org/wiki/North_Island_saddleback
sp.traits.2$mean.value[which(sp.traits.2$species == "Philesturnus_rufusater" & 
                               sp.traits.2$trait == "body.mass")] <- 75 
# source: https://en.wikipedia.org/wiki/K%C4%81k%C4%81p%C5%8D
sp.traits.2$mean.value[which(sp.traits.2$species == "Strigops_habroptila" & 
                               sp.traits.2$trait == "body.mass")] <- 2000

# source: https://en.wikipedia.org/wiki/Common_myna
sp.traits.2$mean.value[which(sp.traits.2$species == "Acridotheres_tristis" & 
                               sp.traits.2$trait == "body.mass")] <- 110

# source: https://en.wikipedia.org/wiki/Silvereye
sp.traits.2$mean.value[which(sp.traits.2$species == "Zosterops_lateralis" & 
                               sp.traits.2$trait == "body.mass")] <- 10

# source: body mass data, where it appears as "Carduelis chloris"
sp.traits.2$mean.value[which(sp.traits.2$species == "Chloris_chloris" & 
                               sp.traits.2$trait == "body.mass")] <- 27.5

# source: https://en.wikipedia.org/wiki/Common_redpoll
sp.traits.2$mean.value[which(sp.traits.2$species == "Acanthis_flammea" & 
                               sp.traits.2$trait == "body.mass")] <- 14

# source: https://www.birdsnz.org.nz/wp-content/uploads/2021/12/Notornis_45_1_49.pdf
sp.traits.2$mean.value[which(sp.traits.2$species == "Cyanoramphus_malherbi" & 
                               sp.traits.2$trait == "body.mass")] <- 45.3

# source: file:///home/david/Downloads/moorebattleyNotornis2003.pdf
sp.traits.2$mean.value[which(sp.traits.2$species == "Anas_chlorotis" & 
                               sp.traits.2$trait == "body.mass")] <- 500

# measurements of callaeas wilsoni, a sister species
# source: https://www.birdsnz.org.nz/wp-content/uploads/2021/12/Notornis_48_4_217.pdf
sp.traits.2$mean.value[which(sp.traits.2$species == "Callaeas_cinereus" & 
                               sp.traits.2$trait == "body.mass")] <- 225


sp.traits.2$sd.value <- sp.traits.2$n <- NULL
sp.traits.2$trait[which(sp.traits.2$trait == "Hand-wing.Index")] <- "HWI"
sp.traits.wide <- pivot_wider(sp.traits.2,names_from = trait,values_from = mean.value)
sp.list <- sort(unique(sp.traits$species))
sp.status <- sp.traits.wide %>% select(species,guild,status)
# -------------------------------------------------------------------------

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

occurrence.species <- data.frame(data.source = "model_occurrence_prob", 
                                 species = sort(unique(occurrence.prob.df$species)))
occurrence.species <- left_join(occurrence.species,sp.status)

presence.species <- data.frame(data.source = "model_presences", 
                               species = sort(unique(presence.df$species[which(presence.df$presence == 1)])))
presence.species <- left_join(presence.species,sp.status)

# -------------------------------------------------------------------------

sp.comm <- read.csv2(paste("results/species_level_communicability_",grid.size,"km.csv",sep=""))

# merge data and maintain only the communicability metric that we use as response
sp.data <- left_join(sp.traits.wide,sp.comm) %>%
  rename(comm = scaled.weighted.communicability) %>%
  subset(!is.na(status))

comm.species <- sp.data %>% 
  filter(!is.na(comm)) %>%
  select(species,guild,status) %>%
  unique() %>%
  mutate(data.source = "communicability")

# -------------------------------------------------------------------------

table(occurrence.species$guild)
table(presence.species$guild)
table(comm.species$guild)

sp.data$comm[which(sp.data$comm == 0)] <- 1e-5
sp.data$comm[which(sp.data$comm == 1)] <- 0.99999

plant.data <- subset(sp.data,guild == "plants") %>%
  select(species,comm,status,
         FRUIT_DIAMETER_mm,
         FRUIT_LENGTH_mm,
         PLANT_MAX_MEAN_VEGETATIVE_HEIGHT_m) %>%
  drop_na()
bird.data <- subset(sp.data, guild == "birds") %>%
  select(species,comm,status,
         # Beak.Depth, Beak.Length_Culmen,
         # Beak.Length_Nares,
         # Beak.Width,
         BILL_LENGTH_mm,
         # BILL_WIDTH_mm,
         body.mass,HWI,
         # Kipps.Distance,
         # Tail.Length,Tarsus.Length,Wing.Length
  ) %>%
  drop_na()
