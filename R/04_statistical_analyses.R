
# Relate species communicability with traits
# -------------------------------------------------------------------------

library(tidyverse)
library(betareg)
library(DHARMa)
library(performance)

# -------------------------------------------------------------------------
grid.size <- 10

sp.comm <- read.csv2(paste("results/species_level_communicability_",grid.size,"km.csv",sep=""))

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

# -------------------------------------------------------------------------

# merge data and maintain only the communicability metric that we use as response
sp.data <- left_join(sp.traits.wide,sp.comm) %>%
  rename(comm = scaled.weighted.communicability) %>%
  subset(!is.na(status))

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

# -------------------------------------------------------------------------

# boxplot
boxplot.status.guild <- ggplot(sp.data,aes(x = guild, y = comm)) + 
  geom_boxplot(aes(fill = status)) + 
  NULL
boxplot.status.guild

# -------------------------------------------------------------------------

# cor(plant.data$FRUIT_DIAMETER_mm,plant.data$FRUIT_LENGTH_mm) # ~0.85
# cor(plant.data$FRUIT_DIAMETER_mm,plant.data$PLANT_MAX_MEAN_VEGETATIVE_HEIGHT_m) # 0.04
plant.glm <- betareg(comm ~ status + FRUIT_DIAMETER_mm + PLANT_MAX_MEAN_VEGETATIVE_HEIGHT_m,data = plant.data)
summary(plant.glm)
# DHARMa::testResiduals(DHARMa::simulateResiduals(plant.glm))
performance::check_model(plant.glm)

# -------------------------------------------------------------------------
# cor(bird.data$BILL_LENGTH_mm,bird.data$BILL_WIDTH_mm) # 0.85
# cor(bird.data$BILL_LENGTH_mm,bird.data$body.mass) # 0.66
# cor(bird.data$body.mass,bird.data$HWI) # 0.18

bird.glm <- betareg(comm ~ status + BILL_LENGTH_mm + body.mass + HWI,data = bird.data)
summary(bird.glm)

performance::check_model(bird.glm)
