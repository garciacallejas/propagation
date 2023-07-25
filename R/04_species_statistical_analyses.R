
# Relate species communicability with traits
# -------------------------------------------------------------------------

library(glmmTMB)
library(tidyverse)
library(betareg)
library(DHARMa)
library(performance)
library(corrplot)
library(scales)
library(broom.mixed)

# -------------------------------------------------------------------------
grid.size <- 10

pop.comm <- read.csv2(paste("results/population_level_communicability_",grid.size,"km.csv",sep=""))
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
sp.status <- unique(sp.traits[,c("species","guild","status")])

# which source to use?
observation.source <- "occurrence_prob"

if(observation.source == "observed_records"){
  sp.obs <- read.csv2(paste("data/sp_observations_long_",grid.size,"km.csv",
                            sep=""))
}else if(observation.source == "occurrence_prob"){
  sp.obs <- read.csv2(paste("results/model_occurrences_",grid.size,"km.csv",
                            sep=""))
  # this is to harmonise datasets
  names(sp.obs)[which(names(sp.obs) == "presence")] <- "observations"
  # keep presences only
  sp.obs <- subset(sp.obs, observations == 1)
  
}

all.sp <- sort(unique(pop.comm$species))

# species presences
sp.presence.absences <- expand.grid(cell_id = unique(sp.obs$cell_id),
                                    species = all.sp) %>%
  left_join(sp.obs) %>% 
  replace_na(list(observations = 0)) %>%
  rename(presence = observations) %>%
  arrange(cell_id)

pop.averages <- pop.comm %>% 
  left_join(sp.presence.absences) %>%
  filter(presence == 1) %>%
  group_by(species,guild) %>%
  summarise(avg.comm = mean(population.weighted.communicability),
            sd.comm = sd(population.weighted.communicability)) %>%
  ungroup() %>%
  mutate(scaled.avg.comm = scales::rescale(avg.comm))

# -------------------------------------------------------------------------
# check population-level and species-level communicability

pop.data <- left_join(sp.traits.wide,pop.averages) %>%
  rename(comm = scaled.avg.comm) %>%
  subset(!is.na(status))

plant.pop.data <- subset(pop.data,guild == "plants") %>% 
  select(species,comm,status,
         FRUIT_DIAMETER_mm,
         FRUIT_LENGTH_mm,
         PLANT_MAX_MEAN_VEGETATIVE_HEIGHT_m) %>%
  drop_na()
bird.pop.data <- subset(pop.data, guild == "birds") %>% 
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

# merge data and maintain only the communicability metric that we use as response
sp.data <- left_join(sp.traits.wide,sp.comm) %>%
  rename(comm = scaled.weighted.communicability,comm2 = weighted.communicability) %>%
  subset(!is.na(status))

plant.data <- subset(sp.data,guild == "plants") %>% 
  select(species,comm,comm2,status,
         FRUIT_DIAMETER_mm,
         FRUIT_LENGTH_mm,
         PLANT_MAX_MEAN_VEGETATIVE_HEIGHT_m) %>%
  drop_na()
bird.data <- subset(sp.data, guild == "birds") %>% 
  select(species,comm,comm2,status,
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
# visualizations
# correlations among response variables
all.bird.traits <- sp.data %>% select(Beak.Depth:WING_LENGTH_cm) %>% drop_na()
all.plant.traits <- sp.data %>% 
  select(FRUIT_DIAMETER_mm:PLANT_MAX_MEAN_VEGETATIVE_HEIGHT_m) %>% drop_na()
names(all.plant.traits)[3] <- "PLANT_HEIGHT_m"

bird.traits.cor <- cor(all.bird.traits)
# corrplot(bird.traits.cor)

plant.traits.cor <- cor(all.plant.traits)
# corrplot(plant.traits.cor)

# boxplot - communicability by guild and status
# boxplot.status.guild <- ggplot(sp.data,aes(x = guild, y = comm)) + 
#   geom_boxplot(aes(fill = status)) + 
#   NULL
# boxplot.status.guild

# -------------------------------------------------------------------------
# population-level glmm

plant.pop.glm <- glmmTMB(comm ~ status + scale(FRUIT_DIAMETER_mm) + scale(PLANT_MAX_MEAN_VEGETATIVE_HEIGHT_m),
                        family = tweedie,
                        data = plant.pop.data)

# summary(plant.pop.glm)
# DHARMa::testResiduals(DHARMa::simulateResiduals(plant.pop.glm))

print(xtable::xtable(as.data.frame(broom::tidy(plant.pop.glm)),
                     floating=FALSE,
                     digits = 3,
                     latex.environments=NULL,
                     booktabs=FALSE))

# -------------------------------------------------------------------------
# cor(bird.data$BILL_LENGTH_mm,bird.data$BILL_WIDTH_mm) # 0.85
# cor(bird.data$BILL_LENGTH_mm,bird.data$body.mass) # 0.66
# cor(bird.data$body.mass,bird.data$HWI) # 0.18

# bird.glm <- betareg(comm ~ status + scale(BILL_LENGTH_mm) + scale(body.mass) + scale(HWI),data = bird.data)
bird.pop.glm <- glmmTMB(comm ~ status + scale(body.mass) + 
                         scale(BILL_LENGTH_mm),  
                       # scale(HWI),
                       family = tweedie(link = "log"), data = bird.pop.data)

# summary(bird.pop.glm)
# DHARMa::testResiduals(DHARMa::simulateResiduals(bird.pop.glm))

print(xtable::xtable(as.data.frame(broom::tidy(bird.pop.glm)),
                     floating=FALSE,
                     digits = 3,
                     latex.environments=NULL,
                     booktabs=FALSE))

# -------------------------------------------------------------------------
# species-level glmm

# cor(plant.data$FRUIT_DIAMETER_mm,plant.data$FRUIT_LENGTH_mm) # ~0.85
# cor(plant.data$FRUIT_DIAMETER_mm,plant.data$PLANT_MAX_MEAN_VEGETATIVE_HEIGHT_m) # 0.04
plant.sp.glm <- glmmTMB(comm2 ~ status + scale(FRUIT_DIAMETER_mm) + scale(PLANT_MAX_MEAN_VEGETATIVE_HEIGHT_m),
                      family = tweedie,
                      data = plant.data)

# summary(plant.sp.glm)
# DHARMa::testResiduals(DHARMa::simulateResiduals(plant.sp.glm))

print(xtable::xtable(as.data.frame(broom::tidy(plant.sp.glm)),
                     floating=FALSE,
                     digits = 3,
                     latex.environments=NULL,
                     booktabs=FALSE))

# -------------------------------------------------------------------------
# cor(bird.data$BILL_LENGTH_mm,bird.data$BILL_WIDTH_mm) # 0.85
# cor(bird.data$BILL_LENGTH_mm,bird.data$body.mass) # 0.66
# cor(bird.data$body.mass,bird.data$HWI) # 0.18

# bird.glm <- betareg(comm ~ status + scale(BILL_LENGTH_mm) + scale(body.mass) + scale(HWI),data = bird.data)
bird.sp.glm <- glmmTMB(comm2 ~ status + scale(body.mass) + 
                       scale(BILL_LENGTH_mm),  
                       # scale(HWI),
                     family = tweedie(link = "log"), data = bird.data)

# summary(bird.sp.glm)
# DHARMa::testResiduals(DHARMa::simulateResiduals(bird.sp.glm))

print(xtable::xtable(as.data.frame(broom::tidy(bird.sp.glm)),
                     floating=FALSE,
                     digits = 3,
                     latex.environments=NULL,
                     booktabs=FALSE))


