# clean up trait datasets

# for birds, Peralta and Avonet, as well as Lislevand for body sizes
# for plants, Peralta for now

library(tidyverse)

# -------------------------------------------------------------------------

# independent list of plants/birds
bird.sp.cells.wide <- read.csv2("data/bird_cell_presences.csv")
plant.sp.cells.wide <- read.csv2("data/plant_cell_presences.csv")

bird.sp <- sort(unique(bird.sp.cells.wide$species))
plant.sp <- sort(unique(plant.sp.cells.wide$species))
all.sp <- sort(unique(c(bird.sp,plant.sp)))
num.sp <- length(all.sp)

# status
status <- read.csv("../datasets/plant-bird interactions and traits/plant_bird_interactions.csv")
plant.status <- unique(status[,c("PLANTSPECIES","PLANTSTATUS")])
bird.status <- unique(status[,c("BIRDSPECIES","BIRDSTATUS")])

# -------------------------------------------------------------------------

# trait datasets
plants <- read.csv("../datasets/plant-bird interactions and traits/plant_traits.csv")
birds.peralta <- read.csv("../datasets/plant-bird interactions and traits/bird_traits.csv")
birds.avonet <- read_csv("../datasets/bird traits/AVONET/AVONET_Raw_Data.csv")
birds.avonet.lifestyle <- read_csv("../datasets/bird traits/AVONET/AVONET1_BirdLife.csv")
birds.body.sizes <- read_tsv("../datasets/bird traits/body_sizes/avian_ssd_jan07.txt")

# NOTE: this is taken from 00_create_spatial_matrix
# in order to have consistent species names
# if I modify this, I should modify both
plants$PLANTSPECIES[which(plants$PLANTSPECIES == "Actinidia_deliciosa")] <- "Actinidia_chinensis"
plants$PLANTSPECIES[which(plants$PLANTSPECIES == "Androstoma_empetrifolium")] <- "Androstoma_empetrifolia"
plants$PLANTSPECIES[which(plants$PLANTSPECIES == "Dendrobenthamia_capitata")] <- "Cornus_capitata"
plants$PLANTSPECIES[which(plants$PLANTSPECIES == "Eriobotrya_japonica")] <- "Rhaphiolepis_loquata"
plants$PLANTSPECIES[which(plants$PLANTSPECIES == "Leucopogon_fraseri")] <- "Styphelia_nesophila"
plants$PLANTSPECIES[which(plants$PLANTSPECIES == "Phormium_cookianum")] <- "Phormium_colensoi"
plants$PLANTSPECIES[which(plants$PLANTSPECIES == "Phyllocladus_alpinus")] <- "Phyllocladus_trichomanoides"
plants$PLANTSPECIES[which(plants$PLANTSPECIES == "Phytolacca_octandra")] <- "Phytolacca_icosandra"
plants$PLANTSPECIES[which(plants$PLANTSPECIES == "Piper_excelsum")] <- "Macropiper_excelsum"
plants$PLANTSPECIES[which(plants$PLANTSPECIES == "Pseudopanax_colensoi")] <- "Neopanax_colensoi"
plants$PLANTSPECIES[which(plants$PLANTSPECIES == "Pseudopanax_colensoi_var_colensoi")] <- "Neopanax_colensoi"
plants$PLANTSPECIES[which(plants$PLANTSPECIES == "Pseudopanax_colensoi_var_ternatus")] <- "Neopanax_colensoi"
plants$PLANTSPECIES[which(plants$PLANTSPECIES == "Solanum_nodiflorum")] <- "Solanum_americanum"
plants$PLANTSPECIES[which(plants$PLANTSPECIES == "Streblus_heterophyllus")] <- "Paratrophis_microphylla"

# maybe more checks needed, e.g. for genus?
plants <- subset(plants, PLANTSPECIES %in% plant.sp)

# -------------------------------------------------------------------------

# peralta
birds.peralta$BIRDSPECIES[which(birds.peralta$BIRDSPECIES == "Callaeas cinerea")] <- "Callaeas cinereus"
birds.peralta$BIRDSPECIES[which(birds.peralta$BIRDSPECIES == "Callaeas wilsoni")] <- "Callaeas cinereus"
birds.peralta$BIRDSPECIES[which(birds.peralta$BIRDSPECIES == "Carduelis chloris")] <- "Chloris chloris"
birds.peralta$BIRDSPECIES[which(birds.peralta$BIRDSPECIES == "Carduelis flammea")] <- "Acanthis flammea"
birds.peralta$BIRDSPECIES[which(birds.peralta$BIRDSPECIES == "Larus novaehollandiae")] <- "Chroicocephalus novaehollandiae"
birds.peralta$BIRDSPECIES[which(birds.peralta$BIRDSPECIES == "Hemiphaga novaeseelandiae chathamensis")] <- "Hemiphaga novaeseelandiae"
birds.peralta$BIRDSPECIES[which(birds.peralta$BIRDSPECIES == "Mohoua novaeseelandiae")] <- "Finschia novaeseelandiae"
birds.peralta$BIRDSPECIES[which(birds.peralta$BIRDSPECIES == "Strigops habroptilus")] <- "Strigops habroptila"

birds.peralta$BIRDSPECIES <- str_replace(birds.peralta$BIRDSPECIES," ","_")

birds.peralta <- subset(birds.peralta, BIRDSPECIES %in% bird.sp)

# avonet
birds.avonet <- birds.avonet[,c("Species1_BirdLife","Sex","Country_WRI","Beak.Length_Culmen","Beak.Length_Nares",
                                "Beak.Width","Beak.Depth","Tarsus.Length","Wing.Length","Kipps.Distance","Secondary1",
                                "Hand-wing.Index","Tail.Length")]

birds.avonet$name <- str_replace(birds.avonet$Species1_BirdLife," ","_")
birds.avonet.2 <- subset(birds.avonet, name %in% bird.sp)

# body sizes
birds.body.sizes <- birds.body.sizes[,c("Species_name","M_mass","F_mass")]
birds.body.sizes$name <- str_replace(birds.body.sizes$Species_name," ","_")
birds.body.sizes.2 <- subset(birds.body.sizes, name %in% bird.sp)

# -------------------------------------------------------------------------

plant.traits <- plants %>%
  select(PLANTSPECIES,FRUIT_LENGTH_mm:PLANT_MAX_MEAN_VEGETATIVE_HEIGHT_m) %>%
  pivot_longer(FRUIT_LENGTH_mm:PLANT_MAX_MEAN_VEGETATIVE_HEIGHT_m,names_to = "trait",values_to = "value") %>%
  group_by(PLANTSPECIES,trait) %>%
  summarise(n = n(),mean.value = mean(value,na.rm = TRUE),sd.value = sd(value,na.rm = TRUE)) %>%
  mutate(guild = "plants") %>%
  rename(species = PLANTSPECIES)

birds.peralta.traits <- birds.peralta %>%
  select(BIRDSPECIES,BILL_LENGTH_mm:WING_LENGTH_cm) %>%
  pivot_longer(BILL_LENGTH_mm:WING_LENGTH_cm,names_to = "trait",values_to = "value") %>%
  group_by(BIRDSPECIES,trait) %>%
  summarise(n = n(),mean.value = mean(value,na.rm = TRUE),sd.value = sd(value,na.rm = TRUE)) %>%
  mutate(guild = "birds") %>%
  rename(species = BIRDSPECIES)

birds.avonet.traits <- birds.avonet.2 %>%
  select(name,Beak.Length_Culmen:Tail.Length) %>%
  pivot_longer(Beak.Length_Culmen:Tail.Length,names_to = "trait",values_to = "value") %>%
  group_by(name,trait) %>%
  summarise(n = n(),mean.value = mean(value,na.rm = TRUE),sd.value = sd(value,na.rm = TRUE)) %>%
  mutate(guild = "birds") %>%
  rename(species = name)

birds.bs.trait <- birds.body.sizes.2 %>%
  select(name,M_mass,F_mass) %>%
  group_by(name) %>%
  summarise(body.mass = (M_mass + F_mass)/2) %>%
  pivot_longer(body.mass,names_to = "trait",values_to = "value") %>%
  rename(species = name) %>%
  mutate(mean.value = na_if(value,-999),sd.value = NA,guild = "birds") %>%
  select(species,guild,trait,mean.value,sd.value)

# -------------------------------------------------------------------------

bird.trait.data <- bind_rows(birds.peralta.traits,birds.avonet.traits, birds.bs.trait)
bird.trait.data.2 <- complete(ungroup(bird.trait.data),species,trait)
bird.trait.data.2$guild <- "birds"

trait.data <- bind_rows(bird.trait.data.2,plant.traits)

# -------------------------------------------------------------------------
# add species status

names(plant.status) <- c("species","status")
names(bird.status) <- c("species","status")

plant.status$species[which(plant.status$species == "Actinidia_deliciosa")] <- "Actinidia_chinensis"
plant.status$species[which(plant.status$species == "Androstoma_empetrifolium")] <- "Androstoma_empetrifolia"
plant.status$species[which(plant.status$species == "Dendrobenthamia_capitata")] <- "Cornus_capitata"
plant.status$species[which(plant.status$species == "Eriobotrya_japonica")] <- "Rhaphiolepis_loquata"
plant.status$species[which(plant.status$species == "Leucopogon_fraseri")] <- "Styphelia_nesophila"
plant.status$species[which(plant.status$species == "Phormium_cookianum")] <- "Phormium_colensoi"
plant.status$species[which(plant.status$species == "Phyllocladus_alpinus")] <- "Phyllocladus_trichomanoides"
plant.status$species[which(plant.status$species == "Phytolacca_octandra")] <- "Phytolacca_icosandra"
plant.status$species[which(plant.status$species == "Piper_excelsum")] <- "Macropiper_excelsum"
plant.status$species[which(plant.status$species == "Pseudopanax_colensoi")] <- "Neopanax_colensoi"
plant.status$species[which(plant.status$species == "Pseudopanax_colensoi_var_colensoi")] <- "Neopanax_colensoi"
plant.status$species[which(plant.status$species == "Pseudopanax_colensoi_var_ternatus")] <- "Neopanax_colensoi"
plant.status$species[which(plant.status$species == "Solanum_nodiflorum")] <- "Solanum_americanum"
plant.status$species[which(plant.status$species == "Streblus_heterophyllus")] <- "Paratrophis_microphylla"

bird.status$species[which(bird.status$species == "Callaeas_cinerea")] <- "Callaeas_cinereus"
bird.status$species[which(bird.status$species == "Callaeas_wilsoni")] <- "Callaeas_cinereus"
bird.status$species[which(bird.status$species == "Carduelis_chloris")] <- "Chloris_chloris"
bird.status$species[which(bird.status$species == "Carduelis_flammea")] <- "Acanthis_flammea"
bird.status$species[which(bird.status$species == "Larus_novaehollandiae")] <- "Chroicocephalus_novaehollandiae"
bird.status$species[which(bird.status$species == "Hemiphaga_novaeseelandiae_chathamensis")] <- "Hemiphaga_novaeseelandiae"
bird.status$species[which(bird.status$species == "Mohoua_novaeseelandiae")] <- "Finschia_novaeseelandiae"
bird.status$species[which(bird.status$species == "Strigops_habroptilus")] <- "Strigops_habroptila"

status.clean <- bind_rows(plant.status,bird.status)

# -------------------------------------------------------------------------

trait.data$status <- status.clean$status[match(trait.data$species,status.clean$species)]
trait.data <- arrange(trait.data, guild, species,trait)
trait.data <- trait.data[,c("species","guild","status","n","trait","mean.value","sd.value")]

# -------------------------------------------------------------------------
# bird lifestyle info

birds.avonet.lifestyle$species <- str_replace(birds.avonet.lifestyle$Species1," ","_")
birds.avonet.lifestyle <- birds.avonet.lifestyle[,c("species","Habitat","Trophic.Level","Trophic.Niche","Primary.Lifestyle")]
birds.avonet.lifestyle.2 <- subset(birds.avonet.lifestyle, species %in% bird.sp)

# -------------------------------------------------------------------------

write.csv2(trait.data,"results/trait_data.csv",row.names = FALSE)
write.csv2(birds.avonet.lifestyle.2,"results/bird_lifestyle.csv",row.names = FALSE)


