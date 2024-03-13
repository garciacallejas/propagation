# clean up trait datasets

# for birds, Peralta and Avonet, as well as Lislevand for body sizes
# for plants, Peralta for now

library(tidyverse)
library(taxize)

# -------------------------------------------------------------------------

# independent list of plants/birds
sp.list <- read.csv2("data/species_list.csv")

bird.sp <- sort(unique(sp.list$species[which(sp.list$guild == "birds")]))
plant.sp <- sort(unique(sp.list$species[which(sp.list$guild == "plants")]))
all.sp <- sort(unique(c(bird.sp,plant.sp)))
num.sp <- length(all.sp)

# status
plant.status <- subset(sp.list,guild == "plants")
bird.status <- subset(sp.list,guild == "birds")

# -------------------------------------------------------------------------

# trait datasets
# these are all publicly available
plants <- read.csv("data/plant-bird interactions and traits/plant_traits.csv")
birds.peralta <- read.csv("data/plant-bird interactions and traits/bird_traits.csv")
birds.avonet <- read_csv("data/AVONET/AVONET_Raw_Data.csv")
birds.avonet.lifestyle <- read_csv("data/AVONET/AVONET1_BirdLife.csv")
birds.body.sizes <- read_tsv("data/body_sizes/avian_ssd_jan07.txt")

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

plants$PLANTSPECIES[which(plants$PLANTSPECIES == "Olea_europaea_subsp_cuspidata")] <- "Olea_europaea"
plants$PLANTSPECIES[which(plants$PLANTSPECIES == "Passiflora_edulis_f_edulis")] <- "Passiflora_edulis"
plants$PLANTSPECIES[which(plants$PLANTSPECIES == "Passiflora_tripartita_var_mollissima")] <- "Passiflora_tripartita"

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

status.clean <- bind_rows(plant.status,bird.status)

# -------------------------------------------------------------------------
# get species IUCN status, using the token from the examples webpage

API <- "9bb4facb6d23f48efbf424bb05c0c1ef1cf6f468393bc745d42179ac4aca5fee"
all.sp.iucn <- str_replace(all.sp,"_"," ")
my.sp.iucn.status <- iucn_summary(all.sp.iucn,key = API)
my.sp.iucn.status.df <- iucn_status(my.sp.iucn.status) %>%
  as.data.frame() %>%
  rownames_to_column(var = "species") %>% set_names(c("species","IUCN_status"))

my.sp.iucn.status.df$species <- str_replace(my.sp.iucn.status.df$species," ","_")
status.clean.2 <- left_join(status.clean,my.sp.iucn.status.df)
# threatened.sp <- my.sp.iucn.status.df %>% filter(IUCN_status %in% c("CR","EN","LR/cd","LR/nt","NT","VU"))

# NOTE: this returns many species without information (NA)
# another source is the NZ Department of Conservation, in the report:
# https://www.doc.govt.nz/globalassets/documents/science-and-technical/nztcs36entire.pdf
# being "only" 22 bird species, here I input the status of each bird manually
clean.bird.sp <- c("Acanthis_flammea","Acanthisitta_chloris","Acridotheres_tristis","Anthornis_melanura",
                   "Anthus_novaeseelandiae","Carduelis_carduelis","Fringilla_coelebs","Gallirallus_australis",
                   "Gerygone_igata","Hemiphaga_novaeseelandiae","Mohoua_albicilla","Nestor_meridionalis",
                   "Nestor_notabilis","Petroica_australis","Petroica_macrocephala","Platycercus_eximius",
                   "Prosthemadera_novaeseelandiae","Rhipidura_fuliginosa","Sturnus_vulgaris","Turdus_merula",
                   "Turdus_philomelos","Zosterops_lateralis" )
NZ_doc_status <- data.frame(species = clean.bird.sp, NZ_doc_status = NA, NZ_doc_qualifier = NA)
# NZ_doc_status <- data.frame(species = clean.bird.sp, NZ_doc_status = NA, NZ_doc_qualifier = NA)

NZ_doc_status$NZ_doc_status[NZ_doc_status$species == "Acanthis_flammea"] <- "not threatened"
NZ_doc_status$NZ_doc_status[NZ_doc_status$species == "Acanthisitta_chloris"] <- "at risk"
NZ_doc_status$NZ_doc_status[NZ_doc_status$species == "Acridotheres_tristis"] <- "not threatened"
NZ_doc_status$NZ_doc_status[NZ_doc_status$species == "Anthornis_melanura"] <- "not threatened"
NZ_doc_status$NZ_doc_status[NZ_doc_status$species == "Anthus_novaeseelandiae"] <- "at risk"
NZ_doc_status$NZ_doc_status[NZ_doc_status$species == "Carduelis_carduelis"] <- "not threatened"
NZ_doc_status$NZ_doc_status[NZ_doc_status$species == "Fringilla_coelebs"] <- "not threatened"
NZ_doc_status$NZ_doc_status[NZ_doc_status$species == "Gallirallus_australis"] <- "relict"
NZ_doc_status$NZ_doc_status[NZ_doc_status$species == "Gerygone_igata"] <- "not threatened"
NZ_doc_status$NZ_doc_status[NZ_doc_status$species == "Hemiphaga_novaeseelandiae"] <- "not threatened"
NZ_doc_status$NZ_doc_status[NZ_doc_status$species == "Mohoua_albicilla"] <- "not threatened"
NZ_doc_status$NZ_doc_status[NZ_doc_status$species == "Nestor_meridionalis"] <- "at risk"
NZ_doc_status$NZ_doc_status[NZ_doc_status$species == "Nestor_notabilis"] <- "endangered"
NZ_doc_status$NZ_doc_status[NZ_doc_status$species == "Petroica_australis"] <- "at risk"
NZ_doc_status$NZ_doc_status[NZ_doc_status$species == "Petroica_macrocephala"] <- "not threatened"
NZ_doc_status$NZ_doc_status[NZ_doc_status$species == "Platycercus_eximius"] <- "not threatened"
NZ_doc_status$NZ_doc_status[NZ_doc_status$species == "Prosthemadera_novaeseelandiae"] <- "not threatened"
NZ_doc_status$NZ_doc_status[NZ_doc_status$species == "Rhipidura_fuliginosa"] <- "not threatened"
NZ_doc_status$NZ_doc_status[NZ_doc_status$species == "Sturnus_vulgaris"] <- "not threatened"
NZ_doc_status$NZ_doc_status[NZ_doc_status$species == "Turdus_merula"] <- "not threatened"
NZ_doc_status$NZ_doc_status[NZ_doc_status$species == "Turdus_philomelos"] <- "not threatened"
NZ_doc_status$NZ_doc_status[NZ_doc_status$species == "Zosterops_lateralis"] <- "not threatened"

# in case it's needed
# NZ_doc_status$NZ_doc_qualifier[NZ_doc_status$species == "Acanthis_flammea"] <- "SO"
# NZ_doc_status$NZ_doc_qualifier[NZ_doc_status$species == "Acanthisitta_chloris"] <- "not threatened"
# NZ_doc_status$NZ_doc_qualifier[NZ_doc_status$species == "Acridotheres_tristis"] <- "SO"
# NZ_doc_status$NZ_doc_qualifier[NZ_doc_status$species == "Anthornis_melanura"] <- "not threatened"
# NZ_doc_status$NZ_doc_qualifier[NZ_doc_status$species == "Anthus_novaeseelandiae"] <- "CI,CR"
# NZ_doc_status$NZ_doc_qualifier[NZ_doc_status$species == "Carduelis_carduelis"] <- "SO"
# NZ_doc_status$NZ_doc_qualifier[NZ_doc_status$species == "Fringilla_coelebs"] <- "SO"
# NZ_doc_status$NZ_doc_qualifier[NZ_doc_status$species == "Gallirallus_australis"] <- "CI"
# NZ_doc_status$NZ_doc_qualifier[NZ_doc_status$species == "Gerygone_igata"] <- "not threatened"
# NZ_doc_status$NZ_doc_qualifier[NZ_doc_status$species == "Hemiphaga_novaeseelandiae"] <- "CD"
# NZ_doc_status$NZ_doc_qualifier[NZ_doc_status$species == "Mohoua_albicilla"] <- "CD,CI,PD"
# NZ_doc_status$NZ_doc_qualifier[NZ_doc_status$species == "Nestor_meridionalis"] <- "CD,CR,PD,RF,PF"
# NZ_doc_status$NZ_doc_qualifier[NZ_doc_status$species == "Nestor_notabilis"] <- "CD,CI,CR"
# NZ_doc_status$NZ_doc_qualifier[NZ_doc_status$species == "Petroica_australis"] <- "CD,CI,DPT,RR"
# NZ_doc_status$NZ_doc_qualifier[NZ_doc_status$species == "Petroica_macrocephala"] <- "not threatened"
# NZ_doc_status$NZ_doc_qualifier[NZ_doc_status$species == "Platycercus_eximius"] <- "SO"
# NZ_doc_status$NZ_doc_qualifier[NZ_doc_status$species == "Prosthemadera_novaeseelandiae"] <- "Inc"
# NZ_doc_status$NZ_doc_qualifier[NZ_doc_status$species == "Rhipidura_fuliginosa"] <- "EF"
# NZ_doc_status$NZ_doc_qualifier[NZ_doc_status$species == "Sturnus_vulgaris"] <- "SO"
# NZ_doc_status$NZ_doc_qualifier[NZ_doc_status$species == "Turdus_merula"] <- "SO"
# NZ_doc_status$NZ_doc_qualifier[NZ_doc_status$species == "Turdus_philomelos"] <- "SO"
# NZ_doc_status$NZ_doc_qualifier[NZ_doc_status$species == "Zosterops_lateralis"] <- "SO"
status.clean.3 <- left_join(status.clean.2,NZ_doc_status)

# -------------------------------------------------------------------------

trait.data$status <- status.clean.3$status[match(trait.data$species,status.clean.3$species)]
trait.data$IUCN_status <- status.clean.3$IUCN_status[match(trait.data$species,status.clean.3$species)]
trait.data$NZ_doc_status <- status.clean.3$NZ_doc_status[match(trait.data$species,status.clean.3$species)]

trait.data <- arrange(trait.data, guild, species,trait)
trait.data <- trait.data[,c("species","guild","status","IUCN_status","NZ_doc_status","n","trait","mean.value","sd.value")]

# -------------------------------------------------------------------------
# bird lifestyle info

birds.avonet.lifestyle$species <- str_replace(birds.avonet.lifestyle$Species1," ","_")
birds.avonet.lifestyle <- birds.avonet.lifestyle[,c("species","Habitat","Trophic.Level","Trophic.Niche","Primary.Lifestyle")]
birds.avonet.lifestyle.2 <- subset(birds.avonet.lifestyle, species %in% bird.sp)

# -------------------------------------------------------------------------

write.csv2(trait.data,"data/trait_data.csv",row.names = FALSE)
write.csv2(birds.avonet.lifestyle.2,"data/bird_lifestyle.csv",row.names = FALSE)


