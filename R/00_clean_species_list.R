
library(tidyverse)

# -------------------------------------------------------------------------
grid.size <- 10 # km

# this dataset is already filtered to the two main islands, 
# so only species that appear here will be considered
sp.obs <- read.csv2(paste("results/sp_observations_long_",grid.size,"km.csv",sep=""))

all.sp.obs <- sort(unique(sp.obs$species))

orig.data <- read.csv("../datasets/plant-bird interactions and traits/plant_bird_interactions.csv")
clean.data <- orig.data

# small issues

# these names were updated in the gbif occ_search
# Actinidia deliciosa is a variety of Actinidia chinensis
clean.data$PLANTSPECIES[which(clean.data$PLANTSPECIES == "Actinidia_deliciosa")] <- "Actinidia_chinensis"
# Androstoma empetrifolium -> a
clean.data$PLANTSPECIES[which(clean.data$PLANTSPECIES == "Androstoma_empetrifolium")] <- "Androstoma_empetrifolia"
# citrus paradisi is not retrieved
# citrus sinensis is not retrieved
# dacrydium cupressinum... downloaded manually
# Dendrobenthamia_capitata is cornus capitata
clean.data$PLANTSPECIES[which(clean.data$PLANTSPECIES == "Dendrobenthamia_capitata")] <- "Cornus_capitata"
# Dysoxylum_spectabile... downloaded manually
# Eriobotrya japonica is thapiolepis loquata
clean.data$PLANTSPECIES[which(clean.data$PLANTSPECIES == "Eriobotrya_japonica")] <- "Rhaphiolepis_loquata"
# Leucopogon_fraseri is Styphelia_nesophila
clean.data$PLANTSPECIES[which(clean.data$PLANTSPECIES == "Leucopogon_fraseri")] <- "Styphelia_nesophila"
# Phormium_cookianum is Phormium_colensoi
clean.data$PLANTSPECIES[which(clean.data$PLANTSPECIES == "Phormium_cookianum")] <- "Phormium_colensoi"
# Phyllocladus_alpinus is Phyllocladus_trichomanoides
clean.data$PLANTSPECIES[which(clean.data$PLANTSPECIES == "Phyllocladus_alpinus")] <- "Phyllocladus_trichomanoides"
# Phytolacca_octandra is Phytolacca_icosandra
clean.data$PLANTSPECIES[which(clean.data$PLANTSPECIES == "Phytolacca_octandra")] <- "Phytolacca_icosandra"
# Piper_excelsum is Macropiper_excelsum
clean.data$PLANTSPECIES[which(clean.data$PLANTSPECIES == "Piper_excelsum")] <- "Macropiper_excelsum"
# Pseudopanax_colensoi is Neopanax_colensoi
clean.data$PLANTSPECIES[which(clean.data$PLANTSPECIES == "Pseudopanax_colensoi")] <- "Neopanax_colensoi"
clean.data$PLANTSPECIES[which(clean.data$PLANTSPECIES == "Pseudopanax_colensoi_var_colensoi")] <- "Neopanax_colensoi"
clean.data$PLANTSPECIES[which(clean.data$PLANTSPECIES == "Pseudopanax_colensoi_var_ternatus")] <- "Neopanax_colensoi"
# Pseudopanax chathamicus only appears in small islands, very rare
# Psidium_guajava is not so common, but should be there
# Solanum_nodiflorum is Solanum_americanum
clean.data$PLANTSPECIES[which(clean.data$PLANTSPECIES == "Solanum_nodiflorum")] <- "Solanum_americanum"
# Streblus_heterophyllus is Paratrophis_microphylla
clean.data$PLANTSPECIES[which(clean.data$PLANTSPECIES == "Streblus_heterophyllus")] <- "Paratrophis_microphylla"
# Podocarpus cunninghamii is Podocarpus hallii
clean.data$PLANTSPECIES[which(clean.data$PLANTSPECIES == "Podocarpus cunninghamii")] <- "Podocarpus hallii"
# Tropaeolum_speciosum downloaded manually

plant.sp <- sort(unique(clean.data$PLANTSPECIES))

# simplify subspecies/varieties
plant.sp <- str_replace(plant.sp,"_"," ")
plant.sp <- sort(unique(gsub("\\_.*","",plant.sp)))
plant.sp <- str_replace(plant.sp," ","_")

# -------------------------------------------------------------------------

# these names were automatically updated in the occ_search
clean.data$BIRDSPECIES[which(clean.data$BIRDSPECIES == "Callaeas_cinerea")] <- "Callaeas_cinereus"
# Callaeas_wilsoni is a subsp of Callaeas_cinereus
clean.data$BIRDSPECIES[which(clean.data$BIRDSPECIES == "Callaeas_wilsoni")] <- "Callaeas_cinereus"
# Carduelis chloris is Chloris chloris
clean.data$BIRDSPECIES[which(clean.data$BIRDSPECIES == "Carduelis_chloris")] <- "Chloris_chloris"
# Carduelis flammea is Acanthis flammea
clean.data$BIRDSPECIES[which(clean.data$BIRDSPECIES == "Carduelis_flammea")] <- "Acanthis_flammea"
# cyanoramphus forbesi only appears in Chatham Island
# Larus novaehollandiae is Chroicocephalus novaehollandiae
clean.data$BIRDSPECIES[which(clean.data$BIRDSPECIES == "Larus_novaehollandiae")] <- "Chroicocephalus_novaehollandiae"
# Mohoua novaeseelandiae is Finschia novaeseelandiae
clean.data$BIRDSPECIES[which(clean.data$BIRDSPECIES == "Mohoua_novaeseelandiae")] <- "Finschia_novaeseelandiae"
# Philesturnus rufusater is found only in offshore islands and sanctuaries
# Strygops habroptilus is Strygops habroptila
clean.data$BIRDSPECIES[which(clean.data$BIRDSPECIES == "Strigops_habroptilus")] <- "Strigops_habroptila"
# simplify this subsp
clean.data$BIRDSPECIES[which(clean.data$BIRDSPECIES == "Hemiphaga_novaeseelandiae_chathamensis")] <- "Hemiphaga_novaeseelandiae"
# simplify this subsp
clean.data$BIRDSPECIES[which(clean.data$BIRDSPECIES == "Cyanoramphus_novaezelandiae_X_forbesi")] <- "Cyanoramphus_novaezelandiae"

bird.sp <- sort(unique(clean.data$BIRDSPECIES))

# simplify subspecies/varieties
bird.sp <- str_replace(bird.sp,"_"," ")
bird.sp <- sort(unique(gsub("\\_.*","",bird.sp)))
bird.sp <- str_replace(bird.sp," ","_")

# -------------------------------------------------------------------------

bird.df <- data.frame(species = bird.sp,guild = "birds")
bird.df$status <- clean.data$BIRDSTATUS[match(bird.df$species,clean.data$BIRDSPECIES)]

plant.df <- data.frame(species = plant.sp,guild = "plants")
plant.df$status <- clean.data$PLANTSTATUS[match(plant.df$species,clean.data$PLANTSPECIES)]

all.sp <- bind_rows(plant.df,bird.df)

# -------------------------------------------------------------------------
write.csv2(all.sp,"data/species_list.csv",row.names = F)
write.csv2(clean.data,"results/plant_bird_interactions_clean.csv",row.names = F)


# -------------------------------------------------------------------------
# sp.files <- list.files(path = "results/sp_observations",
#                        pattern = "*.csv",
#                        full.names = T)
# sp.obs <- sp.files %>% map_dfr(read_csv2)
# sp.obs <- subset(sp.obs, !is.na(year))
# 
# # -------------------------------------------------------------------------
# 
# # issues:
# # set a temporal threshold for observations?
# # min(sp.obs$year,na.rm = T)
# 
# min.year <- 2010
# sp.obs.2 <- subset(sp.obs, year >= min.year)
# 
# # assume that bird populations from adjacent cells are linked
# # in next iterations, this can be improved by accounting for differences in dispersal ability
# # possibly inferred by looking at morphological traits
# 
# # -------------------------------------------------------------------------
# # keep the set of all species
# sp.obs.2$species <- str_replace(sp.obs.2$species," ","_")
# all.sp <- sort(unique(sp.obs.2$species))
# 
# # now, a few bird species are not present in the two main islands, which is 
# # the territory I am considering. 
# # so, remove them
# bird.sp <- bird.sp[which(bird.sp %in% all.sp & bird.sp %in% all.sp.obs)]
# 
# # for plants, I discard taxa identified at the genus level,
# # and a couple of very rare species (Pseudopanax chathamicus,
# # Jasminum officinale, citrus paradisi, citrus sinensis)
# plant.sp <- plant.sp[which(plant.sp %in% all.sp & plant.sp %in% all.sp.obs)]
# # discarded.plants <- plant.sp[which(!(plant.sp %in% all.sp))]
# 
# # -------------------------------------------------------------------------
# # with this clean set of plant and bird species,
# # I need to update the list of interactions to only account for these sp
# clean.int.data <- subset(orig.data, BIRDSPECIES %in% bird.sp & PLANTSPECIES %in% plant.sp)
# 
# # -------------------------------------------------------------------------
# 
# write.csv2(clean.int.data,"results/plant_bird_interactions_clean.csv",row.names = F)
# write.csv2(sp.obs.2,"results/plant_bird_observations_clean.csv",row.names = F)

