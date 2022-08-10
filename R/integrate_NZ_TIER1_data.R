# integrate Tier 1 DOC observations and NVS ones into GBIF

# note from NVS data: use this reference in publications
# GBIF.org (09 July 2022) GBIF Occurrence Download https://doi.org/10.15468/dl.qeptyv

# -------------------------------------------------------------------------
library(tidyverse)

# -------------------------------------------------------------------------

plot.coords <- read.csv2("/home/david/Work/datasets/NZ/Tier1_DOC/nz_plot_locations.csv")
names(plot.coords) <- c("plot","x","y")
# Im using 
bird.names <- read.csv2("/home/david/Work/datasets/NZ/Tier1_DOC/nz_bird_names_simplified.csv")
unresolved.names <- read.csv2("data/unresolved_bird_names.csv")
bird.names <- bird.names[,c("SpeciesName","ScientificName","ScientificNameSpeciesLevel")]

# -------------------------------------------------------------------------
# birds

b1 <- read.csv2("/home/david/Work/datasets/NZ/Tier1_DOC/bird_records_1.csv")

invalid.sp <- c("None_identified","Unidentifiable","Unknown","No species recorded")

b1.1 <- subset(b1,!species %in% invalid.sp)
b1.1$year <- substr(b1.1$date,nchar(b1.1$date)-3,nchar(b1.1$date))

b1.2 <- b1.1[,c("sampling_location","year","species")]
names(b1.2)[1] <- "plot"

b1.3 <- b1.2 %>% group_by(plot,year,species) %>%
  summarise(number = n()) %>%
  mutate(species_group = "birds") 

# -------------------------------------------------------------------------
b2 <- read.csv2("/home/david/Work/datasets/NZ/Tier1_DOC/bird_records_2.csv")

b2.1 <- subset(b2, !SpeciesName %in% invalid.sp & !is.na(Near) & !is.na(Far) & !is.na(VeryFar))
b2.1$year <- substr(b2.1$DateStarted,nchar(b2.1$DateStarted)-3,nchar(b2.1$DateStarted))
b2.1$total <- b2.1$Near + as.numeric(b2.1$Far) + as.numeric(b2.1$VeryFar)

b2.2 <- b2.1 %>%
  filter(!is.na(total)) %>%
  dplyr::select(Place,year,SpeciesName,total) %>%
  group_by(Place,year,SpeciesName) %>%
  summarise(number = n()) %>%
  mutate(species_group = "birds")
names(b2.2) <- names(b1.3)

# -------------------------------------------------------------------------
b3 <- read.csv2("/home/david/Work/datasets/NZ/Tier1_DOC/bird_records_3.csv")

b3.1 <- subset(b3,!SpeciesName %in% invalid.sp & NoSpeciesRecorded == "N")
b3.1$year <- substr(b3.1$DateStarted,nchar(b3.1$DateStarted)-3,nchar(b3.1$DateStarted))

b3.2 <- b3.1 %>% group_by(Place,year,SpeciesName) %>%
  summarise(number = n()) %>%
  mutate(species_group = "birds") 

names(b3.2) <- names(b1.3)

# -------------------------------------------------------------------------
b4 <- read.csv2("/home/david/Work/datasets/NZ/Tier1_DOC/bird_records_4.csv")

b4.1 <- subset(b4,!SpeciesName %in% invalid.sp)
b4.1$year <- substr(b4.1$DateStarted,nchar(b4.1$DateStarted)-3,nchar(b4.1$DateStarted))

b4.2 <- b4.1 %>% group_by(Place,year,SpeciesName) %>%
  summarise(number = n()) %>%
  mutate(species_group = "birds") 
names(b4.2) <- names(b1.3)

# -------------------------------------------------------------------------
# b4 is from bird flocks in the sky - it does not assign a number of individuals
# and birds could be moving, etc, so perhaps it's safer not to include them
# b.all <- bind_rows(b1.3,b2.2,b3.2,b4.2) %>%
b.all <- bind_rows(b1.3,b2.2,b3.2) %>%
  filter(!is.na(species)) %>%
  group_by(plot,year,species_group,species) %>%
  summarise(records = sum(number))

# resolved names for those not in the thesaurus
b.all$resolved.name <- unresolved.names$equivalent.resolved[match(b.all$species,unresolved.names$unresolved.name)]
# modify the original column (species) only for those with a resolved name (i.e. with a wrong original name)
b.all$species[which(!is.na(b.all$resolved.name))] <- b.all$resolved.name[which(!is.na(b.all$resolved.name))]
# now, safely assign scientific names - at the species level, to be consistent
# with the interaction data
b.all$scientific_name <- bird.names$ScientificNameSpeciesLevel[match(b.all$species,bird.names$SpeciesName)]
# check that there are only the NAs contained in the "nz_bird_names" file: finch, dove, shag
# unres <- subset(b.all,is.na(scientific_name))
b.all.2 <- subset(b.all, !is.na(scientific_name))
# tidy columns
b.all.2$species <- b.all.2$scientific_name
b.all.2$scientific_name <- NULL
b.all.2$resolved.name <- NULL

b.all.3 <- left_join(b.all.2,plot.coords) %>% 
  dplyr::select(x,y,year,species_group,species,records)

write.csv2(b.all.3,"data/birds_tier1.csv",row.names = F)

# -------------------------------------------------------------------------
# mammals
# TODO: add scientific names

m1 <- read.csv2("/home/david/Work/datasets/NZ/Tier1_DOC/mammal_records_1.csv")

m1.1 <- subset(m1,!SpeciesName %in% invalid.sp)
m1.1$year <- substr(m1.1$DateStarted,nchar(m1.1$DateStarted)-3,nchar(m1.1$DateStarted))
m1.1$NumberObserved[m1.1$NumberObserved == "Present"] <- 1
m1.1$NumberObserved <- as.numeric(m1.1$NumberObserved)

m1.2 <- m1.1 %>%
  dplyr::select(Place,year,SpeciesName,NumberObserved) %>%
  mutate(species_group = "mammals")
names(m1.2) <- names(b1.3)

# -------------------------------------------------------------------------

m2 <- read.csv2("/home/david/Work/datasets/NZ/Tier1_DOC/mammal_records_2.csv")

m2.1 <- subset(m2,!SpeciesName %in% invalid.sp)
m2.1$year <- substr(m2.1$DateStarted,nchar(m2.1$DateStarted)-3,nchar(m2.1$DateStarted))

m2.2 <- m2.1 %>%
  dplyr::select(Place,year,SpeciesName) %>%
  mutate(number = 1,
         species_group = "mammals")
names(m2.2) <- names(b1.3)


# -------------------------------------------------------------------------
# possums
p1 <- read.csv2("/home/david/Work/datasets/NZ/Tier1_DOC/possum_records_1.csv")
# these codes all represent at least the presence of possums
valid.results <- c("BR","FA","FI","MA","MI","P, NT","P","P, U","ESC")

p1.1 <- p1 %>%
  mutate(year = substr(DateStarted,nchar(DateStarted)-3,nchar(DateStarted))) %>%
  dplyr::select(Place,year,Night1Result,Night2Result) %>%
  filter(Night1Result %in% valid.results | Night2Result %in% valid.results) %>%
  mutate(species = "Possum",
         number = 1,
         species_group = "mammals") %>%
  select(Place,year,species,number,species_group)
names(p1.1)[1] <- "plot"

# -------------------------------------------------------------------------
# ungulates
u1 <- read.csv2("/home/david/Work/datasets/NZ/Tier1_DOC/ungulate_records_1.csv")

u1.1 <- u1 %>%
  mutate(year = substr(DateStarted,nchar(DateStarted)-3,nchar(DateStarted))) %>%
  dplyr::select(Place,year,Non_intactUngulatePellets:WallabyPellets)
names(u1.1) <- c("plot","year","Unknown ungulate","Possum","Rabbit","Hare","Pig1","Pig2","Wallaby")

u1.1$Pig <- ifelse(u1.1$Pig1 == "Y" | u1.1$Pig2 == "Y","Y","N")
u1.1$Pig1 <- NULL
u1.1$Pig2 <- NULL

u1.2 <- u1.1 %>% pivot_longer("Unknown ungulate":"Pig",names_to = "species",values_to = "value") %>%
  filter(value == "Y") %>%
  mutate(number = 1,
         species_group = "mammals") %>%
  select(-value)

# -------------------------------------------------------------------------
# mammal, possum, and ungulate records

m.all <- bind_rows(m1.2,m2.2,p1.1,u1.2) %>% 
  filter(!is.na(species)) %>%
  group_by(plot,year,species_group,species) %>%
  summarise(records = sum(number))

m.all.2 <- left_join(m.all,plot.coords) %>% 
  dplyr::select(x,y,year,species_group,species,records)

write.csv2(m.all.2,"data/mammals_tier1.csv",row.names = F)

