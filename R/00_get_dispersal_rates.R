
# obtain dispersal kernels for bird species based on 
# its relationship with hand-wing index

# INPUTS
# - trait data: trait.data.csv
# - avonet data: /home/david/Work/datasets/NZ/bird traits/AVONET/AVONET1_BirdLife
# - grid size: numeric constant

# OUTPUTS
# dataframe with dispersal "resistance" for every species and distance

# -------------------------------------------------------------------------

library(tidyverse)

# -------------------------------------------------------------------------

trait.data <- read.csv2("data/trait_data.csv")
avonet.data.nz <- read.csv("/home/david/Work/datasets/NZ/bird traits/AVONET/AVONET1_BirdLife.csv")

grid.size <- 10

# -------------------------------------------------------------------------

flightless <- c("Porphyrio_hochstetteri","Apteryx_mantelli","Apteryx_owenii",
                "Apteryx_haastii","Strigops_habroptila","Gallirallus_australis",
                "Anas_aucklandica","Anas_nesiotis")

trait.data$flight <- TRUE
trait.data$flight[trait.data$species %in% flightless] <- FALSE

trait.data$flight[trait.data$guild != "birds"] <- NA

# -------------------------------------------------------------------------

hwi <- subset(trait.data, trait == "Hand-wing.Index")
na.hwi <- subset(hwi, is.na(mean.value))
na.hwi.sp <- na.hwi$species

na.hwi$mean.value <- avonet.data.nz$Hand.Wing.Index[match(na.hwi$species,
                                                          avonet.data.nz$species)]

# from Sheard et al. 2020
# https://doi.org/10.1038/s41467-020-16313-6
na.hwi$mean.value[na.hwi$species == "Chroicocephalus_novaehollandiae"] <- 53.01
na.hwi$mean.value[na.hwi$species == "Finschia_novaeseelandiae"] <- 20.16

hwi$mean.value[which(is.na(hwi$mean.value))] <- na.hwi$mean.value

hwi <- hwi[,c("species","flight","status","mean.value")]
names(hwi)[4] <- "hand.wing.index"

# -------------------------------------------------------------------------
# an excetption is "Acanthis flammea", that has no trait info in the trait_data
# table, but does appear in the list of birds, and also in avonet

acanthis.hwi <- avonet.data.nz$Hand.Wing.Index[which(avonet.data.nz$Species1 == "Acanthis flammea")]
hwi.full <- rbind(hwi, data.frame(species = "Acanthis_flammea",
                                  flight = TRUE,
                                  status = "Exotic",
                                  hand.wing.index = acanthis.hwi))
hwi.full <- arrange(hwi.full,species)
# -------------------------------------------------------------------------

write.csv2(hwi.full,"results/bird_hand_wing_index.csv",row.names = F)

# -------------------------------------------------------------------------

# test distances

# 300 is large enough to assume that further distances are zero
# dist <- seq(0,300,by = grid.size)
# 
# kern.fun <- function(dist,hwi){
#   my.kern <- exp(-dist*(1/hwi))
# }
# 
# kern.df <- expand.grid(dist = dist, hwi = hwi$hand.wing.index)
# kern.df$disp.rate <- kern.fun(kern.df$dist,kern.df$hwi)
# kern.df$species <- hwi$species[match(kern.df$hwi,hwi$hand.wing.index)]

# write.csv2(kern.df, paste("results/species_dispersal_rates_",grid.size,"km.csv",sep=""),row.names = F)

# ggplot(kern.df,aes(x = dist, y = disp.rate, group = species)) + 
#   geom_line(aes(color = species)) +
#   NULL


