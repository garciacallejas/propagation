
# obtain dispersal distances for bird species based on 
# its relationship with hand-wing index

# INPUTS
# - trait data: trait.data.csv
# - avonet data: /home/david/Work/datasets/NZ/bird traits/AVONET/AVONET1_BirdLife
# - man/min dispersal distances (constants)

# OUTPUTS

# -------------------------------------------------------------------------

library(tidyverse)

# -------------------------------------------------------------------------

trait.data <- read.csv2("results/trait_data.csv")
avonet.data.nz <- read.csv("/home/david/Work/datasets/NZ/bird traits/AVONET/AVONET1_BirdLife.csv")

# -------------------------------------------------------------------------

flightless <- c("Porphyrio_hochstetteri","Apteryx_mantelli","Apteryx_owenii",
                "Apteryx_haastii","Strigops_habroptila","Gallirallus_australis",
                "Anas_aucklandica","Anas_nesiotis")

trait.data$flight <- TRUE
trait.data$flight[trait.data$species %in% flightless] <- FALSE

trait.data$flight[trait.data$guild != "birds"] <- NA

# -------------------------------------------------------------------------

max.dispersal <- 20
min.dispersal <- 0

hwi <- subset(trait.data, trait == "Hand-wing.Index")
na.hwi <- subset(hwi, is.na(mean.value))
na.hwi.sp <- na.hwi$species

na.hwi$mean.value <- avonet.data.nz$Hand.Wing.Index[match(na.hwi$species,avonet.data.nz$species)]

# from Sheard et al. 2020
# https://doi.org/10.1038/s41467-020-16313-6
na.hwi$mean.value[na.hwi$species == "Chroicocephalus_novaehollandiae"] <- 53.01
na.hwi$mean.value[na.hwi$species == "Finschia_novaeseelandiae"] <- 20.16

hwi$mean.value[which(is.na(hwi$mean.value))] <- na.hwi$mean.value

my.slope <- (max.dispersal - min.dispersal)/(max(hwi$mean.value,na.rm = T) - min(hwi$mean.value,na.rm = T))

my.intercept <- -my.slope*min(hwi$mean.value,na.rm = T)

hwi$dispersal.distance <- my.slope*hwi$mean.value + my.intercept
hwi$dispersal.distance[which(is.na(hwi$dispersal.distance))] <- 0
hwi$dispersal.distance[which(hwi$flight == F)] <- 0

hwi <- hwi[,c("species","flight","status","mean.value","dispersal.distance")]
names(hwi)[4] <- "hand.wing.index"

# -------------------------------------------------------------------------

write.csv2(hwi,paste("results/bird_dispersal_distances_",max.dispersal,"km.csv",sep=""),row.names = F)
