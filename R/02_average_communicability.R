
# Script to average communicability values over grid cells and guilds

# INPUTS: 
# - pairwise communicability dataframes: "/results/communicability/NZ_networks/*.csv"

# OUTPUTS:
# - average communicability by species/cell df: "results/*level_communicability.csv"

# -------------------------------------------------------------------------

library(tidyverse)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# -------------------------------------------------------------------------

grid.size <- 10

# TODO filter the files read to only those from grid.size km
# so far it works because I only have 10km files, but CAREFUL!
pair.comm <- list.files("results/communicability/NZ_networks/",
                        full.names = T) %>% 
  map_dfr(read.csv2)

# test
all.sp <- unique(pair.comm[,c("sp1","guild.sp1")])

sp.comm <- pair.comm %>% 
  group_by(sp1,guild.sp1) %>%
  summarise(bin.communicability = sum(population.bin.communicability),
            weighted.communicability = sum(population.weighted.communicability)) 
sp.comm$scaled.bin.communicability <- scales::rescale(sp.comm$bin.communicability)
sp.comm$scaled.weighted.communicability <- scales::rescale(sp.comm$weighted.communicability)

names(sp.comm)[1:2] <- c("species","guild")

cell.comm <- pair.comm %>%
  group_by(cell.id.sp1) %>%
  summarise(bin.communicability = sum(population.bin.communicability),
            weighted.communicability = sum(population.weighted.communicability)) 
cell.comm$scaled.bin.communicability <- scales::rescale(cell.comm$bin.communicability)
cell.comm$scaled.weighted.communicability <- scales::rescale(cell.comm$weighted.communicability)

names(cell.comm)[1] <- "cell_id"

# -------------------------------------------------------------------------
# clean "pair.comm" 
pair.comm.clean <- pair.comm %>% 
  select(-X) %>%
  rename(species = sp1,guild = guild.sp1,cell_id = cell.id.sp1)


# -------------------------------------------------------------------------

write.csv2(pair.comm.clean,paste("results/population_level_communicability_",grid.size,"km.csv",sep=""),row.names = FALSE)
write.csv2(sp.comm,paste("results/species_level_communicability_",grid.size,"km.csv",sep=""),row.names = FALSE)
write.csv2(cell.comm,paste("results/cell_level_communicability_",grid.size,"km.csv",sep=""),row.names = FALSE)


