
# Script to average communicability values over grid cells and guilds
# also put together local degrees, that are used later in the statistical analyses

# -------------------------------------------------------------------------
# IMPORTANT NOTE:
# the inputs for this script are not uploaded as they involve a large number of files
# so this will not run as it is. The outputs of this script are already uploaded,
# so to reproduce the results, I recommend using the uploaded dataframes.
# To run this script, one needs to generate the raw results from script NZ08,
# which is computationally demanding.
# -------------------------------------------------------------------------

# INPUTS: 
# - pairwise communicability dataframes: "/results/communicability/NZ_networks/*.csv"
# - local node degrees: "results/sp_degrees/*.csv

# OUTPUTS:
# - average communicability by species/cell df: "results/*level_communicability.csv"
# - single dataframe of local degrees: "results/population_level_degrees_".

# -------------------------------------------------------------------------

library(tidyverse)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# -------------------------------------------------------------------------

grid.size <- 10

pair.comm <- list.files("results/communicability/NZ_networks/",
                        full.names = T) %>% 
  map_dfr(read.csv2)

# put together local degrees in a single dataframe
all.local.degrees <- list.files("results/sp_degrees/",
                                full.names = T) %>% 
  map_dfr(read.csv2)

# -------------------------------------------------------------------------

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
# this just to show the final number of bird and plant species considered
# in the statistical analyses
pair.comm.clean %>% select(species,guild) %>%
  unique() %>%
  count(guild)

# -------------------------------------------------------------------------

write.csv2(pair.comm.clean,paste("results/population_level_communicability_",grid.size,"km.csv",sep=""),row.names = FALSE)
write.csv2(sp.comm,paste("results/species_level_communicability_",grid.size,"km.csv",sep=""),row.names = FALSE)
write.csv2(cell.comm,paste("results/cell_level_communicability_",grid.size,"km.csv",sep=""),row.names = FALSE)
write.csv2(all.local.degrees,paste("results/population_level_degrees_",grid.size,"km.csv",sep=""),row.names = F)

