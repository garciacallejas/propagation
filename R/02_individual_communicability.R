
# Script to average communicability values over individual species

# INPUTS: 
# - pairwise communicability dataframe: "{external_path}/results/communicability_pairwise.Rdata"
# - bird presences in cells: "data/bird_cell_presences.csv"
# - plant presences in cells: "data/plant_cell_presences.csv"

# OUTPUTS:
# - individual communicability df: "results/individual_aggregated_communicability.csv"

# NOTE:
# some inputs are too big for git/github. load them externally and keep 
# the path {external_path} always the same

external_path <- "/home/david/Work/datasets/NZ/"

# -------------------------------------------------------------------------

library(tidyverse)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# -------------------------------------------------------------------------

# this is a huge file (~12GB on memory)
load(paste(external_path,"results/communicability_pairwise.Rdata",sep=""))

bird.sp.cells.wide <- read.csv2("data/bird_cell_presences.csv")
plant.sp.cells.wide <- read.csv2("data/plant_cell_presences.csv")

bird.sp <- sort(unique(bird.sp.cells.wide$species))
plant.sp <- sort(unique(plant.sp.cells.wide$species))
all.sp <- sort(unique(c(bird.sp,plant.sp)))
num.sp <- length(all.sp)

represented.cells <- unique(c(names(bird.sp.cells.wide),names(plant.sp.cells.wide)))
represented.cells <- as.numeric(substr(represented.cells[-1],
                                       start = 2,
                                       stop = nchar(represented.cells[-1])))

# -------------------------------------------------------------------------

df1$guild.sp1 <- ifelse(df1$sp1 %in% bird.sp, "birds", "plants")
df1$guild.sp2 <- ifelse(df1$sp2 %in% bird.sp, "birds", "plants")

df.ind <- df1 %>% 
  group_by(sp1,guild.sp1,grid.id.sp1) %>%
  summarise(aggregated.communicability = sum(binary.communicability)) %>%
  mutate(scaled.aggregated.communicability = range01(aggregated.communicability)) %>%
  rename(species = sp1,
         guild = guild.sp1,
         cell = grid.id.sp1)

# -------------------------------------------------------------------------

write.csv2(df.ind,"results/individual_aggregated_communicability.csv",row.names = FALSE)


