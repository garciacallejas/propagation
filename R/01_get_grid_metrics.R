

library(tidyverse)

# -------------------------------------------------------------------------

# these are for building the list of interactions
clean.int.data <- read.csv2("data/plant_bird_clean_interaction_data.csv")
bird.sp.cells.wide <- read.csv2("data/bird_cell_presences.csv")
plant.sp.cells.wide <- read.csv2("data/plant_cell_presences.csv")

non.matching <- names(bird.sp.cells.wide)[which(!(names(bird.sp.cells.wide) %in% names(plant.sp.cells.wide)))]
plant.sp.cells.wide[,non.matching] <- 0

bird.sp.cells.wide$guild <- "birds"
plant.sp.cells.wide$guild <- "plants"

all.sp.cells.wide <- bind_rows(bird.sp.cells.wide,plant.sp.cells.wide)

all.sp.cells <- pivot_longer(all.sp.cells.wide,2:(ncol(all.sp.cells.wide)-1),
                             names_to = "cell",
                             values_to = "occurrences")
all.sp.cells$cell <- substr(all.sp.cells$cell,
                            start = 2,
                            stop = nchar(all.sp.cells$cell))

cell.richness <- all.sp.cells %>% 
  filter(occurrences > 0) %>%
  group_by(cell,guild) %>%
  summarise(richness = n())

# -------------------------------------------------------------------------

write.csv2(cell.richness,file = "results/cell_richness.csv",row.names = FALSE)
