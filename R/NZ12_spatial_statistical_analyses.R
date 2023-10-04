
# Relate cell-level communicability to different drivers

# INPUTS
# - spatial grid: "data/NZ_grid_"
# - species occurrences: "results/model_occurrences_"
# - cell-level communicability: "results/cell_level_communicability_"
# - cell land uses: "data/land_use_frequencies_"
# - species traits: "data/trait_data.csv"

# OUTPUTS
# - cell statistical models (Table 3)
# - correlations between cell-level variables (Fig. S16)
# -------------------------------------------------------------------------
library(mgcv)
library(tidyverse)
library(betareg)
library(DHARMa)
library(performance)
library(corrplot)
library(scales)
library(vegan)
library(broom)
library(broom.mixed)
library(glmmTMB)
library(sf)
nt <- min(parallel::detectCores(),5)

# -------------------------------------------------------------------------
grid.size <- 10
max.dist <- 3*grid.size

cell.comm <- read.csv2(paste("results/cell_level_communicability_",grid.size,"km.csv",sep=""))
cell.land.uses <- read.csv2(paste("data/land_use_frequencies_",grid.size,"km.csv",sep=""))
sp.occurrences <- read.csv2(paste("results/model_occurrences_",grid.size,"km.csv",sep=""))
sp.traits <- read.csv2("data/trait_data.csv")
sp.status <- unique(sp.traits[,c("species","guild","status")])

NZ.grid <- st_read(paste("data/NZ_grid_",grid.size,"km.shp",sep=""))
cell.id <- sort(unique(NZ.grid$cell_id))

NZ.grid.coords <- st_centroid(NZ.grid) %>% 
  cbind(st_coordinates(.)) %>% 
  as.data.frame() %>%
  select(cell_id,X,Y)

# centroid of the cells
NZ.centroid <- st_centroid(NZ.grid) %>% cbind(st_coordinates(.))

# distance matrix between all centroids
NZ.distances <- st_distance(NZ.centroid, NZ.centroid)
NZ.distances <- units::drop_units(NZ.distances)
rownames(NZ.distances) <- cell.id
colnames(NZ.distances) <- cell.id

# convert to KM
NZ.distances <- NZ.distances/1e3

# cell area
cell.area <- grid.size^2

# -------------------------------------------------------------------------
# obtain several response variables
# 1 - plant and bird richness per cell
cell.richness <- sp.occurrences %>% 
  left_join(sp.status) %>% 
  group_by(cell_id,guild) %>%
  summarise(richness = sum(presence))
cell.richness.wide <- pivot_wider(cell.richness, names_from = guild,values_from = richness)
names(cell.richness.wide)[2:3] <- c("bird_richness","plant_richness")

# 2 - percentage of forest and percentage of shrubland
cell.land.uses.wide <- cell.land.uses %>% 
  filter(habitat_type %in% c("forest","shrubland")) %>%
  pivot_wider(names_from = habitat_type,values_from = frequency)
names(cell.land.uses.wide)[1] <- "cell_id"
names(cell.land.uses.wide)[2:3] <- c("forest_percentage","shrubland_percentage")

# 3 - shannon diversity of habitat types per cell
cell.shannon.div <- cell.land.uses %>%
  group_by(cell.id) %>%
  summarise(shannon_diversity_habitats = diversity(frequency))
names(cell.shannon.div)[1] <- "cell_id"

# 4 - proportion of native bird/plant species
native.sp <- sp.occurrences %>%
  left_join(sp.status) %>%
  group_by(cell_id,guild,status) %>%
  summarise(richness = sum(presence)) %>%
  pivot_wider(names_from = guild,values_from = richness) %>%
  left_join(cell.richness.wide) %>%
  group_by(cell_id,status) %>%
  summarise(birds.proportion = birds/bird_richness,
            plants.proportion = plants/plant_richness) %>%
  filter(status == "Native") %>%
  rename(native_bird_freq = birds.proportion,
         native_plant_freq = plants.proportion) %>%
  select(-status)
native.sp$native_bird_freq[which(is.nan(native.sp$native_bird_freq))] <- 0
native.sp$native_plant_freq[which(is.nan(native.sp$native_plant_freq))] <- 0

# 5 land-uses of neighbouring cells
neigh.land.use <- list()
cell.land.use.area <- cell.land.uses %>%
  mutate(land.use.area = frequency * cell.area)

for(i.id in 1:length(cell.id)){
  my.cell.distances <- NZ.distances[cell.id[i.id],]
  my.neigh.ids <- which(my.cell.distances <= max.dist)
  my.neigh.ids <- my.neigh.ids[-which(my.neigh.ids == cell.id[i.id])]
  my.neigh.land.uses <- subset(cell.land.use.area,cell.id %in% my.neigh.ids)
  total.neigh.area <- length(my.neigh.ids) * cell.area
  neigh.land.use.area <- my.neigh.land.uses %>%
    group_by(habitat_type) %>%
    summarise(neigh.area = sum(land.use.area)/total.neigh.area)
  
  neigh.land.use[[length(neigh.land.use)+1]] <- data.frame(cell_id = cell.id[i.id],
                                                           neigh_forest_area = neigh.land.use.area$neigh.area[neigh.land.use.area$habitat_type == "forest"],
                                                           neigh_shrubland_area = neigh.land.use.area$neigh.area[neigh.land.use.area$habitat_type == "shrubland"],
                                                           neigh_diversity_habitats = diversity(neigh.land.use.area$neigh.area))
  
  
}

neigh.land.use.df <- bind_rows(neigh.land.use)

# -------------------------------------------------------------------------

cell.response.data <- left_join(cell.richness.wide,cell.land.uses.wide)
cell.response.data <- left_join(cell.response.data,cell.shannon.div)
cell.response.data <- left_join(cell.response.data,native.sp)

# in the last version, neighbour variables are not used
# because they are highly correlated with cell-level variables
# cell.response.data <- left_join(cell.response.data,neigh.land.use.df)
# -------------------------------------------------------------------------

# merge data and maintain only the communicability metric that we use as response
cell.data <- left_join(cell.response.data,cell.comm) %>%
  rename(comm = scaled.weighted.communicability, comm2 = weighted.communicability) %>%
  select(!c(bin.communicability,scaled.bin.communicability))

# cell.data$comm[which(cell.data$comm == 0)] <- 1e-5
# cell.data$comm[which(cell.data$comm == 1)] <- 0.99999

cell.data <- left_join(cell.data,NZ.grid.coords)

# -------------------------------------------------------------------------
# visualizations
# correlations among response variables
var.matrix <- as.matrix(cell.response.data[,2:ncol(cell.response.data)])

testRes = cor.mtest(var.matrix, conf.level = 0.95)

var.cor <- cor(var.matrix,method = "spearman")

# png(height=600, width=600, file="results/images/cell_variables_correlations.png")
# 
# corrplot(var.cor, p.mat = testRes$p, method = 'circle', 
#          # type = 'lower', 
#          insig='blank',
#          order = 'AOE', diag = FALSE)$corrPos -> p1
# text(p1$x, p1$y, round(p1$corr, 2))
# 
# dev.off()

# -------------------------------------------------------------------------
# relationship between body mass distributions and land-use categories
# (figure not included in the main text or supp mat)

cell.prevalent.land.use <- cell.land.uses %>% 
  rename(cell_id = cell.id) %>%
  filter(!habitat_type %in% c("marine_benthic","marine_intertidal","marine_neritic")) %>%
  group_by(cell_id) %>%
  filter(frequency == max(frequency)) %>%
  select(-frequency)
cell.prevalent.land.use$habitat_type[which(cell.prevalent.land.use$habitat_type %in% 
                                             c("anthropic_habitats","grassland"))] <- "crops and pasture"
cell.prevalent.land.use$habitat_type[which(cell.prevalent.land.use$habitat_type %in% 
                                             c("rocky_area"))] <- "rocky areas"
cell.prevalent.land.use$habitat_type[which(cell.prevalent.land.use$habitat_type %in% 
                                             c("alpine_areas"))] <- "alpine areas"

cell.lu.species <- left_join(cell.prevalent.land.use,sp.occurrences) %>%
  filter(presence == 1) %>%
  left_join(subset(sp.traits, trait == "body.mass" & guild == "birds")) %>%
  drop_na(guild) %>%
  select(cell_id,habitat_type,species,mean.value) %>%
  rename(mean_body_mass = mean.value) %>%
  drop_na(mean_body_mass)

max.bm <- 1800

land.use.body.mass.plot <- ggplot(subset(cell.lu.species, mean_body_mass < max.bm)) +
  geom_histogram(aes(x = mean_body_mass, fill = habitat_type), color = "black",bins = 100) +
  facet_grid(habitat_type~., scales = "free_y") +
  scale_fill_manual(name = "Land use",values = c("grey80","darkkhaki",
                                                 "darkolivegreen",
                                                 "grey20",
                                                 "darkolivegreen3",
                                                 "deepskyblue4")) +
  theme_bw() +
  theme(strip.background = element_blank(), legend.position = "none")+
  xlab("mean bird body mass (g)") +
  NULL
# land.use.body.mass.plot

# -------------------------------------------------------------------------

cell.gam.spatial <- gam(comm2 ~ 
                          # scale(bird_richness) +
                          # scale(plant_richness) +
                          scale(forest_percentage) +
                          scale(shrubland_percentage) +
                          scale(shannon_diversity_habitats) +
                          # scale(neigh_forest_area) +
                          # scale(neigh_shrubland_area) +
                          # scale(neigh_diversity_habitats) +
                          s(X, Y),
                        family = tw(link='log'), 
                        data = cell.data)
summary(cell.gam.spatial)
# DHARMa::testResiduals(DHARMa::simulateResiduals(cell.gam.spatial))
# performance::check_model(cell.gam.spatial)

# print(xtable::xtable(as.data.frame(broom::tidy(cell.gam.spatial)),
#                      floating=FALSE,
#                      digits = 3,
#                      latex.environments=NULL,
#                      booktabs=FALSE))

# -------------------------------------------------------------------------
# neighbor-level variables

cell.gam.spatial.neigh <- gam(comm2 ~ 
                          # scale(bird_richness) +
                          # scale(plant_richness) +
                          # scale(forest_percentage) +
                          # scale(shrubland_percentage) +
                          # scale(shannon_diversity_habitats) +
                          scale(neigh_forest_area) +
                          scale(neigh_shrubland_area) +
                          scale(neigh_diversity_habitats) +
                          s(X, Y),
                        family = tw(link='log'), 
                        data = cell.data)
summary(cell.gam.spatial.neigh)
# DHARMa::testResiduals(DHARMa::simulateResiduals(cell.gam.spatial))
# performance::check_model(cell.gam.spatial)

# print(xtable::xtable(as.data.frame(broom::tidy(cell.gam.spatial.neigh)),
#                      floating=FALSE,
#                      digits = 3,
#                      latex.environments=NULL,
#                      booktabs=FALSE))

