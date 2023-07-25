
# Relate cell-level communicability to different drivers
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

cell.comm <- read.csv2(paste("results/cell_level_communicability_",grid.size,"km.csv",sep=""))
cell.land.uses <- read.csv2(paste("data/land_use_frequencies_",grid.size,"km.csv",sep=""))
sp.occurrences <- read.csv2(paste("results/model_occurrences_",grid.size,"km.csv",sep=""))
sp.traits <- read.csv2("data/trait_data.csv")
sp.status <- unique(sp.traits[,c("species","guild","status")])

NZ.grid <- st_read(paste("data/NZ_grid_",grid.size,"km.shp",sep=""))

NZ.grid.coords <- st_centroid(NZ.grid) %>% 
  cbind(st_coordinates(.)) %>% 
  as.data.frame() %>%
  select(cell_id,X,Y)

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

cell.response.data <- left_join(cell.richness.wide,cell.land.uses.wide)
cell.response.data <- left_join(cell.response.data,cell.shannon.div)

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

# -------------------------------------------------------------------------

# cell.glm <- betareg(comm ~ scale(bird_richness) + 
#                       scale(plant_richness) + 
#                       scale(forest_percentage) + 
#                       scale(shrubland_percentage) + 
#                       scale(shannon_diversity_habitats),data = cell.data)

# cell.glm <- glmmTMB(comm ~ scale(bird_richness) + 
#                       scale(plant_richness) +
#                       scale(forest_percentage) +
#                       scale(shrubland_percentage) +
#                       scale(shannon_diversity_habitats), 
#                     family = beta_family,
#                     control = glmmTMBControl(parallel = nt),
#                     data = cell.data)

# summary(cell.glm)
# DHARMa::testResiduals(DHARMa::simulateResiduals(cell.glm))
# sims <- simulateResiduals(cell.glm)
# testSpatialAutocorrelation(sims, x = cell.data$X, y = cell.data$Y, plot = FALSE)
# this data clearly has spatial autocorrelation

# second model
# as per https://cran.r-project.org/web/packages/glmmTMB/vignettes/covstruct.html
# a "position" and a dummy group need to be created
cell.data$pos <- numFactor(cell.data$X, cell.data$Y)
cell.data$group <- factor(rep(1, nrow(cell.data)))

cell.data.subset <- cell.data[sample(1:nrow(cell.data),200,replace = F),]

# I assume exponentially decaying spatial autocorrelation - if I'm not mistaken, "exp"
# not working, it does not converge
# cell.glm.spatial <- glmmTMB(comm2 ~ 
#                               # scale(bird_richness) +
#                               # scale(plant_richness) +
#                               scale(forest_percentage) +
#                               scale(shrubland_percentage) +
#                               scale(shannon_diversity_habitats) +
#                               exp(pos + 0 | group)  ,
#                             control = glmmTMBControl(parallel = nt),
#                             family = tweedie(link = "log"),
#                             data = cell.data.subset)
# summary(cell.glm.spatial)
# DHARMa::testResiduals(DHARMa::simulateResiduals(cell.glm.spatial))

cell.gam.spatial <- gam(comm2 ~ 
                          # scale(bird_richness) +
                          # scale(plant_richness) +
                          scale(forest_percentage) +
                          scale(shrubland_percentage) +
                          scale(shannon_diversity_habitats) + 
                          s(X, Y),
                        family = tw(link='log'), 
                        data = cell.data)
summary(cell.gam.spatial)
DHARMa::testResiduals(DHARMa::simulateResiduals(cell.gam.spatial))
# performance::check_model(cell.gam.spatial)

print(xtable::xtable(as.data.frame(broom::tidy(cell.gam.spatial)),
                     floating=FALSE,
                     digits = 3,
                     latex.environments=NULL,
                     booktabs=FALSE))

# -------------------------------------------------------------------------
