
# Relate species communicability with traits

# INPUTS
# - spatial grid: "data/NZ_grid_"
# - species occurrences: "results/model_occurrences_"
# - population-level degrees: "results/population_level_degrees_"
# - population-level communicability: "results/population_level_communicability_"
# - species-level communicability: "results/species_level_communicability_"
# - cell land uses: "data/land_use_frequencies_"
# - species traits: "data/trait_data.csv"
# - species interactions: "results/plant_bird_interactions_clean.csv"

# OUTPUTS
# statistical models (Tables 1,2, S1,S2)
# correlations between species-level variables (Figs. S14,S15)

# -------------------------------------------------------------------------

library(glmmTMB)
library(tidyverse)
library(betareg)
library(DHARMa)
library(performance)
library(corrplot)
library(scales)
library(broom.mixed)
library(effects)
# library(sf)

# -------------------------------------------------------------------------
grid.size <- 10

sp.int.orig <- read.csv2("results/plant_bird_interactions_clean.csv")
all.local.degrees <- read.csv2(paste("results/population_level_degrees_",grid.size,"km.csv",sep=""))
pop.comm <- read.csv2(paste("results/population_level_communicability_",grid.size,"km.csv",sep=""))
sp.comm <- read.csv2(paste("results/species_level_communicability_",grid.size,"km.csv",sep=""))
sp.traits <- read.csv2("data/trait_data.csv")
sp.traits.2 <- sp.traits

# fill by hand some body mass values that are not present in the gathered dataset
# source: https://en.wikipedia.org/wiki/North_Island_saddleback
sp.traits.2$mean.value[which(sp.traits.2$species == "Philesturnus_rufusater" & 
                               sp.traits.2$trait == "body.mass")] <- 75 
# source: https://en.wikipedia.org/wiki/K%C4%81k%C4%81p%C5%8D
sp.traits.2$mean.value[which(sp.traits.2$species == "Strigops_habroptila" & 
                               sp.traits.2$trait == "body.mass")] <- 2000

# source: https://en.wikipedia.org/wiki/Common_myna
sp.traits.2$mean.value[which(sp.traits.2$species == "Acridotheres_tristis" & 
                               sp.traits.2$trait == "body.mass")] <- 110

# source: https://en.wikipedia.org/wiki/Silvereye
sp.traits.2$mean.value[which(sp.traits.2$species == "Zosterops_lateralis" & 
                               sp.traits.2$trait == "body.mass")] <- 10

# source: body mass data, where it appears as "Carduelis chloris"
sp.traits.2$mean.value[which(sp.traits.2$species == "Chloris_chloris" & 
                               sp.traits.2$trait == "body.mass")] <- 27.5

# source: https://en.wikipedia.org/wiki/Common_redpoll
sp.traits.2$mean.value[which(sp.traits.2$species == "Acanthis_flammea" & 
                               sp.traits.2$trait == "body.mass")] <- 14

# source: https://www.birdsnz.org.nz/wp-content/uploads/2021/12/Notornis_45_1_49.pdf
sp.traits.2$mean.value[which(sp.traits.2$species == "Cyanoramphus_malherbi" & 
                               sp.traits.2$trait == "body.mass")] <- 45.3

# source: file:///home/david/Downloads/moorebattleyNotornis2003.pdf
sp.traits.2$mean.value[which(sp.traits.2$species == "Anas_chlorotis" & 
                               sp.traits.2$trait == "body.mass")] <- 500

# measurements of callaeas wilsoni, a sister species
# source: https://www.birdsnz.org.nz/wp-content/uploads/2021/12/Notornis_48_4_217.pdf
sp.traits.2$mean.value[which(sp.traits.2$species == "Callaeas_cinereus" & 
                               sp.traits.2$trait == "body.mass")] <- 225


sp.traits.2$sd.value <- sp.traits.2$n <- NULL
sp.traits.2$trait[which(sp.traits.2$trait == "Hand-wing.Index")] <- "HWI"
sp.traits.wide <- pivot_wider(sp.traits.2,names_from = trait,values_from = mean.value)
sp.list <- sort(unique(sp.traits$species))
sp.status <- unique(sp.traits[,c("species","guild","status")])

# which source to use?
observation.source <- "occurrence_prob"

if(observation.source == "observed_records"){
  sp.obs <- read.csv2(paste("data/sp_observations_long_",grid.size,"km.csv",
                            sep=""))
}else if(observation.source == "occurrence_prob"){
  sp.obs <- read.csv2(paste("results/model_occurrences_",grid.size,"km.csv",
                            sep=""))
  # this is to harmonise datasets
  names(sp.obs)[which(names(sp.obs) == "presence")] <- "observations"
  # keep presences only
  sp.obs <- subset(sp.obs, observations == 1)
  
}

all.sp <- sort(unique(pop.comm$species))

# -------------------------------------------------------------------------
# how many native and exotic plants/birds?
sp.traits.valid.sp <- subset(sp.traits.wide,species %in% all.sp)
table(sp.traits.valid.sp$guild,sp.traits.valid.sp$status)

# -------------------------------------------------------------------------
# average local degree
# this is the average of the realised degrees of every species in their local
# communities
avg.local.deg <- all.local.degrees %>%
  group_by(species) %>%
  summarise(avg.local.degree = mean(deg))

# -------------------------------------------------------------------------
# species degree in the metaweb

sp.int <- subset(sp.int.orig, PLANTSPECIES %in% all.sp & 
                   BIRDSPECIES %in% all.sp) %>%
  select(BIRDSPECIES,PLANTSPECIES) %>%
  unique()

# connectance?
# nrow(sp.int)/(length(unique(sp.int$PLANTSPECIES)) * length(unique(sp.int$BIRDSPECIES)))

plant.deg <- sp.int %>% 
  group_by(PLANTSPECIES) %>%
  summarise(deg = n()) %>%
  rename(species = PLANTSPECIES)

bird.deg <- sp.int %>%
  group_by(BIRDSPECIES) %>%
  summarise(deg = n()) %>%
  rename(species = BIRDSPECIES)

sp.degree <- bind_rows(plant.deg,bird.deg)

# species presences
sp.presence.absences <- expand.grid(cell_id = unique(sp.obs$cell_id),
                                    species = all.sp) %>%
  left_join(sp.obs) %>% 
  replace_na(list(observations = 0)) %>%
  rename(presence = observations) %>%
  arrange(cell_id)

# prevalence per species
sp.presences <- sp.presence.absences %>%
  group_by(species) %>%
  summarise(prevalence = sum(presence))

# population-level averages
pop.averages <- pop.comm %>% 
  left_join(sp.presence.absences) %>%
  filter(presence == 1) %>%
  group_by(species,guild) %>%
  summarise(avg.comm = mean(population.weighted.communicability),
            sd.comm = sd(population.weighted.communicability)) %>%
  ungroup() %>%
  mutate(scaled.avg.comm = scales::rescale(avg.comm))

# -------------------------------------------------------------------------
# check population-level and species-level communicability

pop.data <- left_join(sp.traits.wide,pop.averages) %>%
  rename(comm = scaled.avg.comm) %>%
  subset(!is.na(status)) %>%
  left_join(avg.local.deg)

plant.pop.data <- subset(pop.data,guild == "plants") %>% 
  select(species,comm,
         status,
         avg.local.degree,
         FRUIT_DIAMETER_mm,
         # FRUIT_LENGTH_mm,
         PLANT_MAX_MEAN_VEGETATIVE_HEIGHT_m) %>%
  drop_na()
bird.pop.data <- subset(pop.data, guild == "birds") %>% 
  select(species,comm,
         status,
         avg.local.degree,
         # Beak.Depth, Beak.Length_Culmen,
         # Beak.Length_Nares, 
         # Beak.Width,
         BILL_LENGTH_mm,
         # BILL_WIDTH_mm,
         body.mass,
         HWI,
         # Kipps.Distance,
         # Tail.Length,Tarsus.Length,Wing.Length
  ) %>%
  drop_na()

# merge data and maintain only the communicability metric that we use as response
sp.data <- left_join(sp.traits.wide,sp.comm) %>%
  rename(comm = scaled.weighted.communicability,comm2 = weighted.communicability) %>%
  subset(!is.na(status)) %>%
  left_join(sp.presences) %>%
  left_join(sp.degree)

plant.data <- subset(sp.data,guild == "plants") %>% 
  select(species,comm,comm2,status,
         prevalence,
         deg,
         FRUIT_DIAMETER_mm,
         FRUIT_LENGTH_mm,
         PLANT_MAX_MEAN_VEGETATIVE_HEIGHT_m) %>%
  drop_na()
bird.data <- subset(sp.data, guild == "birds") %>% 
  select(species,comm,comm2,status,
         prevalence,
         deg,
         # Beak.Depth, Beak.Length_Culmen,
         # Beak.Length_Nares, 
         # Beak.Width,
         BILL_LENGTH_mm,
         # BILL_WIDTH_mm,
         body.mass,
         HWI,
         # Kipps.Distance,
         # Tail.Length,Tarsus.Length,Wing.Length
         ) %>%
  drop_na()

# -------------------------------------------------------------------------
# visualizations
# correlations among response variables
all.bird.traits <- sp.data %>% select(Beak.Depth:WING_LENGTH_cm) %>% drop_na()
all.plant.traits <- sp.data %>% 
  select(FRUIT_DIAMETER_mm:PLANT_MAX_MEAN_VEGETATIVE_HEIGHT_m) %>% drop_na()
names(all.plant.traits)[3] <- "PLANT_HEIGHT_m"

bird.traits.cor <- cor(all.bird.traits)
# corrplot(bird.traits.cor)

plant.traits.cor <- cor(all.plant.traits)
# corrplot(plant.traits.cor)

# boxplot - communicability by guild and status
# boxplot.status.guild <- ggplot(sp.data,aes(x = guild, y = comm)) + 
#   geom_boxplot(aes(fill = status)) + 
#   NULL
# boxplot.status.guild

# -------------------------------------------------------------------------
# correlations among the interesting variables

names(bird.data)[which(names(bird.data) == "deg")] <- "metaweb_degree"
names(bird.data)[which(names(bird.data) == "body.mass")] <- "body_mass"
names(bird.data)[which(names(bird.data) == "BILL_LENGTH_mm")] <- "bill_length"

# birds
bird.matrix <- as.matrix(bird.data[,c("prevalence","metaweb_degree","bill_length","body_mass","HWI")])
bird.cor.test <- cor.mtest(bird.data[,c("prevalence","metaweb_degree","bill_length","body_mass","HWI")], conf.level = 0.95)
bird.cor <- cor(bird.matrix, method = "spearman")

# 1
# corrplot(bird.cor)

# 2 - only significant correlations, with p-values
# grDevices::cairo_pdf(height=800, width=800, file="results/images/bird_correlations.pdf")
# png(height=400, width=400, file="results/images/bird_correlations.png")
# 
# corrplot(bird.cor, p.mat = bird.cor.test$p,
#          method = 'circle',
#          # type = 'upper',
#          insig='blank',
#          order = 'AOE', diag = FALSE)$corrPos -> p1
# text(p1$x, p1$y, round(p1$corr, 2))
# 
# dev.off()

# plants
plant.cor.data <- plant.data[,c("prevalence","deg","FRUIT_DIAMETER_mm","PLANT_MAX_MEAN_VEGETATIVE_HEIGHT_m")]
names(plant.cor.data) <- c("prevalence","metaweb_degree","fruit_diam","max_mean_height")

plant.matrix <- as.matrix(plant.cor.data)
plant.cor.test <- cor.mtest(plant.matrix, conf.level = 0.95)
plant.cor <- cor(plant.matrix, method = "spearman")

# 1
# corrplot(plant.cor)

# 2 - only significant correlations, with p-values
# png(height=400, width=400, file="results/images/plant_correlations.png")
# 
# corrplot(plant.cor, p.mat = plant.cor.test$p,
#          method = 'circle',
#          # type = 'upper',
#          insig='blank',
#          order = 'AOE', diag = FALSE)$corrPos -> p1
# text(p1$x, p1$y, round(p1$corr, 2))
# dev.off()
# -------------------------------------------------------------------------
# population-level glmm
plant.pop.data.scaled <- plant.pop.data %>%
  mutate(avg.local.degree.scaled = scale(avg.local.degree),
         fruit.diam.scaled = scale(FRUIT_DIAMETER_mm),
         max.height.scaled = scale(PLANT_MAX_MEAN_VEGETATIVE_HEIGHT_m))

plant.pop.glm.scaled <- glmmTMB(comm ~ status*fruit.diam.scaled + 
                                  avg.local.degree.scaled + 
                                  max.height.scaled,
                                family = Gamma(link="log"), data = plant.pop.data.scaled)

# plant.pop.glm.scaled.simple <- glmmTMB(comm ~ status + fruit.diam.scaled + 
#                                   avg.local.degree.scaled + 
#                                   max.height.scaled,
#                                 family = Gamma(link="log"), data = plant.pop.data.scaled)

# plant.pop.variable.importance <- glmm.hp::glmm.hp(plant.pop.glm.scaled.simple)

# summary(plant.pop.glm.scaled)
# DHARMa::testResiduals(DHARMa::simulateResiduals(plant.pop.glm.scaled))
# plot(effects::allEffects(plant.pop.glm.scaled))
# plot(effect(term = "status:fruit.diam.scaled",mod = plant.pop.glm.scaled))

# print(xtable::xtable(as.data.frame(broom::tidy(plant.pop.glm.scaled)),
#                      floating=FALSE,
#                      digits = 3,
#                      latex.environments=NULL,
#                      booktabs=FALSE))

# -------------------------------------------------------------------------
# cor(bird.data$BILL_LENGTH_mm,bird.data$BILL_WIDTH_mm) # 0.85
# cor(bird.data$BILL_LENGTH_mm,bird.data$body.mass) # 0.66
# cor(bird.data$body.mass,bird.data$HWI) # 0.18

bird.pop.data.scaled <- bird.pop.data %>%
  mutate(body.mass.scaled = scale(body.mass),
         avg.local.degree.scaled = scale(avg.local.degree))

bird.pop.glm.scaled <- glmmTMB(comm ~ status*body.mass.scaled +
                                 avg.local.degree.scaled, 
                               family = Gamma(link = "log"), data = bird.pop.data.scaled)

bird.pop.glm.scaled.simple <- glmmTMB(comm ~ status + body.mass.scaled +
                                 avg.local.degree.scaled, 
                               family = Gamma(link = "log"), data = bird.pop.data.scaled)

# bird.pop.variable.importance <- glmm.hp::glmm.hp(bird.pop.glm.scaled.simple)

# summary(bird.pop.glm.scaled)
# DHARMa::testResiduals(DHARMa::simulateResiduals(bird.pop.glm.scaled))

# print(xtable::xtable(as.data.frame(broom::tidy(bird.pop.glm.scaled)),
#                      floating=FALSE,
#                      digits = 3,
#                      latex.environments=NULL,
#                      booktabs=FALSE))

# -------------------------------------------------------------------------
# species-level glmm

# cor(plant.data$FRUIT_DIAMETER_mm,plant.data$FRUIT_LENGTH_mm) # ~0.85
# cor(plant.data$FRUIT_DIAMETER_mm,plant.data$PLANT_MAX_MEAN_VEGETATIVE_HEIGHT_m) # 0.04

plant.data.scaled <- plant.data %>% 
  mutate(degree.scaled = scale(deg),
         prevalence.scaled = scale(prevalence),
         fruit.diam.scaled = scale(FRUIT_DIAMETER_mm),
         max.height.scaled = scale(PLANT_MAX_MEAN_VEGETATIVE_HEIGHT_m))

# to avoid a single zero, input a very small value, much smaller than the next one
plant.data.scaled$comm[plant.data.scaled$comm == 0] <- 1e-8

plant.sp.glm.scaled <- glmmTMB(comm ~ status*fruit.diam.scaled + 
                                 degree.scaled +
                                 prevalence.scaled +
                                 max.height.scaled,
                               family = Gamma(link="log"),
                               data = plant.data.scaled)


plant.sp.glm.scaled.simple <- glmmTMB(comm ~ status + fruit.diam.scaled + 
                                 degree.scaled +
                                 prevalence.scaled +
                                 max.height.scaled,
                                 family = Gamma(link="log"),
                                 data = plant.data.scaled)

# plant.sp.variable.importance <- glmm.hp::glmm.hp(plant.sp.glm.scaled.simple, 
#                                                  type = "adjR2", commonality = FALSE)

# visual inspection
# ggplot(plant.data.scaled) +
#   geom_point(aes(y = max.height.scaled, x = comm2))
# plot(effects::allEffects(plant.sp.glm.scaled))

# summary(plant.sp.glm.scaled)
# DHARMa::testResiduals(DHARMa::simulateResiduals(plant.sp.glm.scaled))
# plot(effect(term = "status:fruit.diam.scaled",mod = plant.sp.glm.scaled))

# print(xtable::xtable(as.data.frame(broom::tidy(plant.sp.glm.scaled)),
#                      floating=FALSE,
#                      digits = 3,
#                      latex.environments=NULL,
#                      booktabs=FALSE))

# -------------------------------------------------------------------------
# cor(bird.data$BILL_LENGTH_mm,bird.data$BILL_WIDTH_mm) # 0.85
# cor(bird.data$BILL_LENGTH_mm,bird.data$body.mass) # 0.66
# cor(bird.data$body.mass,bird.data$HWI) # 0.18

bird.data.scaled <- bird.data %>%
  mutate(degree.scaled = scale(metaweb_degree),
         prevalence.scaled = scale(prevalence),
         body.mass.scaled = scale(body_mass))

bird.sp.glm.scaled <- glmmTMB(comm ~ status*body.mass.scaled + prevalence.scaled +
                                degree.scaled,
                              family = tweedie(link = "log"), data = bird.data.scaled)

bird.sp.glm.scaled.simple <- glmmTMB(comm ~ status + body.mass.scaled + prevalence.scaled +
                                degree.scaled,
                              family = Gamma(link = "log"), data = bird.data.scaled)

# visual inspection
# ggplot(bird.data.scaled) + 
#   geom_point(aes(y = prevalence.scaled, x = comm2))
# plot(effects::allEffects(bird.sp.glm.scaled))

# bird.sp.variable.importance <- glmm.hp::glmm.hp(bird.sp.glm.scaled.simple, type = "adjR2", commonality = FALSE)

# summary(bird.sp.glm.scaled)
# DHARMa::testResiduals(DHARMa::simulateResiduals(bird.sp.glm.scaled))
# DHARMa::testDispersion(DHARMa::simulateResiduals(bird.sp.glm.scaled))
# plot(effect(term = "status:body.mass.scaled",mod = bird.sp.glm.scaled))

# residuals vs. predictors plot, following Hartig's advice:
# https://stats.stackexchange.com/questions/490680/significant-dispersion-test

# bird.sp.res <- data.frame(res = resid(bird.sp.glm.scaled))
# bird.sp.fit <- data.frame(fit = fitted(bird.sp.glm.scaled))
# bird.sp.data.model <- cbind(bird.data.scaled,bird.sp.res,bird.sp.fit)
# 
# ggplot(bird.sp.data.model, aes(y = res, x = body.mass.scaled)) +
#   geom_point()
# 
# ggplot(bird.sp.data.model, aes(y = res, x = prevalence.scaled)) +
#   geom_point()
# 
# ggplot(bird.sp.data.model, aes(y = res, x = degree.scaled)) +
#   geom_point()
# 
# ggplot(bird.sp.data.model, aes(y = res, x = status)) +
#   geom_point()

# print(xtable::xtable(as.data.frame(broom::tidy(bird.sp.glm.scaled)),
#                      floating=FALSE,
#                      digits = 3,
#                      latex.environments=NULL,
#                      booktabs=FALSE))

# -------------------------------------------------------------------------
# plot effects, in particular the interaction bodymass*status/fruitdiameter*status

# plant population-level communicability
plant.pop.status.pred <- as.data.frame(ggeffects::predict_response(plant.pop.glm.scaled,terms = c("status")))
plant.pop.status.pred$x <- factor(plant.pop.status.pred$x, levels = c("Native","Exotic"))

plant.pop.diam.pred <- as.data.frame(ggeffects::predict_response(plant.pop.glm.scaled,
                                                                terms = c("fruit.diam.scaled [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99]")))
plant.pop.pred <- as.data.frame(ggeffects::predict_response(plant.pop.glm.scaled, 
                                                            c("fruit.diam.scaled [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99]","status")))
plant.pop.pred$group <- factor(plant.pop.pred$group, levels = c("Native","Exotic"))

# plant species-level communicability
plant.sp.status.pred <- as.data.frame(ggeffects::predict_response(plant.sp.glm.scaled,terms = c("status")))
plant.sp.status.pred$x <- factor(plant.sp.status.pred$x, levels = c("Native","Exotic"))

plant.sp.diam.pred <- as.data.frame(ggeffects::predict_response(plant.sp.glm.scaled,
                                                                terms = c("fruit.diam.scaled [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99]")))
plant.sp.pred <- as.data.frame(ggeffects::predict_response(plant.sp.glm.scaled,terms = c("fruit.diam.scaled [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99]","status")))
plant.sp.pred$group <- factor(plant.sp.pred$group, levels = c("Native","Exotic"))

# bird population-level communicability
bird.pop.status.pred <- as.data.frame(ggeffects::predict_response(bird.pop.glm.scaled,terms = c("status")))
bird.pop.status.pred$x <- factor(bird.pop.status.pred$x, levels = c("Native","Exotic"))

bird.pop.size.pred <- as.data.frame(ggeffects::predict_response(bird.pop.glm.scaled,
                                                                 terms = c("body.mass.scaled [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99]")))

bird.pop.pred <- as.data.frame(ggeffects::predict_response(bird.pop.glm.scaled, c("body.mass.scaled [0.1,0.3,0.5,0.7,0.9]","status")))
bird.pop.pred$group <- factor(bird.pop.pred$group, levels = c("Native","Exotic"))

# bird species-level communicability
bird.sp.status.pred <- as.data.frame(ggeffects::predict_response(bird.sp.glm.scaled,terms = c("status")))
bird.sp.status.pred$x <- factor(bird.sp.status.pred$x, levels = c("Native","Exotic"))

bird.sp.size.pred <- as.data.frame(ggeffects::predict_response(bird.sp.glm.scaled,
                                                                terms = c("body.mass.scaled [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99]")))

bird.sp.pred <- as.data.frame(ggeffects::predict_response(bird.sp.glm.scaled, c("body.mass.scaled [0.1,0.3,0.5,0.7,0.9]","status")))
bird.sp.pred$group <- factor(bird.sp.pred$group, levels = c("Native","Exotic"))
# -------------------------------------------------------------------------
# plots

# plant population-level

plant.pop.status.plot <- ggplot(plant.pop.status.pred,aes(y = predicted, x = x)) + 
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, color = x)) +
  geom_point(aes(color = x)) +
  scale_color_OkabeIto(name = "")+#,order = c(2,3,1,4,5))+#c(2:8,1)) +
  theme_bw() +
  xlab("") +
  ylab("population-level communicability") +
  theme(legend.position="none") +
  NULL
# plant.pop.status.plot

plant.pop.diam.plot <- ggplot(plant.pop.diam.pred,aes(y = predicted, x = x)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2) +
  geom_line() +
  # scale_color_OkabeIto(name = "")+#,order = c(2,3,1,4,5))+#c(2:8,1)) +
  theme_bw() +
  xlab("fruit diameter (scaled)") +
  ylab("population-level communicability") +
  NULL
# plant.pop.diam.plot

plant.pop.plot <- ggplot(plant.pop.pred,aes(y = predicted, x = x)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = .2) +
  geom_line(aes(color = group), linewidth = 1.3) +
  scale_color_OkabeIto(name = "")+#,order = c(2,3,1,4,5))+#c(2:8,1)) +
  scale_fill_OkabeIto(guide = "none") +
  theme_bw() +
  xlab("fruit diameter (scaled)") + 
  ylab("population-level communicability") +
  NULL
# plant.pop.plot

# -------------------------------------------------------------------------
# plant species-level

plant.sp.status.plot <- ggplot(plant.sp.status.pred,aes(y = predicted, x = x)) + 
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, color = x)) +
  geom_point(aes(color = x)) +
  scale_color_OkabeIto(name = "")+#,order = c(2,3,1,4,5))+#c(2:8,1)) +
  theme_bw() +
  xlab("") +
  ylab("species-level communicability") +
  theme(legend.position="none") +
  NULL
# plant.sp.status.plot

plant.sp.diam.plot <- ggplot(plant.sp.diam.pred,aes(y = predicted, x = x)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2) +
  geom_line() +
  # scale_color_OkabeIto(name = "")+#,order = c(2,3,1,4,5))+#c(2:8,1)) +
  theme_bw() +
  xlab("fruit diameter (scaled)") +
  ylab("species-level communicability") +
  NULL
# plant.sp.diam.plot

plant.sp.plot <- ggplot(plant.sp.pred,aes(y = predicted, x = x)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = .2) +
  geom_line(aes(color = group), linewidth = 1.3) +
  scale_color_OkabeIto(name = "")+#,order = c(2,3,1,4,5))+#c(2:8,1)) +
  scale_fill_OkabeIto(guide = "none") +
  theme_bw() +
  xlab("fruit diameter (scaled)") + 
  ylab("species-level communicability") +
  NULL
# plant.sp.plot

# -------------------------------------------------------------------------
# bird population-level

bird.pop.status.plot <- ggplot(bird.pop.status.pred,aes(y = predicted, x = x)) + 
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, color = x)) +
  geom_point(aes(color = x)) +
  scale_color_OkabeIto(name = "")+#,order = c(2,3,1,4,5))+#c(2:8,1)) +
  theme_bw() +
  xlab("") +
  ylab("population-level communicability") +
  theme(legend.position="none") +
  NULL
# bird.pop.status.plot

bird.pop.size.plot <- ggplot(bird.pop.size.pred,aes(y = predicted, x = x)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2) +
  geom_line() +
  # scale_color_OkabeIto(name = "")+#,order = c(2,3,1,4,5))+#c(2:8,1)) +
  theme_bw() +
  xlab("body mass (scaled)") +
  ylab("population-level communicability") +
  NULL
# bird.pop.size.plot

bird.pop.plot <- ggplot(bird.pop.pred,aes(y = predicted, x = x)) + 
  # geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = .2) +
  geom_line(aes(color = group), linewidth = 1.3) +
  scale_color_OkabeIto(name = "")+#,order = c(2,3,1,4,5))+#c(2:8,1)) +
  scale_fill_OkabeIto(guide = "none") +
  theme_bw() +
  xlab("body mass (scaled)") + 
  ylab("population-level communicability") +
  # ylim(c(0,0.3)) +
  NULL
# bird.pop.plot

# -------------------------------------------------------------------------
# bird species-level

bird.sp.status.plot <- ggplot(bird.sp.status.pred,aes(y = predicted, x = x)) + 
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, color = x)) +
  geom_point(aes(color = x)) +
  scale_color_OkabeIto(name = "")+#,order = c(2,3,1,4,5))+#c(2:8,1)) +
  theme_bw() +
  xlab("") +
  ylab("species-level communicability") +
  theme(legend.position="none") +
  NULL
# bird.sp.status.plot

bird.sp.size.plot <- ggplot(bird.sp.size.pred,aes(y = predicted, x = x)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2) +
  geom_line() +
  # scale_color_OkabeIto(name = "")+#,order = c(2,3,1,4,5))+#c(2:8,1)) +
  theme_bw() +
  xlab("body mass (scaled)") +
  ylab("species-level communicability") +
  NULL
# bird.sp.size.plot

bird.sp.plot <- ggplot(bird.sp.pred,aes(y = predicted, x = x)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = .2) +
  geom_line(aes(color = group), linewidth = 1.3) +
  scale_color_OkabeIto(name = "")+#,order = c(2,3,1,4,5))+#c(2:8,1)) +
  scale_fill_OkabeIto(guide = "none") +
  theme_bw() +
  xlab("body mass (scaled)") + 
  ylab("species-level communicability") +
  NULL
# bird.sp.plot

# -------------------------------------------------------------------------
ggsave("results/images/pop_plant_interaction_effects.png",
       plot = plant.pop.plot, width = 7, height = 4)

ggsave("results/images/sp_plant_interaction_effects.png",
       plot = plant.sp.plot, width = 7, height = 4)

ggsave("results/images/pop_bird_interaction_effects.png",
       plot = bird.pop.plot, width = 7, height = 4)

ggsave("results/images/sp_bird_interaction_effects.png",
       plot = bird.sp.plot, width = 7, height = 4)
