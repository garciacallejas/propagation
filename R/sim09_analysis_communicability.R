
# statistical analyses on communicability of simulated landscapes

# INPUTS
# - individual communicability dataframes: "results/communicability/.."
# - categories of the different factors

# OUTPUTS
# Fig 2, 3
# numerical correlations

# -------------------------------------------------------------------------
library(tidyverse)
library(ggridges)
library(patchwork)
library(lmerTest)
library(glmmTMB)
library(ggpointdensity)
# range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# -------------------------------------------------------------------------
# which data to analyse? correlated or uncorrelated degree-dispersal?
correlated.deg.disp <- F

# -------------------------------------------------------------------------
# read general data
param <- read.csv2("results/sim_landscape_matrices/parameters_v3.csv")

network.categories <- read.csv2("results/network_gradient_categories.csv")
landscape.categories <- read.csv2("results/spatial_autocorrelation_categories.csv")
dispersal.categories <-  read.csv2("results/dispersal_kernels.csv")
cell.distances <- read.csv2("results/cell_distances.csv")
cell.distances$cell_from <- as.character(cell.distances$cell_from)
cell.distances$cell_to <- as.character(cell.distances$cell_to)

if(correlated.deg.disp){
  netcom.df <- read.csv2("results/sim_correlated_network_level_communicability.csv") 
}else{
  netcom.df <- read.csv2("results/sim_network_level_communicability.csv")
}

species.data.df <- read.csv2("results/sim_species_level_communicability.csv")

# to assign species-level metrics: degree and dispersal distances
local.networks <- read.csv2("results/sim_degree_dist_networks.csv")
local.dispersal <- read.csv2("results/dispersal_kernels.csv")

# species presences/absences
load("results/presence_dataframes.RData")

# -------------------------------------------------------------------------
# recover factors

netcom.df <- left_join(netcom.df,network.categories)
netcom.df <- left_join(netcom.df,landscape.categories)

dispersal.exp.rates <- unique(dispersal.categories[,c(2:3)])
netcom.df <- left_join(netcom.df,dispersal.exp.rates)

netcom.df$landscape.category <- as.factor(netcom.df$spatial.autocorr.value)
netcom.df$network.category <- as.factor(netcom.df$poisson.lambda)
netcom.df$dispersal.category <- factor(netcom.df$exponential.rate,
                                       levels = sort(unique(netcom.df$exponential.rate),decreasing = T))

# network.categories <- network.categories$network.category
# landscape.categories <- landscape.categories$landscape.category
# dispersal.categories <- unique(dispersal.categories$dispersal.category)
# network.categories <- network.categories$network.category
# landscape.categories <- landscape.categories$landscape.category
# dispersal.categories <- unique(dispersal.categories$dispersal.category)

replicates <- param$num.category.replicates

richness <- param$richness
cells <- param$ncol * param$nrow

# -------------------------------------------------------------------------
netcom.df$scaled.communicability <- scales::rescale(netcom.df$raw.communicability)

# test
netcom.df$scaled.communicability[netcom.df$scaled.communicability == 0] <- 1e-6

# land.dist <- ggplot(netcom.df, aes(x = scaled.communicability, y = landscape.category)) + 
#   geom_density_ridges()
# net.dist <- ggplot(netcom.df, aes(x = scaled.communicability, y = network.category)) + 
#   geom_density_ridges()
# disp.dist <- ggplot(netcom.df, aes(x = scaled.communicability, y = dispersal.category)) + 
#   geom_density_ridges()

netcom.df.long <- pivot_longer(netcom.df, dispersal.category:landscape.category,
                               names_to = "category.type",
                               values_to = "category.level")

netcom.df.long$category.type <- as.factor(netcom.df.long$category.type)
netcom.df.long$category.level <- as.factor(netcom.df.long$category.level)

# for some reason this does not work
# comm.boxplot <- ggplot(netcom.df.long, aes(x = category.level, y = scaled.communicability)) + 
#   geom_boxplot() + 
#   facet_grid(.~category.type, drop = TRUE) + 
#   NULL
# comm.boxplot

net.plot <- ggplot(subset(netcom.df.long, category.type == "network.category"), 
                   aes(x = category.level, y = scaled.communicability)) + 
  geom_boxplot() +
  theme_bw() +
  labs(x = "network connectedness (\u03bb)", 
       y = "") +
  ggtitle("a)") +
  NULL
# net.plot

disp.plot <- ggplot(subset(netcom.df.long, category.type == "dispersal.category"), 
                    aes(x = category.level, y = scaled.communicability)) + 
  geom_boxplot() +
  theme_bw() +
  labs(x = "dispersal rate (d)", y = "landscape communicability") +
  ggtitle("b)") +
  NULL
# disp.plot

landscape.plot <- ggplot(subset(netcom.df.long, category.type == "landscape.category"), 
                    aes(x = category.level, y = scaled.communicability)) + 
  geom_boxplot() +
  theme_bw() +
  labs(x = "landscape configuration (c)", y = "") +
  ggtitle("c)") +
  NULL
# landscape.plot


combined.plot <- net.plot + disp.plot + landscape.plot

# -------------------------------------------------------------------------
# ggsave("results/images/simulation_categories.pdf",
#        plot = combined.plot,
#        device = cairo_pdf,
#        width = 8, height = 4,dpi = 300)

# -------------------------------------------------------------------------

# simulation results should not be analysed with formal statistical tests, 
# as they can be tricked by increasing power. Correlations are appropriate.

# In addition, I was not able to find a well-specified model for the species-level
# analysis, after trying several ideas. But the qualitative trends seem robust.

# land.dist/net.dist/disp.dist

# m1 <- lm(log(scaled.communicability) ~ landscape.category + network.category + dispersal.category, 
#          data = netcom.df)

# kolmogorov-smirnov test is significant, but the visual patterns are ok
# DHARMa::testResiduals(m1)
# summary(m1)

# -------------------------------------------------------------------------

sp.data.scaled <- species.data.df %>%
  mutate(across(degree:presences,scales::rescale)) %>%
  mutate(id = paste(network.category,landscape.category,dispersal.category,sep="")) %>%
  select(id,sp,wc,bc,degree:presences)

sp.data.scaled.2 <- sp.data.scaled %>% 
  group_by(id) %>%
  mutate(scaled.wc = scales::rescale(wc))


# -------------------------------------------------------------------------

correlations.df <- data.frame(id = unique(sp.data.scaled.2$id), 
                               dispersal.distance = NA, 
                               degree = NA,
                               presences = NA)

for(i in 1:nrow(correlations.df)){
  my.data <- subset(sp.data.scaled.2,id == correlations.df$id[i])
  
  correlations.df$dispersal.distance[i] <- cor.test(my.data$scaled.wc,my.data$dispersal.distance,method = "spearman")$estimate
  correlations.df$degree[i] <- cor.test(my.data$scaled.wc,my.data$degree,method = "spearman")$estimate
  correlations.df$presences[i] <- cor.test(my.data$scaled.wc,my.data$presences,method = "spearman")$estimate
  
}
cor.df.long <- pivot_longer(correlations.df,dispersal.distance:presences,names_to = "factor",values_to = "rho")

cor.df.long$factor[cor.df.long$factor == "dispersal.distance"] <- "dispersal\n distance"
cor.df.long$factor[cor.df.long$factor == "presences"] <- "number of\npresences"

corr.dist.plots <- ggplot(cor.df.long, aes(x = factor, y = rho)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "", y = expression(rho)) +
  NULL
# corr.dist.plots

sp.data.long <- sp.data.scaled.2 %>% pivot_longer(degree:presences,names_to = "factor",values_to = "value")

sp.data.long$factor[sp.data.long$factor == "dispersal.distance"] <- "dispersal\n distance"
sp.data.long$factor[sp.data.long$factor == "presences"] <- "number of\n presences"

corr.raw.plots <- ggplot(sp.data.long, aes(x = value, y = scaled.wc)) +
  geom_pointdensity() +
  scale_color_continuous(type = "viridis",name = "number of\nobservations") +
  # geom_hex(bins = 1000) +
  # geom_point(alpha = .4) +
  facet_grid(factor~.,scales = "free_y") +
  # scale_fill_continuous(type = "viridis") +
  ylab("scaled weighted communicability") +
  xlab("scaled value") +
  theme_bw() +
  theme(strip.background = element_blank())+
  NULL
# corr.raw.plots

# -------------------------------------------------------------------------

cor.complete <- corr.dist.plots + corr.raw.plots + plot_layout(widths = c(1, 2))

ggsave("results/images/simulations/species_level_communicability_simulations.pdf",
       plot = cor.complete,
       device = cairo_pdf,
       width = 9, height = 4,dpi = 300)

# sp.data.scaled.2$scaled.wc[sp.data.scaled.2$scaled.wc == 0] <- 1e-7
# sp.data.scaled.2$presences[sp.data.scaled.2$presences == 0] <- 1e-7
# 
# t1 <- subset(sp.data.scaled.2,id == sp.data.scaled.2$id[86])

# -------------------------------------------------------------------------

# m.glm <- glm(scaled.wc ~ degree + dispersal.distance + log(presences),
#              data = t1,family = Gamma(link = "log"))
# 
# m.glmm <- glmmTMB(scaled.wc ~ log(presences) + degree + dispersal.distance + (1|id),
#                   data = sp.data.scaled.2,family = Gamma(link = "log"))
# 
# # mrand <- lmer(wc ~ presences + degree + dispersal.distance + (1|id),
# #               data = sp.data.scaled)
# 
# simulationOutput <- DHARMa::simulateResiduals(fittedModel = m.glmm)
# plot(simulationOutput)
# DHARMa::testResiduals(m.glmm)

# sp2 <- subset(sp.data.scaled,degree < 0.6 & presences < 0.75)

# ggplot(t1,aes(y = scaled.wc, x = presences)) +
#   geom_point() +
#   geom_smooth() +
#   NULL
# 
# mg <- glm(wc ~ presences,
#   data = sp2, family = Gamma(link = "log"))
# 
# simulationOutput <- DHARMa::simulateResiduals(fittedModel = mg)
# plot(simulationOutput)
# DHARMa::plotResiduals(simulationOutput, species.data.df$degree)
# 
# DHARMa::testResiduals(mg)
# summary(mg)
